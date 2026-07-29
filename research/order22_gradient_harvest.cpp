#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int kOrder = 22;
constexpr __int128 kTarget =
    static_cast<__int128>(409600000000000LL);
using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using Inverse =
    std::array<std::array<long double, kOrder>, kOrder>;
using MatrixKey = std::array<std::uint64_t, 8>;
using Wide = __int128;

std::string wide_string(Wide value) {
  if (value == 0) return "0";
  bool negative = value < 0;
  if (negative) value = -value;
  std::string result;
  while (value) {
    result.push_back(static_cast<char>('0' + value % 10));
    value /= 10;
  }
  if (negative) result.push_back('-');
  std::reverse(result.begin(), result.end());
  return result;
}

Wide absolute(Wide value) { return value < 0 ? -value : value; }

Wide determinant(const Matrix& source) {
  std::array<std::array<Wide, kOrder>, kOrder> work{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      work[row][column] = source[row][column];
    }
  }
  Wide previous = 1;
  int sign = 1;
  for (int column = 0; column < kOrder - 1; ++column) {
    int pivot_row = column;
    while (
        pivot_row < kOrder && work[pivot_row][column] == 0) {
      ++pivot_row;
    }
    if (pivot_row == kOrder) return 0;
    if (pivot_row != column) {
      std::swap(work[pivot_row], work[column]);
      sign = -sign;
    }
    const Wide pivot = work[column][column];
    for (int row = column + 1; row < kOrder; ++row) {
      for (int inner = column + 1; inner < kOrder; ++inner) {
        const Wide numerator =
            work[row][inner] * pivot -
            work[row][column] * work[column][inner];
        if (column && numerator % previous) {
          throw std::runtime_error(
              "non-exact Bareiss division");
        }
        work[row][inner] =
            column ? numerator / previous : numerator;
      }
      work[row][column] = 0;
    }
    previous = pivot;
  }
  return sign * work[kOrder - 1][kOrder - 1];
}

bool inverse(const Matrix& matrix, Inverse& result) {
  std::array<
      std::array<long double, 2 * kOrder>,
      kOrder>
      work{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      work[row][column] = matrix[row][column];
      work[row][column + kOrder] =
          row == column ? 1.0L : 0.0L;
    }
  }
  for (int column = 0; column < kOrder; ++column) {
    int pivot_row = column;
    for (int row = column + 1; row < kOrder; ++row) {
      if (std::fabs(work[row][column]) >
          std::fabs(work[pivot_row][column])) {
        pivot_row = row;
      }
    }
    if (std::fabs(work[pivot_row][column]) < 1e-24L) {
      return false;
    }
    if (pivot_row != column) {
      std::swap(work[pivot_row], work[column]);
    }
    const long double pivot = work[column][column];
    for (int inner = 0; inner < 2 * kOrder; ++inner) {
      work[column][inner] /= pivot;
    }
    for (int row = 0; row < kOrder; ++row) {
      if (row == column) continue;
      const long double factor = work[row][column];
      if (factor == 0.0L) continue;
      for (int inner = 0; inner < 2 * kOrder; ++inner) {
        work[row][inner] -= factor * work[column][inner];
      }
    }
  }
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      result[row][column] =
          work[row][column + kOrder];
    }
  }
  return true;
}

Matrix read_matrix(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error(
        "cannot read seed " + path.string());
  }
  Matrix matrix{};
  for (auto& row : matrix) {
    for (int& value : row) {
      if (!(input >> value) || (value != -1 && value != 1)) {
        throw std::runtime_error(
            "invalid seed " + path.string());
      }
    }
  }
  std::string extra;
  if (input >> extra) {
    throw std::runtime_error("extra seed data");
  }
  return matrix;
}

void write_matrix(
    const std::filesystem::path& path,
    const Matrix& matrix) {
  std::ofstream output(path);
  if (!output) {
    throw std::runtime_error(
        "cannot write " + path.string());
  }
  for (const auto& row : matrix) {
    for (int column = 0; column < kOrder; ++column) {
      if (column) output << ' ';
      output << row[column];
    }
    output << '\n';
  }
  output.flush();
  if (!output) {
    throw std::runtime_error(
        "cannot flush " + path.string());
  }
}

MatrixKey matrix_key(const Matrix& matrix) {
  MatrixKey key{};
  int bit = 0;
  for (const auto& row : matrix) {
    for (const int value : row) {
      if (value > 0) {
        key[bit / 64] |= std::uint64_t{1} << (bit % 64);
      }
      ++bit;
    }
  }
  return key;
}

struct Option {
  long double ratio = 0.0L;
  bool column = false;
  int index = -1;
  std::array<int, kOrder> values{};
};

bool ascend(
    Matrix& matrix,
    std::mt19937_64& randomizer,
    bool random_improvement,
    const std::string& axis_mode,
    int& steps) {
  steps = 0;
  for (; steps < 500; ++steps) {
    Inverse inv{};
    if (!inverse(matrix, inv)) return false;
    std::vector<Option> improving;
    Option best;
    if (axis_mode != "columns") {
      for (int row = 0; row < kOrder; ++row) {
        Option option;
        option.index = row;
        for (int column = 0; column < kOrder; ++column) {
          const long double value = inv[column][row];
          option.ratio += std::fabs(value);
          option.values[column] = value >= 0 ? 1 : -1;
        }
        bool same = true;
        bool opposite = true;
        for (int column = 0; column < kOrder; ++column) {
          same &=
              option.values[column] == matrix[row][column];
          opposite &=
              option.values[column] == -matrix[row][column];
        }
        if (
            !same && !opposite &&
            option.ratio > 1.0L + 1e-12L) {
          improving.push_back(option);
          if (option.ratio > best.ratio) best = option;
        }
      }
    }
    if (axis_mode != "rows") {
      for (int column = 0; column < kOrder; ++column) {
        Option option;
        option.column = true;
        option.index = column;
        for (int row = 0; row < kOrder; ++row) {
          const long double value = inv[column][row];
          option.ratio += std::fabs(value);
          option.values[row] = value >= 0 ? 1 : -1;
        }
        bool same = true;
        bool opposite = true;
        for (int row = 0; row < kOrder; ++row) {
          same &=
              option.values[row] == matrix[row][column];
          opposite &=
              option.values[row] == -matrix[row][column];
        }
        if (
            !same && !opposite &&
            option.ratio > 1.0L + 1e-12L) {
          improving.push_back(option);
          if (option.ratio > best.ratio) best = option;
        }
      }
    }
    if (improving.empty()) return true;
    Option chosen = best;
    if (random_improvement) {
      chosen = improving[randomizer() % improving.size()];
    }
    if (chosen.column) {
      for (int row = 0; row < kOrder; ++row) {
        matrix[row][chosen.index] = chosen.values[row];
      }
    } else {
      matrix[chosen.index] = chosen.values;
    }
  }
  return false;
}

struct Arguments {
  std::vector<std::filesystem::path> seeds;
  std::filesystem::path output;
  std::uint64_t random_seed = 1;
  double seconds = 60.0;
  std::uint64_t max_attempts = 0;
  int minimum_kick = 2;
  int maximum_kick = 242;
  int random_start_percent = 10;
  bool random_improvement = false;
  std::string axis_mode = "both";
  bool walk = false;
  std::uint64_t max_saved = 20000;
};

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    auto value = [&]() {
      if (++index == argc) {
        throw std::runtime_error("missing option value");
      }
      return std::string(argv[index]);
    };
    if (option == "--seed") {
      result.seeds.emplace_back(value());
    } else if (option == "--output") {
      result.output = value();
    } else if (option == "--random-seed") {
      result.random_seed = std::stoull(value());
    } else if (option == "--seconds") {
      result.seconds = std::stod(value());
    } else if (option == "--max-attempts") {
      result.max_attempts = std::stoull(value());
    } else if (option == "--minimum-kick") {
      result.minimum_kick = std::stoi(value());
    } else if (option == "--maximum-kick") {
      result.maximum_kick = std::stoi(value());
    } else if (option == "--random-start-percent") {
      result.random_start_percent = std::stoi(value());
    } else if (option == "--random-improvement") {
      result.random_improvement = std::stoi(value()) != 0;
    } else if (option == "--axis") {
      result.axis_mode = value();
    } else if (option == "--walk") {
      result.walk = std::stoi(value()) != 0;
    } else if (option == "--max-saved") {
      result.max_saved = std::stoull(value());
    } else {
      throw std::runtime_error("unknown option " + option);
    }
  }
  if (result.seeds.empty() || result.output.empty()) {
    throw std::runtime_error("provide --seed and --output");
  }
  if (!std::isfinite(result.seconds) || result.seconds <= 0.0) {
    throw std::runtime_error("--seconds must be positive");
  }
  if (
      result.minimum_kick < 0 ||
      result.maximum_kick < result.minimum_kick ||
      result.maximum_kick > kOrder * kOrder) {
    throw std::runtime_error("invalid kick range");
  }
  if (
      result.random_start_percent < 0 ||
      result.random_start_percent > 100) {
    throw std::runtime_error(
        "--random-start-percent must be in [0,100]");
  }
  if (
      result.axis_mode != "both" &&
      result.axis_mode != "rows" &&
      result.axis_mode != "columns") {
    throw std::runtime_error(
        "--axis must be both, rows, or columns");
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    if (
        std::filesystem::exists(arguments.output) &&
        !std::filesystem::is_empty(arguments.output)) {
      throw std::runtime_error(
          "--output must be absent or empty");
    }
    std::filesystem::create_directories(arguments.output);
    std::vector<Matrix> seeds;
    for (const auto& path : arguments.seeds) {
      Matrix seed = read_matrix(path);
      const Wide score = absolute(determinant(seed));
      if (score != kTarget) {
        throw std::runtime_error(
            "seed is off target: " + path.string() +
            " |det|=" + wide_string(score));
      }
      seeds.push_back(seed);
    }
    std::mt19937_64 randomizer(arguments.random_seed);
    std::array<int, kOrder * kOrder> coordinates{};
    for (int index = 0; index < kOrder * kOrder; ++index) {
      coordinates[index] = index;
    }
    const auto started = std::chrono::steady_clock::now();
    const auto deadline =
        started +
        std::chrono::duration<double>(arguments.seconds);
    std::uint64_t attempts = 0;
    std::uint64_t completed = 0;
    std::uint64_t target_hits = 0;
    std::uint64_t saved = 0;
    std::uint64_t total_steps = 0;
    std::set<MatrixKey> raw_matrices;
    std::map<Wide, std::uint64_t> score_counts;
    Wide best = 0;
    Matrix walk_matrix = seeds.front();
    bool have_walk_matrix = false;
    while (
        std::chrono::steady_clock::now() < deadline &&
        (!arguments.max_attempts ||
         attempts < arguments.max_attempts)) {
      ++attempts;
      Matrix matrix;
      const bool random_start =
          static_cast<int>(randomizer() % 100) <
          arguments.random_start_percent;
      if (random_start) {
        for (auto& row : matrix) {
          for (int& value : row) {
            value = randomizer() & 1 ? 1 : -1;
          }
        }
      } else {
        matrix =
            arguments.walk && have_walk_matrix
                ? walk_matrix
                : seeds[randomizer() % seeds.size()];
        std::shuffle(
            coordinates.begin(),
            coordinates.end(),
            randomizer);
        const int width =
            arguments.maximum_kick -
            arguments.minimum_kick + 1;
        const int kick =
            arguments.minimum_kick +
            static_cast<int>(randomizer() % width);
        for (int index = 0; index < kick; ++index) {
          const int coordinate = coordinates[index];
          matrix[coordinate / kOrder][coordinate % kOrder] *= -1;
        }
      }
      int steps = 0;
      if (!ascend(
              matrix,
              randomizer,
              arguments.random_improvement,
              arguments.axis_mode,
              steps)) {
        continue;
      }
      if (arguments.walk) {
        walk_matrix = matrix;
        have_walk_matrix = true;
      }
      ++completed;
      total_steps += steps;
      const Wide score = absolute(determinant(matrix));
      ++score_counts[score];
      if (score > best) {
        best = score;
        std::cout
            << "best=" << wide_string(best)
            << " attempts=" << attempts
            << " completed=" << completed
            << " target_hits=" << target_hits << '\n'
            << std::flush;
      }
      if (score != kTarget) continue;
      ++target_hits;
      if (
          saved >= arguments.max_saved ||
          !raw_matrices.insert(matrix_key(matrix)).second) {
        continue;
      }
      std::ostringstream name;
      name
          << "target-" << std::setw(6) << std::setfill('0')
          << saved++ << ".matrix.txt";
      write_matrix(arguments.output / name.str(), matrix);
    }
    const double elapsed =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started)
            .count();
    std::ofstream summary(arguments.output / "summary.txt");
    summary
        << "engine order22-gradient-harvest-v1\n"
        << "prng std::mt19937_64\n"
        << "random_seed " << arguments.random_seed << '\n'
        << "seconds_limit " << arguments.seconds << '\n'
        << "max_attempts " << arguments.max_attempts << '\n'
        << "minimum_kick " << arguments.minimum_kick << '\n'
        << "maximum_kick " << arguments.maximum_kick << '\n'
        << "random_start_percent "
        << arguments.random_start_percent << '\n'
        << "random_improvement "
        << arguments.random_improvement << '\n'
        << "axis " << arguments.axis_mode << '\n'
        << "walk " << arguments.walk << '\n'
        << "max_saved " << arguments.max_saved << '\n';
    for (const auto& path : arguments.seeds) {
      summary << "seed " << path << '\n';
    }
    summary
        << "elapsed_seconds " << elapsed << '\n'
        << "attempts " << attempts << '\n'
        << "completed " << completed << '\n'
        << "target_hits " << target_hits << '\n'
        << "saved_raw " << saved << '\n'
        << "average_ascent_steps "
        << (
               completed
                   ? static_cast<double>(total_steps) / completed
                   : 0.0)
        << '\n'
        << "best " << wide_string(best) << '\n'
        << "score_histogram\n";
    std::size_t histogram_rows = 0;
    for (
        auto iterator = score_counts.rbegin();
        iterator != score_counts.rend() && histogram_rows < 20;
        ++iterator, ++histogram_rows) {
      summary
          << wide_string(iterator->first)
          << ' ' << iterator->second << '\n';
    }
    summary.flush();
    if (!summary) {
      throw std::runtime_error("cannot flush summary");
    }
    std::cout
        << "finished attempts=" << attempts
        << " completed=" << completed
        << " hits=" << target_hits
        << " saved=" << saved
        << " best=" << wide_string(best)
        << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "order22_gradient_harvest: "
              << error.what() << '\n';
    return 2;
  }
}
