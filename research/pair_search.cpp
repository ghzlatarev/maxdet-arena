#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>

namespace {

constexpr int kOrder = 23;
using Wide = __int128_t;
using Matrix = std::array<std::array<Wide, kOrder>, kOrder>;

Wide absolute(Wide value) { return value < 0 ? -value : value; }

std::string wide_to_string(Wide value) {
  if (value == 0) return "0";
  const bool negative = value < 0;
  if (negative) value = -value;
  std::string result;
  while (value != 0) {
    result.push_back(static_cast<char>('0' + value % 10));
    value /= 10;
  }
  if (negative) result.push_back('-');
  std::reverse(result.begin(), result.end());
  return result;
}

Wide exact_determinant(Matrix matrix) {
  Wide previous_pivot = 1;
  int sign = 1;
  for (int column = 0; column < kOrder - 1; ++column) {
    int pivot_row = column;
    while (pivot_row < kOrder && matrix[pivot_row][column] == 0) ++pivot_row;
    if (pivot_row == kOrder) return 0;
    if (pivot_row != column) {
      std::swap(matrix[pivot_row], matrix[column]);
      sign = -sign;
    }

    const Wide pivot = matrix[column][column];
    for (int row = column + 1; row < kOrder; ++row) {
      for (int inner = column + 1; inner < kOrder; ++inner) {
        const Wide numerator =
            matrix[row][inner] * pivot -
            matrix[row][column] * matrix[column][inner];
        if (column != 0 && numerator % previous_pivot != 0) {
          throw std::runtime_error("exact Bareiss division failed");
        }
        matrix[row][inner] =
            column == 0 ? numerator : numerator / previous_pivot;
      }
      matrix[row][column] = 0;
    }
    previous_pivot = pivot;
  }
  return sign * matrix[kOrder - 1][kOrder - 1];
}

Matrix read_matrix(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open start matrix");
  Matrix matrix{};
  for (auto& row : matrix) {
    for (Wide& value : row) {
      int parsed = 0;
      if (!(input >> parsed) || (parsed != -1 && parsed != 1)) {
        throw std::runtime_error("invalid start matrix");
      }
      value = parsed;
    }
  }
  std::string extra;
  if (input >> extra) throw std::runtime_error("extra start matrix data");
  return matrix;
}

void write_matrix(const std::filesystem::path& path, const Matrix& matrix) {
  if (!path.parent_path().empty()) {
    std::filesystem::create_directories(path.parent_path());
  }
  const std::filesystem::path temporary = path.string() + ".tmp";
  {
    std::ofstream output(temporary, std::ios::trunc);
    if (!output) throw std::runtime_error("cannot create output matrix");
    for (const auto& row : matrix) {
      for (int column = 0; column < kOrder; ++column) {
        if (column != 0) output << ' ';
        output << wide_to_string(row[column]);
      }
      output << '\n';
    }
    output.flush();
    if (!output) throw std::runtime_error("cannot flush output matrix");
  }
  std::filesystem::rename(temporary, path);
}

Matrix transpose(const Matrix& matrix) {
  Matrix result{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      result[row][column] = matrix[column][row];
    }
  }
  return result;
}

struct PairResult {
  Matrix matrix{};
  Wide score = 0;
  std::uint64_t assignments = 0;
};

PairResult optimize_row_pair(const Matrix& matrix, int first, int second) {
  if (first == second) throw std::runtime_error("row pair must be distinct");
  const Wide determinant = exact_determinant(matrix);
  if (determinant == 0) {
    throw std::runtime_error("pair optimization requires a nonsingular matrix");
  }

  std::array<Wide, kOrder> first_cofactors{};
  std::array<Wide, kOrder> second_cofactors{};
  for (int column = 0; column < kOrder; ++column) {
    Matrix basis = matrix;
    basis[first].fill(0);
    basis[first][column] = 1;
    first_cofactors[column] = exact_determinant(basis);

    basis = matrix;
    basis[second].fill(0);
    basis[second][column] = 1;
    second_cofactors[column] = exact_determinant(basis);
  }

  std::array<std::array<Wide, kOrder>, kOrder> coefficients{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      const Wide numerator =
          first_cofactors[row] * second_cofactors[column] -
          second_cofactors[row] * first_cofactors[column];
      if (numerator % determinant != 0) {
        throw std::runtime_error("second-cofactor division failed");
      }
      coefficients[row][column] = numerator / determinant;
    }
  }

  // Negating an entire replacement row preserves the absolute determinant, so
  // fix its first sign and enumerate only the remaining 22 signs in Gray order.
  std::array<Wide, kOrder> first_row{};
  first_row.fill(-1);
  first_row[0] = 1;
  std::array<Wide, kOrder> weighted_columns{};
  for (int column = 0; column < kOrder; ++column) {
    for (int row = 0; row < kOrder; ++row) {
      weighted_columns[column] += first_row[row] * coefficients[row][column];
    }
  }

  PairResult result{matrix, absolute(determinant), 0};
  const std::uint64_t assignment_count = std::uint64_t{1} << (kOrder - 1);
  for (std::uint64_t assignment = 0; assignment < assignment_count;
       ++assignment) {
    if (assignment != 0) {
      const int bit = static_cast<int>(std::countr_zero(assignment)) + 1;
      first_row[bit] = -first_row[bit];
      const Wide delta = 2 * first_row[bit];
      for (int column = 0; column < kOrder; ++column) {
        weighted_columns[column] += delta * coefficients[bit][column];
      }
    }

    Wide score = 0;
    for (const Wide value : weighted_columns) score += absolute(value);
    ++result.assignments;
    if (score <= result.score) continue;

    result.score = score;
    for (int column = 0; column < kOrder; ++column) {
      result.matrix[first][column] = first_row[column];
      result.matrix[second][column] =
          weighted_columns[column] < 0 ? -1 : 1;
    }
  }

  const Wide checked = absolute(exact_determinant(result.matrix));
  if (checked != result.score) {
    throw std::runtime_error(
        "joint row score disagrees with exact determinant");
  }
  return result;
}

struct Arguments {
  std::filesystem::path start;
  std::filesystem::path output;
  std::filesystem::path research_output;
  std::filesystem::path log;
  std::uint64_t seed = 23;
  double seconds = 3600.0;
  std::uint64_t passes = 0;
  int kick_size = 0;
};

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    auto value = [&]() -> std::string {
      if (++index >= argc) throw std::runtime_error("missing option value");
      return argv[index];
    };
    if (option == "--start") arguments.start = value();
    else if (option == "--output") arguments.output = value();
    else if (option == "--research-output") arguments.research_output = value();
    else if (option == "--log") arguments.log = value();
    else if (option == "--seed") arguments.seed = std::stoull(value());
    else if (option == "--seconds") arguments.seconds = std::stod(value());
    else if (option == "--passes") arguments.passes = std::stoull(value());
    else if (option == "--kick-size") arguments.kick_size = std::stoi(value());
    else throw std::runtime_error("unknown option: " + option);
  }
  if (arguments.start.empty()) throw std::runtime_error("--start is required");
  if (arguments.output.empty()) throw std::runtime_error("--output is required");
  if (arguments.log.empty()) throw std::runtime_error("--log is required");
  auto normalized_path = [](const std::filesystem::path& path) {
    return std::filesystem::absolute(path).lexically_normal();
  };
  const auto start = normalized_path(arguments.start);
  const auto output = normalized_path(arguments.output);
  const auto log = normalized_path(arguments.log);
  if (output == log || output == start || log == start) {
    throw std::runtime_error("--start, --output, and --log must be distinct");
  }
  if (!arguments.research_output.empty()) {
    const auto research = normalized_path(arguments.research_output);
    if (research == output || research == log || research == start) {
      throw std::runtime_error(
          "--research-output must not alias start, output, or log");
    }
  }
  if (!std::isfinite(arguments.seconds) ||
      (arguments.passes == 0 && arguments.seconds <= 0)) {
    throw std::runtime_error("--seconds must be finite and positive");
  }
  if (arguments.kick_size < 0 ||
      arguments.kick_size > kOrder * kOrder) {
    throw std::runtime_error("--kick-size must be between 0 and 529");
  }
  return arguments;
}

void log_record(std::ofstream& log, Wide best_score, Wide state_score,
                Wide pair_score, double elapsed, std::uint64_t assignments,
                int first, int second, bool columns, const char* event,
                std::uint64_t seed) {
  log << "{\"absolute_determinant\":\"" << wide_to_string(best_score)
      << "\",\"assignments\":" << assignments
      << ",\"elapsed_seconds\":" << std::fixed << std::setprecision(3)
      << elapsed
      << ",\"event\":\"" << event
      << "\",\"first\":" << first
      << ",\"orientation\":\"" << (columns ? "columns" : "rows")
      << "\",\"pair_determinant\":\"" << wide_to_string(pair_score)
      << "\",\"second\":" << second
      << ",\"seed\":" << seed
      << ",\"state_determinant\":\"" << wide_to_string(state_score)
      << "\"}\n";
  log.flush();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    std::mt19937_64 randomizer(arguments.seed);
    if (!arguments.log.parent_path().empty()) {
      std::filesystem::create_directories(arguments.log.parent_path());
    }
    std::ofstream log(arguments.log, std::ios::app);
    if (!log) throw std::runtime_error("cannot open research log");

    Matrix best_matrix = read_matrix(arguments.start);
    Wide best_score = absolute(exact_determinant(best_matrix));
    write_matrix(arguments.output, best_matrix);
    Matrix state_matrix = best_matrix;
    Wide state_score = best_score;
    Wide research_score = 0;
    auto checkpoint_research_state = [&]() {
      if (!arguments.research_output.empty() &&
          state_score < best_score &&
          state_score > research_score) {
        research_score = state_score;
        write_matrix(arguments.research_output, state_matrix);
      }
    };

    std::array<std::pair<int, int>, kOrder * (kOrder - 1)> pairs{};
    int pair_count = 0;
    for (int orientation = 0; orientation < 2; ++orientation) {
      for (int first = 0; first < kOrder; ++first) {
        for (int second = first + 1; second < kOrder; ++second) {
          pairs[pair_count++] = {
              orientation * kOrder + first,
              orientation * kOrder + second,
          };
        }
      }
    }
    std::shuffle(pairs.begin(), pairs.begin() + pair_count, randomizer);

    const auto started = std::chrono::steady_clock::now();
    const auto deadline =
        started + std::chrono::duration<double>(arguments.seconds);
    std::uint64_t total_assignments = 0;
    std::uint64_t completed_passes = 0;
    int examined = 0;
    log_record(
        log, best_score, state_score, state_score, 0.0, 0, -1, -1, false,
        "start", arguments.seed);

    std::array<int, kOrder * kOrder> coordinates{};
    std::iota(coordinates.begin(), coordinates.end(), 0);
    auto kick_from_best = [&]() {
      for (int attempt = 0; attempt < 32; ++attempt) {
        state_matrix = best_matrix;
        std::shuffle(coordinates.begin(), coordinates.end(), randomizer);
        for (int index = 0; index < arguments.kick_size; ++index) {
          const int encoded = coordinates[index];
          state_matrix[encoded / kOrder][encoded % kOrder] *= -1;
        }
        state_score = absolute(exact_determinant(state_matrix));
        if (state_score != 0) {
          checkpoint_research_state();
          return;
        }
      }
      throw std::runtime_error("could not produce a nonsingular kicked state");
    };

    while (arguments.passes != 0
               ? completed_passes < arguments.passes
               : std::chrono::steady_clock::now() < deadline) {
      if (examined == pair_count) {
        ++completed_passes;
        if (arguments.passes != 0 &&
            completed_passes >= arguments.passes) {
          break;
        }
        std::shuffle(pairs.begin(), pairs.begin() + pair_count, randomizer);
        examined = 0;
      }
      if (examined == 0 && arguments.kick_size != 0) {
        kick_from_best();
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        log_record(
            log, best_score, state_score, state_score, elapsed,
            total_assignments, -1, -1, false, "kick", arguments.seed);
      }
      const auto encoded = pairs[examined++];
      const bool columns = encoded.first >= kOrder;
      const int first = encoded.first % kOrder;
      const int second = encoded.second % kOrder;
      const Matrix oriented = columns ? transpose(state_matrix) : state_matrix;
      PairResult result = optimize_row_pair(oriented, first, second);
      total_assignments += result.assignments;
      Matrix candidate = columns ? transpose(result.matrix) : result.matrix;
      const Wide candidate_score = absolute(exact_determinant(candidate));
      if (candidate_score != result.score) {
        throw std::runtime_error(
            "orientation score disagrees with exact determinant");
      }
      const bool moved = candidate_score > state_score;
      bool improved = false;
      if (moved) {
        state_score = candidate_score;
        state_matrix = candidate;
        if (state_score > best_score) {
          improved = true;
          best_score = state_score;
          best_matrix = state_matrix;
          write_matrix(arguments.output, best_matrix);
        }
        checkpoint_research_state();
      }
      const double elapsed = std::chrono::duration<double>(
          std::chrono::steady_clock::now() - started).count();
      log_record(
          log, best_score, state_score, candidate_score, elapsed,
          total_assignments, first, second, columns,
          improved ? "new_best" : moved ? "move" : "pair_finished",
          arguments.seed);
      std::cout << (columns ? "columns " : "rows ")
                << first << ',' << second
                << " pair |det|=" << wide_to_string(candidate_score)
                << " state=" << wide_to_string(state_score)
                << " best=" << wide_to_string(best_score)
                << " elapsed=" << std::fixed << std::setprecision(1)
                << elapsed << "s\n" << std::flush;
    }

    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    log_record(
        log, best_score, state_score, best_score, elapsed, total_assignments,
        -1, -1, false, "finished", arguments.seed);
    std::cout << "finished |det|=" << wide_to_string(best_score)
              << " research |det|=" << wide_to_string(research_score)
              << " assignments=" << total_assignments << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "pair_search: " << error.what() << '\n';
    return 2;
  }
}
