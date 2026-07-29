#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <system_error>

namespace {

constexpr int kCoreOrder = 22;
constexpr int kOutputOrder = 23;
constexpr std::uint64_t kAssignments = std::uint64_t{1} << (kCoreOrder - 1);
using Wide = __int128_t;

template <std::size_t Order>
using Matrix = std::array<std::array<int, Order>, Order>;

using Core = Matrix<kCoreOrder>;
using Candidate = Matrix<kOutputOrder>;
using Adjugate =
    std::array<std::array<Wide, kCoreOrder>, kCoreOrder>;

struct Options {
  std::filesystem::path start;
  std::filesystem::path output;
  std::filesystem::path log;
};

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

Wide absolute(Wide value) { return value < 0 ? -value : value; }

template <std::size_t Order>
Wide exact_determinant(Matrix<Order> matrix) {
  constexpr int size = static_cast<int>(Order);
  std::array<std::array<Wide, Order>, Order> work{};
  for (int row = 0; row < size; ++row) {
    for (int column = 0; column < size; ++column) {
      work[row][column] = matrix[row][column];
    }
  }

  Wide previous_pivot = 1;
  int sign = 1;
  for (int column = 0; column < size - 1; ++column) {
    int pivot_row = column;
    while (pivot_row < size && work[pivot_row][column] == 0) {
      ++pivot_row;
    }
    if (pivot_row == size) return 0;
    if (pivot_row != column) {
      std::swap(work[pivot_row], work[column]);
      sign = -sign;
    }

    const Wide pivot = work[column][column];
    for (int row = column + 1; row < size; ++row) {
      for (int inner = column + 1; inner < size; ++inner) {
        const Wide numerator =
            work[row][inner] * pivot -
            work[row][column] * work[column][inner];
        if (column != 0 && numerator % previous_pivot != 0) {
          throw std::runtime_error("exact Bareiss division failed");
        }
        work[row][inner] =
            column == 0 ? numerator : numerator / previous_pivot;
      }
      work[row][column] = 0;
    }
    previous_pivot = pivot;
  }
  return sign * work[size - 1][size - 1];
}

Core read_core(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open core matrix");
  Core matrix{};
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      if (!(input >> matrix[row][column]) ||
          (matrix[row][column] != -1 && matrix[row][column] != 1)) {
        throw std::runtime_error("invalid 22 by 22 sign core");
      }
    }
  }
  std::string extra;
  if (input >> extra) throw std::runtime_error("extra core matrix data");
  return matrix;
}

Adjugate exact_adjugate(const Core& core) {
  Adjugate adjugate{};
  for (int excluded_row = 0; excluded_row < kCoreOrder; ++excluded_row) {
    for (int excluded_column = 0; excluded_column < kCoreOrder;
         ++excluded_column) {
      Matrix<kCoreOrder - 1> minor{};
      int minor_row = 0;
      for (int row = 0; row < kCoreOrder; ++row) {
        if (row == excluded_row) continue;
        int minor_column = 0;
        for (int column = 0; column < kCoreOrder; ++column) {
          if (column == excluded_column) continue;
          minor[minor_row][minor_column] = core[row][column];
          ++minor_column;
        }
        ++minor_row;
      }
      Wide cofactor = exact_determinant(minor);
      if ((excluded_row + excluded_column) % 2 != 0) {
        cofactor = -cofactor;
      }
      adjugate[excluded_column][excluded_row] = cofactor;
    }
  }
  return adjugate;
}

void validate_adjugate(
    const Core& core,
    const Adjugate& adjugate,
    Wide determinant) {
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      Wide product = 0;
      for (int inner = 0; inner < kCoreOrder; ++inner) {
        product +=
            static_cast<Wide>(core[row][inner]) *
            adjugate[inner][column];
      }
      const Wide expected = row == column ? determinant : 0;
      if (product != expected) {
        throw std::runtime_error("exact adjugate identity check failed");
      }
    }
  }
}

Candidate bordered_candidate(
    const Core& core,
    const Adjugate& adjugate,
    const std::array<int, kCoreOrder>& border_column,
    Wide core_determinant) {
  Candidate candidate{};
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      candidate[row][column] = core[row][column];
    }
    candidate[row][kCoreOrder] = border_column[row];
  }

  for (int row = 0; row < kCoreOrder; ++row) {
    Wide value = 0;
    for (int column = 0; column < kCoreOrder; ++column) {
      value +=
          adjugate[row][column] * border_column[column];
    }
    candidate[kCoreOrder][row] = value > 0 ? -1 : 1;
  }
  candidate[kCoreOrder][kCoreOrder] = core_determinant > 0 ? 1 : -1;
  return candidate;
}

void write_matrix(
    const std::filesystem::path& path,
    const Candidate& matrix) {
  if (!path.parent_path().empty()) {
    std::filesystem::create_directories(path.parent_path());
  }
  if (std::filesystem::exists(path) &&
      !std::filesystem::is_regular_file(path)) {
    throw std::runtime_error("output path is not a regular file");
  }
  const std::filesystem::path temporary = path.string() + ".tmp";
  {
    std::ofstream output(temporary, std::ios::trunc);
    if (!output) throw std::runtime_error("cannot create output matrix");
    for (const auto& row : matrix) {
      for (int column = 0; column < kOutputOrder; ++column) {
        if (column != 0) output << ' ';
        output << row[column];
      }
      output << '\n';
    }
    output.flush();
    if (!output) throw std::runtime_error("cannot flush output matrix");
  }
  std::error_code error;
  std::filesystem::rename(temporary, path, error);
  if (error) {
    std::filesystem::remove(temporary);
    throw std::runtime_error("cannot replace output matrix");
  }
}

std::string json_escape(const std::string& value) {
  std::string result;
  for (const char character : value) {
    switch (character) {
      case '\\':
        result += "\\\\";
        break;
      case '"':
        result += "\\\"";
        break;
      case '\n':
        result += "\\n";
        break;
      case '\r':
        result += "\\r";
        break;
      case '\t':
        result += "\\t";
        break;
      default:
        if (static_cast<unsigned char>(character) < 0x20U) {
          throw std::runtime_error("control character in log path");
        }
        result += character;
    }
  }
  return result;
}

Options parse_options(int argc, char** argv) {
  Options options;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--help") {
      std::cout
          << "Usage: order22_border --start CORE --output MATRIX --log JSONL\n"
          << "Exhaust all 2^21 border columns of a fixed 22x22 sign core.\n";
      std::exit(0);
    }
    if (index + 1 >= argc) {
      throw std::runtime_error("missing value for " + argument);
    }
    const std::filesystem::path value = argv[++index];
    if (argument == "--start") {
      options.start = value;
    } else if (argument == "--output") {
      options.output = value;
    } else if (argument == "--log") {
      options.log = value;
    } else {
      throw std::runtime_error("unknown argument: " + argument);
    }
  }
  if (options.start.empty() || options.output.empty() || options.log.empty()) {
    throw std::runtime_error("--start, --output, and --log are required");
  }

  const auto start = std::filesystem::absolute(options.start).lexically_normal();
  const auto output =
      std::filesystem::absolute(options.output).lexically_normal();
  const auto log = std::filesystem::absolute(options.log).lexically_normal();
  if (start == output || start == log || output == log) {
    throw std::runtime_error("start, output, and log paths must differ");
  }
  options.start = start;
  options.output = output;
  options.log = log;
  return options;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse_options(argc, argv);
    const auto start_time = std::chrono::steady_clock::now();
    const Core core = read_core(options.start);
    const Wide core_determinant = exact_determinant(core);
    if (core_determinant == 0) {
      throw std::runtime_error("core matrix is singular");
    }
    const Adjugate adjugate = exact_adjugate(core);
    validate_adjugate(core, adjugate, core_determinant);

    if (!options.log.parent_path().empty()) {
      std::filesystem::create_directories(options.log.parent_path());
    }
    std::ofstream log(options.log, std::ios::trunc);
    if (!log) throw std::runtime_error("cannot create log");
    log << "{\"absolute_core_determinant\":\""
        << wide_to_string(absolute(core_determinant))
        << "\",\"assignments_total\":" << kAssignments
        << ",\"event\":\"start\",\"start\":\""
        << json_escape(options.start.string()) << "\"}\n";

    std::array<int, kCoreOrder> border_column{};
    border_column.fill(1);
    std::array<Wide, kCoreOrder> products{};
    for (int row = 0; row < kCoreOrder; ++row) {
      for (int column = 0; column < kCoreOrder; ++column) {
        products[row] += adjugate[row][column];
      }
    }

    Wide best_score = absolute(core_determinant);
    for (const Wide value : products) best_score += absolute(value);
    Candidate best =
        bordered_candidate(core, adjugate, border_column, core_determinant);
    if (absolute(exact_determinant(best)) != best_score) {
      throw std::runtime_error("initial border identity check failed");
    }
    write_matrix(options.output, best);
    std::uint64_t improvements = 0;
    std::uint64_t best_multiplicity = 1;

    for (std::uint64_t assignment = 1; assignment < kAssignments;
         ++assignment) {
      const unsigned bit = std::countr_zero(assignment) + 1U;
      const int old_sign = border_column[bit];
      border_column[bit] = -old_sign;
      for (int row = 0; row < kCoreOrder; ++row) {
        products[row] -=
            static_cast<Wide>(2 * old_sign) * adjugate[row][bit];
      }

      Wide score = absolute(core_determinant);
      for (const Wide value : products) score += absolute(value);
      if (score < best_score) continue;
      if (score == best_score) {
        ++best_multiplicity;
        continue;
      }

      Candidate candidate = bordered_candidate(
          core,
          adjugate,
          border_column,
          core_determinant);
      const Wide checked_score = absolute(exact_determinant(candidate));
      if (checked_score != score) {
        throw std::runtime_error("border determinant identity check failed");
      }
      best_score = score;
      best = candidate;
      ++improvements;
      best_multiplicity = 1;
      write_matrix(options.output, best);
      log << "{\"absolute_determinant\":\""
          << wide_to_string(best_score)
          << "\",\"assignment\":" << assignment
          << ",\"event\":\"new_best\",\"improvements\":" << improvements
          << "}\n";
    }

    const Wide final_score = absolute(exact_determinant(best));
    if (final_score != best_score) {
      throw std::runtime_error("final exact determinant check failed");
    }
    const double elapsed = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - start_time)
                               .count();
    log << "{\"absolute_determinant\":\""
        << wide_to_string(best_score)
        << "\",\"assignments_completed\":" << kAssignments
        << ",\"best_border_columns_up_to_global_sign\":"
        << best_multiplicity
        << ",\"complete\":true,\"elapsed_seconds\":" << elapsed
        << ",\"event\":\"finished\",\"improvements\":" << improvements
        << "}\n";
    log.flush();
    if (!log) throw std::runtime_error("cannot flush log");

    std::cout << "completed " << kAssignments
              << " exact border scores; best=" << wide_to_string(best_score)
              << " elapsed_seconds=" << elapsed << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
