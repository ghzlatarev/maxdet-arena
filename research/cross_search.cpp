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
#include <string_view>

namespace {

constexpr int kOrder = 23;
constexpr int kFree = kOrder - 1;
constexpr int kEntries = kOrder * kOrder;
using Wide = __int128_t;
using Matrix = std::array<std::array<Wide, kOrder>, kOrder>;
using Clock = std::chrono::steady_clock;

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
      const Wide left = matrix[row][column];
      for (int inner = column + 1; inner < kOrder; ++inner) {
        const Wide numerator =
            matrix[row][inner] * pivot -
            left * matrix[column][inner];
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

Matrix sign_normalized(Matrix matrix) {
  for (int column = 0; column < kOrder; ++column) {
    if (matrix[0][column] == -1) {
      for (int row = 0; row < kOrder; ++row) {
        matrix[row][column] *= -1;
      }
    }
  }
  for (int row = 0; row < kOrder; ++row) {
    if (matrix[row][0] == -1) {
      for (int column = 0; column < kOrder; ++column) {
        matrix[row][column] *= -1;
      }
    }
  }
  return matrix;
}

bool sign_equivalent(const Matrix& left, const Matrix& right) {
  return sign_normalized(left) == sign_normalized(right);
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

struct Arguments {
  std::filesystem::path start;
  std::filesystem::path output;
  std::filesystem::path log;
  std::uint64_t seed = 23;
  double seconds = 3600.0;
  double heartbeat_seconds = 30.0;
  std::uint64_t passes = 0;
  std::uint64_t settle_passes = 2;
  int kick_size = 0;
};

std::uint64_t parse_unsigned(
    std::string_view text, std::string_view option) {
  std::size_t parsed = 0;
  const std::string copy(text);
  const unsigned long long value = std::stoull(copy, &parsed);
  if (parsed != copy.size()) {
    throw std::runtime_error("invalid integer for " + std::string(option));
  }
  return static_cast<std::uint64_t>(value);
}

int parse_integer(std::string_view text, std::string_view option) {
  std::size_t parsed = 0;
  const std::string copy(text);
  const long value = std::stol(copy, &parsed);
  if (parsed != copy.size() ||
      value < std::numeric_limits<int>::min() ||
      value > std::numeric_limits<int>::max()) {
    throw std::runtime_error("invalid integer for " + std::string(option));
  }
  return static_cast<int>(value);
}

double parse_double(std::string_view text, std::string_view option) {
  std::size_t parsed = 0;
  const std::string copy(text);
  const double value = std::stod(copy, &parsed);
  if (parsed != copy.size() || !std::isfinite(value)) {
    throw std::runtime_error("invalid number for " + std::string(option));
  }
  return value;
}

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
    else if (option == "--log") arguments.log = value();
    else if (option == "--seed") {
      arguments.seed = parse_unsigned(value(), option);
    } else if (option == "--seconds") {
      arguments.seconds = parse_double(value(), option);
    } else if (option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds = parse_double(value(), option);
    } else if (option == "--passes") {
      arguments.passes = parse_unsigned(value(), option);
    } else if (option == "--settle-passes") {
      arguments.settle_passes = parse_unsigned(value(), option);
    } else if (option == "--kick-size") {
      arguments.kick_size = parse_integer(value(), option);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }

  if (arguments.start.empty()) throw std::runtime_error("--start is required");
  if (arguments.output.empty()) {
    throw std::runtime_error("--output is required");
  }
  if (arguments.log.empty()) throw std::runtime_error("--log is required");
  if (arguments.seconds <= 0 && arguments.passes == 0) {
    throw std::runtime_error(
        "--seconds must be positive when --passes is zero");
  }
  if (arguments.heartbeat_seconds < 0) {
    throw std::runtime_error("--heartbeat-seconds must be non-negative");
  }
  if (arguments.settle_passes == 0) {
    throw std::runtime_error("--settle-passes must be positive");
  }
  if (arguments.kick_size < 0 || arguments.kick_size > kEntries) {
    throw std::runtime_error("--kick-size must be between 0 and 529");
  }

  auto normalized = [](const std::filesystem::path& path) {
    return std::filesystem::absolute(path).lexically_normal();
  };
  const auto start = normalized(arguments.start);
  const auto output = normalized(arguments.output);
  const auto log = normalized(arguments.log);
  if (start == output || start == log || output == log) {
    throw std::runtime_error(
        "--start, --output, and --log must be distinct");
  }
  return arguments;
}

struct CrossResult {
  Matrix matrix{};
  Wide score = 0;
  std::uint64_t assignments = 0;
  bool complete = true;
};

CrossResult optimize_cross(
    const Matrix& matrix,
    int selected_row,
    int selected_column,
    Clock::time_point deadline) {
  CrossResult result{matrix, absolute(exact_determinant(matrix)), 0, true};

  std::array<int, kFree> row_columns{};
  std::array<int, kFree> column_rows{};
  int next_column = 0;
  int next_row = 0;
  for (int coordinate = 0; coordinate < kOrder; ++coordinate) {
    if (coordinate != selected_column) row_columns[next_column++] = coordinate;
    if (coordinate != selected_row) column_rows[next_row++] = coordinate;
  }

  // With the selected row and column zeroed, the determinant is
  //
  //   center * d + x^T H y.
  //
  // Each entry of H is obtained independently by an exact basis determinant.
  // This is setup work for a 2^21 neighborhood, so the 485 Bareiss calls are
  // small compared with enumerating the row signs.
  Matrix basis = matrix;
  basis[selected_row].fill(0);
  for (int row : column_rows) basis[row][selected_column] = 0;

  Matrix central_basis = basis;
  central_basis[selected_row][selected_column] = 1;
  const Wide central_coefficient = exact_determinant(central_basis);

  std::array<std::array<Wide, kFree>, kFree> coefficients{};
  Matrix coefficient_basis = basis;
  for (int row_sign = 0; row_sign < kFree; ++row_sign) {
    if (Clock::now() >= deadline) {
      result.complete = false;
      return result;
    }
    const int column = row_columns[row_sign];
    coefficient_basis[selected_row][column] = 1;
    for (int column_sign = 0; column_sign < kFree; ++column_sign) {
      const int row = column_rows[column_sign];
      coefficient_basis[row][selected_column] = 1;
      coefficients[row_sign][column_sign] =
          exact_determinant(coefficient_basis);
      coefficient_basis[row][selected_column] = 0;
    }
    coefficient_basis[selected_row][column] = 0;
  }

  // Negating x together with the center negates the determinant, so fix the
  // first row sign and enumerate only 2^21 assignments.  For each x, the best
  // center and all 22 column signs are independent sign choices:
  //
  //   max |center*d + sum_i (x^T H)_i y_i|
  //       = |d| + sum_i |(x^T H)_i|.
  std::array<int, kFree> row_signs{};
  row_signs.fill(-1);
  row_signs[0] = 1;
  std::array<Wide, kFree> weighted_columns{};
  for (int column_sign = 0; column_sign < kFree; ++column_sign) {
    for (int row_sign = 0; row_sign < kFree; ++row_sign) {
      weighted_columns[column_sign] +=
          row_signs[row_sign] * coefficients[row_sign][column_sign];
    }
  }

  std::array<int, kFree> best_row_signs{};
  std::array<int, kFree> best_column_signs{};
  int best_center = 1;
  bool found_better = false;
  const std::uint64_t assignment_count = std::uint64_t{1} << (kFree - 1);
  for (std::uint64_t assignment = 0; assignment < assignment_count;
       ++assignment) {
    if (assignment != 0) {
      const int bit =
          static_cast<int>(std::countr_zero(assignment)) + 1;
      row_signs[bit] = -row_signs[bit];
      const int delta = 2 * row_signs[bit];
      for (int column_sign = 0; column_sign < kFree; ++column_sign) {
        weighted_columns[column_sign] +=
            delta * coefficients[bit][column_sign];
      }
    }

    Wide score = absolute(central_coefficient);
    for (Wide value : weighted_columns) score += absolute(value);
    ++result.assignments;
    if (score > result.score) {
      result.score = score;
      best_row_signs = row_signs;
      best_center = central_coefficient < 0 ? -1 : 1;
      for (int index = 0; index < kFree; ++index) {
        best_column_signs[index] =
            weighted_columns[index] < 0 ? -1 : 1;
      }
      found_better = true;
    }

    if ((assignment & 65535U) == 65535U &&
        Clock::now() >= deadline) {
      result.complete = false;
      break;
    }
  }

  if (found_better) {
    for (int index = 0; index < kFree; ++index) {
      result.matrix[selected_row][row_columns[index]] =
          best_row_signs[index];
      result.matrix[column_rows[index]][selected_column] =
          best_column_signs[index];
    }
    result.matrix[selected_row][selected_column] = best_center;
  }

  const Wide checked = absolute(exact_determinant(result.matrix));
  if (checked != result.score) {
    throw std::runtime_error(
        "cross score disagrees with exact determinant");
  }
  return result;
}

void log_record(
    std::ofstream& log,
    const Arguments& arguments,
    Wide best_score,
    Wide state_score,
    Wide cross_score,
    double elapsed,
    std::uint64_t assignments,
    std::uint64_t completed_passes,
    int row,
    int column,
    const char* event) {
  log << "{\"absolute_determinant\":\"" << wide_to_string(best_score)
      << "\",\"assignments\":" << assignments
      << ",\"column\":" << column
      << ",\"completed_passes\":" << completed_passes
      << ",\"cross_determinant\":\"" << wide_to_string(cross_score)
      << "\",\"elapsed_seconds\":" << std::fixed << std::setprecision(3)
      << elapsed
      << ",\"event\":\"" << event
      << "\",\"row\":" << row
      << ",\"seed\":" << arguments.seed
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
    if (best_score == 0) {
      throw std::runtime_error("start matrix must be nonsingular");
    }
    Matrix state_matrix = best_matrix;
    Wide state_score = best_score;
    write_matrix(arguments.output, best_matrix);

    std::array<std::pair<int, int>, kEntries> crosses{};
    int cross_count = 0;
    for (int row = 0; row < kOrder; ++row) {
      for (int column = 0; column < kOrder; ++column) {
        crosses[cross_count++] = {row, column};
      }
    }
    std::array<int, kEntries> coordinates{};
    std::iota(coordinates.begin(), coordinates.end(), 0);

    const auto started = Clock::now();
    const auto deadline =
        arguments.passes == 0
            ? started + std::chrono::duration_cast<Clock::duration>(
                            std::chrono::duration<double>(arguments.seconds))
            : Clock::time_point::max();
    auto next_heartbeat =
        started + std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(
                          arguments.heartbeat_seconds));
    auto elapsed = [&]() {
      return std::chrono::duration<double>(Clock::now() - started).count();
    };
    auto out_of_time = [&]() {
      return arguments.passes == 0 && Clock::now() >= deadline;
    };

    std::uint64_t total_assignments = 0;
    std::uint64_t completed_passes = 0;
    std::uint64_t basin_passes = 0;
    bool kick_pending = arguments.kick_size != 0;
    bool timed_out = false;
    log_record(
        log, arguments, best_score, state_score, state_score, 0.0, 0, 0,
        -1, -1, "start");

    auto kick_from_best = [&]() {
      for (int attempt = 0; attempt < 64; ++attempt) {
        state_matrix = best_matrix;
        std::shuffle(coordinates.begin(), coordinates.end(), randomizer);
        for (int index = 0; index < arguments.kick_size; ++index) {
          const int encoded = coordinates[index];
          state_matrix[encoded / kOrder][encoded % kOrder] *= -1;
        }
        state_score = absolute(exact_determinant(state_matrix));
        if (state_score != 0) return;
      }
      throw std::runtime_error("could not produce a nonsingular kick");
    };

    while (!out_of_time() &&
           (arguments.passes == 0 ||
            completed_passes < arguments.passes)) {
      if (kick_pending) {
        kick_from_best();
        kick_pending = false;
        log_record(
            log, arguments, best_score, state_score, state_score, elapsed(),
            total_assignments, completed_passes, -1, -1, "kick");
      }

      std::shuffle(crosses.begin(), crosses.end(), randomizer);
      std::uint64_t moves_in_pass = 0;
      bool pass_complete = true;
      for (const auto& [row, column] : crosses) {
        if (out_of_time()) {
          pass_complete = false;
          timed_out = true;
          break;
        }
        CrossResult result =
            optimize_cross(state_matrix, row, column, deadline);
        total_assignments += result.assignments;
        if (result.score > state_score) {
          state_matrix = result.matrix;
          state_score = result.score;
          ++moves_in_pass;
          const bool new_best = state_score > best_score;
          const bool new_elite =
              !new_best && state_score == best_score &&
              !sign_equivalent(state_matrix, best_matrix);
          if (new_best) {
            best_score = state_score;
            best_matrix = state_matrix;
            write_matrix(arguments.output, best_matrix);
          } else if (new_elite) {
            // Equal-score factorizations can have different cross
            // neighborhoods.  Retain the newest exact elite so the next
            // kicked basin does not always restart from the original one.
            best_matrix = state_matrix;
            write_matrix(arguments.output, best_matrix);
          }
          log_record(
              log, arguments, best_score, state_score, result.score,
              elapsed(), total_assignments, completed_passes, row, column,
              new_best ? "new_best"
                       : (new_elite ? "new_elite" : "move"));
          std::cout << (new_best ? "new best"
                                : (new_elite ? "new elite" : "move"))
                    << " |det|=" << wide_to_string(state_score)
                    << " cross=" << row << ',' << column
                    << " assignments=" << total_assignments << '\n'
                    << std::flush;
        }
        if (!result.complete) {
          pass_complete = false;
          timed_out = true;
          break;
        }
        if (arguments.heartbeat_seconds > 0 &&
            Clock::now() >= next_heartbeat) {
          log_record(
              log, arguments, best_score, state_score, result.score,
              elapsed(), total_assignments, completed_passes, row, column,
              "heartbeat");
          std::cout << "heartbeat best=" << wide_to_string(best_score)
                    << " state=" << wide_to_string(state_score)
                    << " passes=" << completed_passes
                    << " assignments=" << total_assignments << '\n'
                    << std::flush;
          next_heartbeat =
              Clock::now() +
              std::chrono::duration_cast<Clock::duration>(
                  std::chrono::duration<double>(
                      arguments.heartbeat_seconds));
        }
      }

      if (!pass_complete) break;
      ++completed_passes;
      ++basin_passes;
      log_record(
          log, arguments, best_score, state_score, state_score, elapsed(),
          total_assignments, completed_passes, -1, -1,
          moves_in_pass == 0 ? "basin_settled" : "pass_finished");

      if (arguments.kick_size == 0) {
        if (moves_in_pass == 0) break;
      } else if (
          moves_in_pass == 0 ||
          basin_passes >= arguments.settle_passes) {
        kick_pending = true;
        basin_passes = 0;
      }
    }

    log_record(
        log, arguments, best_score, state_score, best_score, elapsed(),
        total_assignments, completed_passes, -1, -1,
        timed_out || out_of_time() ? "time_limit" : "finished");
    std::cout << "finished best |det|=" << wide_to_string(best_score)
              << " state |det|=" << wide_to_string(state_score)
              << " assignments=" << total_assignments
              << " completed_passes=" << completed_passes
              << " elapsed=" << std::fixed << std::setprecision(3)
              << elapsed() << "s\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "cross_search: " << error.what() << '\n';
    return 2;
  }
}
