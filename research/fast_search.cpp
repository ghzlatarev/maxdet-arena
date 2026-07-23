#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int kOrder = 23;
using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using Wide = __int128_t;

struct State {
  Matrix matrix{};
  std::array<std::array<long double, kOrder>, kOrder> inverse{};
  long double log_abs_det = -std::numeric_limits<long double>::infinity();
  bool nonsingular = false;
  std::uint64_t accepted_since_rebuild = 0;
};

std::string wide_to_string(Wide value) {
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

Wide exact_determinant(Matrix matrix) {
  std::array<std::array<Wide, kOrder>, kOrder> work{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      work[row][column] = matrix[row][column];
    }
  }

  Wide previous_pivot = 1;
  int sign = 1;
  for (int column = 0; column < kOrder - 1; ++column) {
    int pivot_row = column;
    while (pivot_row < kOrder && work[pivot_row][column] == 0) ++pivot_row;
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
  return sign * work[kOrder - 1][kOrder - 1];
}

Wide absolute(Wide value) { return value < 0 ? -value : value; }

Matrix read_matrix(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open start matrix");
  Matrix matrix{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      if (!(input >> matrix[row][column]) ||
          (matrix[row][column] != -1 && matrix[row][column] != 1)) {
        throw std::runtime_error("invalid start matrix");
      }
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
  const auto temporary = path.string() + ".tmp";
  {
    std::ofstream output(temporary, std::ios::trunc);
    if (!output) throw std::runtime_error("cannot create output matrix");
    for (const auto& row : matrix) {
      for (int column = 0; column < kOrder; ++column) {
        if (column) output << ' ';
        output << row[column];
      }
      output << '\n';
    }
    output.flush();
    if (!output) throw std::runtime_error("cannot flush output matrix");
  }
  std::filesystem::rename(temporary, path);
}

bool rebuild(State& state) {
  std::array<std::array<long double, 2 * kOrder>, kOrder> augmented{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      augmented[row][column] = state.matrix[row][column];
      augmented[row][column + kOrder] = row == column ? 1.0L : 0.0L;
    }
  }

  long double log_abs_det = 0.0L;
  for (int column = 0; column < kOrder; ++column) {
    int pivot_row = column;
    for (int row = column + 1; row < kOrder; ++row) {
      if (std::fabs(augmented[row][column]) >
          std::fabs(augmented[pivot_row][column])) {
        pivot_row = row;
      }
    }
    const long double pivot_candidate = augmented[pivot_row][column];
    if (std::fabs(pivot_candidate) < 1e-24L) {
      state.nonsingular = false;
      state.log_abs_det = -std::numeric_limits<long double>::infinity();
      return false;
    }
    if (pivot_row != column) {
      std::swap(augmented[pivot_row], augmented[column]);
    }
    const long double pivot = augmented[column][column];
    log_abs_det += std::log(std::fabs(pivot));
    for (int inner = 0; inner < 2 * kOrder; ++inner) {
      augmented[column][inner] /= pivot;
    }
    for (int row = 0; row < kOrder; ++row) {
      if (row == column) continue;
      const long double factor = augmented[row][column];
      if (factor == 0.0L) continue;
      for (int inner = 0; inner < 2 * kOrder; ++inner) {
        augmented[row][inner] -= factor * augmented[column][inner];
      }
    }
  }

  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      state.inverse[row][column] = augmented[row][column + kOrder];
    }
  }
  state.nonsingular = true;
  state.log_abs_det = log_abs_det;
  state.accepted_since_rebuild = 0;
  return true;
}

long double flip_log_delta(const State& state, int row, int column) {
  if (!state.nonsingular) {
    return std::numeric_limits<long double>::infinity();
  }
  const long double delta = -2.0L * state.matrix[row][column];
  const long double ratio = 1.0L + delta * state.inverse[column][row];
  if (std::fabs(ratio) < 1e-24L) {
    return -std::numeric_limits<long double>::infinity();
  }
  return std::log(std::fabs(ratio));
}

bool apply_flip(State& state, int row, int column) {
  if (!state.nonsingular) {
    state.matrix[row][column] *= -1;
    return rebuild(state);
  }

  const long double delta = -2.0L * state.matrix[row][column];
  const long double ratio = 1.0L + delta * state.inverse[column][row];
  if (std::fabs(ratio) < 1e-20L) {
    state.matrix[row][column] *= -1;
    return rebuild(state);
  }

  std::array<long double, kOrder> inverse_column{};
  std::array<long double, kOrder> inverse_row{};
  for (int index = 0; index < kOrder; ++index) {
    inverse_column[index] = state.inverse[index][row];
    inverse_row[index] = state.inverse[column][index];
  }
  const long double factor = delta / ratio;
  for (int inner_row = 0; inner_row < kOrder; ++inner_row) {
    for (int inner_column = 0; inner_column < kOrder; ++inner_column) {
      state.inverse[inner_row][inner_column] -=
          factor * inverse_column[inner_row] * inverse_row[inner_column];
    }
  }
  state.matrix[row][column] *= -1;
  state.log_abs_det += std::log(std::fabs(ratio));
  ++state.accepted_since_rebuild;
  if (state.accepted_since_rebuild >= 128) return rebuild(state);
  return true;
}

Matrix random_matrix(std::mt19937_64& randomizer) {
  Matrix matrix{};
  for (auto& row : matrix) {
    for (int& value : row) value = (randomizer() & 1U) ? 1 : -1;
  }
  return matrix;
}

struct Arguments {
  std::filesystem::path start;
  std::filesystem::path output;
  std::filesystem::path log;
  std::string mode = "hybrid";
  std::uint64_t seed = 23;
  double seconds = 3600.0;
  double heartbeat_seconds = 60.0;
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
    else if (option == "--log") arguments.log = value();
    else if (option == "--mode") arguments.mode = value();
    else if (option == "--seed") arguments.seed = std::stoull(value());
    else if (option == "--seconds") arguments.seconds = std::stod(value());
    else if (option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds = std::stod(value());
    }
    else throw std::runtime_error("unknown option: " + option);
  }
  if (arguments.output.empty()) throw std::runtime_error("--output is required");
  if (arguments.log.empty()) throw std::runtime_error("--log is required");
  if (!std::isfinite(arguments.seconds) || arguments.seconds <= 0) {
    throw std::runtime_error("--seconds must be finite and positive");
  }
  if (!std::isfinite(arguments.heartbeat_seconds) ||
      arguments.heartbeat_seconds < 0) {
    throw std::runtime_error(
        "--heartbeat-seconds must be finite and non-negative");
  }
  if (arguments.mode != "hill" && arguments.mode != "anneal" &&
      arguments.mode != "hybrid") {
    throw std::runtime_error("--mode must be hill, anneal, or hybrid");
  }
  return arguments;
}

void log_record(std::ofstream& log, const Arguments& arguments,
                std::uint64_t epoch, Wide score, const char* event,
                double elapsed_seconds, std::uint64_t accepted) {
  log << "{\"absolute_determinant\":\"" << wide_to_string(score)
      << "\",\"accepted\":" << accepted
      << ",\"elapsed_seconds\":" << std::fixed << std::setprecision(3)
      << elapsed_seconds
      << ",\"epoch\":" << epoch
      << ",\"event\":\"" << event
      << "\",\"mode\":\"" << arguments.mode
      << "\",\"seed\":" << arguments.seed << "}\n";
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

    State state;
    state.matrix = arguments.start.empty()
                       ? random_matrix(randomizer)
                       : read_matrix(arguments.start);
    rebuild(state);

    Matrix best_matrix = state.matrix;
    Wide best_score = absolute(exact_determinant(best_matrix));
    write_matrix(arguments.output, best_matrix);

    const auto started = std::chrono::steady_clock::now();
    const auto deadline =
        started + std::chrono::duration<double>(arguments.seconds);
    auto next_heartbeat =
        started + std::chrono::duration<double>(arguments.heartbeat_seconds);
    std::uniform_real_distribution<long double> unit(0.0L, 1.0L);
    std::uniform_int_distribution<int> coordinate(0, kOrder - 1);
    std::uint64_t epoch = 0;
    std::uint64_t accepted = 0;
    log_record(log, arguments, 0, best_score, "start", 0.0, accepted);

    while (std::chrono::steady_clock::now() < deadline) {
      ++epoch;
      if (!state.nonsingular) {
        state.matrix = random_matrix(randomizer);
        rebuild(state);
        continue;
      }

      if (arguments.mode == "anneal") {
        const int row = coordinate(randomizer);
        const int column = coordinate(randomizer);
        const long double change = flip_log_delta(state, row, column);
        const long double elapsed = std::chrono::duration<long double>(
            std::chrono::steady_clock::now() - started).count();
        const long double phase = std::fmod(elapsed, 8.0L) / 8.0L;
        const long double temperature = 0.25L * (1.0L - phase) + 0.005L;
        if (change > 0.0L || std::log(unit(randomizer)) < change / temperature) {
          apply_flip(state, row, column);
          ++accepted;
        }
      } else {
        long double best_change = 1e-16L;
        int best_row = -1;
        int best_column = -1;
        for (int row = 0; row < kOrder; ++row) {
          for (int column = 0; column < kOrder; ++column) {
            const long double change = flip_log_delta(state, row, column);
            if (change > best_change) {
              best_change = change;
              best_row = row;
              best_column = column;
            }
          }
        }
        if (best_row >= 0) {
          apply_flip(state, best_row, best_column);
          ++accepted;
        } else {
          const Wide local_score = absolute(exact_determinant(state.matrix));
          if (local_score > best_score) {
            best_score = local_score;
            best_matrix = state.matrix;
            write_matrix(arguments.output, best_matrix);
            const double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
            log_record(
                log, arguments, epoch, best_score, "new_best", elapsed, accepted);
            std::cout << "new best |det|=" << wide_to_string(best_score)
                      << " epoch=" << epoch << '\n' << std::flush;
          }

          if (arguments.mode == "hill") {
            state.matrix = random_matrix(randomizer);
            rebuild(state);
          } else {
            const int kick_size = 2 + static_cast<int>(randomizer() % 7);
            for (int kick = 0; kick < kick_size; ++kick) {
              apply_flip(state, coordinate(randomizer), coordinate(randomizer));
            }
            rebuild(state);
          }
        }
      }

      if ((epoch & 4095U) == 0U) {
        rebuild(state);
        const Wide current_score = absolute(exact_determinant(state.matrix));
        if (current_score > best_score) {
          best_score = current_score;
          best_matrix = state.matrix;
          write_matrix(arguments.output, best_matrix);
          const double elapsed = std::chrono::duration<double>(
              std::chrono::steady_clock::now() - started).count();
          log_record(
              log, arguments, epoch, best_score, "new_best", elapsed, accepted);
          std::cout << "new best |det|=" << wide_to_string(best_score)
                    << " epoch=" << epoch << '\n' << std::flush;
        }
      }

      const auto now = std::chrono::steady_clock::now();
      if (arguments.heartbeat_seconds > 0 && now >= next_heartbeat) {
        const double elapsed =
            std::chrono::duration<double>(now - started).count();
        log_record(
            log, arguments, epoch, best_score, "heartbeat", elapsed, accepted);
        next_heartbeat =
            now + std::chrono::duration<double>(arguments.heartbeat_seconds);
      }
    }

    rebuild(state);
    const Wide final_score = absolute(exact_determinant(state.matrix));
    if (final_score > best_score) {
      best_score = final_score;
      best_matrix = state.matrix;
      write_matrix(arguments.output, best_matrix);
      const double elapsed = std::chrono::duration<double>(
          std::chrono::steady_clock::now() - started).count();
      log_record(
          log, arguments, epoch, best_score, "new_best", elapsed, accepted);
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    log_record(
        log, arguments, epoch, best_score, "finished", elapsed, accepted);
    std::cout << "finished |det|=" << wide_to_string(best_score)
              << " epochs=" << epoch << " accepted=" << accepted << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "fast_search: " << error.what() << '\n';
    return 2;
  }
}
