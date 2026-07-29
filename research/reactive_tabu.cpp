#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

namespace {

constexpr int kOrder = 23;
constexpr int kEntries = kOrder * kOrder;
constexpr std::uint64_t kRebuildInterval = 64;
constexpr int kMinimumTenure = 7;
constexpr int kMaximumTenure = 96;
constexpr std::size_t kVisitTableSize = 1U << 18U;
constexpr long double kSingularRatio = 1.0e-18L;
constexpr long double kAspirationProbeSlack = 1.0e-9L;

using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using Inverse =
    std::array<std::array<long double, kOrder>, kOrder>;
using Wide = __int128_t;
using Clock = std::chrono::steady_clock;

volatile std::sig_atomic_t stop_requested = 0;

void request_stop(int) { stop_requested = 1; }

struct State {
  Matrix matrix{};
  Inverse inverse{};
  long double log_abs_determinant =
      -std::numeric_limits<long double>::infinity();
  bool nonsingular = false;
  std::uint64_t moves_since_rebuild = 0;
};

struct Arguments {
  std::filesystem::path start;
  std::filesystem::path output;
  std::filesystem::path log;
  std::uint64_t seed = 23;
  double seconds = 3600.0;
  double heartbeat_seconds = 60.0;
};

struct Move {
  int row = -1;
  int column = -1;
  long double projected_log_abs_determinant =
      -std::numeric_limits<long double>::infinity();
  bool aspiration = false;
};

struct Visit {
  std::uint64_t hash = 0;
  std::uint64_t iteration = 0;
  bool occupied = false;
};

struct Statistics {
  std::uint64_t iterations = 0;
  std::uint64_t moves = 0;
  std::uint64_t inverse_rebuilds = 0;
  std::uint64_t exact_checks = 0;
  std::uint64_t aspiration_probes = 0;
  std::uint64_t cycles = 0;
  std::uint64_t singular_moves_rejected = 0;
};

Wide absolute(Wide value) { return value < 0 ? -value : value; }

std::string wide_to_string(Wide value) {
  if (value == 0) return "0";
  const bool negative = value < 0;
  if (negative) value = -value;
  std::string result;
  while (value != 0) {
    result.push_back(
        static_cast<char>('0' + static_cast<int>(value % 10)));
    value /= 10;
  }
  if (negative) result.push_back('-');
  std::reverse(result.begin(), result.end());
  return result;
}

Wide exact_determinant(const Matrix& matrix) {
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
    while (pivot_row < kOrder && work[pivot_row][column] == 0) {
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
  return static_cast<Wide>(sign) * work[kOrder - 1][kOrder - 1];
}

Matrix read_matrix(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("cannot open start matrix: " + path.string());
  }
  Matrix matrix{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      if (!(input >> matrix[row][column]) ||
          (matrix[row][column] != -1 && matrix[row][column] != 1)) {
        throw std::runtime_error(
            "start matrix must contain exactly 23x23 entries in {-1,+1}");
      }
    }
  }
  std::string extra;
  if (input >> extra) {
    throw std::runtime_error("start matrix contains extra data");
  }
  return matrix;
}

std::string matrix_bytes(const Matrix& matrix) {
  std::string bytes;
  bytes.reserve(1700);
  for (const auto& row : matrix) {
    for (int column = 0; column < kOrder; ++column) {
      if (column != 0) bytes.push_back(' ');
      bytes += row[column] == 1 ? "1" : "-1";
    }
    bytes.push_back('\n');
  }
  return bytes;
}

void write_all(int descriptor, std::string_view bytes) {
  std::size_t offset = 0;
  while (offset < bytes.size()) {
    const ssize_t written =
        ::write(descriptor, bytes.data() + offset, bytes.size() - offset);
    if (written < 0) {
      if (errno == EINTR) continue;
      throw std::runtime_error(
          "cannot write checkpoint: " + std::string(std::strerror(errno)));
    }
    if (written == 0) {
      throw std::runtime_error("short write while writing checkpoint");
    }
    offset += static_cast<std::size_t>(written);
  }
}

void sync_directory(const std::filesystem::path& directory) {
  int flags = O_RDONLY;
#ifdef O_DIRECTORY
  flags |= O_DIRECTORY;
#endif
#ifdef O_CLOEXEC
  flags |= O_CLOEXEC;
#endif
  const int descriptor = ::open(directory.c_str(), flags);
  if (descriptor < 0) {
    throw std::runtime_error(
        "cannot open checkpoint directory for sync: " +
        std::string(std::strerror(errno)));
  }
  const int result = ::fsync(descriptor);
  const int saved_errno = errno;
  ::close(descriptor);
  if (result != 0) {
    throw std::runtime_error(
        "cannot sync checkpoint directory: " +
        std::string(std::strerror(saved_errno)));
  }
}

void atomic_write_matrix(const std::filesystem::path& path,
                         const Matrix& matrix, std::uint64_t nonce) {
  const std::filesystem::path directory =
      path.parent_path().empty() ? std::filesystem::path(".")
                                 : path.parent_path();
  std::filesystem::create_directories(directory);
  const std::string bytes = matrix_bytes(matrix);

  std::filesystem::path temporary;
  int descriptor = -1;
  for (std::uint64_t attempt = 0; attempt < 64; ++attempt) {
    temporary =
        directory /
        ("." + path.filename().string() + ".reactive-tabu-" +
         std::to_string(static_cast<long long>(::getpid())) + "-" +
         std::to_string(nonce) + "-" + std::to_string(attempt) + ".tmp");
    int flags = O_WRONLY | O_CREAT | O_EXCL;
#ifdef O_CLOEXEC
    flags |= O_CLOEXEC;
#endif
    descriptor = ::open(temporary.c_str(), flags, 0600);
    if (descriptor >= 0) break;
    if (errno != EEXIST) {
      throw std::runtime_error(
          "cannot create checkpoint temporary file: " +
          std::string(std::strerror(errno)));
    }
  }
  if (descriptor < 0) {
    throw std::runtime_error(
        "cannot allocate a unique checkpoint temporary file");
  }

  bool renamed = false;
  try {
    write_all(descriptor, bytes);
    if (::fsync(descriptor) != 0) {
      throw std::runtime_error(
          "cannot sync checkpoint: " + std::string(std::strerror(errno)));
    }
    if (::close(descriptor) != 0) {
      descriptor = -1;
      throw std::runtime_error(
          "cannot close checkpoint: " + std::string(std::strerror(errno)));
    }
    descriptor = -1;
    if (::rename(temporary.c_str(), path.c_str()) != 0) {
      throw std::runtime_error(
          "cannot atomically install checkpoint: " +
          std::string(std::strerror(errno)));
    }
    renamed = true;
    sync_directory(directory);
  } catch (...) {
    if (descriptor >= 0) ::close(descriptor);
    if (!renamed) {
      std::error_code ignored;
      std::filesystem::remove(temporary, ignored);
    }
    throw;
  }
}

bool rebuild_inverse(State& state) {
  std::array<std::array<long double, 2 * kOrder>, kOrder> augmented{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      augmented[row][column] =
          static_cast<long double>(state.matrix[row][column]);
      augmented[row][column + kOrder] =
          row == column ? 1.0L : 0.0L;
    }
  }

  long double log_abs_determinant = 0.0L;
  for (int column = 0; column < kOrder; ++column) {
    int pivot_row = column;
    for (int row = column + 1; row < kOrder; ++row) {
      if (std::fabs(augmented[row][column]) >
          std::fabs(augmented[pivot_row][column])) {
        pivot_row = row;
      }
    }
    if (std::fabs(augmented[pivot_row][column]) < kSingularRatio) {
      state.nonsingular = false;
      state.log_abs_determinant =
          -std::numeric_limits<long double>::infinity();
      return false;
    }
    if (pivot_row != column) {
      std::swap(augmented[pivot_row], augmented[column]);
    }
    const long double pivot = augmented[column][column];
    log_abs_determinant += std::log(std::fabs(pivot));
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
      state.inverse[row][column] =
          augmented[row][column + kOrder];
    }
  }
  state.log_abs_determinant = log_abs_determinant;
  state.nonsingular = true;
  state.moves_since_rebuild = 0;
  return true;
}

long double flip_ratio(const State& state, int row, int column) {
  const long double delta =
      -2.0L * static_cast<long double>(state.matrix[row][column]);
  return 1.0L + delta * state.inverse[column][row];
}

bool apply_flip(State& state, int row, int column) {
  const long double delta =
      -2.0L * static_cast<long double>(state.matrix[row][column]);
  const long double ratio =
      1.0L + delta * state.inverse[column][row];
  if (std::fabs(ratio) < kSingularRatio) {
    state.matrix[row][column] *= -1;
    return rebuild_inverse(state);
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
          factor * inverse_column[inner_row] *
          inverse_row[inner_column];
    }
  }
  state.matrix[row][column] *= -1;
  state.log_abs_determinant += std::log(std::fabs(ratio));
  ++state.moves_since_rebuild;
  return true;
}

std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

std::array<std::uint64_t, kEntries> make_zobrist(std::uint64_t seed) {
  std::array<std::uint64_t, kEntries> values{};
  std::uint64_t state = seed ^ 0x4d61784465743233ULL;
  for (std::uint64_t& value : values) {
    state = splitmix64(state);
    value = state;
  }
  return values;
}

std::uint64_t matrix_hash(
    const Matrix& matrix,
    const std::array<std::uint64_t, kEntries>& zobrist) {
  std::uint64_t hash = 0;
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      if (matrix[row][column] == 1) {
        hash ^= zobrist[row * kOrder + column];
      }
    }
  }
  return hash;
}

std::uint64_t strict_unsigned(std::string_view text,
                              std::string_view option) {
  if (text.empty() ||
      !std::all_of(text.begin(), text.end(),
                   [](unsigned char character) {
                     return character >= '0' && character <= '9';
                   })) {
    throw std::runtime_error(
        std::string(option) + " must be a non-negative integer");
  }
  std::size_t consumed = 0;
  const std::uint64_t result =
      std::stoull(std::string(text), &consumed);
  if (consumed != text.size()) {
    throw std::runtime_error(
        std::string(option) + " must be a non-negative integer");
  }
  return result;
}

double strict_double(std::string_view text, std::string_view option,
                     bool allow_zero) {
  std::size_t consumed = 0;
  const double result = std::stod(std::string(text), &consumed);
  if (consumed != text.size() || !std::isfinite(result) ||
      result < 0.0 || (!allow_zero && result == 0.0)) {
    throw std::runtime_error(
        std::string(option) +
        (allow_zero ? " must be finite and non-negative"
                    : " must be finite and positive"));
  }
  return result;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    auto value = [&]() -> std::string_view {
      ++index;
      if (index >= argc) {
        throw std::runtime_error("missing value for " + option);
      }
      return argv[index];
    };
    if (option == "--start") {
      arguments.start = value();
    } else if (option == "--output") {
      arguments.output = value();
    } else if (option == "--log") {
      arguments.log = value();
    } else if (option == "--seed") {
      arguments.seed = strict_unsigned(value(), option);
    } else if (option == "--seconds") {
      arguments.seconds = strict_double(value(), option, false);
    } else if (option == "--heartbeat" ||
               option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds =
          strict_double(value(), option, true);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (arguments.start.empty()) {
    throw std::runtime_error("--start is required");
  }
  if (arguments.output.empty()) {
    throw std::runtime_error("--output is required");
  }
  if (arguments.log.empty()) {
    throw std::runtime_error("--log is required");
  }
  return arguments;
}

void log_record(std::ofstream& log, const Arguments& arguments,
                const Statistics& statistics, const char* event,
                double elapsed_seconds, Wide best_score,
                Wide last_exact_score, int tenure) {
  log << "{\"absolute_determinant\":\""
      << wide_to_string(best_score)
      << "\",\"aspiration_probes\":"
      << statistics.aspiration_probes
      << ",\"cycles\":" << statistics.cycles
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(3) << elapsed_seconds
      << ",\"event\":\"" << event
      << "\",\"exact_checks\":" << statistics.exact_checks
      << ",\"inverse_rebuilds\":"
      << statistics.inverse_rebuilds
      << ",\"iterations\":" << statistics.iterations
      << ",\"last_exact_determinant\":\""
      << wide_to_string(last_exact_score)
      << "\",\"moves\":" << statistics.moves
      << ",\"seed\":" << arguments.seed
      << ",\"singular_moves_rejected\":"
      << statistics.singular_moves_rejected
      << ",\"tenure\":" << tenure << "}\n";
  log.flush();
  if (!log) throw std::runtime_error("cannot append research log");
}

bool exact_check_and_promote(
    const State& state, Matrix& best_matrix, Wide& best_score,
    Wide& last_exact_score, const Arguments& arguments,
    Statistics& statistics, std::ofstream& log, const Clock::time_point& started,
    int tenure, std::uint64_t checkpoint_nonce) {
  const Wide exact_score = absolute(exact_determinant(state.matrix));
  ++statistics.exact_checks;
  last_exact_score = exact_score;
  if (exact_score <= best_score) return false;

  // The exact Bareiss result above is the sole checkpoint ranking gate.
  best_score = exact_score;
  best_matrix = state.matrix;
  atomic_write_matrix(arguments.output, best_matrix, checkpoint_nonce);
  const double elapsed =
      std::chrono::duration<double>(Clock::now() - started).count();
  log_record(log, arguments, statistics, "new_best", elapsed, best_score,
             last_exact_score, tenure);
  std::cout << "new best |det|=" << wide_to_string(best_score)
            << " iteration=" << statistics.iterations << '\n'
            << std::flush;
  return true;
}

Move choose_move(
    const State& state,
    const std::array<std::uint64_t, kEntries>& tabu_until,
    std::uint64_t iteration, Wide best_score, std::mt19937_64& randomizer,
    Statistics& statistics) {
  Move best;
  const long double best_exact_log =
      std::log(static_cast<long double>(best_score));
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      const int coordinate = row * kOrder + column;
      const long double ratio = flip_ratio(state, row, column);
      const long double magnitude = std::fabs(ratio);
      const long double projected =
          magnitude < kSingularRatio
              ? -std::numeric_limits<long double>::infinity()
              : state.log_abs_determinant + std::log(magnitude);
      const bool tabu = iteration < tabu_until[coordinate];
      bool aspiration = false;
      if (tabu &&
          projected + kAspirationProbeSlack >= best_exact_log) {
        Matrix candidate = state.matrix;
        candidate[row][column] *= -1;
        ++statistics.aspiration_probes;
        ++statistics.exact_checks;
        aspiration =
            absolute(exact_determinant(candidate)) > best_score;
      }
      if (tabu && !aspiration) continue;

      const bool better =
          best.row < 0 ||
          projected > best.projected_log_abs_determinant;
      const bool tied =
          best.row >= 0 &&
          std::fabs(projected -
                    best.projected_log_abs_determinant) <=
              1.0e-18L;
      if (better || (tied && (randomizer() & 1U) != 0U)) {
        best = Move{row, column, projected, aspiration};
      }
    }
  }
  if (best.row < 0) {
    throw std::runtime_error("no admissible tabu move");
  }
  return best;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    if (!arguments.log.parent_path().empty()) {
      std::filesystem::create_directories(arguments.log.parent_path());
    }
    std::ofstream log(arguments.log, std::ios::app);
    if (!log) throw std::runtime_error("cannot open research log");

    State state;
    state.matrix = read_matrix(arguments.start);
    if (!rebuild_inverse(state)) {
      throw std::runtime_error("start matrix must be nonsingular");
    }

    Statistics statistics;
    statistics.inverse_rebuilds = 1;
    Matrix best_matrix = state.matrix;
    Wide best_score = absolute(exact_determinant(state.matrix));
    ++statistics.exact_checks;
    if (best_score == 0) {
      throw std::runtime_error("start matrix must have nonzero determinant");
    }
    Wide last_exact_score = best_score;
    std::uint64_t checkpoint_nonce = 0;
    atomic_write_matrix(
        arguments.output, best_matrix, checkpoint_nonce++);

    std::mt19937_64 randomizer(arguments.seed);
    const auto zobrist = make_zobrist(arguments.seed);
    std::uint64_t current_hash = matrix_hash(state.matrix, zobrist);
    std::vector<Visit> visits(kVisitTableSize);
    visits[current_hash & (kVisitTableSize - 1U)] =
        Visit{current_hash, 0, true};
    std::array<std::uint64_t, kEntries> tabu_until{};
    const int baseline_tenure =
        kMinimumTenure + static_cast<int>(arguments.seed % 5U);
    int tenure = baseline_tenure;
    std::uint64_t last_cycle_iteration = 0;

    const auto started = Clock::now();
    const auto deadline =
        started + std::chrono::duration<double>(arguments.seconds);
    auto next_heartbeat =
        started +
        std::chrono::duration<double>(arguments.heartbeat_seconds);

    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);
    log_record(log, arguments, statistics, "start", 0.0, best_score,
               last_exact_score, tenure);

    while (!stop_requested && Clock::now() < deadline) {
      ++statistics.iterations;
      const Move move = choose_move(
          state, tabu_until, statistics.iterations, best_score,
          randomizer, statistics);
      const int coordinate = move.row * kOrder + move.column;
      const Matrix before = state.matrix;

      if (!apply_flip(state, move.row, move.column)) {
        state.matrix = before;
        if (!rebuild_inverse(state)) {
          throw std::runtime_error(
              "cannot recover inverse after a singular move");
        }
        ++statistics.inverse_rebuilds;
        ++statistics.singular_moves_rejected;
        tabu_until[coordinate] =
            statistics.iterations +
            static_cast<std::uint64_t>(kMaximumTenure);
        continue;
      }
      ++statistics.moves;
      current_hash ^= zobrist[coordinate];
      const int jitter = static_cast<int>(randomizer() % 5U);
      tabu_until[coordinate] =
          statistics.iterations +
          static_cast<std::uint64_t>(tenure + jitter + 1);

      Visit& visit =
          visits[current_hash & (kVisitTableSize - 1U)];
      if (visit.occupied && visit.hash == current_hash &&
          statistics.iterations > visit.iteration) {
        const std::uint64_t cycle_length =
            statistics.iterations - visit.iteration;
        if (cycle_length <=
            static_cast<std::uint64_t>(4 * kMaximumTenure)) {
          ++statistics.cycles;
          last_cycle_iteration = statistics.iterations;
          tenure = std::min(
              kMaximumTenure,
              tenure + 2 +
                  static_cast<int>(
                      std::min<std::uint64_t>(cycle_length / 8U, 8U)));
        }
      }
      visit = Visit{current_hash, statistics.iterations, true};
      if (statistics.iterations - last_cycle_iteration >= 512U &&
          (statistics.iterations & 127U) == 0U &&
          tenure > baseline_tenure) {
        --tenure;
      }

      const long double best_exact_log =
          std::log(static_cast<long double>(best_score));
      const bool near_best =
          move.aspiration ||
          state.log_abs_determinant + kAspirationProbeSlack >=
              best_exact_log;
      if (near_best) {
        exact_check_and_promote(
            state, best_matrix, best_score, last_exact_score, arguments,
            statistics, log, started, tenure, checkpoint_nonce++);
      }

      if (state.moves_since_rebuild >= kRebuildInterval) {
        if (!rebuild_inverse(state)) {
          throw std::runtime_error(
              "current matrix became singular during inverse rebuild");
        }
        ++statistics.inverse_rebuilds;
        exact_check_and_promote(
            state, best_matrix, best_score, last_exact_score, arguments,
            statistics, log, started, tenure, checkpoint_nonce++);
      }

      const auto now = Clock::now();
      if (arguments.heartbeat_seconds > 0.0 &&
          now >= next_heartbeat) {
        const double elapsed =
            std::chrono::duration<double>(now - started).count();
        log_record(log, arguments, statistics, "heartbeat", elapsed,
                   best_score, last_exact_score, tenure);
        next_heartbeat =
            now +
            std::chrono::duration<double>(
                arguments.heartbeat_seconds);
      }
    }

    if (!rebuild_inverse(state)) {
      throw std::runtime_error(
          "current matrix is singular at final inverse rebuild");
    }
    ++statistics.inverse_rebuilds;
    exact_check_and_promote(
        state, best_matrix, best_score, last_exact_score, arguments,
        statistics, log, started, tenure, checkpoint_nonce++);
    const double elapsed =
        std::chrono::duration<double>(Clock::now() - started).count();
    log_record(
        log, arguments, statistics,
        stop_requested ? "stopped" : "finished", elapsed, best_score,
        last_exact_score, tenure);
    std::cout << "finished |det|=" << wide_to_string(best_score)
              << " iterations=" << statistics.iterations
              << " moves=" << statistics.moves
              << " tenure=" << tenure << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "reactive_tabu: " << error.what() << '\n';
    return 2;
  }
}
