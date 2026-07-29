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
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

// Exact reactive tabu search in the dephased 22 by 22 binary core.
//
// Every nonsingular 23 by 23 sign matrix is switching-equivalent to
//
//       [ 1   1 ]
//   H = [ 1  J-2B ],
//
// where B is a 22 by 22 {0,1}-matrix and
//
//   det(H) = 2^22 det(B).
//
// Thus the order-23 determinant quotient is det(B) itself.  Hadamard's
// inequality bounds |det(B)| by 1,089,457,290, and every entry of adj(B)
// (a 21 by 21 binary determinant) by 278,624,678.  The engine maintains
// det(B) and adj(B) exactly in int64_t.  A rank-one update B' = B + u v^T
// has
//
//   d'   = d + v^T adj(B) u,
//   adj(B') = (d' adj(B) - adj(B) u v^T adj(B)) / d.
//
// The numerator divisions are exact.  __int128 intermediates remain far
// below their range.  Besides the 484 core-bit flips, the neighborhood
// includes 22 full-row complements, 22 full-column complements, and one
// whole-core complement.  Those complement moves are the dephased images of
// flipping entries in the sign-matrix border, including the top-left entry,
// so the gauge quotient retains all 529 original single-entry directions.

namespace {

namespace fs = std::filesystem;

constexpr int kSignOrder = 23;
constexpr int kCoreOrder = 22;
constexpr int kCoreEntries = kCoreOrder * kCoreOrder;
constexpr int kMoveCount = kCoreEntries + 2 * kCoreOrder + 1;
constexpr std::int64_t kDeterminantBound = INT64_C(1089457290);
constexpr std::int64_t kCofactorBound = INT64_C(278624678);
constexpr std::int64_t kScale = INT64_C(1) << 22;
constexpr int kMinimumTenure = 7;
constexpr int kMaximumTenure = 96;
constexpr std::size_t kVisitTableSize = 1U << 18U;
constexpr std::uint64_t kIdentityCheckInterval = UINT64_C(4096);
constexpr std::uint64_t kDeterminantCheckInterval = UINT64_C(65536);

using SignMatrix =
    std::array<std::array<int, kSignOrder>, kSignOrder>;
using CoreMatrix =
    std::array<std::array<std::uint8_t, kCoreOrder>, kCoreOrder>;
using Adjugate =
    std::array<std::array<std::int64_t, kCoreOrder>, kCoreOrder>;
using Vector = std::array<std::int64_t, kCoreOrder>;
using Wide = __int128_t;
using Clock = std::chrono::steady_clock;

volatile std::sig_atomic_t stop_requested = 0;

void request_stop(int) { stop_requested = 1; }

struct Arguments {
  fs::path start;
  fs::path output;
  fs::path log;
  std::uint64_t seed = 23;
  double seconds = 3600.0;
  double heartbeat_seconds = 30.0;
  std::uint64_t max_iterations = 0;
  std::uint64_t breakout_interval = 0;
  int breakout_flips = 12;
  int breakout_attempts = 64;
};

struct State {
  CoreMatrix core{};
  Adjugate adjugate{};
  std::int64_t determinant = 0;
};

enum class MoveKind : std::uint8_t {
  kBit,
  kRowComplement,
  kColumnComplement,
  kWholeComplement,
};

struct Move {
  MoveKind kind = MoveKind::kBit;
  int first = -1;
  int second = -1;
  int id = -1;
  std::int64_t determinant = 0;
  bool aspiration = false;
};

struct Visit {
  std::uint64_t hash = 0;
  std::uint64_t iteration = 0;
  bool occupied = false;
};

struct Statistics {
  std::uint64_t iterations = 0;
  std::uint64_t candidate_evaluations = 0;
  std::uint64_t singular_candidates = 0;
  std::uint64_t bit_moves = 0;
  std::uint64_t row_complements = 0;
  std::uint64_t column_complements = 0;
  std::uint64_t whole_complements = 0;
  std::uint64_t cycles = 0;
  std::uint64_t identity_checks = 0;
  std::uint64_t determinant_checks = 0;
  std::uint64_t promotions = 0;
  std::uint64_t breakouts = 0;
  std::uint64_t breakout_attempts = 0;
  std::uint64_t breakout_singular = 0;
};

std::uint64_t magnitude(std::int64_t value) {
  return value < 0 ? static_cast<std::uint64_t>(-value)
                   : static_cast<std::uint64_t>(value);
}

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

std::string json_escape(std::string_view input) {
  std::ostringstream output;
  for (char character : input) {
    switch (character) {
      case '"':
        output << "\\\"";
        break;
      case '\\':
        output << "\\\\";
        break;
      case '\b':
        output << "\\b";
        break;
      case '\f':
        output << "\\f";
        break;
      case '\n':
        output << "\\n";
        break;
      case '\r':
        output << "\\r";
        break;
      case '\t':
        output << "\\t";
        break;
      default:
        if (static_cast<unsigned char>(character) < 0x20U) {
          output << "\\u" << std::hex << std::setw(4)
                 << std::setfill('0')
                 << static_cast<unsigned>(
                        static_cast<unsigned char>(character))
                 << std::dec;
        } else {
          output << character;
        }
    }
  }
  return output.str();
}

std::uint64_t strict_unsigned(std::string_view text,
                              std::string_view option) {
  if (text.empty() ||
      !std::all_of(
          text.begin(), text.end(),
          [](unsigned char character) {
            return character >= static_cast<unsigned char>('0') &&
                   character <= static_cast<unsigned char>('9');
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
    } else if (option == "--max-iterations") {
      arguments.max_iterations = strict_unsigned(value(), option);
    } else if (option == "--breakout-interval") {
      arguments.breakout_interval =
          strict_unsigned(value(), option);
    } else if (option == "--breakout-flips") {
      const std::uint64_t parsed = strict_unsigned(value(), option);
      if (parsed == 0 ||
          parsed > static_cast<std::uint64_t>(kCoreEntries)) {
        throw std::runtime_error(
            "--breakout-flips must be in [1,484]");
      }
      arguments.breakout_flips = static_cast<int>(parsed);
    } else if (option == "--breakout-attempts") {
      const std::uint64_t parsed = strict_unsigned(value(), option);
      if (parsed == 0 ||
          parsed >
              static_cast<std::uint64_t>(
                  std::numeric_limits<int>::max())) {
        throw std::runtime_error(
            "--breakout-attempts must be in [1,INT_MAX]");
      }
      arguments.breakout_attempts = static_cast<int>(parsed);
    } else if (option == "--help") {
      std::cout
          << "usage: core_adjugate_tabu --start MATRIX --output MATRIX "
             "--log JSONL [--seed N] [--seconds S] "
             "[--heartbeat-seconds S] [--max-iterations N] "
             "[--breakout-interval N] [--breakout-flips N] "
             "[--breakout-attempts N]\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (arguments.start.empty() || arguments.output.empty() ||
      arguments.log.empty()) {
    throw std::runtime_error(
        "--start, --output, and --log are required");
  }
  if (arguments.start == arguments.output ||
      arguments.start == arguments.log ||
      arguments.output == arguments.log) {
    throw std::runtime_error("input and output paths must be distinct");
  }
  return arguments;
}

SignMatrix read_sign_matrix(const fs::path& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("cannot open matrix: " + path.string());
  }
  SignMatrix matrix{};
  std::string line;
  for (int row = 0; row < kSignOrder; ++row) {
    if (!std::getline(input, line)) {
      throw std::runtime_error("matrix has fewer than 23 rows");
    }
    std::istringstream row_input(line);
    for (int column = 0; column < kSignOrder; ++column) {
      if (!(row_input >> matrix[row][column]) ||
          (matrix[row][column] != -1 &&
           matrix[row][column] != 1)) {
        throw std::runtime_error("invalid sign-matrix entry");
      }
    }
    std::string extra;
    if (row_input >> extra) {
      throw std::runtime_error(
          "matrix row has more than 23 entries");
    }
  }
  while (std::getline(input, line)) {
    if (line.find_first_not_of(" \t\r\n") != std::string::npos) {
      throw std::runtime_error("matrix has more than 23 rows");
    }
  }
  return matrix;
}

CoreMatrix dephase_to_core(const SignMatrix& matrix) {
  CoreMatrix core{};
  for (int row = 1; row < kSignOrder; ++row) {
    for (int column = 1; column < kSignOrder; ++column) {
      const int dephased =
          matrix[row][column] * matrix[row][0] *
          matrix[0][column] * matrix[0][0];
      core[row - 1][column - 1] =
          static_cast<std::uint8_t>(dephased == -1 ? 1 : 0);
    }
  }
  return core;
}

SignMatrix core_to_sign(const CoreMatrix& core) {
  SignMatrix matrix{};
  for (auto& row : matrix) row.fill(1);
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      matrix[row + 1][column + 1] =
          core[row][column] == 0U ? 1 : -1;
    }
  }
  return matrix;
}

Wide bareiss(const std::vector<std::vector<Wide>>& source) {
  if (source.empty()) return 1;
  std::vector<std::vector<Wide>> work = source;
  const int order = static_cast<int>(work.size());
  Wide previous = 1;
  int sign = 1;
  for (int column = 0; column < order - 1; ++column) {
    int pivot = column;
    while (pivot < order && work[pivot][column] == 0) ++pivot;
    if (pivot == order) return 0;
    if (pivot != column) {
      std::swap(work[pivot], work[column]);
      sign = -sign;
    }
    const Wide pivot_value = work[column][column];
    for (int row = column + 1; row < order; ++row) {
      for (int inner = column + 1; inner < order; ++inner) {
        const Wide numerator =
            work[row][inner] * pivot_value -
            work[row][column] * work[column][inner];
        if (numerator % previous != 0) {
          throw std::runtime_error("non-exact Bareiss division");
        }
        work[row][inner] = numerator / previous;
      }
      work[row][column] = 0;
    }
    previous = pivot_value;
  }
  return static_cast<Wide>(sign) * work.back().back();
}

std::vector<std::vector<Wide>> core_as_wide(
    const CoreMatrix& core, int omitted_row = -1,
    int omitted_column = -1) {
  std::vector<std::vector<Wide>> result;
  result.reserve(
      static_cast<std::size_t>(
          kCoreOrder - (omitted_row >= 0 ? 1 : 0)));
  for (int row = 0; row < kCoreOrder; ++row) {
    if (row == omitted_row) continue;
    std::vector<Wide> output_row;
    output_row.reserve(
        static_cast<std::size_t>(
            kCoreOrder - (omitted_column >= 0 ? 1 : 0)));
    for (int column = 0; column < kCoreOrder; ++column) {
      if (column == omitted_column) continue;
      output_row.push_back(core[row][column]);
    }
    result.push_back(std::move(output_row));
  }
  return result;
}

std::int64_t checked_narrow(Wide value, std::int64_t bound,
                            std::string_view name) {
  if (value < -static_cast<Wide>(bound) ||
      value > static_cast<Wide>(bound)) {
    throw std::runtime_error(
        std::string(name) + " exceeded proved bound: " +
        wide_to_string(value));
  }
  return static_cast<std::int64_t>(value);
}

std::int64_t exact_core_determinant(const CoreMatrix& core) {
  return checked_narrow(
      bareiss(core_as_wide(core)), kDeterminantBound,
      "core determinant");
}

Adjugate exact_adjugate(const CoreMatrix& core) {
  Adjugate result{};
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      Wide cofactor =
          bareiss(core_as_wide(core, row, column));
      if (((row + column) & 1) != 0) cofactor = -cofactor;
      result[column][row] = checked_narrow(
          cofactor, kCofactorBound, "core cofactor");
    }
  }
  return result;
}

void check_adjugate_identity(const State& state) {
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      Wide value = 0;
      for (int inner = 0; inner < kCoreOrder; ++inner) {
        value += static_cast<Wide>(state.core[row][inner]) *
                 state.adjugate[inner][column];
      }
      const Wide expected =
          row == column ? state.determinant : 0;
      if (value != expected) {
        throw std::runtime_error(
            "B*adj(B) invariant failed at (" +
            std::to_string(row) + "," +
            std::to_string(column) + ")");
      }
    }
  }
}

std::string sign_matrix_text(const SignMatrix& matrix) {
  std::ostringstream output;
  for (int row = 0; row < kSignOrder; ++row) {
    for (int column = 0; column < kSignOrder; ++column) {
      if (column != 0) output << ' ';
      output << matrix[row][column];
    }
    output << '\n';
  }
  return output.str();
}

void sync_directory(const fs::path& directory) {
  const int descriptor =
      ::open(directory.c_str(), O_RDONLY | O_DIRECTORY);
  if (descriptor < 0) {
    throw std::runtime_error(
        "cannot open output directory for sync: " +
        std::string(std::strerror(errno)));
  }
  if (::fsync(descriptor) != 0) {
    const int saved = errno;
    ::close(descriptor);
    throw std::runtime_error(
        "cannot sync output directory: " +
        std::string(std::strerror(saved)));
  }
  if (::close(descriptor) != 0) {
    throw std::runtime_error(
        "cannot close output directory: " +
        std::string(std::strerror(errno)));
  }
}

void atomic_write(const fs::path& path, std::string_view content,
                  std::uint64_t nonce) {
  const fs::path directory =
      path.parent_path().empty() ? fs::path(".") : path.parent_path();
  fs::create_directories(directory);
  const fs::path temporary =
      directory /
      (path.filename().string() + ".tmp." +
       std::to_string(static_cast<unsigned long>(::getpid())) +
       "." + std::to_string(nonce));
  int descriptor = -1;
  bool renamed = false;
  try {
    descriptor = ::open(
        temporary.c_str(), O_WRONLY | O_CREAT | O_EXCL,
        S_IRUSR | S_IWUSR);
    if (descriptor < 0) {
      throw std::runtime_error(
          "cannot create atomic output: " +
          std::string(std::strerror(errno)));
    }
    std::size_t offset = 0;
    while (offset < content.size()) {
      const ssize_t written = ::write(
          descriptor, content.data() + offset,
          content.size() - offset);
      if (written < 0) {
        if (errno == EINTR) continue;
        throw std::runtime_error(
            "cannot write atomic output: " +
            std::string(std::strerror(errno)));
      }
      offset += static_cast<std::size_t>(written);
    }
    if (::fsync(descriptor) != 0) {
      throw std::runtime_error(
          "cannot sync atomic output: " +
          std::string(std::strerror(errno)));
    }
    if (::close(descriptor) != 0) {
      descriptor = -1;
      throw std::runtime_error(
          "cannot close atomic output: " +
          std::string(std::strerror(errno)));
    }
    descriptor = -1;
    if (::rename(temporary.c_str(), path.c_str()) != 0) {
      throw std::runtime_error(
          "cannot install atomic output: " +
          std::string(std::strerror(errno)));
    }
    renamed = true;
    sync_directory(directory);
  } catch (...) {
    if (descriptor >= 0) ::close(descriptor);
    if (!renamed) {
      std::error_code ignored;
      fs::remove(temporary, ignored);
    }
    throw;
  }
}

std::uint64_t splitmix64(std::uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30U)) *
          UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27U)) *
          UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31U);
}

std::array<std::uint64_t, kCoreEntries> make_zobrist(
    std::uint64_t seed) {
  std::array<std::uint64_t, kCoreEntries> values{};
  std::uint64_t state =
      seed ^ UINT64_C(0x436f726541646a32);
  for (std::uint64_t& value : values) {
    state = splitmix64(state);
    value = state;
  }
  return values;
}

std::uint64_t core_hash(
    const CoreMatrix& core,
    const std::array<std::uint64_t, kCoreEntries>& zobrist) {
  std::uint64_t hash = 0;
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      if (core[row][column] != 0U) {
        hash ^= zobrist[row * kCoreOrder + column];
      }
    }
  }
  return hash;
}

std::int64_t bit_candidate_determinant(
    const State& state, int row, int column) {
  const std::int64_t delta =
      state.core[row][column] == 0U ? 1 : -1;
  return state.determinant +
         delta * state.adjugate[column][row];
}

std::int64_t row_complement_determinant(
    const State& state, int row) {
  Wide result = state.determinant;
  for (int column = 0; column < kCoreOrder; ++column) {
    const std::int64_t delta =
        state.core[row][column] == 0U ? 1 : -1;
    result += static_cast<Wide>(delta) *
              state.adjugate[column][row];
  }
  return checked_narrow(
      result, kDeterminantBound,
      "row-complement determinant");
}

std::int64_t column_complement_determinant(
    const State& state, int column) {
  Wide result = state.determinant;
  for (int row = 0; row < kCoreOrder; ++row) {
    const std::int64_t delta =
        state.core[row][column] == 0U ? 1 : -1;
    result += static_cast<Wide>(delta) *
              state.adjugate[column][row];
  }
  return checked_narrow(
      result, kDeterminantBound,
      "column-complement determinant");
}

std::int64_t whole_complement_determinant(const State& state) {
  Wide result = state.determinant;
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      result -= state.adjugate[row][column];
    }
  }
  return checked_narrow(
      result, kDeterminantBound,
      "whole-complement determinant");
}

Move choose_move(
    const State& state,
    const std::array<std::uint64_t, kMoveCount>& tabu_until,
    std::uint64_t iteration, std::uint64_t best_magnitude,
    std::mt19937_64& randomizer, Statistics& statistics) {
  Move best;
  std::uint64_t best_candidate_magnitude = 0;
  auto consider = [&](Move candidate) {
    ++statistics.candidate_evaluations;
    if (candidate.determinant == 0) {
      ++statistics.singular_candidates;
      return;
    }
    candidate.aspiration =
        magnitude(candidate.determinant) > best_magnitude;
    if (iteration < tabu_until[candidate.id] &&
        !candidate.aspiration) {
      return;
    }
    const std::uint64_t candidate_magnitude =
        magnitude(candidate.determinant);
    const bool better =
        best.id < 0 ||
        candidate_magnitude > best_candidate_magnitude;
    const bool tied =
        best.id >= 0 &&
        candidate_magnitude == best_candidate_magnitude;
    if (better ||
        (tied && (randomizer() & UINT64_C(1)) != 0U)) {
      best = candidate;
      best_candidate_magnitude = candidate_magnitude;
    }
  };

  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      consider(Move{
          MoveKind::kBit, row, column,
          row * kCoreOrder + column,
          bit_candidate_determinant(state, row, column),
          false});
    }
  }
  for (int row = 0; row < kCoreOrder; ++row) {
    consider(Move{
        MoveKind::kRowComplement, row, -1,
        kCoreEntries + row,
        row_complement_determinant(state, row), false});
  }
  for (int column = 0; column < kCoreOrder; ++column) {
    consider(Move{
        MoveKind::kColumnComplement, column, -1,
        kCoreEntries + kCoreOrder + column,
        column_complement_determinant(state, column), false});
  }
  consider(Move{
      MoveKind::kWholeComplement, -1, -1, kMoveCount - 1,
      whole_complement_determinant(state), false});
  if (best.id < 0) {
    throw std::runtime_error("no nonsingular admissible move");
  }
  return best;
}

void apply_rank_one(State& state, const Vector& u,
                    const Vector& v,
                    std::int64_t new_determinant) {
  if (state.determinant == 0 || new_determinant == 0) {
    throw std::runtime_error(
        "rank-one update requires nonsingular endpoints");
  }
  Vector left{};
  Vector right{};
  for (int row = 0; row < kCoreOrder; ++row) {
    Wide value = 0;
    for (int inner = 0; inner < kCoreOrder; ++inner) {
      value += static_cast<Wide>(state.adjugate[row][inner]) *
               u[inner];
    }
    left[row] = checked_narrow(
        value, INT64_MAX, "adjugate-u product");
  }
  for (int column = 0; column < kCoreOrder; ++column) {
    Wide value = 0;
    for (int inner = 0; inner < kCoreOrder; ++inner) {
      value += static_cast<Wide>(v[inner]) *
               state.adjugate[inner][column];
    }
    right[column] = checked_narrow(
        value, INT64_MAX, "v-adjugate product");
  }

  Adjugate updated{};
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      const Wide numerator =
          static_cast<Wide>(new_determinant) *
              state.adjugate[row][column] -
          static_cast<Wide>(left[row]) * right[column];
      if (numerator % state.determinant != 0) {
        throw std::runtime_error(
            "non-exact adjugate rank-one division");
      }
      updated[row][column] = checked_narrow(
          numerator / state.determinant, kCofactorBound,
          "updated cofactor");
    }
  }
  state.adjugate = updated;
  state.determinant = new_determinant;
}

void apply_move(State& state, const Move& move,
                std::uint64_t& hash,
                const std::array<std::uint64_t, kCoreEntries>&
                    zobrist,
                Statistics& statistics) {
  Vector u{};
  Vector v{};
  if (move.kind == MoveKind::kBit) {
    const int row = move.first;
    const int column = move.second;
    u[row] = state.core[row][column] == 0U ? 1 : -1;
    v[column] = 1;
    apply_rank_one(state, u, v, move.determinant);
    state.core[row][column] ^= 1U;
    hash ^= zobrist[row * kCoreOrder + column];
    ++statistics.bit_moves;
    return;
  }
  if (move.kind == MoveKind::kRowComplement) {
    const int row = move.first;
    u[row] = 1;
    for (int column = 0; column < kCoreOrder; ++column) {
      v[column] =
          state.core[row][column] == 0U ? 1 : -1;
    }
    apply_rank_one(state, u, v, move.determinant);
    for (int column = 0; column < kCoreOrder; ++column) {
      state.core[row][column] ^= 1U;
      hash ^= zobrist[row * kCoreOrder + column];
    }
    ++statistics.row_complements;
    return;
  }
  if (move.kind == MoveKind::kColumnComplement) {
    const int column = move.first;
    v[column] = 1;
    for (int row = 0; row < kCoreOrder; ++row) {
      u[row] =
          state.core[row][column] == 0U ? 1 : -1;
    }
    apply_rank_one(state, u, v, move.determinant);
    for (int row = 0; row < kCoreOrder; ++row) {
      state.core[row][column] ^= 1U;
      hash ^= zobrist[row * kCoreOrder + column];
    }
    ++statistics.column_complements;
    return;
  }

  // For even order, det(J-B) = det(B-J).  Apply the rank-one update
  // B-J = B - 11^T, then negate the matrix.  Negating an even-order
  // matrix preserves its determinant and negates its adjugate.
  u.fill(-1);
  v.fill(1);
  apply_rank_one(state, u, v, move.determinant);
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      state.core[row][column] ^= 1U;
      hash ^= zobrist[row * kCoreOrder + column];
      state.adjugate[row][column] =
          -state.adjugate[row][column];
    }
  }
  ++statistics.whole_complements;
}

bool apply_exact_breakout(
    State& state, const CoreMatrix& best_core,
    std::uint64_t& hash,
    const std::array<std::uint64_t, kCoreEntries>& zobrist,
    const Arguments& arguments, std::mt19937_64& randomizer,
    Statistics& statistics) {
  std::array<int, kCoreEntries> entries{};
  std::iota(entries.begin(), entries.end(), 0);
  for (int attempt = 0; attempt < arguments.breakout_attempts;
       ++attempt) {
    ++statistics.breakout_attempts;
    std::shuffle(entries.begin(), entries.end(), randomizer);
    CoreMatrix candidate = best_core;
    for (int index = 0; index < arguments.breakout_flips; ++index) {
      const int flat = entries[static_cast<std::size_t>(index)];
      candidate[flat / kCoreOrder][flat % kCoreOrder] ^= 1U;
    }
    const std::int64_t determinant =
        exact_core_determinant(candidate);
    if (determinant == 0) {
      ++statistics.breakout_singular;
      continue;
    }
    state.core = candidate;
    state.determinant = determinant;
    state.adjugate = exact_adjugate(state.core);
    check_adjugate_identity(state);
    hash = core_hash(state.core, zobrist);
    ++statistics.breakouts;
    return true;
  }
  return false;
}

void log_record(std::ofstream& log, const Arguments& arguments,
                const Statistics& statistics, const char* event,
                double elapsed, const State& state,
                std::uint64_t best_magnitude, int tenure) {
  log << "{\"best_absolute_determinant\":\""
      << wide_to_string(
             static_cast<Wide>(best_magnitude) * kScale)
      << "\",\"best_core_quotient\":" << best_magnitude
      << ",\"bit_moves\":" << statistics.bit_moves
      << ",\"breakout_attempts\":"
      << statistics.breakout_attempts
      << ",\"breakout_flips\":" << arguments.breakout_flips
      << ",\"breakout_interval\":"
      << arguments.breakout_interval
      << ",\"breakout_singular\":"
      << statistics.breakout_singular
      << ",\"breakouts\":" << statistics.breakouts
      << ",\"candidate_evaluations\":"
      << statistics.candidate_evaluations
      << ",\"column_complements\":"
      << statistics.column_complements
      << ",\"core_determinant\":" << state.determinant
      << ",\"cycles\":" << statistics.cycles
      << ",\"determinant_checks\":"
      << statistics.determinant_checks
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(6) << elapsed
      << ",\"engine\":\""
      << (arguments.breakout_interval == 0
              ? "core-adjugate-reactive-tabu-v1"
              : "core-adjugate-reactive-tabu-v1+exact-breakout")
      << "\""
      << ",\"event\":\"" << event << "\""
      << ",\"identity_checks\":" << statistics.identity_checks
      << ",\"iterations\":" << statistics.iterations
      << ",\"max_iterations\":" << arguments.max_iterations
      << ",\"output_path\":\""
      << json_escape(arguments.output.string()) << "\""
      << ",\"promotions\":" << statistics.promotions
      << ",\"row_complements\":" << statistics.row_complements
      << ",\"seconds\":" << arguments.seconds
      << ",\"seed\":" << arguments.seed
      << ",\"singular_candidates\":"
      << statistics.singular_candidates
      << ",\"tenure\":" << tenure
      << ",\"whole_complements\":"
      << statistics.whole_complements << "}\n";
  log.flush();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    if (!arguments.log.parent_path().empty()) {
      fs::create_directories(arguments.log.parent_path());
    }
    std::ofstream log(arguments.log, std::ios::app);
    if (!log) {
      throw std::runtime_error("cannot open research log");
    }

    const SignMatrix input = read_sign_matrix(arguments.start);
    State state;
    state.core = dephase_to_core(input);
    state.determinant = exact_core_determinant(state.core);
    if (state.determinant == 0) {
      throw std::runtime_error("start matrix must be nonsingular");
    }
    state.adjugate = exact_adjugate(state.core);
    check_adjugate_identity(state);

    const Wide input_determinant =
        bareiss([&]() {
          std::vector<std::vector<Wide>> result(
              kSignOrder,
              std::vector<Wide>(kSignOrder));
          for (int row = 0; row < kSignOrder; ++row) {
            for (int column = 0; column < kSignOrder; ++column) {
              result[row][column] = input[row][column];
            }
          }
          return result;
        }());
    if (input_determinant !=
            static_cast<Wide>(state.determinant) * kScale &&
        input_determinant !=
            -static_cast<Wide>(state.determinant) * kScale) {
      throw std::runtime_error(
          "dephased core determinant disagrees with sign matrix");
    }

    CoreMatrix best_core = state.core;
    std::uint64_t best_magnitude = magnitude(state.determinant);
    std::uint64_t checkpoint_nonce = 0;
    atomic_write(
        arguments.output,
        sign_matrix_text(core_to_sign(best_core)),
        checkpoint_nonce++);

    std::mt19937_64 randomizer(arguments.seed);
    const auto zobrist = make_zobrist(arguments.seed);
    std::uint64_t hash = core_hash(state.core, zobrist);
    std::vector<Visit> visits(kVisitTableSize);
    visits[hash & (kVisitTableSize - 1U)] =
        Visit{hash, 0, true};
    std::array<std::uint64_t, kMoveCount> tabu_until{};
    const int baseline_tenure =
        kMinimumTenure +
        static_cast<int>(arguments.seed % UINT64_C(5));
    int tenure = baseline_tenure;
    std::uint64_t last_cycle_iteration = 0;
    std::uint64_t last_progress_iteration = 0;
    Statistics statistics;

    const auto started = Clock::now();
    const auto deadline =
        started + std::chrono::duration<double>(arguments.seconds);
    auto next_heartbeat =
        started +
        std::chrono::duration<double>(
            arguments.heartbeat_seconds);
    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);
    auto promote_current = [&]() {
      const std::uint64_t current_magnitude =
          magnitude(state.determinant);
      if (current_magnitude <= best_magnitude) {
        return false;
      }
      best_magnitude = current_magnitude;
      best_core = state.core;
      last_progress_iteration = statistics.iterations;
      ++statistics.promotions;
      atomic_write(
          arguments.output,
          sign_matrix_text(core_to_sign(best_core)),
          checkpoint_nonce++);
      const double elapsed =
          std::chrono::duration<double>(
              Clock::now() - started)
              .count();
      log_record(
          log, arguments, statistics, "new_best", elapsed,
          state, best_magnitude, tenure);
      std::cout
          << "new best |det|="
          << wide_to_string(
                 static_cast<Wide>(best_magnitude) * kScale)
          << " quotient=" << best_magnitude
          << " iteration=" << statistics.iterations << '\n'
          << std::flush;
      return true;
    };
    log_record(
        log, arguments, statistics, "start", 0.0, state,
        best_magnitude, tenure);

    while (!stop_requested && Clock::now() < deadline &&
           (arguments.max_iterations == 0 ||
            statistics.iterations <
                arguments.max_iterations)) {
      ++statistics.iterations;
      const Move move = choose_move(
          state, tabu_until, statistics.iterations,
          best_magnitude, randomizer, statistics);
      apply_move(
          state, move, hash, zobrist, statistics);

      const int jitter =
          static_cast<int>(randomizer() % UINT64_C(5));
      tabu_until[move.id] =
          statistics.iterations +
          static_cast<std::uint64_t>(tenure + jitter + 1);

      Visit& visit =
          visits[hash & (kVisitTableSize - 1U)];
      if (visit.occupied && visit.hash == hash &&
          statistics.iterations > visit.iteration) {
        const std::uint64_t cycle_length =
            statistics.iterations - visit.iteration;
        if (cycle_length <=
            static_cast<std::uint64_t>(
                4 * kMaximumTenure)) {
          ++statistics.cycles;
          last_cycle_iteration = statistics.iterations;
          tenure = std::min(
              kMaximumTenure,
              tenure + 2 +
                  static_cast<int>(
                      std::min<std::uint64_t>(
                          cycle_length / UINT64_C(8),
                          UINT64_C(8))));
        }
      }
      visit = Visit{hash, statistics.iterations, true};
      if (statistics.iterations - last_cycle_iteration >=
              UINT64_C(512) &&
          (statistics.iterations & UINT64_C(127)) == 0 &&
          tenure > baseline_tenure) {
        --tenure;
      }

      static_cast<void>(promote_current());

      if (arguments.breakout_interval != 0 &&
          statistics.iterations - last_progress_iteration >=
              arguments.breakout_interval) {
        const bool kicked = apply_exact_breakout(
            state, best_core, hash, zobrist, arguments,
            randomizer, statistics);
        if (!kicked) {
          state.core = best_core;
          state.determinant =
              exact_core_determinant(state.core);
          state.adjugate = exact_adjugate(state.core);
          check_adjugate_identity(state);
          hash = core_hash(state.core, zobrist);
        }
        tabu_until.fill(0);
        std::fill(visits.begin(), visits.end(), Visit{});
        visits[hash & (kVisitTableSize - 1U)] =
            Visit{hash, statistics.iterations, true};
        tenure = baseline_tenure;
        last_cycle_iteration = statistics.iterations;
        last_progress_iteration = statistics.iterations;
        const double elapsed =
            std::chrono::duration<double>(
                Clock::now() - started)
                .count();
        log_record(
            log, arguments, statistics,
            kicked ? "breakout" : "breakout_restart",
            elapsed, state, best_magnitude, tenure);
        static_cast<void>(promote_current());
      }

      if (statistics.iterations %
              kIdentityCheckInterval ==
          0) {
        check_adjugate_identity(state);
        ++statistics.identity_checks;
      }
      if (statistics.iterations %
              kDeterminantCheckInterval ==
          0) {
        const std::int64_t checked =
            exact_core_determinant(state.core);
        if (checked != state.determinant) {
          throw std::runtime_error(
              "incremental determinant invariant failed");
        }
        ++statistics.determinant_checks;
      }

      const auto now = Clock::now();
      if (arguments.heartbeat_seconds > 0.0 &&
          now >= next_heartbeat) {
        const double elapsed =
            std::chrono::duration<double>(
                now - started)
                .count();
        log_record(
            log, arguments, statistics, "heartbeat", elapsed,
            state, best_magnitude, tenure);
        next_heartbeat =
            now +
            std::chrono::duration<double>(
                arguments.heartbeat_seconds);
      }
    }

    check_adjugate_identity(state);
    ++statistics.identity_checks;
    if (exact_core_determinant(state.core) !=
        state.determinant) {
      throw std::runtime_error(
          "final determinant invariant failed");
    }
    ++statistics.determinant_checks;
    atomic_write(
        arguments.output,
        sign_matrix_text(core_to_sign(best_core)),
        checkpoint_nonce++);
    const double elapsed =
        std::chrono::duration<double>(
            Clock::now() - started)
            .count();
    log_record(
        log, arguments, statistics,
        stop_requested ? "stopped" : "finished", elapsed,
        state, best_magnitude, tenure);
    std::cout
        << "finished |det|="
        << wide_to_string(
               static_cast<Wide>(best_magnitude) * kScale)
        << " quotient=" << best_magnitude
        << " iterations=" << statistics.iterations
        << " candidate_evaluations="
        << statistics.candidate_evaluations << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "core_adjugate_tabu: " << error.what() << '\n';
    return 2;
  }
}
