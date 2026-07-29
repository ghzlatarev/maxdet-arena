// Exact-pair QUBO trust-region pilot for order-23 MaxDet matrices.
//
// Around a nonsingular center A, choose a support of M entries and evaluate
// the signed determinant quotient q(A) = det(A) / 2^22 exactly modulo
// p = 2^32 - 5.  The order-23 Hadamard bound makes the centered residue of q
// unique.  Exact singleton and pair values define the quadratic interpolation
//
//   q_hat(x) = q0 + sum_i (qi-q0)x_i
//                  + sum_{i<j}(qij-qi-qj+q0)x_i*x_j.
//
// Fixed-cardinality simulated annealing and swap descent optimize q_hat under
// row/column caps.  The strongest proposals are then scored by the true
// rank-k determinant lemma, not by the surrogate.  The trust radius adapts to
// measured model error and the working center periodically restarts.
//
// This is an experimental search tool.  Any apparent frontier improvement
// still requires ./arena verify.

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr int kOrder = 23;
constexpr int kEntries = kOrder * kOrder;
constexpr std::uint64_t kTwo22 = UINT64_C(1) << 22;
constexpr std::uint64_t kPrime = UINT64_C(4294967291);  // 2^32 - 5.

using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using ModMatrix =
    std::array<std::array<std::uint64_t, kOrder>, kOrder>;
using Wide = __int128_t;
using UWide = __uint128_t;
using Clock = std::chrono::steady_clock;

std::atomic<bool> g_stop{false};

void request_stop(int) {
  g_stop.store(true, std::memory_order_relaxed);
}

std::uint64_t add_mod(std::uint64_t left, std::uint64_t right) {
  const std::uint64_t sum = left + right;
  return sum >= kPrime ? sum - kPrime : sum;
}

std::uint64_t sub_mod(std::uint64_t left, std::uint64_t right) {
  return left >= right ? left - right : left + kPrime - right;
}

// Inputs to mul_mod are below p, so the product fits in uint64_t.  Two
// pseudo-Mersenne folds suffice for p = 2^32 - 5.
std::uint64_t reduce_mod(std::uint64_t value) {
  std::uint64_t reduced =
      static_cast<std::uint32_t>(value) + 5 * (value >> 32);
  reduced =
      static_cast<std::uint32_t>(reduced) + 5 * (reduced >> 32);
  if (reduced >= kPrime) {
    reduced -= kPrime;
  }
  return reduced;
}

std::uint64_t mul_mod(std::uint64_t left, std::uint64_t right) {
  return reduce_mod(left * right);
}

std::uint64_t pow_mod(std::uint64_t base, std::uint64_t exponent) {
  std::uint64_t result = 1;
  while (exponent != 0) {
    if ((exponent & 1U) != 0) {
      result = mul_mod(result, base);
    }
    base = mul_mod(base, base);
    exponent >>= 1;
  }
  return result;
}

std::uint64_t signed_mod(int value) {
  if (value >= 0) {
    return static_cast<std::uint64_t>(value);
  }
  return kPrime - static_cast<std::uint64_t>(-value);
}

UWide power_u128(std::uint64_t base, int exponent) {
  UWide result = 1;
  for (int index = 0; index < exponent; ++index) {
    result *= base;
  }
  return result;
}

std::uint64_t integer_sqrt(UWide value) {
  std::uint64_t low = 0;
  std::uint64_t high = UINT64_C(1) << 63;
  while (low + 1 < high) {
    const std::uint64_t middle = low + (high - low) / 2;
    if (static_cast<UWide>(middle) * middle <= value) {
      low = middle;
    } else {
      high = middle;
    }
  }
  return low;
}

std::uint64_t absolute_i64(std::int64_t value) {
  return value < 0 ? static_cast<std::uint64_t>(-value)
                   : static_cast<std::uint64_t>(value);
}

std::string json_escape(std::string_view input) {
  std::ostringstream output;
  for (const unsigned char character : input) {
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
        if (character < 0x20) {
          output << "\\u" << std::hex << std::setw(4)
                 << std::setfill('0') << static_cast<int>(character)
                 << std::dec;
        } else {
          output << static_cast<char>(character);
        }
    }
  }
  return output.str();
}

Matrix read_matrix(const fs::path& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("cannot open matrix: " + path.string());
  }
  Matrix matrix{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      if (!(input >> matrix[row][column]) ||
          (matrix[row][column] != -1 && matrix[row][column] != 1)) {
        throw std::runtime_error(
            "matrix must contain exactly 23x23 entries in {-1,+1}: " +
            path.string());
      }
    }
  }
  std::string extra;
  if (input >> extra) {
    throw std::runtime_error("matrix contains extra data: " + path.string());
  }
  return matrix;
}

std::string matrix_text(const Matrix& matrix) {
  std::ostringstream output;
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      if (column != 0) {
        output << ' ';
      }
      output << matrix[row][column];
    }
    output << '\n';
  }
  return output.str();
}

void atomic_write(const fs::path& destination, const std::string& bytes,
                  std::uint64_t nonce) {
  const fs::path directory = destination.parent_path().empty()
                                 ? fs::path(".")
                                 : destination.parent_path();
  fs::create_directories(directory);
  const fs::path temporary =
      directory /
      ("." + destination.filename().string() + ".qubo-" +
       std::to_string(nonce) + ".tmp");
  {
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) {
      throw std::runtime_error("cannot create temporary output: " +
                               temporary.string());
    }
    output << bytes;
    output.flush();
    if (!output) {
      throw std::runtime_error("failed writing temporary output: " +
                               temporary.string());
    }
  }
  std::error_code error;
  fs::rename(temporary, destination, error);
  if (error) {
    std::error_code ignored;
    fs::remove(temporary, ignored);
    throw std::runtime_error("cannot install output " + destination.string() +
                             ": " + error.message());
  }
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
    if (pivot_row == kOrder) {
      return 0;
    }
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
          throw std::runtime_error("Bareiss division was not exact");
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

std::int64_t bareiss_quotient(const Matrix& matrix) {
  const Wide determinant = exact_determinant(matrix);
  if (determinant % static_cast<Wide>(kTwo22) != 0) {
    throw std::runtime_error(
        "exact determinant violates universal 2^22 divisibility");
  }
  const Wide quotient = determinant / static_cast<Wide>(kTwo22);
  if (quotient < std::numeric_limits<std::int64_t>::min() ||
      quotient > std::numeric_limits<std::int64_t>::max()) {
    throw std::runtime_error("exact quotient does not fit int64");
  }
  return static_cast<std::int64_t>(quotient);
}

std::uint64_t determinant_mod(Matrix matrix) {
  std::array<std::array<std::uint64_t, kOrder>, kOrder> work{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      work[row][column] =
          matrix[row][column] == 1 ? 1 : kPrime - 1;
    }
  }

  std::uint64_t determinant = 1;
  bool negative = false;
  for (int column = 0; column < kOrder; ++column) {
    int pivot = column;
    while (pivot < kOrder && work[pivot][column] == 0) {
      ++pivot;
    }
    if (pivot == kOrder) {
      return 0;
    }
    if (pivot != column) {
      std::swap(work[pivot], work[column]);
      negative = !negative;
    }
    const std::uint64_t pivot_value = work[column][column];
    determinant = mul_mod(determinant, pivot_value);
    const std::uint64_t inverse_pivot =
        pow_mod(pivot_value, kPrime - 2);
    for (int row = column + 1; row < kOrder; ++row) {
      if (work[row][column] == 0) {
        continue;
      }
      const std::uint64_t factor =
          mul_mod(work[row][column], inverse_pivot);
      for (int inner = column; inner < kOrder; ++inner) {
        work[row][inner] =
            sub_mod(work[row][inner],
                    mul_mod(factor, work[column][inner]));
      }
    }
  }
  if (negative && determinant != 0) {
    determinant = kPrime - determinant;
  }
  return determinant;
}

struct ModAnalysis {
  std::uint64_t determinant = 0;
  ModMatrix inverse{};
  bool nonsingular = false;
};

ModAnalysis analyze_matrix(const Matrix& matrix) {
  std::array<std::array<std::uint64_t, 2 * kOrder>, kOrder> work{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      work[row][column] =
          matrix[row][column] == 1 ? 1 : kPrime - 1;
    }
    work[row][kOrder + row] = 1;
  }

  std::uint64_t determinant = 1;
  bool negative = false;
  for (int column = 0; column < kOrder; ++column) {
    int pivot = column;
    while (pivot < kOrder && work[pivot][column] == 0) {
      ++pivot;
    }
    if (pivot == kOrder) {
      return ModAnalysis{};
    }
    if (pivot != column) {
      std::swap(work[pivot], work[column]);
      negative = !negative;
    }
    const std::uint64_t pivot_value = work[column][column];
    determinant = mul_mod(determinant, pivot_value);
    const std::uint64_t inverse_pivot =
        pow_mod(pivot_value, kPrime - 2);
    for (int inner = 0; inner < 2 * kOrder; ++inner) {
      work[column][inner] =
          mul_mod(work[column][inner], inverse_pivot);
    }
    for (int row = 0; row < kOrder; ++row) {
      if (row == column || work[row][column] == 0) {
        continue;
      }
      const std::uint64_t factor = work[row][column];
      for (int inner = 0; inner < 2 * kOrder; ++inner) {
        work[row][inner] =
            sub_mod(work[row][inner],
                    mul_mod(factor, work[column][inner]));
      }
    }
  }
  if (negative && determinant != 0) {
    determinant = kPrime - determinant;
  }

  ModAnalysis result;
  result.determinant = determinant;
  result.nonsingular = true;
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      result.inverse[row][column] = work[row][kOrder + column];
    }
  }
  return result;
}

std::uint64_t determinant_mod_dynamic(
    std::vector<std::vector<std::uint64_t>> matrix) {
  const std::size_t order = matrix.size();
  std::uint64_t determinant = 1;
  bool negative = false;
  for (std::size_t column = 0; column < order; ++column) {
    std::size_t pivot = column;
    while (pivot < order && matrix[pivot][column] == 0) {
      ++pivot;
    }
    if (pivot == order) {
      return 0;
    }
    if (pivot != column) {
      std::swap(matrix[pivot], matrix[column]);
      negative = !negative;
    }
    const std::uint64_t pivot_value = matrix[column][column];
    determinant = mul_mod(determinant, pivot_value);
    const std::uint64_t inverse_pivot =
        pow_mod(pivot_value, kPrime - 2);
    for (std::size_t row = column + 1; row < order; ++row) {
      if (matrix[row][column] == 0) {
        continue;
      }
      const std::uint64_t factor =
          mul_mod(matrix[row][column], inverse_pivot);
      for (std::size_t inner = column; inner < order; ++inner) {
        matrix[row][inner] =
            sub_mod(matrix[row][inner],
                    mul_mod(factor, matrix[column][inner]));
      }
    }
  }
  if (negative && determinant != 0) {
    determinant = kPrime - determinant;
  }
  return determinant;
}

struct ExactContext {
  std::uint64_t inverse_two22 = 0;
  std::uint64_t quotient_bound = 0;

  std::int64_t recover(std::uint64_t determinant_residue) const {
    const std::uint64_t quotient_residue =
        mul_mod(determinant_residue, inverse_two22);
    const std::int64_t quotient =
        quotient_residue <= kPrime / 2
            ? static_cast<std::int64_t>(quotient_residue)
            : static_cast<std::int64_t>(quotient_residue) -
                  static_cast<std::int64_t>(kPrime);
    if (absolute_i64(quotient) > quotient_bound) {
      throw std::runtime_error(
          "centered determinant quotient violates Hadamard bound");
    }
    return quotient;
  }
};

struct Entry {
  int row = 0;
  int column = 0;
  int flat = 0;
};

std::uint64_t flip_delta_mod(int entry_value) {
  return signed_mod(-2 * entry_value);
}

std::int64_t singleton_quotient(const Matrix& matrix,
                                const ModAnalysis& analysis,
                                const Entry& entry,
                                const ExactContext& exact) {
  const std::uint64_t delta =
      flip_delta_mod(matrix[entry.row][entry.column]);
  const std::uint64_t ratio =
      add_mod(1, mul_mod(delta,
                         analysis.inverse[entry.column][entry.row]));
  return exact.recover(mul_mod(analysis.determinant, ratio));
}

std::int64_t pair_quotient(const Matrix& matrix,
                           const ModAnalysis& analysis,
                           const Entry& first,
                           const Entry& second,
                           const ExactContext& exact) {
  const std::uint64_t first_delta =
      flip_delta_mod(matrix[first.row][first.column]);
  const std::uint64_t second_delta =
      flip_delta_mod(matrix[second.row][second.column]);

  const std::uint64_t m00 =
      add_mod(1, mul_mod(first_delta,
                         analysis.inverse[first.column][first.row]));
  const std::uint64_t m01 =
      mul_mod(second_delta,
              analysis.inverse[first.column][second.row]);
  const std::uint64_t m10 =
      mul_mod(first_delta,
              analysis.inverse[second.column][first.row]);
  const std::uint64_t m11 =
      add_mod(1, mul_mod(second_delta,
                         analysis.inverse[second.column][second.row]));
  const std::uint64_t ratio =
      sub_mod(mul_mod(m00, m11), mul_mod(m01, m10));
  return exact.recover(mul_mod(analysis.determinant, ratio));
}

std::int64_t mask_quotient(
    const Matrix& matrix,
    const ModAnalysis& analysis,
    const std::vector<Entry>& support,
    const std::vector<unsigned char>& selected,
    const ExactContext& exact) {
  std::vector<std::size_t> active;
  active.reserve(selected.size());
  for (std::size_t index = 0; index < selected.size(); ++index) {
    if (selected[index] != 0) {
      active.push_back(index);
    }
  }
  if (active.empty()) {
    return exact.recover(analysis.determinant);
  }

  const std::size_t rank = active.size();
  std::vector<std::vector<std::uint64_t>> lemma(
      rank, std::vector<std::uint64_t>(rank, 0));
  for (std::size_t left = 0; left < rank; ++left) {
    const Entry& observed = support[active[left]];
    for (std::size_t right = 0; right < rank; ++right) {
      const Entry& changed = support[active[right]];
      const std::uint64_t delta =
          flip_delta_mod(matrix[changed.row][changed.column]);
      std::uint64_t value =
          mul_mod(delta,
                  analysis.inverse[observed.column][changed.row]);
      if (left == right) {
        value = add_mod(value, 1);
      }
      lemma[left][right] = value;
    }
  }
  const std::uint64_t ratio = determinant_mod_dynamic(std::move(lemma));
  return exact.recover(mul_mod(analysis.determinant, ratio));
}

Matrix apply_mask(const Matrix& center,
                  const std::vector<Entry>& support,
                  const std::vector<unsigned char>& selected) {
  Matrix result = center;
  for (std::size_t index = 0; index < selected.size(); ++index) {
    if (selected[index] != 0) {
      const Entry& entry = support[index];
      result[entry.row][entry.column] *= -1;
    }
  }
  return result;
}

struct Logger {
  std::ofstream stream;

  explicit Logger(const fs::path& path) {
    const fs::path directory =
        path.parent_path().empty() ? fs::path(".") : path.parent_path();
    fs::create_directories(directory);
    if (fs::exists(path)) {
      throw std::runtime_error("refusing to overwrite log: " +
                               path.string());
    }
    stream.open(path, std::ios::out | std::ios::trunc);
    if (!stream) {
      throw std::runtime_error("cannot create log: " + path.string());
    }
  }

  void line(const std::string& text) {
    stream << text << '\n';
    stream.flush();
    if (!stream) {
      throw std::runtime_error("failed writing JSONL log");
    }
  }
};

struct Options {
  fs::path start;
  fs::path output;
  fs::path research_output;
  fs::path log;
  std::uint64_t seed = 31001;
  std::uint64_t frontier = UINT64_C(2779447296000000);
  double seconds = 300.0;
  double heartbeat_seconds = 15.0;
  int support_size = 96;
  int exploit_percent = 67;
  int min_flips = 8;
  int max_flips = 24;
  int initial_max_flips = 12;
  int max_per_row = 5;
  int max_per_column = 5;
  int optimizer_starts = 48;
  int optimizer_steps = 1800;
  int true_pool = 32;
  int rounds_per_restart = 8;
  int differential_samples = 96;
};

std::uint64_t parse_u64(const std::string& value,
                        std::string_view name) {
  std::size_t consumed = 0;
  const std::uint64_t result = std::stoull(value, &consumed);
  if (consumed != value.size()) {
    throw std::runtime_error("invalid integer for " + std::string(name));
  }
  return result;
}

int parse_int(const std::string& value, std::string_view name) {
  std::size_t consumed = 0;
  const long result = std::stol(value, &consumed);
  if (consumed != value.size() ||
      result < std::numeric_limits<int>::min() ||
      result > std::numeric_limits<int>::max()) {
    throw std::runtime_error("invalid integer for " + std::string(name));
  }
  return static_cast<int>(result);
}

double parse_double(const std::string& value, std::string_view name) {
  std::size_t consumed = 0;
  const double result = std::stod(value, &consumed);
  if (consumed != value.size() || !std::isfinite(result)) {
    throw std::runtime_error("invalid number for " + std::string(name));
  }
  return result;
}

void usage(const char* program) {
  std::cerr
      << "usage: " << program
      << " --start MATRIX --output MATRIX --research-output MATRIX"
      << " --log JSONL [options]\n"
      << "  --seconds N --seed N --frontier N --support N\n"
      << "  --min-flips N --max-flips N --initial-max-flips N\n"
      << "  --max-per-row N --max-per-column N\n"
      << "  --optimizer-starts N --optimizer-steps N --true-pool N\n"
      << "  --rounds-per-restart N --differential-samples N\n"
      << "  --exploit-percent N --heartbeat-seconds N\n";
}

Options parse_options(int argc, char** argv) {
  Options options;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    auto value = [&](std::string_view name) -> std::string {
      if (index + 1 >= argc) {
        throw std::runtime_error("missing value after " +
                                 std::string(name));
      }
      return argv[++index];
    };
    if (argument == "--start") {
      options.start = value(argument);
    } else if (argument == "--output") {
      options.output = value(argument);
    } else if (argument == "--research-output") {
      options.research_output = value(argument);
    } else if (argument == "--log") {
      options.log = value(argument);
    } else if (argument == "--seed") {
      options.seed = parse_u64(value(argument), argument);
    } else if (argument == "--frontier") {
      options.frontier = parse_u64(value(argument), argument);
    } else if (argument == "--seconds") {
      options.seconds = parse_double(value(argument), argument);
    } else if (argument == "--heartbeat-seconds") {
      options.heartbeat_seconds =
          parse_double(value(argument), argument);
    } else if (argument == "--support") {
      options.support_size = parse_int(value(argument), argument);
    } else if (argument == "--exploit-percent") {
      options.exploit_percent = parse_int(value(argument), argument);
    } else if (argument == "--min-flips") {
      options.min_flips = parse_int(value(argument), argument);
    } else if (argument == "--max-flips") {
      options.max_flips = parse_int(value(argument), argument);
    } else if (argument == "--initial-max-flips") {
      options.initial_max_flips = parse_int(value(argument), argument);
    } else if (argument == "--max-per-row") {
      options.max_per_row = parse_int(value(argument), argument);
    } else if (argument == "--max-per-column") {
      options.max_per_column = parse_int(value(argument), argument);
    } else if (argument == "--optimizer-starts") {
      options.optimizer_starts = parse_int(value(argument), argument);
    } else if (argument == "--optimizer-steps") {
      options.optimizer_steps = parse_int(value(argument), argument);
    } else if (argument == "--true-pool") {
      options.true_pool = parse_int(value(argument), argument);
    } else if (argument == "--rounds-per-restart") {
      options.rounds_per_restart =
          parse_int(value(argument), argument);
    } else if (argument == "--differential-samples") {
      options.differential_samples =
          parse_int(value(argument), argument);
    } else if (argument == "--help" || argument == "-h") {
      usage(argv[0]);
      std::exit(0);
    } else {
      throw std::runtime_error("unknown argument: " + argument);
    }
  }

  if (options.start.empty() || options.output.empty() ||
      options.research_output.empty() || options.log.empty()) {
    throw std::runtime_error(
        "--start, --output, --research-output, and --log are required");
  }
  if (options.output == options.research_output ||
      options.output == options.log ||
      options.research_output == options.log) {
    throw std::runtime_error("output and log paths must be distinct");
  }
  if (!(options.seconds > 0.0) ||
      !(options.heartbeat_seconds > 0.0)) {
    throw std::runtime_error("time limits must be positive");
  }
  if (options.support_size < options.min_flips ||
      options.support_size > 128) {
    throw std::runtime_error(
        "--support must be between --min-flips and 128");
  }
  if (options.exploit_percent < 0 || options.exploit_percent > 100) {
    throw std::runtime_error("--exploit-percent must be in [0,100]");
  }
  if (options.min_flips < 3 ||
      options.max_flips < options.min_flips ||
      options.max_flips > 24 ||
      options.initial_max_flips < options.min_flips ||
      options.initial_max_flips > options.max_flips) {
    throw std::runtime_error(
        "require 3 <= min-flips <= initial-max-flips <= max-flips <= 24");
  }
  if (options.max_per_row < 1 || options.max_per_column < 1 ||
      options.max_per_row * kOrder < options.max_flips ||
      options.max_per_column * kOrder < options.max_flips) {
    throw std::runtime_error("row/column caps cannot realize max-flips");
  }
  if (options.optimizer_starts < 1 || options.optimizer_steps < 0 ||
      options.true_pool < 1 || options.rounds_per_restart < 1 ||
      options.differential_samples < 0) {
    throw std::runtime_error("search counts must be nonnegative and useful");
  }
  return options;
}

std::vector<int> random_flat_mask(std::mt19937_64& generator,
                                  int cardinality,
                                  int max_per_row,
                                  int max_per_column) {
  std::vector<int> order(kEntries);
  std::iota(order.begin(), order.end(), 0);
  for (int attempt = 0; attempt < 128; ++attempt) {
    std::shuffle(order.begin(), order.end(), generator);
    std::array<int, kOrder> row_counts{};
    std::array<int, kOrder> column_counts{};
    std::vector<int> result;
    result.reserve(static_cast<std::size_t>(cardinality));
    for (const int flat : order) {
      const int row = flat / kOrder;
      const int column = flat % kOrder;
      if (row_counts[row] >= max_per_row ||
          column_counts[column] >= max_per_column) {
        continue;
      }
      result.push_back(flat);
      ++row_counts[row];
      ++column_counts[column];
      if (static_cast<int>(result.size()) == cardinality) {
        return result;
      }
    }
  }
  throw std::runtime_error("cannot construct constrained random mask");
}

Matrix apply_flat_mask(const Matrix& center,
                       const std::vector<int>& mask) {
  Matrix result = center;
  for (const int flat : mask) {
    result[flat / kOrder][flat % kOrder] *= -1;
  }
  return result;
}

void run_differential_tests(const Matrix& base,
                            const ExactContext& exact,
                            int sample_count,
                            std::mt19937_64& generator,
                            Logger& logger) {
  const auto started = Clock::now();
  const ModAnalysis base_analysis = analyze_matrix(base);
  if (!base_analysis.nonsingular) {
    throw std::runtime_error("differential-test base is singular");
  }
  const std::int64_t base_mod_q = exact.recover(base_analysis.determinant);
  const std::int64_t base_bareiss_q = bareiss_quotient(base);
  if (base_mod_q != base_bareiss_q) {
    throw std::runtime_error(
        "modular/Bareiss mismatch on differential-test base");
  }

  for (int sample = 0; sample < sample_count; ++sample) {
    static constexpr std::array<int, 5> sizes{1, 2, 8, 16, 24};
    const int cardinality =
        sizes[static_cast<std::size_t>(sample) % sizes.size()];
    const std::vector<int> flat =
        random_flat_mask(generator, cardinality, 6, 6);
    const Matrix candidate = apply_flat_mask(base, flat);
    const std::int64_t modular_q =
        exact.recover(determinant_mod(candidate));
    const std::int64_t bareiss_q = bareiss_quotient(candidate);
    if (modular_q != bareiss_q) {
      std::ostringstream message;
      message << "modular/Bareiss mismatch in sample " << sample
              << " cardinality " << cardinality << ": modular "
              << modular_q << " Bareiss " << bareiss_q;
      throw std::runtime_error(message.str());
    }
    if (cardinality == 1) {
      const Entry entry{flat[0] / kOrder, flat[0] % kOrder, flat[0]};
      if (singleton_quotient(base, base_analysis, entry, exact) !=
          bareiss_q) {
        throw std::runtime_error(
            "rank-one formula failed differential check");
      }
    } else if (cardinality == 2) {
      const Entry first{flat[0] / kOrder, flat[0] % kOrder, flat[0]};
      const Entry second{flat[1] / kOrder, flat[1] % kOrder, flat[1]};
      if (pair_quotient(base, base_analysis, first, second, exact) !=
          bareiss_q) {
        throw std::runtime_error(
            "rank-two formula failed differential check");
      }
    }
  }

  const double elapsed =
      std::chrono::duration<double>(Clock::now() - started).count();
  std::ostringstream event;
  event << "{\"event\":\"differential_check\",\"passed\":true"
        << ",\"samples\":" << sample_count + 1
        << ",\"bareiss_base_quotient\":" << base_bareiss_q
        << ",\"elapsed_seconds\":" << std::fixed
        << std::setprecision(6) << elapsed << '}';
  logger.line(event.str());
}

struct RankedEntry {
  Entry entry;
  std::int64_t oriented_singleton = 0;
};

std::vector<Entry> select_support(const Matrix& center,
                                  const ModAnalysis& analysis,
                                  const ExactContext& exact,
                                  int support_size,
                                  int exploit_percent,
                                  std::mt19937_64& generator) {
  const std::int64_t center_q = exact.recover(analysis.determinant);
  const int orientation = center_q >= 0 ? 1 : -1;
  std::vector<RankedEntry> ranked;
  ranked.reserve(kEntries);
  for (int flat = 0; flat < kEntries; ++flat) {
    const Entry entry{flat / kOrder, flat % kOrder, flat};
    const std::int64_t q =
        singleton_quotient(center, analysis, entry, exact);
    ranked.push_back(
        RankedEntry{entry, static_cast<std::int64_t>(orientation) * q});
  }
  std::sort(ranked.begin(), ranked.end(),
            [](const RankedEntry& left, const RankedEntry& right) {
              if (left.oriented_singleton != right.oriented_singleton) {
                return left.oriented_singleton >
                       right.oriented_singleton;
              }
              return left.entry.flat < right.entry.flat;
            });

  const int exploit_count =
      support_size * exploit_percent / 100;
  const int soft_line_cap =
      std::max(4, (support_size + kOrder - 1) / kOrder + 2);
  std::array<int, kOrder> row_counts{};
  std::array<int, kOrder> column_counts{};
  std::array<bool, kEntries> used{};
  std::vector<Entry> support;
  support.reserve(static_cast<std::size_t>(support_size));

  auto add = [&](const Entry& entry, bool enforce_cap) {
    if (used[entry.flat]) {
      return false;
    }
    if (enforce_cap &&
        (row_counts[entry.row] >= soft_line_cap ||
         column_counts[entry.column] >= soft_line_cap)) {
      return false;
    }
    used[entry.flat] = true;
    ++row_counts[entry.row];
    ++column_counts[entry.column];
    support.push_back(entry);
    return true;
  };

  for (const RankedEntry& candidate : ranked) {
    if (static_cast<int>(support.size()) >= exploit_count) {
      break;
    }
    static_cast<void>(add(candidate.entry, true));
  }
  for (const RankedEntry& candidate : ranked) {
    if (static_cast<int>(support.size()) >= exploit_count) {
      break;
    }
    static_cast<void>(add(candidate.entry, false));
  }

  std::vector<int> random_order(kEntries);
  std::iota(random_order.begin(), random_order.end(), 0);
  std::shuffle(random_order.begin(), random_order.end(), generator);
  for (const int rank_index : random_order) {
    if (static_cast<int>(support.size()) >= support_size) {
      break;
    }
    static_cast<void>(add(ranked[rank_index].entry, true));
  }
  for (const int rank_index : random_order) {
    if (static_cast<int>(support.size()) >= support_size) {
      break;
    }
    static_cast<void>(add(ranked[rank_index].entry, false));
  }
  if (static_cast<int>(support.size()) != support_size) {
    throw std::runtime_error("failed to select requested support");
  }
  return support;
}

struct QuadraticModel {
  std::vector<Entry> support;
  std::int64_t base = 0;
  std::vector<std::int64_t> singleton;
  std::vector<std::int64_t> linear;
  std::vector<std::vector<std::int64_t>> interaction;
  std::uint64_t exact_pair_evaluations = 0;
};

QuadraticModel build_model(const Matrix& center,
                           const ModAnalysis& analysis,
                           const ExactContext& exact,
                           std::vector<Entry> support) {
  QuadraticModel model;
  model.support = std::move(support);
  const std::int64_t signed_base = exact.recover(analysis.determinant);
  const int orientation = signed_base >= 0 ? 1 : -1;
  model.base = static_cast<std::int64_t>(absolute_i64(signed_base));
  const std::size_t size = model.support.size();
  model.singleton.resize(size);
  model.linear.resize(size);
  model.interaction.assign(
      size, std::vector<std::int64_t>(size, 0));

  for (std::size_t index = 0; index < size; ++index) {
    const std::int64_t signed_q =
        singleton_quotient(center, analysis, model.support[index], exact);
    model.singleton[index] =
        static_cast<std::int64_t>(orientation) * signed_q;
    model.linear[index] = model.singleton[index] - model.base;
  }
  for (std::size_t first = 0; first < size; ++first) {
    for (std::size_t second = first + 1; second < size; ++second) {
      const std::int64_t signed_q =
          pair_quotient(center, analysis, model.support[first],
                        model.support[second], exact);
      const std::int64_t oriented_q =
          static_cast<std::int64_t>(orientation) * signed_q;
      const std::int64_t coefficient =
          oriented_q - model.singleton[first] -
          model.singleton[second] + model.base;
      model.interaction[first][second] = coefficient;
      model.interaction[second][first] = coefficient;
      ++model.exact_pair_evaluations;
    }
  }
  return model;
}

std::int64_t predict(const QuadraticModel& model,
                     const std::vector<unsigned char>& selected) {
  std::int64_t prediction = model.base;
  for (std::size_t first = 0; first < selected.size(); ++first) {
    if (selected[first] == 0) {
      continue;
    }
    prediction += model.linear[first];
    for (std::size_t second = first + 1; second < selected.size();
         ++second) {
      if (selected[second] != 0) {
        prediction += model.interaction[first][second];
      }
    }
  }
  return prediction;
}

std::string mask_key(const std::vector<unsigned char>& selected) {
  std::string key;
  key.reserve(selected.size());
  for (const unsigned char bit : selected) {
    key.push_back(bit == 0 ? '0' : '1');
  }
  return key;
}

struct Proposal {
  std::vector<unsigned char> selected;
  std::int64_t predicted = std::numeric_limits<std::int64_t>::min();
  int cardinality = 0;
};

bool feasible_add_after_remove(
    const QuadraticModel& model,
    const std::array<int, kOrder>& row_counts,
    const std::array<int, kOrder>& column_counts,
    int remove,
    int add,
    int max_per_row,
    int max_per_column) {
  const Entry& removed = model.support[static_cast<std::size_t>(remove)];
  const Entry& added = model.support[static_cast<std::size_t>(add)];
  const int row_after_remove =
      row_counts[added.row] - (removed.row == added.row ? 1 : 0);
  const int column_after_remove =
      column_counts[added.column] -
      (removed.column == added.column ? 1 : 0);
  return row_after_remove < max_per_row &&
         column_after_remove < max_per_column;
}

std::vector<unsigned char> construct_initial(
    const QuadraticModel& model,
    int cardinality,
    int max_per_row,
    int max_per_column,
    bool deterministic,
    std::mt19937_64& generator) {
  const std::size_t size = model.support.size();
  std::vector<unsigned char> selected(size, 0);
  std::vector<std::int64_t> field = model.linear;
  std::array<int, kOrder> row_counts{};
  std::array<int, kOrder> column_counts{};

  for (int step = 0; step < cardinality; ++step) {
    std::vector<std::pair<std::int64_t, int>> choices;
    choices.reserve(size);
    for (std::size_t index = 0; index < size; ++index) {
      if (selected[index] != 0) {
        continue;
      }
      const Entry& entry = model.support[index];
      if (row_counts[entry.row] >= max_per_row ||
          column_counts[entry.column] >= max_per_column) {
        continue;
      }
      choices.emplace_back(field[index], static_cast<int>(index));
    }
    if (choices.empty()) {
      throw std::runtime_error(
          "support cannot realize requested proposal cardinality");
    }
    std::sort(choices.begin(), choices.end(),
              [](const auto& left, const auto& right) {
                if (left.first != right.first) {
                  return left.first > right.first;
                }
                return left.second < right.second;
              });
    std::size_t chosen_rank = 0;
    if (!deterministic) {
      const std::size_t shortlist =
          std::min<std::size_t>(choices.size(), 12);
      std::uniform_int_distribution<std::size_t> choose(0,
                                                        shortlist - 1);
      chosen_rank = choose(generator);
    }
    const int chosen = choices[chosen_rank].second;
    selected[static_cast<std::size_t>(chosen)] = 1;
    const Entry& entry =
        model.support[static_cast<std::size_t>(chosen)];
    ++row_counts[entry.row];
    ++column_counts[entry.column];
    for (std::size_t index = 0; index < size; ++index) {
      field[index] +=
          model.interaction[index][static_cast<std::size_t>(chosen)];
    }
  }
  return selected;
}

void initialize_fields(
    const QuadraticModel& model,
    const std::vector<unsigned char>& selected,
    std::vector<std::int64_t>& field,
    std::array<int, kOrder>& row_counts,
    std::array<int, kOrder>& column_counts) {
  field = model.linear;
  row_counts.fill(0);
  column_counts.fill(0);
  for (std::size_t chosen = 0; chosen < selected.size(); ++chosen) {
    if (selected[chosen] == 0) {
      continue;
    }
    const Entry& entry = model.support[chosen];
    ++row_counts[entry.row];
    ++column_counts[entry.column];
    for (std::size_t index = 0; index < selected.size(); ++index) {
      field[index] += model.interaction[index][chosen];
    }
  }
}

int random_member(const std::vector<unsigned char>& selected,
                  bool want_selected,
                  std::mt19937_64& generator) {
  std::uniform_int_distribution<int> distribution(
      0, static_cast<int>(selected.size()) - 1);
  for (;;) {
    const int candidate = distribution(generator);
    if ((selected[static_cast<std::size_t>(candidate)] != 0) ==
        want_selected) {
      return candidate;
    }
  }
}

void apply_swap(const QuadraticModel& model,
                int remove,
                int add,
                std::vector<unsigned char>& selected,
                std::vector<std::int64_t>& field,
                std::array<int, kOrder>& row_counts,
                std::array<int, kOrder>& column_counts,
                std::int64_t& prediction) {
  const std::size_t remove_index = static_cast<std::size_t>(remove);
  const std::size_t add_index = static_cast<std::size_t>(add);
  const std::int64_t delta =
      -field[remove_index] + field[add_index] -
      model.interaction[remove_index][add_index];
  prediction += delta;
  selected[remove_index] = 0;
  selected[add_index] = 1;
  const Entry& removed = model.support[remove_index];
  const Entry& added = model.support[add_index];
  --row_counts[removed.row];
  --column_counts[removed.column];
  ++row_counts[added.row];
  ++column_counts[added.column];
  for (std::size_t index = 0; index < selected.size(); ++index) {
    field[index] -= model.interaction[index][remove_index];
    field[index] += model.interaction[index][add_index];
  }
}

Proposal optimize_one(const QuadraticModel& model,
                      int cardinality,
                      int steps,
                      int max_per_row,
                      int max_per_column,
                      bool deterministic_start,
                      std::mt19937_64& generator) {
  std::vector<unsigned char> selected =
      construct_initial(model, cardinality, max_per_row, max_per_column,
                        deterministic_start, generator);
  std::vector<std::int64_t> field;
  std::array<int, kOrder> row_counts{};
  std::array<int, kOrder> column_counts{};
  initialize_fields(model, selected, field, row_counts, column_counts);
  std::int64_t prediction = predict(model, selected);
  Proposal best{selected, prediction, cardinality};

  std::vector<std::int64_t> sampled_scales;
  sampled_scales.reserve(64);
  for (int sample = 0; sample < 64; ++sample) {
    const int remove = random_member(selected, true, generator);
    const int add = random_member(selected, false, generator);
    if (!feasible_add_after_remove(model, row_counts, column_counts,
                                   remove, add, max_per_row,
                                   max_per_column)) {
      continue;
    }
    const std::int64_t delta =
        -field[static_cast<std::size_t>(remove)] +
        field[static_cast<std::size_t>(add)] -
        model.interaction[static_cast<std::size_t>(remove)]
                         [static_cast<std::size_t>(add)];
    sampled_scales.push_back(
        static_cast<std::int64_t>(absolute_i64(delta)));
  }
  double initial_temperature = 1.0;
  if (!sampled_scales.empty()) {
    const auto middle =
        sampled_scales.begin() +
        static_cast<std::ptrdiff_t>(sampled_scales.size() / 2);
    std::nth_element(sampled_scales.begin(), middle,
                     sampled_scales.end());
    initial_temperature =
        std::max(1.0, static_cast<double>(*middle));
  }
  std::uniform_real_distribution<double> unit(0.0, 1.0);

  for (int step = 0; step < steps; ++step) {
    const int remove = random_member(selected, true, generator);
    const int add = random_member(selected, false, generator);
    if (!feasible_add_after_remove(model, row_counts, column_counts,
                                   remove, add, max_per_row,
                                   max_per_column)) {
      continue;
    }
    const std::int64_t delta =
        -field[static_cast<std::size_t>(remove)] +
        field[static_cast<std::size_t>(add)] -
        model.interaction[static_cast<std::size_t>(remove)]
                         [static_cast<std::size_t>(add)];
    const double progress =
        steps <= 1 ? 1.0
                   : static_cast<double>(step) /
                         static_cast<double>(steps - 1);
    const double temperature =
        initial_temperature * std::pow(0.002, progress);
    const bool accept =
        delta >= 0 ||
        unit(generator) <
            std::exp(static_cast<double>(delta) /
                     std::max(temperature, 1.0));
    if (!accept) {
      continue;
    }
    apply_swap(model, remove, add, selected, field, row_counts,
               column_counts, prediction);
    if (prediction > best.predicted) {
      best.selected = selected;
      best.predicted = prediction;
    }
  }

  selected = best.selected;
  prediction = best.predicted;
  initialize_fields(model, selected, field, row_counts, column_counts);
  for (int pass = 0; pass < 64; ++pass) {
    std::int64_t best_delta = 0;
    int best_remove = -1;
    int best_add = -1;
    for (std::size_t remove = 0; remove < selected.size(); ++remove) {
      if (selected[remove] == 0) {
        continue;
      }
      for (std::size_t add = 0; add < selected.size(); ++add) {
        if (selected[add] != 0 ||
            !feasible_add_after_remove(
                model, row_counts, column_counts,
                static_cast<int>(remove), static_cast<int>(add),
                max_per_row, max_per_column)) {
          continue;
        }
        const std::int64_t delta =
            -field[remove] + field[add] -
            model.interaction[remove][add];
        if (delta > best_delta) {
          best_delta = delta;
          best_remove = static_cast<int>(remove);
          best_add = static_cast<int>(add);
        }
      }
    }
    if (best_delta <= 0) {
      break;
    }
    apply_swap(model, best_remove, best_add, selected, field,
               row_counts, column_counts, prediction);
  }
  if (prediction != predict(model, selected)) {
    throw std::runtime_error("incremental QUBO score drifted");
  }
  return Proposal{std::move(selected), prediction, cardinality};
}

std::vector<int> cardinalities_for(int minimum, int maximum) {
  std::vector<int> result;
  for (int cardinality = minimum; cardinality <= maximum;
       cardinality += 4) {
    result.push_back(cardinality);
  }
  if (result.empty() || result.back() != maximum) {
    result.push_back(maximum);
  }
  return result;
}

std::vector<Proposal> optimize_model(const QuadraticModel& model,
                                     const Options& options,
                                     int trust_max,
                                     std::mt19937_64& generator) {
  const std::vector<int> cardinalities =
      cardinalities_for(options.min_flips, trust_max);
  std::vector<Proposal> proposals;
  proposals.reserve(static_cast<std::size_t>(options.optimizer_starts));
  std::unordered_set<std::string> seen;
  for (int start = 0; start < options.optimizer_starts; ++start) {
    if (g_stop.load(std::memory_order_relaxed)) {
      break;
    }
    const int cardinality =
        cardinalities[static_cast<std::size_t>(start) %
                      cardinalities.size()];
    Proposal proposal =
        optimize_one(model, cardinality, options.optimizer_steps,
                     options.max_per_row, options.max_per_column,
                     start < static_cast<int>(cardinalities.size()),
                     generator);
    const std::string key = mask_key(proposal.selected);
    if (seen.insert(key).second) {
      proposals.push_back(std::move(proposal));
    }
  }
  std::sort(proposals.begin(), proposals.end(),
            [](const Proposal& left, const Proposal& right) {
              if (left.predicted != right.predicted) {
                return left.predicted > right.predicted;
              }
              return left.cardinality < right.cardinality;
            });
  return proposals;
}

struct Evaluation {
  Proposal proposal;
  std::int64_t quotient = 0;
  std::int64_t oriented = 0;
  std::uint64_t score = 0;
  Matrix matrix{};
};

double median(std::vector<double> values) {
  if (values.empty()) {
    return 0.0;
  }
  const std::size_t middle = values.size() / 2;
  std::nth_element(values.begin(),
                   values.begin() +
                       static_cast<std::ptrdiff_t>(middle),
                   values.end());
  const double upper = values[middle];
  if ((values.size() & 1U) != 0) {
    return upper;
  }
  std::nth_element(values.begin(),
                   values.begin() +
                       static_cast<std::ptrdiff_t>(middle - 1),
                   values.begin() +
                       static_cast<std::ptrdiff_t>(middle));
  return (values[middle - 1] + upper) / 2.0;
}

double correlation(const std::vector<double>& left,
                   const std::vector<double>& right) {
  if (left.size() != right.size() || left.size() < 2) {
    return 0.0;
  }
  const double left_mean =
      std::accumulate(left.begin(), left.end(), 0.0) /
      static_cast<double>(left.size());
  const double right_mean =
      std::accumulate(right.begin(), right.end(), 0.0) /
      static_cast<double>(right.size());
  double covariance = 0.0;
  double left_variance = 0.0;
  double right_variance = 0.0;
  for (std::size_t index = 0; index < left.size(); ++index) {
    const double left_delta = left[index] - left_mean;
    const double right_delta = right[index] - right_mean;
    covariance += left_delta * right_delta;
    left_variance += left_delta * left_delta;
    right_variance += right_delta * right_delta;
  }
  if (left_variance == 0.0 || right_variance == 0.0) {
    return 0.0;
  }
  return covariance / std::sqrt(left_variance * right_variance);
}

struct Accuracy {
  double median_relative_error = 0.0;
  double pearson = 0.0;
  std::int64_t best_predicted = 0;
  std::int64_t best_actual_oriented = 0;
};

Accuracy assess_accuracy(const QuadraticModel& model,
                         const std::vector<Evaluation>& evaluations) {
  Accuracy result;
  if (evaluations.empty()) {
    return result;
  }
  std::vector<double> relative_errors;
  std::vector<double> predicted_deltas;
  std::vector<double> actual_deltas;
  relative_errors.reserve(evaluations.size());
  predicted_deltas.reserve(evaluations.size());
  actual_deltas.reserve(evaluations.size());
  result.best_predicted = evaluations.front().proposal.predicted;
  result.best_actual_oriented = evaluations.front().oriented;
  for (const Evaluation& evaluation : evaluations) {
    const double predicted_delta =
        static_cast<double>(evaluation.proposal.predicted - model.base);
    const double actual_delta =
        static_cast<double>(evaluation.oriented - model.base);
    const double denominator = std::max(1.0, std::fabs(predicted_delta));
    relative_errors.push_back(
        std::fabs(actual_delta - predicted_delta) / denominator);
    predicted_deltas.push_back(predicted_delta);
    actual_deltas.push_back(actual_delta);
    result.best_predicted =
        std::max(result.best_predicted, evaluation.proposal.predicted);
    result.best_actual_oriented =
        std::max(result.best_actual_oriented, evaluation.oriented);
  }
  result.median_relative_error = median(std::move(relative_errors));
  result.pearson = correlation(predicted_deltas, actual_deltas);
  return result;
}

std::int64_t validate_with_bareiss(const Matrix& matrix,
                                   std::int64_t modular_q,
                                   std::string_view reason) {
  const std::int64_t bareiss_q = bareiss_quotient(matrix);
  if (bareiss_q != modular_q) {
    throw std::runtime_error(
        std::string(reason) + " modular/Bareiss mismatch: modular " +
        std::to_string(modular_q) + " Bareiss " +
        std::to_string(bareiss_q));
  }
  return bareiss_q;
}

struct SearchState {
  Matrix matrix{};
  std::int64_t quotient = 0;
};

SearchState make_random_restart(const SearchState& global,
                                const ExactContext& exact,
                                const Options& options,
                                int trust_max,
                                std::mt19937_64& generator) {
  std::uniform_int_distribution<int> cardinality(
      options.min_flips, trust_max);
  for (int attempt = 0; attempt < 128; ++attempt) {
    const std::vector<int> flat =
        random_flat_mask(generator, cardinality(generator),
                         options.max_per_row, options.max_per_column);
    Matrix candidate = apply_flat_mask(global.matrix, flat);
    const std::int64_t quotient =
        exact.recover(determinant_mod(candidate));
    if (quotient != 0) {
      return SearchState{std::move(candidate), quotient};
    }
  }
  throw std::runtime_error("failed to construct nonsingular restart");
}

std::string support_json(const std::vector<Entry>& support) {
  std::ostringstream output;
  output << '[';
  for (std::size_t index = 0; index < support.size(); ++index) {
    if (index != 0) {
      output << ',';
    }
    output << '[' << support[index].row + 1 << ','
           << support[index].column + 1 << ']';
  }
  output << ']';
  return output.str();
}

std::uint64_t support_fingerprint(const std::vector<Entry>& support) {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  for (const Entry& entry : support) {
    hash ^= static_cast<std::uint64_t>(entry.flat + 1);
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

std::string elapsed_json(double value) {
  std::ostringstream output;
  output << std::fixed << std::setprecision(6) << value;
  return output.str();
}

int run(const Options& options) {
  if (fs::exists(options.log)) {
    throw std::runtime_error("refusing to overwrite log: " +
                             options.log.string());
  }
  if (fs::exists(options.output) ||
      fs::exists(options.research_output)) {
    throw std::runtime_error(
        "output paths must be fresh for a reproducible campaign");
  }

  Logger logger(options.log);
  std::signal(SIGINT, request_stop);
  std::signal(SIGTERM, request_stop);

  const UWide hadamard_squared = power_u128(kOrder, kOrder);
  const std::uint64_t hadamard_floor =
      integer_sqrt(hadamard_squared);
  ExactContext exact;
  exact.quotient_bound = hadamard_floor / kTwo22;
  if (2 * exact.quotient_bound >= kPrime) {
    throw std::runtime_error(
        "single-prime exact recovery is ambiguous");
  }
  exact.inverse_two22 = pow_mod(kTwo22 % kPrime, kPrime - 2);

  const Matrix start = read_matrix(options.start);
  const ModAnalysis start_analysis = analyze_matrix(start);
  if (!start_analysis.nonsingular) {
    throw std::runtime_error("start matrix is singular");
  }
  const std::int64_t start_q =
      exact.recover(start_analysis.determinant);
  validate_with_bareiss(start, start_q, "start");
  const std::uint64_t start_score_q = absolute_i64(start_q);

  std::mt19937_64 generator(options.seed);
  std::ostringstream started_event;
  started_event
      << "{\"event\":\"started\",\"algorithm\":"
      << "\"exact-pair-qubo-trust-region-v1\""
      << ",\"start\":\"" << json_escape(options.start.string())
      << "\",\"output\":\"" << json_escape(options.output.string())
      << "\",\"research_output\":\""
      << json_escape(options.research_output.string())
      << "\",\"seed\":" << options.seed
      << ",\"seconds\":" << elapsed_json(options.seconds)
      << ",\"frontier\":\"" << options.frontier
      << "\",\"prime\":" << kPrime
      << ",\"divisor\":" << kTwo22
      << ",\"hadamard_quotient_bound\":" << exact.quotient_bound
      << ",\"start_quotient\":" << start_q
      << ",\"start_score\":\"" << start_score_q * kTwo22
      << "\",\"support\":" << options.support_size
      << ",\"exploit_percent\":" << options.exploit_percent
      << ",\"min_flips\":" << options.min_flips
      << ",\"max_flips\":" << options.max_flips
      << ",\"initial_max_flips\":" << options.initial_max_flips
      << ",\"max_per_row\":" << options.max_per_row
      << ",\"max_per_column\":" << options.max_per_column
      << ",\"optimizer_starts\":" << options.optimizer_starts
      << ",\"optimizer_steps\":" << options.optimizer_steps
      << ",\"true_pool\":" << options.true_pool
      << ",\"rounds_per_restart\":" << options.rounds_per_restart
      << ",\"differential_samples\":"
      << options.differential_samples << '}';
  logger.line(started_event.str());

  run_differential_tests(start, exact, options.differential_samples,
                         generator, logger);

  std::uint64_t checkpoint_nonce = 0;
  atomic_write(options.output, matrix_text(start), ++checkpoint_nonce);
  SearchState global{start, start_q};
  SearchState working = global;
  SearchState research_best{};
  bool have_research_best = false;

  const auto campaign_started = Clock::now();
  const auto deadline =
      campaign_started + std::chrono::duration<double>(options.seconds);
  auto next_heartbeat =
      campaign_started +
      std::chrono::duration<double>(options.heartbeat_seconds);

  std::uint64_t round = 0;
  std::uint64_t restart = 0;
  std::uint64_t models_built = 0;
  std::uint64_t pairs_built = 0;
  std::uint64_t proposals_generated = 0;
  std::uint64_t proposals_evaluated = 0;
  std::uint64_t accepted_moves = 0;
  std::uint64_t exact_promotions = 0;
  int trust_max = options.initial_max_flips;
  int consecutive_rejections = 0;
  std::uniform_real_distribution<double> unit(0.0, 1.0);

  while (!g_stop.load(std::memory_order_relaxed) &&
         Clock::now() < deadline) {
    if (round != 0 &&
        round % static_cast<std::uint64_t>(
                    options.rounds_per_restart) ==
            0) {
      ++restart;
      trust_max = options.initial_max_flips;
      consecutive_rejections = 0;
      std::string restart_kind = "global";
      if (restart % 3 == 1) {
        working = make_random_restart(global, exact, options,
                                      trust_max, generator);
        restart_kind = "random_kick";
      } else if (restart % 3 == 2 && have_research_best) {
        working = research_best;
        restart_kind = "research_best";
      } else {
        working = global;
      }
      std::ostringstream restart_event;
      restart_event
          << "{\"event\":\"restart\",\"restart\":" << restart
          << ",\"round\":" << round << ",\"kind\":\""
          << restart_kind << "\",\"working_score\":\""
          << absolute_i64(working.quotient) * kTwo22
          << "\",\"global_score\":\""
          << absolute_i64(global.quotient) * kTwo22 << "\"}";
      logger.line(restart_event.str());
    }

    const auto model_started = Clock::now();
    const ModAnalysis analysis = analyze_matrix(working.matrix);
    if (!analysis.nonsingular) {
      working = global;
      continue;
    }
    const std::int64_t rebuilt_q = exact.recover(analysis.determinant);
    if (rebuilt_q != working.quotient) {
      throw std::runtime_error("working quotient drifted");
    }
    std::vector<Entry> support =
        select_support(working.matrix, analysis, exact,
                       options.support_size, options.exploit_percent,
                       generator);
    QuadraticModel model =
        build_model(working.matrix, analysis, exact, std::move(support));
    const double model_seconds =
        std::chrono::duration<double>(Clock::now() - model_started)
            .count();
    ++models_built;
    pairs_built += model.exact_pair_evaluations;

    const auto optimizer_started = Clock::now();
    std::vector<Proposal> proposals =
        optimize_model(model, options, trust_max, generator);
    const double optimizer_seconds =
        std::chrono::duration<double>(Clock::now() -
                                      optimizer_started)
            .count();
    proposals_generated += proposals.size();

    if (proposals.empty()) {
      ++round;
      continue;
    }
    const std::size_t evaluate_count =
        std::min<std::size_t>(
            proposals.size(),
            static_cast<std::size_t>(options.true_pool));
    std::vector<Evaluation> evaluations;
    evaluations.reserve(evaluate_count);
    const int orientation = working.quotient >= 0 ? 1 : -1;
    for (std::size_t index = 0; index < evaluate_count; ++index) {
      const Proposal& proposal = proposals[index];
      const std::int64_t quotient =
          mask_quotient(working.matrix, analysis, model.support,
                        proposal.selected, exact);
      Matrix candidate =
          apply_mask(working.matrix, model.support, proposal.selected);
      evaluations.push_back(
          Evaluation{proposal, quotient,
                     static_cast<std::int64_t>(orientation) * quotient,
                     absolute_i64(quotient), std::move(candidate)});
      ++proposals_evaluated;
    }
    std::sort(evaluations.begin(), evaluations.end(),
              [](const Evaluation& left, const Evaluation& right) {
                if (left.score != right.score) {
                  return left.score > right.score;
                }
                return left.oriented > right.oriented;
              });
    const Accuracy accuracy = assess_accuracy(model, evaluations);

    Evaluation& best = evaluations.front();
    const bool distinct_from_start = best.matrix != start;
    if (distinct_from_start &&
        (!have_research_best ||
         best.score > absolute_i64(research_best.quotient))) {
      validate_with_bareiss(best.matrix, best.quotient,
                            "research promotion");
      research_best = SearchState{best.matrix, best.quotient};
      have_research_best = true;
      atomic_write(options.research_output,
                   matrix_text(research_best.matrix),
                   ++checkpoint_nonce);
      ++exact_promotions;
    }

    bool global_promoted = false;
    if (best.score > absolute_i64(global.quotient)) {
      validate_with_bareiss(best.matrix, best.quotient,
                            "global promotion");
      global = SearchState{best.matrix, best.quotient};
      atomic_write(options.output, matrix_text(global.matrix),
                   ++checkpoint_nonce);
      ++exact_promotions;
      global_promoted = true;
    }

    const std::uint64_t working_score =
        absolute_i64(working.quotient);
    bool accept = best.score >= working_score;
    double acceptance_probability = accept ? 1.0 : 0.0;
    if (!accept) {
      std::vector<double> singleton_losses;
      singleton_losses.reserve(model.singleton.size());
      for (const std::int64_t singleton : model.singleton) {
        if (singleton < model.base) {
          singleton_losses.push_back(
              static_cast<double>(model.base - singleton));
        }
      }
      const double loss_scale =
          std::max(1.0, median(std::move(singleton_losses)));
      const double campaign_progress = std::clamp(
          std::chrono::duration<double>(Clock::now() -
                                        campaign_started)
                  .count() /
              options.seconds,
          0.0, 1.0);
      const double temperature =
          loss_scale *
          (2.0 + static_cast<double>(trust_max) / 4.0) *
          (1.0 - 0.85 * campaign_progress);
      const double loss =
          static_cast<double>(working_score - best.score);
      acceptance_probability =
          std::exp(-loss / std::max(temperature, 1.0));
      accept = unit(generator) < acceptance_probability ||
               consecutive_rejections >= 2;
    }
    if (accept && best.quotient != 0) {
      working = SearchState{best.matrix, best.quotient};
      ++accepted_moves;
      consecutive_rejections = 0;
    } else {
      ++consecutive_rejections;
    }

    const int old_trust_max = trust_max;
    if (accuracy.median_relative_error < 0.35 &&
        accuracy.pearson > 0.30) {
      trust_max = std::min(options.max_flips, trust_max + 4);
    } else if (accuracy.median_relative_error > 0.85 ||
               accuracy.pearson < 0.0) {
      trust_max = std::max(options.min_flips, trust_max - 4);
    }

    const double elapsed =
        std::chrono::duration<double>(Clock::now() -
                                      campaign_started)
            .count();
    std::ostringstream event;
    event << "{\"event\":\"round\",\"round\":" << round
          << ",\"restart\":" << restart
          << ",\"elapsed_seconds\":" << elapsed_json(elapsed)
          << ",\"center_score\":\"" << model.base * kTwo22
          << "\",\"support_size\":" << model.support.size()
          << ",\"support_fingerprint\":"
          << support_fingerprint(model.support);
    if (round %
            static_cast<std::uint64_t>(options.rounds_per_restart) ==
        0) {
      event << ",\"support\":" << support_json(model.support);
    }
    event
          << ",\"pair_evaluations\":"
          << model.exact_pair_evaluations
          << ",\"model_seconds\":" << elapsed_json(model_seconds)
          << ",\"optimizer_seconds\":"
          << elapsed_json(optimizer_seconds)
          << ",\"unique_proposals\":" << proposals.size()
          << ",\"true_evaluations\":" << evaluations.size()
          << ",\"best_predicted_quotient\":"
          << accuracy.best_predicted
          << ",\"best_actual_oriented_quotient\":"
          << accuracy.best_actual_oriented
          << ",\"best_actual_quotient\":" << best.quotient
          << ",\"best_actual_score\":\"" << best.score * kTwo22
          << "\",\"model_median_relative_error\":"
          << std::fixed << std::setprecision(6)
          << accuracy.median_relative_error
          << ",\"model_pearson\":" << accuracy.pearson
          << ",\"trust_max_before\":" << old_trust_max
          << ",\"trust_max_after\":" << trust_max
          << ",\"accepted\":" << (accept ? "true" : "false")
          << ",\"acceptance_probability\":"
          << acceptance_probability
          << ",\"global_promoted\":"
          << (global_promoted ? "true" : "false")
          << ",\"global_score\":\""
          << absolute_i64(global.quotient) * kTwo22
          << "\",\"above_frontier_unverified\":"
          << (absolute_i64(global.quotient) * kTwo22 >
                      options.frontier
                  ? "true"
                  : "false")
          << '}';
    logger.line(event.str());

    ++round;
    if (Clock::now() >= next_heartbeat) {
      const double heartbeat_elapsed =
          std::chrono::duration<double>(Clock::now() -
                                        campaign_started)
              .count();
      std::ostringstream heartbeat;
      heartbeat
          << "{\"event\":\"heartbeat\",\"elapsed_seconds\":"
          << elapsed_json(heartbeat_elapsed) << ",\"rounds\":"
          << round << ",\"models_built\":" << models_built
          << ",\"pairs_built\":" << pairs_built
          << ",\"proposals_evaluated\":" << proposals_evaluated
          << ",\"accepted_moves\":" << accepted_moves
          << ",\"global_score\":\""
          << absolute_i64(global.quotient) * kTwo22 << "\"}";
      logger.line(heartbeat.str());
      next_heartbeat =
          Clock::now() +
          std::chrono::duration<double>(options.heartbeat_seconds);
    }
  }

  const double total_elapsed =
      std::chrono::duration<double>(Clock::now() - campaign_started)
          .count();
  const std::uint64_t global_score =
      absolute_i64(global.quotient) * kTwo22;
  std::ostringstream finished;
  finished << "{\"event\":\"finished\",\"stopped_by_signal\":"
           << (g_stop.load(std::memory_order_relaxed) ? "true"
                                                     : "false")
           << ",\"elapsed_seconds\":" << elapsed_json(total_elapsed)
           << ",\"rounds\":" << round
           << ",\"restarts\":" << restart
           << ",\"models_built\":" << models_built
           << ",\"pairs_built\":" << pairs_built
           << ",\"proposals_generated\":" << proposals_generated
           << ",\"proposals_evaluated\":" << proposals_evaluated
           << ",\"accepted_moves\":" << accepted_moves
           << ",\"exact_promotions\":" << exact_promotions
           << ",\"global_quotient\":" << global.quotient
           << ",\"global_determinant\":\""
           << static_cast<std::int64_t>(global.quotient *
                                        static_cast<std::int64_t>(
                                            kTwo22))
           << "\",\"global_score\":\"" << global_score
           << "\",\"frontier\":\"" << options.frontier
           << "\",\"above_frontier_unverified\":"
           << (global_score > options.frontier ? "true" : "false");
  if (have_research_best) {
    finished << ",\"research_best_score\":\""
             << absolute_i64(research_best.quotient) * kTwo22 << "\"";
  }
  finished << '}';
  logger.line(finished.str());

  std::cout << "QUBO trust-region pilot finished\n"
            << "rounds: " << round << "\n"
            << "models: " << models_built << "\n"
            << "exact pair evaluations: " << pairs_built << "\n"
            << "true block evaluations: " << proposals_evaluated
            << "\n"
            << "best exact score: " << global_score << "\n"
            << "frontier: " << options.frontier << "\n"
            << "output: " << options.output << "\n"
            << "log: " << options.log << "\n";
  if (global_score > options.frontier) {
    std::cout << "UNVERIFIED FRONTIER CANDIDATE: run ./arena verify "
              << options.output << "\n";
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    return run(parse_options(argc, argv));
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
