// Exact fixed-radius tabu search over a fiber of neutral MaxDet centers.
//
// A state is (C_i, S), where C_i is one of several exact equal-score
// 23-by-23 sign matrices and S is a fixed-cardinality set of flipped entries.
// Ordinary moves exchange one entry in S for one outside it.  Periodically the
// same S is transplanted to a different center C_j.  This searches transverse
// perturbations of a determinant-preserving network without paying the deep
// entrywise valley between its centers.
//
// Every swap is ranked exactly modulo p = 2^32 - 5.  For an order-23 sign
// matrix, det(A) is divisible by 2^22, and
//
//   floor(sqrt(23^23))/2^22 = 1,089,457,290 < p/2.
//
// Hence the centered residue of det(A)/2^22 uniquely recovers the signed
// integer determinant inside the Hadamard bound.  A rank-two determinant
// lemma evaluates a one-out/one-in exchange in O(1), while a finite-field
// Woodbury update maintains the inverse in O(23^2).  Retained artifacts still
// cross the trusted boundary only after an independent `./arena verify`.

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

namespace {

constexpr int kOrder = 23;
constexpr int kEntries = kOrder * kOrder;
constexpr int kMinimumRadius = 4;
constexpr int kMaximumRadius = 64;
constexpr std::uint64_t kPrime = UINT64_C(4294967291);
constexpr std::uint64_t kTwo22 = UINT64_C(1) << 22U;
constexpr std::uint64_t kHadamardQuotientBound = UINT64_C(1089457290);
constexpr std::size_t kVisitTableSize = 1U << 18U;
constexpr std::string_view kEngine = "neutral-fiber-modular-tabu-v1";

static_assert(2U * kHadamardQuotientBound < kPrime);

using Clock = std::chrono::steady_clock;
using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using ModMatrix =
    std::array<std::array<std::uint64_t, kOrder>, kOrder>;

volatile std::sig_atomic_t stop_requested = 0;

void request_stop(int) { stop_requested = 1; }

std::uint64_t add_mod(std::uint64_t first, std::uint64_t second) {
  const std::uint64_t sum = first + second;
  return sum >= kPrime ? sum - kPrime : sum;
}

std::uint64_t sub_mod(std::uint64_t first, std::uint64_t second) {
  return first >= second ? first - second : first + kPrime - second;
}

std::uint64_t reduce_mod(std::uint64_t value) {
  std::uint64_t reduced =
      static_cast<std::uint32_t>(value) + 5U * (value >> 32U);
  reduced =
      static_cast<std::uint32_t>(reduced) + 5U * (reduced >> 32U);
  if (reduced >= kPrime) reduced -= kPrime;
  return reduced;
}

std::uint64_t mul_mod(std::uint64_t first, std::uint64_t second) {
  return reduce_mod(first * second);
}

std::uint64_t neg_mod(std::uint64_t value) {
  return value == 0 ? 0 : kPrime - value;
}

std::uint64_t pow_mod(std::uint64_t base, std::uint64_t exponent) {
  std::uint64_t result = 1;
  while (exponent != 0) {
    if ((exponent & 1U) != 0) result = mul_mod(result, base);
    base = mul_mod(base, base);
    exponent >>= 1U;
  }
  return result;
}

std::uint64_t inverse_mod(std::uint64_t value) {
  if (value == 0) {
    throw std::runtime_error("attempted to invert zero modulo p");
  }
  return pow_mod(value, kPrime - 2U);
}

struct Factorization {
  std::uint64_t determinant = 0;
  ModMatrix inverse{};
  bool nonsingular = false;
};

Factorization factorize(const Matrix& matrix) {
  std::array<std::array<std::uint64_t, 2 * kOrder>, kOrder> work{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      work[row][column] =
          matrix[row][column] == 1 ? 1U : kPrime - 1U;
      work[row][kOrder + column] = row == column ? 1U : 0U;
    }
  }

  std::uint64_t determinant = 1;
  bool negative = false;
  for (int column = 0; column < kOrder; ++column) {
    int pivot = column;
    while (pivot < kOrder && work[pivot][column] == 0) ++pivot;
    if (pivot == kOrder) return Factorization{};
    if (pivot != column) {
      std::swap(work[pivot], work[column]);
      negative = !negative;
    }
    const std::uint64_t pivot_value = work[column][column];
    determinant = mul_mod(determinant, pivot_value);
    const std::uint64_t pivot_inverse = inverse_mod(pivot_value);
    for (int inner = 0; inner < 2 * kOrder; ++inner) {
      work[column][inner] =
          mul_mod(work[column][inner], pivot_inverse);
    }
    for (int row = 0; row < kOrder; ++row) {
      if (row == column || work[row][column] == 0) continue;
      const std::uint64_t multiplier = work[row][column];
      for (int inner = 0; inner < 2 * kOrder; ++inner) {
        work[row][inner] =
            sub_mod(
                work[row][inner],
                mul_mod(multiplier, work[column][inner]));
      }
    }
  }
  if (negative) determinant = neg_mod(determinant);

  Factorization result;
  result.determinant = determinant;
  result.nonsingular = true;
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      result.inverse[row][column] =
          work[row][kOrder + column];
    }
  }
  return result;
}

struct ExactDeterminant {
  std::int64_t signed_value = 0;
  std::uint64_t absolute_value = 0;
};

ExactDeterminant recover_determinant(std::uint64_t residue) {
  static const std::uint64_t two22_inverse =
      inverse_mod(kTwo22 % kPrime);
  const std::uint64_t quotient_residue =
      mul_mod(residue, two22_inverse);
  const std::int64_t quotient =
      quotient_residue <= kPrime / 2U
          ? static_cast<std::int64_t>(quotient_residue)
          : static_cast<std::int64_t>(quotient_residue) -
                static_cast<std::int64_t>(kPrime);
  const std::uint64_t absolute_quotient =
      quotient < 0
          ? static_cast<std::uint64_t>(-quotient)
          : static_cast<std::uint64_t>(quotient);
  if (absolute_quotient > kHadamardQuotientBound) {
    throw std::runtime_error(
        "modular determinant violates the order-23 Hadamard bound");
  }
  const std::int64_t signed_value =
      quotient * static_cast<std::int64_t>(kTwo22);
  const std::uint64_t absolute_value =
      signed_value < 0
          ? static_cast<std::uint64_t>(-signed_value)
          : static_cast<std::uint64_t>(signed_value);
  return ExactDeterminant{signed_value, absolute_value};
}

std::string read_file_bytes(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("cannot open matrix: " + path.string());
  }
  std::string bytes{
      std::istreambuf_iterator<char>(input),
      std::istreambuf_iterator<char>()};
  if (input.bad()) {
    throw std::runtime_error("cannot read matrix: " + path.string());
  }
  return bytes;
}

Matrix parse_matrix(
    std::string_view bytes, const std::filesystem::path& path) {
  std::istringstream input{std::string(bytes)};
  Matrix matrix{};
  for (auto& row : matrix) {
    for (int& value : row) {
      if (!(input >> value) || (value != -1 && value != 1)) {
        throw std::runtime_error(
            "matrix must contain exactly 23x23 entries in {-1,+1}: " +
            path.string());
      }
    }
  }
  std::string extra;
  if (input >> extra) {
    throw std::runtime_error(
        "matrix contains extra data: " + path.string());
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

template <typename Predicate>
std::string bits_hex(Predicate predicate) {
  constexpr char digits[] = "0123456789abcdef";
  std::string output;
  output.reserve((kEntries + 3) / 4);
  unsigned nibble = 0;
  int used = 0;
  for (int entry = 0; entry < kEntries; ++entry) {
    nibble =
        (nibble << 1U) |
        static_cast<unsigned>(predicate(entry) ? 1U : 0U);
    ++used;
    if (used == 4) {
      output.push_back(digits[nibble]);
      nibble = 0;
      used = 0;
    }
  }
  if (used != 0) {
    nibble <<= static_cast<unsigned>(4 - used);
    output.push_back(digits[nibble]);
  }
  return output;
}

std::string matrix_sign_bits_hex(const Matrix& matrix) {
  return bits_hex([&](int entry) {
    return matrix[entry / kOrder][entry % kOrder] == 1;
  });
}

std::string json_escape(std::string_view input) {
  constexpr char digits[] = "0123456789abcdef";
  std::string output;
  output.reserve(input.size() + 8);
  for (const unsigned char character : input) {
    switch (character) {
      case '"':
        output += "\\\"";
        break;
      case '\\':
        output += "\\\\";
        break;
      case '\b':
        output += "\\b";
        break;
      case '\f':
        output += "\\f";
        break;
      case '\n':
        output += "\\n";
        break;
      case '\r':
        output += "\\r";
        break;
      case '\t':
        output += "\\t";
        break;
      default:
        if (character < 0x20U) {
          output += "\\u00";
          output.push_back(digits[character >> 4U]);
          output.push_back(digits[character & 0x0fU]);
        } else {
          output.push_back(static_cast<char>(character));
        }
    }
  }
  return output;
}

void write_all(int descriptor, std::string_view bytes) {
  std::size_t offset = 0;
  while (offset < bytes.size()) {
    const ssize_t written =
        ::write(descriptor, bytes.data() + offset, bytes.size() - offset);
    if (written < 0) {
      if (errno == EINTR) continue;
      throw std::runtime_error(
          "cannot write atomic output: " +
          std::string(std::strerror(errno)));
    }
    if (written == 0) {
      throw std::runtime_error("short write to atomic output");
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
        "cannot open output directory for sync: " +
        std::string(std::strerror(errno)));
  }
  const int status = ::fsync(descriptor);
  const int saved_errno = errno;
  ::close(descriptor);
  if (status != 0) {
    throw std::runtime_error(
        "cannot sync output directory: " +
        std::string(std::strerror(saved_errno)));
  }
}

void atomic_write(
    const std::filesystem::path& path, std::string_view bytes,
    std::string_view tag, std::uint64_t nonce) {
  const std::filesystem::path directory =
      path.parent_path().empty() ? std::filesystem::path(".")
                                 : path.parent_path();
  std::filesystem::create_directories(directory);
  std::filesystem::path temporary;
  int descriptor = -1;
  for (std::uint64_t attempt = 0; attempt < 64; ++attempt) {
    temporary =
        directory /
        ("." + path.filename().string() + "." + std::string(tag) + "-" +
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
          "cannot create atomic output: " +
          std::string(std::strerror(errno)));
    }
  }
  if (descriptor < 0) {
    throw std::runtime_error("cannot allocate atomic temporary");
  }

  bool renamed = false;
  try {
    write_all(descriptor, bytes);
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
      std::filesystem::remove(temporary, ignored);
    }
    throw;
  }
}

std::uint64_t strict_unsigned(
    std::string_view text, std::string_view option) {
  if (text.empty() ||
      !std::all_of(
          text.begin(), text.end(),
          [](unsigned char character) {
            return character >= '0' && character <= '9';
          })) {
    throw std::runtime_error(
        std::string(option) + " must be a non-negative integer");
  }
  std::size_t consumed = 0;
  const std::uint64_t value =
      std::stoull(std::string(text), &consumed);
  if (consumed != text.size()) {
    throw std::runtime_error(
        std::string(option) + " must be a non-negative integer");
  }
  return value;
}

int strict_integer(std::string_view text, std::string_view option) {
  const std::uint64_t value = strict_unsigned(text, option);
  if (value >
      static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error(std::string(option) + " is too large");
  }
  return static_cast<int>(value);
}

double strict_double(std::string_view text, std::string_view option) {
  std::size_t consumed = 0;
  const double value = std::stod(std::string(text), &consumed);
  if (consumed != text.size() || !std::isfinite(value) || value < 0.0) {
    throw std::runtime_error(
        std::string(option) + " must be finite and non-negative");
  }
  return value;
}

std::vector<int> parse_radii(std::string_view text) {
  if (text.empty()) throw std::runtime_error("--radii must not be empty");
  std::vector<int> radii;
  std::size_t begin = 0;
  while (begin <= text.size()) {
    const std::size_t comma = text.find(',', begin);
    const std::size_t end =
        comma == std::string_view::npos ? text.size() : comma;
    if (end == begin) {
      throw std::runtime_error(
          "--radii must be a comma-separated integer list");
    }
    const int radius =
        strict_integer(text.substr(begin, end - begin), "--radii");
    if (radius < kMinimumRadius || radius > kMaximumRadius) {
      throw std::runtime_error(
          "--radii entries must be between 4 and 64");
    }
    if (std::find(radii.begin(), radii.end(), radius) != radii.end()) {
      throw std::runtime_error("--radii entries must be unique");
    }
    radii.push_back(radius);
    if (comma == std::string_view::npos) break;
    begin = comma + 1;
  }
  return radii;
}

struct Arguments {
  std::filesystem::path center_dir;
  std::filesystem::path output;
  std::filesystem::path tie_output;
  std::filesystem::path checkpoint;
  std::filesystem::path log;
  std::vector<int> radii{8, 12, 16, 24};
  std::uint64_t seed = 23;
  std::uint64_t maximum_iterations = 0;
  std::uint64_t restart_iterations = 2048;
  std::uint64_t rebuild_interval = 4096;
  std::uint64_t macro_period = 64;
  std::size_t macro_pool = 4;
  std::size_t swap_samples = 0;
  int baseline_tenure = 13;
  int maximum_tenure = 192;
  double seconds = 3600.0;
  double heartbeat_seconds = 30.0;
  std::uint64_t score_floor = 0;
  bool score_floor_was_set = false;
};

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
    if (option == "--center-dir") {
      arguments.center_dir = value();
    } else if (option == "--output") {
      arguments.output = value();
    } else if (option == "--tie-output") {
      arguments.tie_output = value();
    } else if (option == "--checkpoint") {
      arguments.checkpoint = value();
    } else if (option == "--log") {
      arguments.log = value();
    } else if (option == "--radii" || option == "--radius") {
      arguments.radii = parse_radii(value());
    } else if (option == "--seed") {
      arguments.seed = strict_unsigned(value(), option);
    } else if (option == "--max-iterations") {
      arguments.maximum_iterations = strict_unsigned(value(), option);
    } else if (option == "--restart-iterations") {
      arguments.restart_iterations = strict_unsigned(value(), option);
    } else if (option == "--rebuild-interval") {
      arguments.rebuild_interval = strict_unsigned(value(), option);
    } else if (option == "--macro-period") {
      arguments.macro_period = strict_unsigned(value(), option);
    } else if (option == "--macro-pool") {
      arguments.macro_pool =
          static_cast<std::size_t>(strict_unsigned(value(), option));
    } else if (option == "--swap-samples") {
      arguments.swap_samples =
          static_cast<std::size_t>(strict_unsigned(value(), option));
    } else if (option == "--tabu-tenure") {
      arguments.baseline_tenure = strict_integer(value(), option);
    } else if (option == "--max-tabu-tenure") {
      arguments.maximum_tenure = strict_integer(value(), option);
    } else if (option == "--seconds") {
      arguments.seconds = strict_double(value(), option);
    } else if (
        option == "--heartbeat" ||
        option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds = strict_double(value(), option);
    } else if (option == "--score-floor") {
      arguments.score_floor = strict_unsigned(value(), option);
      arguments.score_floor_was_set = true;
    } else if (option == "--help") {
      std::cout
          << "usage: neutral_fiber_tabu --center-dir DIR --output MATRIX "
             "--checkpoint JSON --log JSONL [options]\n\n"
          << "  --tie-output MATRIX       first non-center exact floor tie\n"
          << "  --score-floor INTEGER     exact promotion floor\n"
          << "  --radii 8,12,16,24        fixed transverse radii [4,64]\n"
          << "  --macro-period N          swaps between center transplants\n"
          << "  --macro-pool N            random choice among best N centers\n"
          << "  --swap-samples N          sampled swaps; 0 evaluates all\n"
          << "  --restart-iterations N    swaps per restart; 0 disables\n"
          << "  --rebuild-interval N      modular verification interval\n"
          << "  --seconds S               wall-clock bound\n"
          << "  --max-iterations N        deterministic bound; 0 disables\n"
          << "  --seed N                  deterministic seed\n"
          << "  --heartbeat-seconds S     0 disables heartbeats\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (
      arguments.center_dir.empty() || arguments.output.empty() ||
      arguments.checkpoint.empty() || arguments.log.empty()) {
    throw std::runtime_error(
        "--center-dir, --output, --checkpoint, and --log are required");
  }
  if (arguments.seconds == 0.0 && arguments.maximum_iterations == 0) {
    throw std::runtime_error(
        "at least one of --seconds or --max-iterations must be positive");
  }
  if (
      arguments.baseline_tenure <= 0 ||
      arguments.maximum_tenure < arguments.baseline_tenure) {
    throw std::runtime_error(
        "tabu tenures must satisfy 1 <= baseline <= maximum");
  }
  if (arguments.macro_pool == 0) {
    throw std::runtime_error("--macro-pool must be positive");
  }
  return arguments;
}

std::filesystem::path resolved_target(
    const std::filesystem::path& path) {
  const std::filesystem::path parent =
      path.parent_path().empty() ? std::filesystem::path(".")
                                 : path.parent_path();
  return std::filesystem::weakly_canonical(parent) / path.filename();
}

struct Center {
  std::filesystem::path path;
  Matrix matrix{};
  std::string sign_bits;
  std::uint64_t score = 0;
};

std::vector<Center> load_centers(const Arguments& arguments) {
  if (!std::filesystem::is_directory(arguments.center_dir)) {
    throw std::runtime_error("--center-dir is not a directory");
  }
  std::vector<std::filesystem::path> paths;
  for (const auto& entry :
       std::filesystem::directory_iterator(arguments.center_dir)) {
    if (!entry.is_regular_file()) continue;
    const std::string filename = entry.path().filename().string();
    if (
        filename.starts_with("tie-") &&
        filename.ends_with(".matrix.txt")) {
      paths.push_back(entry.path());
    }
  }
  std::sort(paths.begin(), paths.end());
  if (paths.size() < 2) {
    throw std::runtime_error(
        "--center-dir must contain at least two tie-*.matrix.txt files");
  }

  std::vector<Center> centers;
  std::vector<std::string> seen;
  for (const auto& path : paths) {
    const Matrix matrix = parse_matrix(read_file_bytes(path), path);
    const Factorization factorization = factorize(matrix);
    if (!factorization.nonsingular) {
      throw std::runtime_error("center is singular: " + path.string());
    }
    const ExactDeterminant exact =
        recover_determinant(factorization.determinant);
    const std::string sign_bits = matrix_sign_bits_hex(matrix);
    if (std::find(seen.begin(), seen.end(), sign_bits) != seen.end()) {
      throw std::runtime_error(
          "center directory contains duplicate matrices");
    }
    seen.push_back(sign_bits);
    centers.push_back(Center{path, matrix, sign_bits, exact.absolute_value});
  }
  const std::uint64_t score = centers.front().score;
  if (!std::all_of(
          centers.begin(), centers.end(),
          [&](const Center& center) { return center.score == score; })) {
    throw std::runtime_error(
        "all neutral fiber centers must have the same exact score");
  }
  return centers;
}

void ensure_fresh_paths(
    const Arguments& arguments, const std::vector<Center>& centers) {
  std::vector<std::pair<std::string, std::filesystem::path>> targets = {
      {"output", arguments.output},
      {"checkpoint", arguments.checkpoint},
      {"log", arguments.log},
  };
  if (!arguments.tie_output.empty()) {
    targets.emplace_back("tie output", arguments.tie_output);
  }
  std::vector<std::filesystem::path> normalized_targets;
  for (const auto& [label, target] : targets) {
    if (std::filesystem::exists(target)) {
      throw std::runtime_error(
          "refusing to overwrite existing " + label + ": " +
          target.string());
    }
    const auto normalized = resolved_target(target);
    if (
        std::find(
            normalized_targets.begin(), normalized_targets.end(),
            normalized) != normalized_targets.end()) {
      throw std::runtime_error("output paths must be distinct");
    }
    for (const Center& center : centers) {
      if (normalized == std::filesystem::canonical(center.path)) {
        throw std::runtime_error(label + " aliases a center matrix");
      }
    }
    normalized_targets.push_back(normalized);
  }
}

std::uint64_t splitmix64(std::uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value =
      (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
  value =
      (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31U);
}

struct State {
  Matrix matrix{};
  ModMatrix inverse{};
  std::array<unsigned char, kEntries> selected{};
  std::uint64_t determinant_residue = 0;
  std::uint64_t score = 0;
  std::uint64_t subset_hash = 0;
  std::size_t center = 0;
  int radius = 0;
};

struct Move {
  int removed = -1;
  int added = -1;
  std::uint64_t determinant_ratio = 0;
  std::uint64_t determinant_residue = 0;
  std::uint64_t score = 0;
  bool aspiration = false;
};

struct Visit {
  std::uint64_t hash = 0;
  std::uint64_t iteration = 0;
  bool occupied = false;
};

struct Statistics {
  std::uint64_t iterations = 0;
  std::uint64_t accepted_swaps = 0;
  std::uint64_t downhill_swaps = 0;
  std::uint64_t equal_swaps = 0;
  std::uint64_t exact_swap_candidates = 0;
  std::uint64_t singular_swap_candidates = 0;
  std::uint64_t restarts = 0;
  std::uint64_t restart_retries = 0;
  std::uint64_t modular_rebuilds = 0;
  std::uint64_t invariant_checks = 0;
  std::uint64_t macro_events = 0;
  std::uint64_t macro_center_evaluations = 0;
  std::uint64_t macro_downhill = 0;
  std::uint64_t macro_equal = 0;
  std::uint64_t macro_uphill = 0;
  std::uint64_t strict_promotion_candidates = 0;
  std::uint64_t frontier_ties = 0;
  std::uint64_t cycles = 0;
  std::uint64_t tabu_resets = 0;
  std::vector<std::uint64_t> center_visits;
};

std::uint64_t calculate_subset_hash(
    const std::array<unsigned char, kEntries>& selected,
    const std::array<std::uint64_t, kEntries>& zobrist) {
  std::uint64_t result = 0;
  for (int entry = 0; entry < kEntries; ++entry) {
    if (selected[entry] != 0) result ^= zobrist[entry];
  }
  return result;
}

std::uint64_t fiber_hash(
    const State& state,
    const std::vector<std::uint64_t>& center_zobrist) {
  return state.subset_hash ^ center_zobrist[state.center];
}

Matrix matrix_at_center(
    const Matrix& center,
    const std::array<unsigned char, kEntries>& selected) {
  Matrix result = center;
  for (int entry = 0; entry < kEntries; ++entry) {
    if (selected[entry] != 0) {
      result[entry / kOrder][entry % kOrder] *= -1;
    }
  }
  return result;
}

void verify_state(
    const State& state, const std::vector<Center>& centers,
    const std::array<std::uint64_t, kEntries>& zobrist,
    Statistics& statistics) {
  int count = 0;
  for (int entry = 0; entry < kEntries; ++entry) {
    const int row = entry / kOrder;
    const int column = entry % kOrder;
    const bool differs =
        state.matrix[row][column] !=
        centers[state.center].matrix[row][column];
    if (differs != (state.selected[entry] != 0)) {
      throw std::runtime_error(
          "matrix disagrees with center-relative flip set");
    }
    if (differs) ++count;
  }
  if (count != state.radius) {
    throw std::runtime_error("fixed-radius invariant was violated");
  }
  if (
      calculate_subset_hash(state.selected, zobrist) !=
      state.subset_hash) {
    throw std::runtime_error("subset hash invariant was violated");
  }
  const Factorization rebuilt = factorize(state.matrix);
  ++statistics.modular_rebuilds;
  ++statistics.invariant_checks;
  if (
      !rebuilt.nonsingular ||
      rebuilt.determinant != state.determinant_residue ||
      rebuilt.inverse != state.inverse ||
      recover_determinant(rebuilt.determinant).absolute_value !=
          state.score) {
    throw std::runtime_error(
        "incremental modular determinant/inverse check failed");
  }
}

State random_state(
    const std::vector<Center>& centers, int radius,
    const std::array<std::uint64_t, kEntries>& zobrist,
    std::mt19937_64& randomizer, Statistics& statistics) {
  std::array<int, kEntries> entries{};
  std::iota(entries.begin(), entries.end(), 0);
  for (std::uint64_t attempt = 0; attempt < 4096; ++attempt) {
    State state;
    state.center =
        static_cast<std::size_t>(
            randomizer() %
            static_cast<std::uint64_t>(centers.size()));
    state.radius = radius;
    std::shuffle(entries.begin(), entries.end(), randomizer);
    for (int index = 0; index < radius; ++index) {
      state.selected[entries[index]] = 1;
      state.subset_hash ^= zobrist[entries[index]];
    }
    state.matrix =
        matrix_at_center(
            centers[state.center].matrix, state.selected);
    const Factorization factorization = factorize(state.matrix);
    ++statistics.modular_rebuilds;
    if (!factorization.nonsingular) {
      ++statistics.restart_retries;
      continue;
    }
    state.determinant_residue = factorization.determinant;
    state.inverse = factorization.inverse;
    state.score =
        recover_determinant(state.determinant_residue).absolute_value;
    ++statistics.restarts;
    ++statistics.center_visits[state.center];
    return state;
  }
  throw std::runtime_error(
      "could not construct a nonsingular random fiber state");
}

std::uint64_t flip_delta(int matrix_value) {
  return matrix_value == 1 ? kPrime - 2U : 2U;
}

std::uint64_t rank_two_ratio(
    const State& state, int removed, int added) {
  const int first_row = removed / kOrder;
  const int first_column = removed % kOrder;
  const int second_row = added / kOrder;
  const int second_column = added % kOrder;
  const std::uint64_t first_delta =
      flip_delta(state.matrix[first_row][first_column]);
  const std::uint64_t second_delta =
      flip_delta(state.matrix[second_row][second_column]);
  const std::uint64_t k00 =
      add_mod(
          1U,
          mul_mod(
              first_delta,
              state.inverse[first_column][first_row]));
  const std::uint64_t k01 =
      mul_mod(
          second_delta,
          state.inverse[first_column][second_row]);
  const std::uint64_t k10 =
      mul_mod(
          first_delta,
          state.inverse[second_column][first_row]);
  const std::uint64_t k11 =
      add_mod(
          1U,
          mul_mod(
              second_delta,
              state.inverse[second_column][second_row]));
  return sub_mod(mul_mod(k00, k11), mul_mod(k01, k10));
}

Move choose_move(
    const State& state,
    const std::array<std::uint64_t, kEntries>& add_tabu_until,
    std::uint64_t iteration, std::uint64_t best_score,
    std::size_t swap_samples, std::mt19937_64& randomizer,
    Statistics& statistics, bool ignore_tabu,
    const Clock::time_point& deadline, bool& halted) {
  std::vector<int> selected;
  std::vector<int> unselected;
  selected.reserve(static_cast<std::size_t>(state.radius));
  unselected.reserve(
      static_cast<std::size_t>(kEntries - state.radius));
  for (int entry = 0; entry < kEntries; ++entry) {
    (state.selected[entry] != 0 ? selected : unselected)
        .push_back(entry);
  }
  const std::size_t total = selected.size() * unselected.size();
  const std::size_t evaluated =
      swap_samples == 0 ? total : std::min(swap_samples, total);
  std::vector<std::size_t> indices;
  if (evaluated != total) {
    indices.resize(total);
    std::iota(indices.begin(), indices.end(), 0U);
    for (std::size_t index = 0; index < evaluated; ++index) {
      const std::size_t remaining = total - index;
      const std::size_t chosen =
          index +
          static_cast<std::size_t>(
              randomizer() %
              static_cast<std::uint64_t>(remaining));
      std::swap(indices[index], indices[chosen]);
    }
    indices.resize(evaluated);
  }

  Move best;
  std::uint64_t equal_best = 0;
  auto evaluate = [&](std::size_t flat_index) {
    const int removed =
        selected[flat_index / unselected.size()];
    const int added =
        unselected[flat_index % unselected.size()];
    const std::uint64_t ratio =
        rank_two_ratio(state, removed, added);
    ++statistics.exact_swap_candidates;
    if (ratio == 0) {
      ++statistics.singular_swap_candidates;
      return;
    }
    const std::uint64_t determinant =
        mul_mod(state.determinant_residue, ratio);
    const std::uint64_t score =
        recover_determinant(determinant).absolute_value;
    const bool tabu =
        !ignore_tabu && iteration < add_tabu_until[added];
    const bool aspiration = tabu && score > best_score;
    if (tabu && !aspiration) return;
    if (best.removed < 0 || score > best.score) {
      best =
          Move{removed, added, ratio, determinant, score, aspiration};
      equal_best = 1;
    } else if (score == best.score) {
      ++equal_best;
      if (randomizer() % equal_best == 0) {
        best =
            Move{removed, added, ratio, determinant, score, aspiration};
      }
    }
  };

  if (evaluated == total) {
    for (std::size_t index = 0; index < total; ++index) {
      evaluate(index);
      if (
          (index & 8191U) == 0U &&
          (stop_requested || Clock::now() >= deadline)) {
        halted = true;
        break;
      }
    }
  } else {
    for (std::size_t index = 0; index < indices.size(); ++index) {
      evaluate(indices[index]);
      if (
          (index & 8191U) == 0U &&
          (stop_requested || Clock::now() >= deadline)) {
        halted = true;
        break;
      }
    }
  }
  return best;
}

void apply_move(State& state, const Move& move) {
  const int first_row = move.removed / kOrder;
  const int first_column = move.removed % kOrder;
  const int second_row = move.added / kOrder;
  const int second_column = move.added % kOrder;
  const std::uint64_t first_delta =
      flip_delta(state.matrix[first_row][first_column]);
  const std::uint64_t second_delta =
      flip_delta(state.matrix[second_row][second_column]);

  const std::uint64_t k00 =
      add_mod(
          1U,
          mul_mod(
              first_delta,
              state.inverse[first_column][first_row]));
  const std::uint64_t k01 =
      mul_mod(
          second_delta,
          state.inverse[first_column][second_row]);
  const std::uint64_t k10 =
      mul_mod(
          first_delta,
          state.inverse[second_column][first_row]);
  const std::uint64_t k11 =
      add_mod(
          1U,
          mul_mod(
              second_delta,
              state.inverse[second_column][second_row]));
  const std::uint64_t ratio =
      sub_mod(mul_mod(k00, k11), mul_mod(k01, k10));
  if (ratio == 0 || ratio != move.determinant_ratio) {
    throw std::runtime_error("selected modular move ratio changed");
  }
  const std::uint64_t inverse_ratio = inverse_mod(ratio);
  const std::uint64_t inverse_k00 =
      mul_mod(k11, inverse_ratio);
  const std::uint64_t inverse_k01 =
      mul_mod(neg_mod(k01), inverse_ratio);
  const std::uint64_t inverse_k10 =
      mul_mod(neg_mod(k10), inverse_ratio);
  const std::uint64_t inverse_k11 =
      mul_mod(k00, inverse_ratio);

  ModMatrix next_inverse{};
  for (int row = 0; row < kOrder; ++row) {
    const std::uint64_t first_u =
        mul_mod(first_delta, state.inverse[row][first_row]);
    const std::uint64_t second_u =
        mul_mod(second_delta, state.inverse[row][second_row]);
    for (int column = 0; column < kOrder; ++column) {
      const std::uint64_t first_v =
          state.inverse[first_column][column];
      const std::uint64_t second_v =
          state.inverse[second_column][column];
      const std::uint64_t first_inner =
          add_mod(
              mul_mod(inverse_k00, first_v),
              mul_mod(inverse_k01, second_v));
      const std::uint64_t second_inner =
          add_mod(
              mul_mod(inverse_k10, first_v),
              mul_mod(inverse_k11, second_v));
      const std::uint64_t correction =
          add_mod(
              mul_mod(first_u, first_inner),
              mul_mod(second_u, second_inner));
      next_inverse[row][column] =
          sub_mod(state.inverse[row][column], correction);
    }
  }
  state.inverse = next_inverse;
  state.matrix[first_row][first_column] *= -1;
  state.matrix[second_row][second_column] *= -1;
  state.selected[move.removed] = 0;
  state.selected[move.added] = 1;
  state.determinant_residue = move.determinant_residue;
  state.score = move.score;
}

bool transplant_center(
    State& state, const std::vector<Center>& centers,
    std::optional<std::size_t>& previous_center,
    std::size_t macro_pool, std::mt19937_64& randomizer,
    Statistics& statistics) {
  struct Candidate {
    std::size_t center = 0;
    Matrix matrix{};
    Factorization factorization{};
    std::uint64_t score = 0;
  };
  std::vector<Candidate> candidates;
  candidates.reserve(centers.size() - 1U);
  for (std::size_t index = 0; index < centers.size(); ++index) {
    if (index == state.center) continue;
    if (
        previous_center.has_value() &&
        index == *previous_center && centers.size() > 2U) {
      continue;
    }
    Matrix matrix =
        matrix_at_center(centers[index].matrix, state.selected);
    Factorization factorization = factorize(matrix);
    ++statistics.macro_center_evaluations;
    ++statistics.modular_rebuilds;
    if (!factorization.nonsingular) continue;
    const std::uint64_t score =
        recover_determinant(factorization.determinant).absolute_value;
    candidates.push_back(
        Candidate{
            index, std::move(matrix), std::move(factorization), score});
  }
  if (candidates.empty()) return false;
  std::sort(
      candidates.begin(), candidates.end(),
      [](const Candidate& first, const Candidate& second) {
        if (first.score != second.score) {
          return first.score > second.score;
        }
        return first.center < second.center;
      });
  const std::size_t pool = std::min(macro_pool, candidates.size());
  const std::size_t choice =
      static_cast<std::size_t>(
          randomizer() % static_cast<std::uint64_t>(pool));
  Candidate selected = std::move(candidates[choice]);
  const std::uint64_t previous_score = state.score;
  const std::size_t old_center = state.center;
  state.center = selected.center;
  state.matrix = std::move(selected.matrix);
  state.inverse = std::move(selected.factorization.inverse);
  state.determinant_residue = selected.factorization.determinant;
  state.score = selected.score;
  previous_center = old_center;
  ++statistics.macro_events;
  ++statistics.center_visits[state.center];
  if (state.score > previous_score) {
    ++statistics.macro_uphill;
  } else if (state.score < previous_score) {
    ++statistics.macro_downhill;
  } else {
    ++statistics.macro_equal;
  }
  return true;
}

bool matches_center(
    const Matrix& matrix, const std::vector<Center>& centers) {
  return std::any_of(
      centers.begin(), centers.end(),
      [&](const Center& center) { return matrix == center.matrix; });
}

std::string radii_json(const std::vector<int>& radii) {
  std::ostringstream output;
  output << '[';
  for (std::size_t index = 0; index < radii.size(); ++index) {
    if (index != 0) output << ',';
    output << radii[index];
  }
  output << ']';
  return output.str();
}

std::string center_visits_json(
    const std::vector<std::uint64_t>& visits) {
  std::ostringstream output;
  output << '[';
  for (std::size_t index = 0; index < visits.size(); ++index) {
    if (index != 0) output << ',';
    output << visits[index];
  }
  output << ']';
  return output.str();
}

std::string centers_json(const std::vector<Center>& centers) {
  std::ostringstream output;
  output << '[';
  for (std::size_t index = 0; index < centers.size(); ++index) {
    if (index != 0) output << ',';
    output << "{\"index\":" << index
           << ",\"path\":\""
           << json_escape(centers[index].path.string())
           << "\",\"score\":\"" << centers[index].score
           << "\",\"sign_bits_hex\":\""
           << centers[index].sign_bits << "\"}";
  }
  output << ']';
  return output.str();
}

std::string make_record(
    std::string_view event, const Arguments& arguments,
    const Statistics& statistics, const State& state,
    const std::vector<Center>& centers, const Matrix& best_matrix,
    std::uint64_t best_score, int tenure, double elapsed_seconds,
    bool include_centers) {
  std::ostringstream output;
  output
      << "{\"accepted_swaps\":" << statistics.accepted_swaps
      << ",\"best_score\":\"" << best_score
      << "\",\"best_sign_bits_hex\":\""
      << matrix_sign_bits_hex(best_matrix)
      << "\",\"center_count\":" << centers.size()
      << ",\"center_index\":" << state.center
      << ",\"center_path\":\""
      << json_escape(centers[state.center].path.string())
      << "\",\"center_visits\":"
      << center_visits_json(statistics.center_visits)
      << ",\"current_flip_bits_hex\":\""
      << bits_hex(
             [&](int entry) {
               return state.selected[entry] != 0;
             })
      << "\",\"current_score\":\"" << state.score
      << "\",\"current_sign_bits_hex\":\""
      << matrix_sign_bits_hex(state.matrix)
      << "\",\"cycles\":" << statistics.cycles
      << ",\"downhill_swaps\":" << statistics.downhill_swaps
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(6) << elapsed_seconds
      << ",\"engine\":\"" << kEngine
      << "\",\"equal_swaps\":" << statistics.equal_swaps
      << ",\"event\":\"" << json_escape(event)
      << "\",\"exact_recovery\":"
         "\"centered det/2^22 residue is unique under Hadamard bound\""
      << ",\"exact_swap_candidates\":"
      << statistics.exact_swap_candidates
      << ",\"frontier_ties\":" << statistics.frontier_ties
      << ",\"hadamard_quotient_bound\":"
      << kHadamardQuotientBound
      << ",\"invariant_checks\":" << statistics.invariant_checks
      << ",\"iterations\":" << statistics.iterations
      << ",\"macro_center_evaluations\":"
      << statistics.macro_center_evaluations
      << ",\"macro_downhill\":" << statistics.macro_downhill
      << ",\"macro_equal\":" << statistics.macro_equal
      << ",\"macro_events\":" << statistics.macro_events
      << ",\"macro_period\":" << arguments.macro_period
      << ",\"macro_pool\":" << arguments.macro_pool
      << ",\"macro_uphill\":" << statistics.macro_uphill
      << ",\"max_iterations\":" << arguments.maximum_iterations
      << ",\"modular_prime\":" << kPrime
      << ",\"modular_rebuilds\":" << statistics.modular_rebuilds
      << ",\"output_path\":\""
      << json_escape(arguments.output.string())
      << "\",\"radii\":" << radii_json(arguments.radii)
      << ",\"radius\":" << state.radius
      << ",\"restart_iterations\":"
      << arguments.restart_iterations
      << ",\"restart_retries\":" << statistics.restart_retries
      << ",\"restarts\":" << statistics.restarts
      << ",\"score_floor\":\"" << arguments.score_floor
      << "\",\"seconds\":" << std::fixed
      << std::setprecision(6) << arguments.seconds
      << ",\"seed\":" << arguments.seed
      << ",\"singular_swap_candidates\":"
      << statistics.singular_swap_candidates
      << ",\"strict_promotion_candidates\":"
      << statistics.strict_promotion_candidates
      << ",\"swap_samples\":" << arguments.swap_samples
      << ",\"tabu_resets\":" << statistics.tabu_resets
      << ",\"tenure\":" << tenure;
  if (include_centers) {
    output << ",\"centers\":" << centers_json(centers);
  }
  output << '}';
  return output.str();
}

void emit_record(
    std::ofstream& log, std::string_view event,
    const Arguments& arguments, const Statistics& statistics,
    const State& state, const std::vector<Center>& centers,
    const Matrix& best_matrix, std::uint64_t best_score, int tenure,
    const Clock::time_point& started, std::uint64_t& checkpoint_nonce,
    bool include_centers = false) {
  const double elapsed =
      std::chrono::duration<double>(Clock::now() - started).count();
  const std::string record =
      make_record(
          event, arguments, statistics, state, centers, best_matrix,
          best_score, tenure, elapsed, include_centers);
  log << record << '\n';
  log.flush();
  if (!log) {
    throw std::runtime_error("cannot append research log");
  }
  atomic_write(
      arguments.checkpoint, record + "\n", "neutral-fiber-checkpoint",
      checkpoint_nonce++);
}

void consider_candidate(
    const State& state, const Arguments& arguments,
    const std::vector<Center>& centers, Statistics& statistics,
    Matrix& best_matrix, std::uint64_t& best_score, bool& tie_written,
    std::ofstream& log, int tenure, const Clock::time_point& started,
    std::uint64_t& checkpoint_nonce, std::uint64_t& matrix_nonce) {
  if (state.score > best_score) {
    best_score = state.score;
    best_matrix = state.matrix;
    ++statistics.strict_promotion_candidates;
    atomic_write(
        arguments.output, matrix_bytes(best_matrix),
        "neutral-fiber-best", matrix_nonce++);
    emit_record(
        log, "promotion_candidate", arguments, statistics, state,
        centers, best_matrix, best_score, tenure, started,
        checkpoint_nonce);
    std::cout
        << "unverified promotion candidate |det|=" << best_score
        << " center=" << state.center << " radius=" << state.radius
        << " iteration=" << statistics.iterations << '\n'
        << std::flush;
  }
  if (
      !tie_written && !arguments.tie_output.empty() &&
      state.score == arguments.score_floor &&
      !matches_center(state.matrix, centers)) {
    atomic_write(
        arguments.tie_output, matrix_bytes(state.matrix),
        "neutral-fiber-tie", matrix_nonce++);
    tie_written = true;
    ++statistics.frontier_ties;
    emit_record(
        log, "frontier_tie_candidate", arguments, statistics, state,
        centers, best_matrix, best_score, tenure, started,
        checkpoint_nonce);
    std::cout
        << "retained unverified floor-tie candidate center="
        << state.center << " radius=" << state.radius
        << " iteration=" << statistics.iterations << '\n'
        << std::flush;
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    Arguments arguments = parse_arguments(argc, argv);
    const std::vector<Center> centers = load_centers(arguments);
    ensure_fresh_paths(arguments, centers);
    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);

    const std::uint64_t center_score = centers.front().score;
    if (!arguments.score_floor_was_set) {
      arguments.score_floor = center_score;
    }
    if (center_score < arguments.score_floor) {
      throw std::runtime_error(
          "neutral center score is below --score-floor");
    }

    for (const auto& path :
         {arguments.output, arguments.checkpoint, arguments.log}) {
      if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
      }
    }
    if (
        !arguments.tie_output.empty() &&
        !arguments.tie_output.parent_path().empty()) {
      std::filesystem::create_directories(
          arguments.tie_output.parent_path());
    }
    std::ofstream log(arguments.log, std::ios::out);
    if (!log) {
      throw std::runtime_error(
          "cannot create research log: " + arguments.log.string());
    }

    Statistics statistics;
    statistics.center_visits.assign(centers.size(), 0);
    std::mt19937_64 randomizer(arguments.seed);
    std::array<std::uint64_t, kEntries> zobrist{};
    std::uint64_t random_state_value =
        arguments.seed ^ UINT64_C(0x4e65757472616c46);
    for (std::uint64_t& value : zobrist) {
      random_state_value = splitmix64(random_state_value);
      value = random_state_value;
    }
    std::vector<std::uint64_t> center_zobrist(centers.size());
    for (std::uint64_t& value : center_zobrist) {
      random_state_value = splitmix64(random_state_value);
      value = random_state_value;
    }

    Matrix best_matrix = centers.front().matrix;
    std::uint64_t best_score = center_score;
    std::uint64_t matrix_nonce = 0;
    std::uint64_t checkpoint_nonce = 0;
    atomic_write(
        arguments.output, matrix_bytes(best_matrix),
        "neutral-fiber-best", matrix_nonce++);

    std::size_t radius_index = 0;
    State state =
        random_state(
            centers, arguments.radii[radius_index], zobrist,
            randomizer, statistics);
    verify_state(state, centers, zobrist, statistics);
    std::optional<std::size_t> previous_center;
    std::array<std::uint64_t, kEntries> add_tabu_until{};
    std::vector<Visit> visits(kVisitTableSize);
    const std::uint64_t initial_hash =
        fiber_hash(state, center_zobrist);
    visits[initial_hash & (kVisitTableSize - 1U)] =
        Visit{initial_hash, 0, true};
    int tenure = arguments.baseline_tenure;
    std::uint64_t last_cycle_iteration = 0;
    std::uint64_t restart_start_swaps = 0;
    std::uint64_t last_verified_swaps = 0;
    bool tie_written = false;

    const auto started = Clock::now();
    const Clock::time_point deadline =
        arguments.seconds == 0.0
            ? Clock::time_point::max()
            : started +
                  std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(arguments.seconds));
    Clock::time_point next_heartbeat =
        arguments.heartbeat_seconds == 0.0
            ? Clock::time_point::max()
            : started +
                  std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(
                          arguments.heartbeat_seconds));

    emit_record(
        log, "start", arguments, statistics, state, centers,
        best_matrix, best_score, tenure, started, checkpoint_nonce,
        true);
    consider_candidate(
        state, arguments, centers, statistics, best_matrix, best_score,
        tie_written, log, tenure, started, checkpoint_nonce,
        matrix_nonce);
    std::cout
        << "start centers=" << centers.size()
        << " center_score=" << center_score
        << " score_floor=" << arguments.score_floor
        << " radii=" << radii_json(arguments.radii)
        << " exact_prime=" << kPrime << '\n'
        << std::flush;

    while (!stop_requested) {
      if (Clock::now() >= deadline) break;
      if (
          arguments.maximum_iterations != 0 &&
          statistics.iterations >= arguments.maximum_iterations) {
        break;
      }
      if (
          arguments.restart_iterations != 0 &&
          statistics.accepted_swaps - restart_start_swaps >=
              arguments.restart_iterations) {
        radius_index =
            (radius_index + 1U) % arguments.radii.size();
        state =
            random_state(
                centers, arguments.radii[radius_index], zobrist,
                randomizer, statistics);
        verify_state(state, centers, zobrist, statistics);
        restart_start_swaps = statistics.accepted_swaps;
        last_verified_swaps = statistics.accepted_swaps;
        previous_center.reset();
        add_tabu_until.fill(0);
        std::fill(visits.begin(), visits.end(), Visit{});
        const std::uint64_t hash =
            fiber_hash(state, center_zobrist);
        visits[hash & (kVisitTableSize - 1U)] =
            Visit{hash, statistics.iterations, true};
        tenure = arguments.baseline_tenure;
        last_cycle_iteration = statistics.iterations;
        consider_candidate(
            state, arguments, centers, statistics, best_matrix,
            best_score, tie_written, log, tenure, started,
            checkpoint_nonce, matrix_nonce);
        emit_record(
            log, "restart", arguments, statistics, state, centers,
            best_matrix, best_score, tenure, started,
            checkpoint_nonce);
      }

      ++statistics.iterations;
      bool halted = false;
      Move move =
          choose_move(
              state, add_tabu_until, statistics.iterations, best_score,
              arguments.swap_samples, randomizer, statistics, false,
              deadline, halted);
      if (halted) break;
      if (move.removed < 0) {
        add_tabu_until.fill(0);
        ++statistics.tabu_resets;
        move =
            choose_move(
                state, add_tabu_until, statistics.iterations,
                best_score, arguments.swap_samples, randomizer,
                statistics, true, deadline, halted);
      }
      if (halted) break;
      if (move.removed < 0) {
        throw std::runtime_error(
            "no nonsingular exact exchange was available");
      }

      const std::uint64_t previous_score = state.score;
      const int removed = move.removed;
      apply_move(state, move);
      state.subset_hash ^=
          zobrist[move.removed] ^ zobrist[move.added];
      ++statistics.accepted_swaps;
      if (state.score < previous_score) {
        ++statistics.downhill_swaps;
      } else if (state.score == previous_score) {
        ++statistics.equal_swaps;
      }
      const int jitter =
          static_cast<int>(randomizer() % 7U);
      add_tabu_until[removed] =
          statistics.iterations +
          static_cast<std::uint64_t>(tenure + jitter + 1);

      if (
          arguments.macro_period != 0 &&
          statistics.accepted_swaps % arguments.macro_period == 0) {
        transplant_center(
            state, centers, previous_center, arguments.macro_pool,
            randomizer, statistics);
      }

      if (
          arguments.rebuild_interval != 0 &&
          statistics.accepted_swaps - last_verified_swaps >=
              arguments.rebuild_interval) {
        verify_state(state, centers, zobrist, statistics);
        last_verified_swaps = statistics.accepted_swaps;
      }

      const std::uint64_t hash =
          fiber_hash(state, center_zobrist);
      Visit& visit = visits[hash & (kVisitTableSize - 1U)];
      if (
          visit.occupied && visit.hash == hash &&
          statistics.iterations > visit.iteration) {
        const std::uint64_t cycle_length =
            statistics.iterations - visit.iteration;
        if (cycle_length <= 1024U) {
          ++statistics.cycles;
          last_cycle_iteration = statistics.iterations;
          tenure =
              std::min(
                  arguments.maximum_tenure,
                  tenure + 2 +
                      static_cast<int>(
                          std::min<std::uint64_t>(
                              cycle_length / 16U, 12U)));
        }
      }
      visit = Visit{hash, statistics.iterations, true};
      if (
          statistics.iterations - last_cycle_iteration >= 2048U &&
          (statistics.iterations & 127U) == 0U &&
          tenure > arguments.baseline_tenure) {
        --tenure;
      }

      consider_candidate(
          state, arguments, centers, statistics, best_matrix,
          best_score, tie_written, log, tenure, started,
          checkpoint_nonce, matrix_nonce);

      const auto now = Clock::now();
      if (now >= next_heartbeat) {
        emit_record(
            log, "heartbeat", arguments, statistics, state, centers,
            best_matrix, best_score, tenure, started,
            checkpoint_nonce);
        do {
          next_heartbeat +=
              std::chrono::duration_cast<Clock::duration>(
                  std::chrono::duration<double>(
                      arguments.heartbeat_seconds));
        } while (next_heartbeat <= now);
      }
    }

    verify_state(state, centers, zobrist, statistics);
    emit_record(
        log, stop_requested ? "stopped" : "finished", arguments,
        statistics, state, centers, best_matrix, best_score, tenure,
        started, checkpoint_nonce);
    std::cout
        << (stop_requested ? "stopped" : "finished")
        << " iterations=" << statistics.iterations
        << " exact_swap_candidates="
        << statistics.exact_swap_candidates
        << " accepted_swaps=" << statistics.accepted_swaps
        << " macro_events=" << statistics.macro_events
        << " best=" << best_score << '\n'
        << std::flush;
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
