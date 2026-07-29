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
#include <limits>
#include <queue>
#include <random>
#include <sstream>
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
constexpr std::uint64_t kProbeInterval = 1U << 20U;
constexpr std::string_view kScreeningSemantics =
    "all_combinations_floating_score_exact_margin_gate_"
    "no_rounding_certificate";

using Clock = std::chrono::steady_clock;
using Wide = __int128_t;
using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using Inverse = std::array<std::array<double, kOrder>, kOrder>;

volatile std::sig_atomic_t stop_requested = 0;

void request_stop(int) { stop_requested = 1; }

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

std::string json_escape(std::string_view input) {
  std::string result;
  result.reserve(input.size() + 8);
  for (const unsigned char character : input) {
    switch (character) {
      case '"':
        result += "\\\"";
        break;
      case '\\':
        result += "\\\\";
        break;
      case '\b':
        result += "\\b";
        break;
      case '\f':
        result += "\\f";
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
        if (character < 0x20U) {
          constexpr char hexadecimal[] = "0123456789abcdef";
          result += "\\u00";
          result.push_back(hexadecimal[character >> 4U]);
          result.push_back(hexadecimal[character & 0x0fU]);
        } else {
          result.push_back(static_cast<char>(character));
        }
    }
  }
  return result;
}

Wide parse_wide(const std::string& text, std::string_view option) {
  if (text.empty()) {
    throw std::runtime_error("empty value for " + std::string(option));
  }
  using UnsignedWide = __uint128_t;
  constexpr UnsignedWide maximum =
      static_cast<UnsignedWide>(std::numeric_limits<Wide>::max());
  UnsignedWide value = 0;
  for (const char character : text) {
    if (character < '0' || character > '9') {
      throw std::runtime_error("invalid value for " + std::string(option));
    }
    const unsigned digit = static_cast<unsigned>(character - '0');
    if (value > (maximum - digit) / 10U) {
      throw std::runtime_error("value out of range for " +
                               std::string(option));
    }
    value = value * 10U + digit;
  }
  return static_cast<Wide>(value);
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

std::string matrix_sign_bits_hex(const Matrix& matrix) {
  constexpr char digits[] = "0123456789abcdef";
  std::string result;
  result.reserve((kEntries + 3) / 4);
  unsigned nibble = 0;
  int used = 0;
  for (const auto& row : matrix) {
    for (const int value : row) {
      nibble = (nibble << 1U) | static_cast<unsigned>(value == 1);
      ++used;
      if (used == 4) {
        result.push_back(digits[nibble]);
        nibble = 0;
        used = 0;
      }
    }
  }
  if (used != 0) {
    nibble <<= static_cast<unsigned>(4 - used);
    result.push_back(digits[nibble]);
  }
  return result;
}

void write_all(int descriptor, std::string_view bytes) {
  std::size_t offset = 0;
  while (offset < bytes.size()) {
    const ssize_t written =
        ::write(descriptor, bytes.data() + offset, bytes.size() - offset);
    if (written < 0) {
      if (errno == EINTR) continue;
      throw std::runtime_error("write failed: " +
                               std::string(std::strerror(errno)));
    }
    if (written == 0) throw std::runtime_error("short write");
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
    throw std::runtime_error("cannot open output directory for sync: " +
                             std::string(std::strerror(errno)));
  }
  const int result = ::fsync(descriptor);
  const int saved_errno = errno;
  ::close(descriptor);
  if (result != 0) {
    throw std::runtime_error("cannot sync output directory: " +
                             std::string(std::strerror(saved_errno)));
  }
}

void atomic_write(const std::filesystem::path& path, std::string_view bytes,
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
    descriptor = ::open(temporary.c_str(), flags, 0644);
    if (descriptor >= 0) break;
    if (errno != EEXIST) {
      throw std::runtime_error("cannot create atomic output: " +
                               std::string(std::strerror(errno)));
    }
  }
  if (descriptor < 0) {
    throw std::runtime_error("cannot allocate atomic output temporary");
  }

  bool renamed = false;
  try {
    write_all(descriptor, bytes);
    if (::fsync(descriptor) != 0) {
      throw std::runtime_error("cannot sync atomic output: " +
                               std::string(std::strerror(errno)));
    }
    if (::close(descriptor) != 0) {
      descriptor = -1;
      throw std::runtime_error("cannot close atomic output: " +
                               std::string(std::strerror(errno)));
    }
    descriptor = -1;
    if (::rename(temporary.c_str(), path.c_str()) != 0) {
      throw std::runtime_error("cannot install atomic output: " +
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

void create_log(const std::filesystem::path& path) {
  const std::filesystem::path directory =
      path.parent_path().empty() ? std::filesystem::path(".")
                                 : path.parent_path();
  std::filesystem::create_directories(directory);
  int flags = O_WRONLY | O_CREAT | O_EXCL;
#ifdef O_CLOEXEC
  flags |= O_CLOEXEC;
#endif
  const int descriptor = ::open(path.c_str(), flags, 0644);
  if (descriptor < 0) {
    throw std::runtime_error(
        errno == EEXIST
            ? "refusing to append to an existing log: " + path.string()
            : "cannot create log: " + std::string(std::strerror(errno)));
  }
  if (::fsync(descriptor) != 0) {
    const int saved_errno = errno;
    ::close(descriptor);
    throw std::runtime_error("cannot sync log: " +
                             std::string(std::strerror(saved_errno)));
  }
  ::close(descriptor);
  sync_directory(directory);
}

void append_log(const std::filesystem::path& path, std::string_view record) {
  int flags = O_WRONLY | O_APPEND;
#ifdef O_CLOEXEC
  flags |= O_CLOEXEC;
#endif
  const int descriptor = ::open(path.c_str(), flags);
  if (descriptor < 0) {
    throw std::runtime_error("cannot open log: " +
                             std::string(std::strerror(errno)));
  }
  write_all(descriptor, record);
  write_all(descriptor, "\n");
  if (::fsync(descriptor) != 0) {
    const int saved_errno = errno;
    ::close(descriptor);
    throw std::runtime_error("cannot sync log: " +
                             std::string(std::strerror(saved_errno)));
  }
  if (::close(descriptor) != 0) {
    throw std::runtime_error("cannot close log: " +
                             std::string(std::strerror(errno)));
  }
}

struct Arguments {
  std::filesystem::path start;
  std::filesystem::path output;
  std::filesystem::path research_output;
  std::filesystem::path log;
  std::filesystem::path snapshot;
  std::uint64_t seed = 23;
  std::uint64_t shard_count = 1;
  std::uint64_t shard_index = 0;
  std::size_t top_pool = 256;
  std::size_t calibration_samples = 256;
  Wide exact_margin = 1000000000000LL;
  double seconds = 3600.0;
  double heartbeat_seconds = 30.0;
};

std::uint64_t parse_unsigned(const std::string& text,
                             std::string_view option) {
  std::size_t consumed = 0;
  const unsigned long long value = std::stoull(text, &consumed);
  if (consumed != text.size()) {
    throw std::runtime_error("invalid value for " + std::string(option));
  }
  return static_cast<std::uint64_t>(value);
}

double parse_double(const std::string& text, std::string_view option) {
  std::size_t consumed = 0;
  const double value = std::stod(text, &consumed);
  if (consumed != text.size() || !std::isfinite(value)) {
    throw std::runtime_error("invalid value for " + std::string(option));
  }
  return value;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    auto require_value = [&]() -> std::string {
      if (++index >= argc) {
        throw std::runtime_error("missing value for " + option);
      }
      return argv[index];
    };
    if (option == "--start") {
      arguments.start = require_value();
    } else if (option == "--output") {
      arguments.output = require_value();
    } else if (option == "--research-output") {
      arguments.research_output = require_value();
    } else if (option == "--log") {
      arguments.log = require_value();
    } else if (option == "--snapshot") {
      arguments.snapshot = require_value();
    } else if (option == "--seed") {
      arguments.seed = parse_unsigned(require_value(), option);
    } else if (option == "--shard-count") {
      arguments.shard_count = parse_unsigned(require_value(), option);
    } else if (option == "--shard-index") {
      arguments.shard_index = parse_unsigned(require_value(), option);
    } else if (option == "--top-pool") {
      arguments.top_pool =
          static_cast<std::size_t>(parse_unsigned(require_value(), option));
    } else if (option == "--calibration-samples") {
      arguments.calibration_samples =
          static_cast<std::size_t>(parse_unsigned(require_value(), option));
    } else if (option == "--exact-margin") {
      arguments.exact_margin = parse_wide(require_value(), option);
    } else if (option == "--seconds") {
      arguments.seconds = parse_double(require_value(), option);
    } else if (option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds = parse_double(require_value(), option);
    } else if (option == "--help") {
      std::cout
          << "usage: radius4_screen --start MATRIX --output MATRIX"
          << " --research-output MATRIX --log JSONL [--snapshot JSON]"
          << " [--seed N] [--shard-count N] [--shard-index I]"
          << " [--top-pool 256] [--calibration-samples 256]"
          << " [--exact-margin 1000000000000]"
          << " [--seconds 3600] [--heartbeat-seconds 30]\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }

  if (arguments.start.empty() || arguments.output.empty() ||
      arguments.research_output.empty() || arguments.log.empty()) {
    throw std::runtime_error(
        "--start, --output, --research-output, and --log are required");
  }
  if (arguments.snapshot.empty()) {
    arguments.snapshot = arguments.log.string() + ".snapshot.json";
  }
  if (arguments.shard_count == 0 ||
      arguments.shard_index >= arguments.shard_count) {
    throw std::runtime_error(
        "--shard-index must be smaller than positive --shard-count");
  }
  if (arguments.top_pool == 0 || arguments.calibration_samples == 0) {
    throw std::runtime_error(
        "--top-pool and --calibration-samples must be positive");
  }
  if (arguments.exact_margin <= 0 || !(arguments.seconds > 0.0) ||
      !(arguments.heartbeat_seconds > 0.0)) {
    throw std::runtime_error(
        "margin and time parameters must all be positive");
  }

  auto normalized = [](const std::filesystem::path& path) {
    std::error_code error;
    const std::filesystem::path result =
        std::filesystem::weakly_canonical(
            std::filesystem::absolute(path), error);
    if (!error) return result;
    return std::filesystem::absolute(path).lexically_normal();
  };
  std::array<std::filesystem::path, 5> paths{
      normalized(arguments.start),
      normalized(arguments.output),
      normalized(arguments.research_output),
      normalized(arguments.log),
      normalized(arguments.snapshot)};
  std::sort(paths.begin(), paths.end());
  if (std::adjacent_find(paths.begin(), paths.end()) != paths.end()) {
    throw std::runtime_error("input and output paths must all be distinct");
  }
  return arguments;
}

Inverse exact_cofactor_inverse(const Matrix& matrix, Wide determinant) {
  Inverse inverse{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      Matrix basis = matrix;
      basis[row].fill(0);
      basis[row][column] = 1;
      const Wide cofactor = exact_determinant(basis);
      inverse[column][row] =
          static_cast<double>(cofactor) / static_cast<double>(determinant);
    }
  }
  return inverse;
}

std::vector<double> update_table(const Matrix& matrix,
                                 const Inverse& inverse) {
  std::vector<double> table(
      static_cast<std::size_t>(kEntries) * kEntries);
  for (int left = 0; left < kEntries; ++left) {
    const int left_column = left % kOrder;
    for (int right = 0; right < kEntries; ++right) {
      const int right_row = right / kOrder;
      const int right_column = right % kOrder;
      const double delta =
          -2.0 * static_cast<double>(matrix[right_row][right_column]);
      table[static_cast<std::size_t>(left) * kEntries + right] =
          delta * inverse[left_column][right_row];
    }
  }
  return table;
}

inline double two_by_two(double a, double b, double c, double d) {
  return a * d - b * c;
}

inline double four_flip_ratio(
    const std::vector<double>& table,
    const std::array<std::uint16_t, 4>& entries) {
  double matrix[4][4];
  for (int row = 0; row < 4; ++row) {
    for (int column = 0; column < 4; ++column) {
      matrix[row][column] =
          table[static_cast<std::size_t>(entries[row]) * kEntries +
                entries[column]] +
          (row == column ? 1.0 : 0.0);
    }
  }

  const double first01 =
      two_by_two(matrix[0][0], matrix[0][1],
                 matrix[1][0], matrix[1][1]);
  const double first02 =
      two_by_two(matrix[0][0], matrix[0][2],
                 matrix[1][0], matrix[1][2]);
  const double first03 =
      two_by_two(matrix[0][0], matrix[0][3],
                 matrix[1][0], matrix[1][3]);
  const double first12 =
      two_by_two(matrix[0][1], matrix[0][2],
                 matrix[1][1], matrix[1][2]);
  const double first13 =
      two_by_two(matrix[0][1], matrix[0][3],
                 matrix[1][1], matrix[1][3]);
  const double first23 =
      two_by_two(matrix[0][2], matrix[0][3],
                 matrix[1][2], matrix[1][3]);

  const double second23 =
      two_by_two(matrix[2][2], matrix[2][3],
                 matrix[3][2], matrix[3][3]);
  const double second13 =
      two_by_two(matrix[2][1], matrix[2][3],
                 matrix[3][1], matrix[3][3]);
  const double second12 =
      two_by_two(matrix[2][1], matrix[2][2],
                 matrix[3][1], matrix[3][2]);
  const double second03 =
      two_by_two(matrix[2][0], matrix[2][3],
                 matrix[3][0], matrix[3][3]);
  const double second02 =
      two_by_two(matrix[2][0], matrix[2][2],
                 matrix[3][0], matrix[3][2]);
  const double second01 =
      two_by_two(matrix[2][0], matrix[2][1],
                 matrix[3][0], matrix[3][1]);

  return first01 * second23 - first02 * second13 +
         first03 * second12 + first12 * second03 -
         first13 * second02 + first23 * second01;
}

Matrix materialize(const Matrix& baseline,
                   const std::array<std::uint16_t, 4>& entries) {
  Matrix result = baseline;
  for (const std::uint16_t entry : entries) {
    result[entry / kOrder][entry % kOrder] *= -1;
  }
  return result;
}

std::uint64_t choose_three(std::uint64_t count) {
  if (count < 3) return 0;
  return count * (count - 1) * (count - 2) / 6;
}

std::uint64_t expected_in_shard(std::uint64_t shard_count,
                                std::uint64_t shard_index) {
  std::uint64_t result = 0;
  for (int first = 0; first < kEntries - 3; ++first) {
    if (static_cast<std::uint64_t>(first) % shard_count == shard_index) {
      result += choose_three(
          static_cast<std::uint64_t>(kEntries - first - 1));
    }
  }
  return result;
}

struct ApproximateCandidate {
  double predicted_score = 0.0;
  std::array<std::uint16_t, 4> entries{};
};

struct ApproximateMinHeap {
  bool operator()(const ApproximateCandidate& left,
                  const ApproximateCandidate& right) const {
    return left.predicted_score > right.predicted_score;
  }
};

struct Statistics {
  std::uint64_t screened = 0;
  std::uint64_t near_gate_candidates = 0;
  std::uint64_t near_gate_exact_checks = 0;
  std::uint64_t final_pool_exact_checks = 0;
  std::uint64_t promotions = 0;
  double approximate_best_score = 0.0;
  double calibration_maximum_absolute_error = 0.0;
  Wide research_best_score = 0;
};

struct Runner {
  const Arguments& arguments;
  const Matrix& baseline;
  const std::vector<double>& table;
  Wide baseline_determinant;
  Wide baseline_score;
  std::uint64_t expected;
  Statistics statistics;
  Matrix best_matrix;
  Wide best_score;
  Matrix research_best_matrix;
  Clock::time_point started;
  Clock::time_point deadline;
  Clock::time_point next_heartbeat;
  std::uint64_t snapshot_nonce = 0;

  double elapsed_seconds() const {
    return std::chrono::duration<double>(Clock::now() - started).count();
  }

  bool should_stop() const {
    return stop_requested || Clock::now() >= deadline;
  }

  std::string record(std::string_view event, bool complete,
                     std::string_view termination) const {
    std::ostringstream output;
    output << "{\"event\":\"" << event << "\""
           << ",\"screening_semantics\":\"" << kScreeningSemantics << "\""
           << ",\"start\":{\"path\":\""
           << json_escape(arguments.start.string())
           << "\",\"row_major_sign_bits_hex\":\""
           << matrix_sign_bits_hex(baseline) << "\"}"
           << ",\"complete\":" << (complete ? "true" : "false")
           << ",\"termination\":\"" << termination << "\""
           << ",\"seed\":" << arguments.seed
           << ",\"shard_count\":" << arguments.shard_count
           << ",\"shard_index\":" << arguments.shard_index
           << ",\"expected_combinations\":" << expected
           << ",\"screened_combinations\":" << statistics.screened
           << ",\"top_pool\":" << arguments.top_pool
           << ",\"calibration_samples\":"
           << arguments.calibration_samples
           << ",\"calibration_maximum_absolute_error\":"
           << std::fixed << std::setprecision(6)
           << statistics.calibration_maximum_absolute_error
           << ",\"exact_margin\":\""
           << wide_to_string(arguments.exact_margin) << "\""
           << ",\"near_gate_candidates\":"
           << statistics.near_gate_candidates
           << ",\"near_gate_exact_checks\":"
           << statistics.near_gate_exact_checks
           << ",\"final_pool_exact_checks\":"
           << statistics.final_pool_exact_checks
           << ",\"promotions\":" << statistics.promotions
           << ",\"approximate_best_score\":" << std::fixed
           << std::setprecision(3) << statistics.approximate_best_score
           << ",\"baseline_score\":\"" << wide_to_string(baseline_score)
           << "\",\"best_score\":\"" << wide_to_string(best_score)
           << "\",\"research_best_score\":\""
           << wide_to_string(statistics.research_best_score) << "\""
           << ",\"elapsed_seconds\":" << std::fixed << std::setprecision(6)
           << elapsed_seconds() << "}";
    return output.str();
  }

  void persist(std::string_view event, bool complete,
               std::string_view termination) {
    const std::string json = record(event, complete, termination);
    append_log(arguments.log, json);
    atomic_write(arguments.snapshot, json + "\n", "radius4-snapshot",
                 ++snapshot_nonce);
  }

  void exact_gate(const std::array<std::uint16_t, 4>& entries) {
    const Matrix candidate = materialize(baseline, entries);
    const Wide score = absolute(exact_determinant(candidate));
    ++statistics.near_gate_exact_checks;
    if (score > statistics.research_best_score) {
      statistics.research_best_score = score;
      research_best_matrix = candidate;
    }
    if (score <= best_score) return;
    best_score = score;
    best_matrix = candidate;
    ++statistics.promotions;
    atomic_write(arguments.output, matrix_bytes(best_matrix),
                 "radius4-promotion", statistics.near_gate_exact_checks);
    persist("promotion", false, "running");
    std::cout << "promotion |det|=" << wide_to_string(score)
              << " shard=" << arguments.shard_index << '\n'
              << std::flush;
  }
};

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);

    for (const auto& path :
         {arguments.output, arguments.research_output,
          arguments.log, arguments.snapshot}) {
      std::error_code error;
      const auto status = std::filesystem::symlink_status(path, error);
      const bool missing =
          error == std::errc::no_such_file_or_directory;
      if (missing) {
        error.clear();
      }
      if (error) {
        throw std::runtime_error(
            "cannot inspect output path before starting: " + path.string());
      }
      if (!missing &&
          status.type() != std::filesystem::file_type::not_found) {
        throw std::runtime_error(
            "refusing to overwrite existing output: " + path.string());
      }
    }

    const Matrix baseline = read_matrix(arguments.start);
    const Wide baseline_determinant = exact_determinant(baseline);
    if (baseline_determinant == 0) {
      throw std::runtime_error("start matrix must be nonsingular");
    }
    const Wide baseline_score = absolute(baseline_determinant);
    if (arguments.exact_margin >= baseline_score) {
      throw std::runtime_error(
          "--exact-margin must be smaller than the start score");
    }
    const Inverse inverse =
        exact_cofactor_inverse(baseline, baseline_determinant);
    const std::vector<double> table = update_table(baseline, inverse);
    const std::uint64_t expected =
        expected_in_shard(arguments.shard_count, arguments.shard_index);

    create_log(arguments.log);
    atomic_write(arguments.output, matrix_bytes(baseline), "radius4-best", 0);

    const auto started = Clock::now();
    Runner runner{
        arguments,
        baseline,
        table,
        baseline_determinant,
        baseline_score,
        expected,
        {},
        baseline,
        baseline_score,
        baseline,
        started,
        started + std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(arguments.seconds)),
        started + std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(
                          arguments.heartbeat_seconds)),
        0};

    std::mt19937_64 randomizer(arguments.seed);
    std::array<std::uint16_t, kEntries> permutation{};
    for (int entry = 0; entry < kEntries; ++entry) {
      permutation[entry] = static_cast<std::uint16_t>(entry);
    }
    for (std::size_t sample = 0; sample < arguments.calibration_samples;
         ++sample) {
      std::shuffle(permutation.begin(), permutation.end(), randomizer);
      std::array<std::uint16_t, 4> entries{
          permutation[0], permutation[1], permutation[2], permutation[3]};
      std::sort(entries.begin(), entries.end());
      const double predicted =
          static_cast<double>(baseline_determinant) *
          four_flip_ratio(table, entries);
      const Wide exact =
          exact_determinant(materialize(baseline, entries));
      const double error =
          std::fabs(predicted - static_cast<double>(exact));
      runner.statistics.calibration_maximum_absolute_error =
          std::max(runner.statistics.calibration_maximum_absolute_error,
                   error);
    }
    if (runner.statistics.calibration_maximum_absolute_error >
        static_cast<double>(arguments.exact_margin) * 0.01) {
      throw std::runtime_error(
          "calibration error exceeds one percent of the exact gate margin");
    }

    runner.persist("start", false, "running");
    std::cout << "start shard=" << arguments.shard_index << '/'
              << arguments.shard_count << " combinations=" << expected
              << " exact_margin=" << wide_to_string(arguments.exact_margin)
              << " calibration_error="
              << runner.statistics.calibration_maximum_absolute_error << '\n'
              << std::flush;

    std::priority_queue<
        ApproximateCandidate, std::vector<ApproximateCandidate>,
        ApproximateMinHeap>
        top;
    const double exact_gate =
        static_cast<double>(baseline_score - arguments.exact_margin);
    bool complete = true;

    for (int first = 0; first < kEntries - 3 && complete; ++first) {
      if (static_cast<std::uint64_t>(first) % arguments.shard_count !=
          arguments.shard_index) {
        continue;
      }
      for (int second = first + 1; second < kEntries - 2 && complete;
           ++second) {
        for (int third = second + 1; third < kEntries - 1 && complete;
             ++third) {
          for (int fourth = third + 1; fourth < kEntries; ++fourth) {
            const std::array<std::uint16_t, 4> entries{
                static_cast<std::uint16_t>(first),
                static_cast<std::uint16_t>(second),
                static_cast<std::uint16_t>(third),
                static_cast<std::uint16_t>(fourth)};
            const double ratio = four_flip_ratio(table, entries);
            const double predicted_score =
                std::fabs(static_cast<double>(baseline_determinant) * ratio);
            ++runner.statistics.screened;
            runner.statistics.approximate_best_score =
                std::max(runner.statistics.approximate_best_score,
                         predicted_score);

            ApproximateCandidate candidate{predicted_score, entries};
            if (top.size() < arguments.top_pool) {
              top.push(candidate);
            } else if (predicted_score > top.top().predicted_score) {
              top.pop();
              top.push(candidate);
            }

            if (predicted_score >= exact_gate) {
              ++runner.statistics.near_gate_candidates;
              runner.exact_gate(entries);
            }

            if ((runner.statistics.screened & (kProbeInterval - 1U)) == 0U) {
              const auto now = Clock::now();
              if (stop_requested || now >= runner.deadline) {
                complete = false;
                break;
              }
              if (now >= runner.next_heartbeat) {
                runner.persist("heartbeat", false, "running");
                runner.next_heartbeat =
                    now + std::chrono::duration_cast<Clock::duration>(
                              std::chrono::duration<double>(
                                  arguments.heartbeat_seconds));
              }
            }
          }
        }
      }
    }
    if (runner.statistics.screened != expected) complete = false;

    std::vector<ApproximateCandidate> finalists;
    finalists.reserve(top.size());
    while (!top.empty()) {
      finalists.push_back(top.top());
      top.pop();
    }
    std::sort(
        finalists.begin(), finalists.end(),
        [](const ApproximateCandidate& left,
           const ApproximateCandidate& right) {
          return left.predicted_score > right.predicted_score;
        });
    for (const ApproximateCandidate& finalist : finalists) {
      const Matrix candidate = materialize(baseline, finalist.entries);
      const Wide score = absolute(exact_determinant(candidate));
      ++runner.statistics.final_pool_exact_checks;
      if (score > runner.statistics.research_best_score) {
        runner.statistics.research_best_score = score;
        runner.research_best_matrix = candidate;
      }
      if (score > runner.best_score) {
        runner.best_score = score;
        runner.best_matrix = candidate;
        ++runner.statistics.promotions;
        atomic_write(arguments.output, matrix_bytes(runner.best_matrix),
                     "radius4-final-promotion",
                     runner.statistics.final_pool_exact_checks);
      }
    }
    if (runner.statistics.research_best_score != 0) {
      atomic_write(arguments.research_output,
                   matrix_bytes(runner.research_best_matrix),
                   "radius4-research-best",
                   runner.statistics.final_pool_exact_checks);
    }

    const std::string termination =
        complete ? "shard_completed"
                 : (stop_requested ? "signal" : "time_limit");
    runner.persist("finished", complete, termination);
    std::cout << "finished complete=" << (complete ? "true" : "false")
              << " screened=" << runner.statistics.screened << '/'
              << expected
              << " exact_checks="
              << (runner.statistics.near_gate_exact_checks +
                  runner.statistics.final_pool_exact_checks)
              << " best=" << wide_to_string(runner.best_score)
              << " research_best="
              << wide_to_string(runner.statistics.research_best_score)
              << '\n';
    return complete ? 0 : 3;
  } catch (const std::exception& error) {
    std::cerr << "radius4_screen: " << error.what() << '\n';
    return 2;
  }
}
