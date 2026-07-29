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
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

namespace {

constexpr int kOrder = 23;
constexpr int kEntries = kOrder * kOrder;
constexpr int kMaskWords = (kEntries + 63) / 64;
constexpr long double kSingularTolerance = 1.0e-24L;

using Clock = std::chrono::steady_clock;
using Wide = __int128_t;
using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using Inverse = std::array<std::array<long double, kOrder>, kOrder>;
using Mask = std::array<std::uint64_t, kMaskWords>;

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
    throw std::runtime_error("cannot open elite matrix: " + path.string());
  }
  Matrix matrix{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      if (!(input >> matrix[row][column]) ||
          (matrix[row][column] != -1 && matrix[row][column] != 1)) {
        throw std::runtime_error(
            "elite matrix must contain exactly 23x23 entries in {-1,+1}");
      }
    }
  }
  std::string extra;
  if (input >> extra) {
    throw std::runtime_error("elite matrix contains extra data: " +
                             path.string());
  }
  return matrix;
}

bool invert_matrix(const Matrix& matrix, Inverse& inverse) {
  std::array<std::array<long double, 2 * kOrder>, kOrder> work{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      work[row][column] = static_cast<long double>(matrix[row][column]);
      work[row][column + kOrder] = row == column ? 1.0L : 0.0L;
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
    if (std::fabs(work[pivot_row][column]) < kSingularTolerance) {
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
      inverse[row][column] = work[row][column + kOrder];
    }
  }
  return true;
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
    if (written == 0) {
      throw std::runtime_error("short write");
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
            ? "refusing to append to an existing telemetry log: " +
                  path.string()
            : "cannot create telemetry log: " +
                  std::string(std::strerror(errno)));
  }
  if (::fsync(descriptor) != 0) {
    const int saved_errno = errno;
    ::close(descriptor);
    throw std::runtime_error("cannot sync telemetry log: " +
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
    throw std::runtime_error("cannot open telemetry log: " +
                             std::string(std::strerror(errno)));
  }
  try {
    write_all(descriptor, record);
    write_all(descriptor, "\n");
    if (::fsync(descriptor) != 0) {
      throw std::runtime_error("cannot sync telemetry log: " +
                               std::string(std::strerror(errno)));
    }
    if (::close(descriptor) != 0) {
      throw std::runtime_error("cannot close telemetry log: " +
                               std::string(std::strerror(errno)));
    }
  } catch (...) {
    ::close(descriptor);
    throw;
  }
}

struct Arguments {
  std::vector<std::filesystem::path> elites;
  std::filesystem::path output;
  std::filesystem::path log;
  std::filesystem::path snapshot;
  std::uint64_t seed = 23;
  std::size_t beam_width = 1200;
  std::size_t exact_pool = 4800;
  std::size_t maximum_pairs = 0;
  double random_fraction = 0.15;
  double seconds = 3600.0;
  bool bidirectional = true;
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
    if (option == "--elite") {
      arguments.elites.emplace_back(require_value());
    } else if (option == "--output") {
      arguments.output = require_value();
    } else if (option == "--log") {
      arguments.log = require_value();
    } else if (option == "--snapshot") {
      arguments.snapshot = require_value();
    } else if (option == "--seed") {
      arguments.seed = parse_unsigned(require_value(), option);
    } else if (option == "--beam-width") {
      arguments.beam_width =
          static_cast<std::size_t>(parse_unsigned(require_value(), option));
    } else if (option == "--exact-pool") {
      arguments.exact_pool =
          static_cast<std::size_t>(parse_unsigned(require_value(), option));
    } else if (option == "--max-pairs") {
      arguments.maximum_pairs =
          static_cast<std::size_t>(parse_unsigned(require_value(), option));
    } else if (option == "--random-fraction") {
      arguments.random_fraction = parse_double(require_value(), option);
    } else if (option == "--seconds") {
      arguments.seconds = parse_double(require_value(), option);
    } else if (option == "--one-way") {
      arguments.bidirectional = false;
    } else if (option == "--help") {
      std::cout
          << "usage: path_relink --elite MATRIX --elite MATRIX [...]"
          << " --output MATRIX --log JSONL [--snapshot JSON]"
          << " [--seed N] [--beam-width 1200] [--exact-pool 4800]"
          << " [--random-fraction 0.15] [--max-pairs 0]"
          << " [--seconds 3600] [--one-way]\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }

  if (arguments.elites.size() < 2) {
    throw std::runtime_error("at least two --elite matrices are required");
  }
  if (arguments.output.empty() || arguments.log.empty()) {
    throw std::runtime_error("--output and --log are required");
  }
  if (arguments.snapshot.empty()) {
    arguments.snapshot = arguments.log.string() + ".snapshot.json";
  }
  if (arguments.beam_width == 0 ||
      arguments.exact_pool < arguments.beam_width) {
    throw std::runtime_error(
        "--exact-pool must be at least the positive --beam-width");
  }
  if (!(arguments.random_fraction >= 0.0 &&
        arguments.random_fraction <= 0.5)) {
    throw std::runtime_error("--random-fraction must be in [0,0.5]");
  }
  if (!(arguments.seconds > 0.0)) {
    throw std::runtime_error("--seconds must be positive");
  }

  auto normalized = [](const std::filesystem::path& path) {
    std::error_code error;
    const std::filesystem::path result =
        std::filesystem::weakly_canonical(
            std::filesystem::absolute(path), error);
    if (!error) return result;
    return std::filesystem::absolute(path).lexically_normal();
  };
  const auto output = normalized(arguments.output);
  const auto log = normalized(arguments.log);
  const auto snapshot = normalized(arguments.snapshot);
  if (output == log || output == snapshot || log == snapshot) {
    throw std::runtime_error(
        "--output, --log, and --snapshot must be distinct");
  }
  for (const auto& elite : arguments.elites) {
    const auto input = normalized(elite);
    if (input == output || input == log || input == snapshot) {
      throw std::runtime_error("an output path aliases an elite input");
    }
  }
  return arguments;
}

struct MaskHash {
  std::size_t operator()(const Mask& mask) const noexcept {
    std::uint64_t result = 0x9e3779b97f4a7c15ULL;
    for (const std::uint64_t word : mask) {
      const std::uint64_t mixed =
          word + 0x9e3779b97f4a7c15ULL + (result << 6U) + (result >> 2U);
      result ^= mixed;
    }
    return static_cast<std::size_t>(result);
  }
};

bool mask_contains(const Mask& mask, int index) {
  return (mask[static_cast<std::size_t>(index) / 64U] &
          (std::uint64_t{1} << (static_cast<unsigned>(index) & 63U))) != 0;
}

void mask_insert(Mask& mask, int index) {
  mask[static_cast<std::size_t>(index) / 64U] |=
      std::uint64_t{1} << (static_cast<unsigned>(index) & 63U);
}

struct State {
  Mask mask{};
  Matrix matrix{};
  Inverse inverse{};
  Wide score = 0;
  long double projected_log_score =
      -std::numeric_limits<long double>::infinity();
  std::uint64_t diversity_key = 0;
};

struct Candidate {
  Mask mask{};
  std::size_t parent = 0;
  std::uint16_t difference_index = 0;
  long double projected_log_score =
      -std::numeric_limits<long double>::infinity();
  std::uint64_t diversity_key = 0;
};

struct Statistics {
  std::uint64_t generated = 0;
  std::uint64_t unique_candidates = 0;
  std::uint64_t exact_checks = 0;
  std::uint64_t singular_candidates = 0;
  std::uint64_t retained_states = 0;
  std::uint64_t completed_depths = 0;
  std::uint64_t completed_pairs = 0;
  std::uint64_t promotions = 0;
};

struct PairJob {
  std::size_t source = 0;
  std::size_t target = 0;
};

struct Campaign {
  const Arguments& arguments;
  const std::vector<Matrix>& elites;
  const std::vector<Wide>& elite_scores;
  std::mt19937_64 randomizer;
  Statistics statistics;
  Matrix best_matrix{};
  Wide best_score = 0;
  Wide depth_best_score = 0;
  Clock::time_point started;
  Clock::time_point deadline;
  std::uint64_t snapshot_nonce = 0;

  double elapsed_seconds() const {
    return std::chrono::duration<double>(Clock::now() - started).count();
  }

  bool should_stop() const {
    return stop_requested || Clock::now() >= deadline;
  }

  void persist_record(const std::string& record) {
    append_log(arguments.log, record);
    atomic_write(arguments.snapshot, record + "\n", "path-relink-snapshot",
                 ++snapshot_nonce);
  }

  std::string event_record(std::string_view event, const PairJob& job,
                           int depth, int difference_count,
                           std::size_t beam_count,
                           std::uint64_t depth_generated,
                           std::uint64_t depth_unique,
                           std::uint64_t depth_exact,
                           Wide pair_best, bool complete,
                           std::string_view termination) const {
    std::ostringstream output;
    output << "{\"event\":\"" << event << "\""
           << ",\"seed\":" << arguments.seed
           << ",\"beam_width\":" << arguments.beam_width
           << ",\"exact_pool\":" << arguments.exact_pool
           << ",\"random_fraction\":" << arguments.random_fraction
           << ",\"bidirectional\":"
           << (arguments.bidirectional ? "true" : "false")
           << ",\"maximum_pairs\":" << arguments.maximum_pairs
           << ",\"elites\":[";
    for (std::size_t index = 0; index < elites.size(); ++index) {
      if (index != 0) output << ',';
      output << "{\"path\":\""
             << json_escape(arguments.elites[index].string())
             << "\",\"absolute_determinant\":\""
             << wide_to_string(elite_scores[index])
             << "\",\"row_major_sign_bits_hex\":\""
             << matrix_sign_bits_hex(elites[index]) << "\"}";
    }
    output << ']'
           << ",\"source_index\":" << job.source
           << ",\"target_index\":" << job.target
           << ",\"source\":\""
           << json_escape(arguments.elites[job.source].string()) << "\""
           << ",\"target\":\""
           << json_escape(arguments.elites[job.target].string()) << "\""
           << ",\"depth\":" << depth
           << ",\"difference_entries\":" << difference_count
           << ",\"beam_states\":" << beam_count
           << ",\"depth_generated\":" << depth_generated
           << ",\"depth_unique\":" << depth_unique
           << ",\"depth_exact_checks\":" << depth_exact
           << ",\"generated\":" << statistics.generated
           << ",\"unique_candidates\":" << statistics.unique_candidates
           << ",\"exact_checks\":" << statistics.exact_checks
           << ",\"singular_candidates\":"
           << statistics.singular_candidates
           << ",\"retained_states\":" << statistics.retained_states
           << ",\"completed_depths\":" << statistics.completed_depths
           << ",\"completed_pairs\":" << statistics.completed_pairs
           << ",\"promotions\":" << statistics.promotions
           << ",\"pair_best_score\":\"" << wide_to_string(pair_best) << "\""
           << ",\"depth_best_score\":\""
           << wide_to_string(depth_best_score) << "\""
           << ",\"best_score\":\"" << wide_to_string(best_score) << "\""
           << ",\"elapsed_seconds\":" << std::fixed << std::setprecision(6)
           << elapsed_seconds()
           << ",\"complete\":" << (complete ? "true" : "false")
           << ",\"termination\":\"" << termination << "\"}";
    return output.str();
  }

  void promote_if_better(const Matrix& matrix, Wide score,
                         const PairJob& job, int depth,
                         int difference_count) {
    if (score <= best_score) return;
    best_score = score;
    best_matrix = matrix;
    ++statistics.promotions;
    atomic_write(arguments.output, matrix_bytes(best_matrix),
                 "path-relink-candidate", statistics.exact_checks);
    persist_record(event_record("promotion", job, depth, difference_count, 0,
                                0, 0, 0, score, false, "running"));
    std::cout << "promotion |det|=" << wide_to_string(score)
              << " source=" << job.source << " target=" << job.target
              << " depth=" << depth << '/' << difference_count << '\n'
              << std::flush;
  }

  bool run_pair(const PairJob& job) {
    std::vector<int> differences;
    differences.reserve(kEntries);
    for (int entry = 0; entry < kEntries; ++entry) {
      const int row = entry / kOrder;
      const int column = entry % kOrder;
      if (elites[job.source][row][column] !=
          elites[job.target][row][column]) {
        differences.push_back(entry);
      }
    }

    const int difference_count = static_cast<int>(differences.size());
    if (difference_count == 0) {
      depth_best_score = elite_scores[job.source];
      ++statistics.completed_pairs;
      persist_record(event_record("pair_finished", job, 0, 0, 0, 0, 0, 0,
                                  elite_scores[job.source], true,
                                  "identical_endpoints"));
      return true;
    }

    State initial;
    initial.matrix = elites[job.source];
    initial.score = elite_scores[job.source];
    initial.projected_log_score =
        std::log(static_cast<long double>(initial.score));
    initial.diversity_key = randomizer();
    if (!invert_matrix(initial.matrix, initial.inverse)) {
      throw std::runtime_error("nonsingular elite failed numerical inversion");
    }

    Wide pair_best =
        std::max(elite_scores[job.source], elite_scores[job.target]);
    depth_best_score = initial.score;
    std::vector<State> beam;
    beam.push_back(std::move(initial));
    persist_record(event_record("pair_started", job, 0, difference_count,
                                beam.size(), 0, 0, 0, pair_best, false,
                                "running"));

    for (int depth = 1; depth <= difference_count; ++depth) {
      depth_best_score = 0;
      if (should_stop()) {
        persist_record(event_record(
            "pair_interrupted", job, depth - 1, difference_count, beam.size(),
            0, 0, 0, pair_best, false,
            stop_requested ? "signal" : "time_limit"));
        return false;
      }

      const std::uint64_t generated_before = statistics.generated;
      const std::uint64_t unique_before = statistics.unique_candidates;
      const std::uint64_t exact_before = statistics.exact_checks;
      std::unordered_map<Mask, Candidate, MaskHash> unique;
      const std::size_t estimated =
          beam.size() *
          static_cast<std::size_t>(difference_count - depth + 1);
      unique.reserve(std::min<std::size_t>(estimated, 1000000));

      std::uint64_t stop_probe = 0;
      for (std::size_t parent_index = 0; parent_index < beam.size();
           ++parent_index) {
        const State& parent = beam[parent_index];
        for (int difference_index = 0;
             difference_index < difference_count; ++difference_index) {
          if (mask_contains(parent.mask, difference_index)) continue;
          ++statistics.generated;
          if ((++stop_probe & 4095U) == 0U && should_stop()) break;

          const int entry = differences[difference_index];
          const int row = entry / kOrder;
          const int column = entry % kOrder;
          const long double delta =
              -2.0L * static_cast<long double>(parent.matrix[row][column]);
          const long double ratio =
              1.0L + delta * parent.inverse[column][row];
          if (std::fabs(ratio) < kSingularTolerance) continue;

          Candidate candidate;
          candidate.mask = parent.mask;
          mask_insert(candidate.mask, difference_index);
          candidate.parent = parent_index;
          candidate.difference_index =
              static_cast<std::uint16_t>(difference_index);
          candidate.projected_log_score =
              std::log(static_cast<long double>(parent.score)) +
              std::log(std::fabs(ratio));
          candidate.diversity_key = randomizer();

          const auto [position, inserted] =
              unique.try_emplace(candidate.mask, candidate);
          if (!inserted &&
              candidate.projected_log_score >
                  position->second.projected_log_score) {
            position->second = candidate;
          }
        }
        if (should_stop()) break;
      }
      if (should_stop()) {
        persist_record(event_record(
            "pair_interrupted", job, depth - 1, difference_count, beam.size(),
            statistics.generated - generated_before,
            statistics.unique_candidates - unique_before,
            statistics.exact_checks - exact_before, pair_best, false,
            stop_requested ? "signal" : "time_limit"));
        return false;
      }

      statistics.unique_candidates += unique.size();
      std::vector<Candidate> ranked;
      ranked.reserve(unique.size());
      for (auto& [mask, candidate] : unique) {
        static_cast<void>(mask);
        ranked.push_back(std::move(candidate));
      }
      std::sort(
          ranked.begin(), ranked.end(),
          [](const Candidate& left, const Candidate& right) {
            if (left.projected_log_score != right.projected_log_score) {
              return left.projected_log_score >
                     right.projected_log_score;
            }
            return left.diversity_key < right.diversity_key;
          });
      if (ranked.size() > arguments.exact_pool) {
        ranked.resize(arguments.exact_pool);
      }

      std::vector<State> evaluated;
      evaluated.reserve(ranked.size());
      for (const Candidate& candidate : ranked) {
        State child;
        child.mask = candidate.mask;
        child.matrix = beam[candidate.parent].matrix;
        const int entry = differences[candidate.difference_index];
        child.matrix[entry / kOrder][entry % kOrder] *= -1;
        child.score = absolute(exact_determinant(child.matrix));
        ++statistics.exact_checks;
        if (child.score == 0 || !invert_matrix(child.matrix, child.inverse)) {
          ++statistics.singular_candidates;
          continue;
        }
        child.projected_log_score =
            std::log(static_cast<long double>(child.score));
        child.diversity_key = candidate.diversity_key;
        depth_best_score = std::max(depth_best_score, child.score);
        pair_best = std::max(pair_best, child.score);
        promote_if_better(child.matrix, child.score, job, depth,
                          difference_count);
        evaluated.push_back(std::move(child));
        if ((statistics.exact_checks & 1023U) == 0U && should_stop()) break;
      }
      if (should_stop()) {
        persist_record(event_record(
            "pair_interrupted", job, depth - 1, difference_count, beam.size(),
            statistics.generated - generated_before,
            statistics.unique_candidates - unique_before,
            statistics.exact_checks - exact_before, pair_best, false,
            stop_requested ? "signal" : "time_limit"));
        return false;
      }
      if (evaluated.empty()) {
        ++statistics.completed_pairs;
        persist_record(event_record(
            "pair_finished", job, depth, difference_count, 0,
            statistics.generated - generated_before,
            statistics.unique_candidates - unique_before,
            statistics.exact_checks - exact_before, pair_best, true,
            "no_nonsingular_states"));
        return true;
      }

      std::sort(evaluated.begin(), evaluated.end(),
                [](const State& left, const State& right) {
                  if (left.score != right.score) {
                    return left.score > right.score;
                  }
                  return left.diversity_key < right.diversity_key;
                });
      const std::size_t keep =
          std::min(arguments.beam_width, evaluated.size());
      const std::size_t random_count = std::min(
          keep,
          static_cast<std::size_t>(std::llround(
              static_cast<double>(keep) * arguments.random_fraction)));
      const std::size_t top_count = keep - random_count;

      std::vector<State> next;
      next.reserve(keep);
      for (std::size_t index = 0; index < top_count; ++index) {
        next.push_back(std::move(evaluated[index]));
      }
      std::vector<std::size_t> tail_indices;
      tail_indices.reserve(evaluated.size() - top_count);
      for (std::size_t index = top_count; index < evaluated.size(); ++index) {
        tail_indices.push_back(index);
      }
      std::shuffle(tail_indices.begin(), tail_indices.end(), randomizer);
      for (std::size_t index = 0;
           index < random_count && index < tail_indices.size(); ++index) {
        next.push_back(std::move(evaluated[tail_indices[index]]));
      }
      beam = std::move(next);
      statistics.retained_states += beam.size();
      ++statistics.completed_depths;

      const bool pair_complete = depth == difference_count;
      if (pair_complete) ++statistics.completed_pairs;
      persist_record(event_record(
          pair_complete ? "pair_finished" : "depth_finished",
          job, depth, difference_count, beam.size(),
          statistics.generated - generated_before,
          statistics.unique_candidates - unique_before,
          statistics.exact_checks - exact_before, pair_best,
          pair_complete,
          pair_complete ? "endpoint_reached" : "running"));
    }
    return true;
  }
};

std::vector<PairJob> pair_schedule(std::size_t elite_count,
                                   bool bidirectional) {
  std::vector<PairJob> jobs;
  auto append = [&](std::size_t first, std::size_t second) {
    jobs.push_back({first, second});
    if (bidirectional) jobs.push_back({second, first});
  };

  // Visit every nontrivial switch from the parent elite before considering
  // paths between switched representatives.
  for (std::size_t target = 1; target < elite_count; ++target) {
    append(0, target);
  }
  for (std::size_t source = 1; source < elite_count; ++source) {
    for (std::size_t target = source + 1; target < elite_count; ++target) {
      append(source, target);
    }
  }
  return jobs;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);

    for (const auto& path :
         {arguments.output, arguments.log, arguments.snapshot}) {
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

    std::vector<Matrix> elites;
    std::vector<Wide> elite_scores;
    elites.reserve(arguments.elites.size());
    elite_scores.reserve(arguments.elites.size());
    Matrix best_matrix{};
    Wide best_score = 0;
    for (const auto& path : arguments.elites) {
      Matrix matrix = read_matrix(path);
      const Wide score = absolute(exact_determinant(matrix));
      if (score == 0) {
        throw std::runtime_error("elite matrix is singular: " + path.string());
      }
      if (score > best_score) {
        best_score = score;
        best_matrix = matrix;
      }
      elites.push_back(std::move(matrix));
      elite_scores.push_back(score);
    }

    create_log(arguments.log);
    atomic_write(arguments.output, matrix_bytes(best_matrix),
                 "path-relink-candidate", 0);

    Campaign campaign{
        arguments,
        elites,
        elite_scores,
        std::mt19937_64(arguments.seed),
        {},
        best_matrix,
        best_score,
        best_score,
        Clock::now(),
        {},
        0};
    campaign.deadline =
        campaign.started + std::chrono::duration_cast<Clock::duration>(
                               std::chrono::duration<double>(
                                   arguments.seconds));

    std::ostringstream start;
    start << "{\"event\":\"start\",\"seed\":" << arguments.seed
          << ",\"beam_width\":" << arguments.beam_width
          << ",\"exact_pool\":" << arguments.exact_pool
          << ",\"random_fraction\":" << arguments.random_fraction
          << ",\"bidirectional\":"
          << (arguments.bidirectional ? "true" : "false")
          << ",\"elite_count\":" << arguments.elites.size()
          << ",\"elites\":[";
    for (std::size_t index = 0; index < arguments.elites.size(); ++index) {
      if (index != 0) start << ',';
      start << "{\"path\":\""
            << json_escape(arguments.elites[index].string()) << "\""
            << ",\"absolute_determinant\":\""
            << wide_to_string(elite_scores[index]) << "\"}";
    }
    start << "],\"best_score\":\"" << wide_to_string(best_score)
          << "\",\"complete\":false,\"termination\":\"running\"}";
    campaign.persist_record(start.str());

    std::vector<PairJob> jobs =
        pair_schedule(elites.size(), arguments.bidirectional);
    if (arguments.maximum_pairs != 0 &&
        jobs.size() > arguments.maximum_pairs) {
      jobs.resize(arguments.maximum_pairs);
    }

    std::cout << "start elites=" << elites.size()
              << " pairs=" << jobs.size()
              << " beam_width=" << arguments.beam_width
              << " exact_pool=" << arguments.exact_pool
              << " |det|=" << wide_to_string(best_score) << '\n'
              << std::flush;

    bool complete = true;
    PairJob last_job{};
    for (const PairJob& job : jobs) {
      last_job = job;
      if (!campaign.run_pair(job)) {
        complete = false;
        break;
      }
    }
    if (campaign.statistics.completed_pairs != jobs.size()) complete = false;

    const std::string termination =
        complete ? "all_pairs_completed"
                 : (stop_requested ? "signal" : "time_limit");
    campaign.persist_record(campaign.event_record(
        "finished", last_job, 0, 0, 0, 0, 0, 0, campaign.best_score,
        complete, termination));

    std::cout << "finished complete=" << (complete ? "true" : "false")
              << " pairs=" << campaign.statistics.completed_pairs << '/'
              << jobs.size()
              << " exact_checks=" << campaign.statistics.exact_checks
              << " |det|=" << wide_to_string(campaign.best_score) << '\n';
    return complete ? 0 : 3;
  } catch (const std::exception& error) {
    std::cerr << "path_relink: " << error.what() << '\n';
    return 2;
  }
}
