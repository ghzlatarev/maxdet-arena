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
constexpr int kMaximumFlips = 24;
constexpr long double kSingularTolerance = 1.0e-22L;

// For a flip set S, write A_S = A + U V^T.  The matrix determinant
// lemma gives det(A_S) / det(A) = det(K), where
// K = I + V^T A^{-1} U has order |S| <= 24.  A beam state stores K^{-1};
// adjoining one arbitrary cell then costs O(|S|^2) via its Schur complement.
// Floating point only ranks candidates.  Every retained beam state and every
// promoted swap is reconstructed as a {-1,+1} matrix and checked with exact
// fraction-free Bareiss elimination before it can update the checkpoint.

using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using Inverse =
    std::array<std::array<long double, kOrder>, kOrder>;
using Wide = __int128_t;
using Clock = std::chrono::steady_clock;

volatile std::sig_atomic_t stop_requested = 0;

void request_stop(int) { stop_requested = 1; }

struct Arguments {
  std::filesystem::path start;
  std::filesystem::path output;
  std::filesystem::path tie_output;
  std::filesystem::path log;
  std::uint64_t seed = 23;
  std::uint64_t runs = 1;
  std::size_t beam_width = 50000;
  int minimum_flips = 4;
  int maximum_flips = 8;
  int maximum_per_row = 3;
  int maximum_per_column = 3;
  std::size_t refine_states = 0;
  std::size_t random_refinements = 0;
  int swap_passes = 0;
  double seconds = 3600.0;
  double heartbeat_seconds = 30.0;
};

struct State {
  std::array<std::uint16_t, kMaximumFlips> entries{};
  std::array<long double, kMaximumFlips * kMaximumFlips> small_inverse{};
  std::array<std::uint8_t, kOrder> row_counts{};
  std::array<std::uint8_t, kOrder> column_counts{};
  std::uint16_t last_rank = 0;
  std::uint8_t size = 0;
  long double determinant_ratio = 1.0L;
  long double log_abs_ratio = 0.0L;
  Wide exact_score = 0;
};

struct Candidate {
  std::size_t parent = 0;
  std::uint16_t entry = 0;
  std::uint16_t rank = 0;
  long double log_abs_ratio =
      -std::numeric_limits<long double>::infinity();
};

struct CandidateMinHeap {
  bool operator()(const Candidate& left, const Candidate& right) const {
    if (left.log_abs_ratio != right.log_abs_ratio) {
      return left.log_abs_ratio > right.log_abs_ratio;
    }
    if (left.parent != right.parent) return left.parent > right.parent;
    return left.rank > right.rank;
  }
};

struct Statistics {
  std::uint64_t generated = 0;
  std::uint64_t exact_checks = 0;
  std::uint64_t completed_depths = 0;
  std::uint64_t completed_runs = 0;
  std::uint64_t swap_neighbors = 0;
  std::uint64_t swap_improvements = 0;
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
        ("." + path.filename().string() + ".multiflip-" +
         std::to_string(static_cast<long long>(::getpid())) + "-" +
         std::to_string(nonce) + "-" + std::to_string(attempt) + ".tmp");
    int flags = O_WRONLY | O_CREAT | O_EXCL;
#ifdef O_CLOEXEC
    flags |= O_CLOEXEC;
#endif
    descriptor = ::open(temporary.c_str(), flags, 0644);
    if (descriptor >= 0) break;
    if (errno != EEXIST) {
      throw std::runtime_error(
          "cannot create checkpoint: " +
          std::string(std::strerror(errno)));
    }
  }
  if (descriptor < 0) {
    throw std::runtime_error("cannot allocate checkpoint temporary file");
  }

  try {
    write_all(descriptor, bytes);
    if (::fsync(descriptor) != 0) {
      throw std::runtime_error(
          "cannot sync checkpoint: " +
          std::string(std::strerror(errno)));
    }
    if (::close(descriptor) != 0) {
      descriptor = -1;
      throw std::runtime_error(
          "cannot close checkpoint: " +
          std::string(std::strerror(errno)));
    }
    descriptor = -1;
    std::filesystem::rename(temporary, path);
    sync_directory(directory);
  } catch (...) {
    const int saved_errno = errno;
    if (descriptor >= 0) ::close(descriptor);
    std::error_code ignored;
    std::filesystem::remove(temporary, ignored);
    errno = saved_errno;
    throw;
  }
}

Inverse invert_matrix(const Matrix& matrix) {
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
      throw std::runtime_error("start matrix is numerically singular");
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

  Inverse inverse{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      inverse[row][column] = work[row][column + kOrder];
    }
  }
  return inverse;
}

long double small_at(const State& state, int row, int column) {
  return state.small_inverse[
      static_cast<std::size_t>(row) * kMaximumFlips + column];
}

long double extension_schur(const State& state, std::uint16_t entry,
                             const Matrix& matrix, const Inverse& inverse) {
  const int new_row = entry / kOrder;
  const int new_column = entry % kOrder;
  const long double new_delta =
      -2.0L * static_cast<long double>(matrix[new_row][new_column]);
  const int size = state.size;

  std::array<long double, kMaximumFlips> right{};
  std::array<long double, kMaximumFlips> lower{};
  std::array<long double, kMaximumFlips> transformed{};
  for (int index = 0; index < size; ++index) {
    const int old_entry = state.entries[index];
    const int old_row = old_entry / kOrder;
    const int old_column = old_entry % kOrder;
    const long double old_delta =
        -2.0L * static_cast<long double>(matrix[old_row][old_column]);
    right[index] = inverse[old_column][new_row] * new_delta;
    lower[index] = inverse[new_column][old_row] * old_delta;
  }
  for (int row = 0; row < size; ++row) {
    for (int column = 0; column < size; ++column) {
      transformed[row] +=
          small_at(state, row, column) * right[column];
    }
  }

  long double schur =
      1.0L + inverse[new_column][new_row] * new_delta;
  for (int index = 0; index < size; ++index) {
    schur -= lower[index] * transformed[index];
  }
  return schur;
}

State extend_state(const State& state, const Candidate& candidate,
                   const Matrix& matrix, const Inverse& inverse) {
  State result = state;
  const int size = state.size;
  const int new_row = candidate.entry / kOrder;
  const int new_column = candidate.entry % kOrder;
  const long double new_delta =
      -2.0L * static_cast<long double>(matrix[new_row][new_column]);

  std::array<long double, kMaximumFlips> right{};
  std::array<long double, kMaximumFlips> lower{};
  std::array<long double, kMaximumFlips> transformed_right{};
  std::array<long double, kMaximumFlips> transformed_lower{};
  for (int index = 0; index < size; ++index) {
    const int old_entry = state.entries[index];
    const int old_row = old_entry / kOrder;
    const int old_column = old_entry % kOrder;
    const long double old_delta =
        -2.0L * static_cast<long double>(matrix[old_row][old_column]);
    right[index] = inverse[old_column][new_row] * new_delta;
    lower[index] = inverse[new_column][old_row] * old_delta;
  }
  for (int row = 0; row < size; ++row) {
    for (int column = 0; column < size; ++column) {
      transformed_right[row] +=
          small_at(state, row, column) * right[column];
      transformed_lower[column] +=
          lower[row] * small_at(state, row, column);
    }
  }

  long double schur =
      1.0L + inverse[new_column][new_row] * new_delta;
  for (int index = 0; index < size; ++index) {
    schur -= lower[index] * transformed_right[index];
  }
  if (std::fabs(schur) < kSingularTolerance) {
    throw std::runtime_error("retained extension became numerically singular");
  }

  for (int row = 0; row < size; ++row) {
    for (int column = 0; column < size; ++column) {
      result.small_inverse[
          static_cast<std::size_t>(row) * kMaximumFlips + column] =
          small_at(state, row, column) +
          transformed_right[row] * transformed_lower[column] / schur;
    }
    result.small_inverse[
        static_cast<std::size_t>(row) * kMaximumFlips + size] =
        -transformed_right[row] / schur;
    result.small_inverse[
        static_cast<std::size_t>(size) * kMaximumFlips + row] =
        -transformed_lower[row] / schur;
  }
  result.small_inverse[
      static_cast<std::size_t>(size) * kMaximumFlips + size] =
      1.0L / schur;
  result.entries[size] = candidate.entry;
  result.last_rank = candidate.rank;
  result.size = static_cast<std::uint8_t>(size + 1);
  ++result.row_counts[new_row];
  ++result.column_counts[new_column];
  result.determinant_ratio *= schur;
  result.log_abs_ratio = candidate.log_abs_ratio;
  result.exact_score = 0;
  return result;
}

Matrix materialize(const Matrix& baseline, const State& state) {
  Matrix result = baseline;
  for (int index = 0; index < state.size; ++index) {
    const int entry = state.entries[index];
    result[entry / kOrder][entry % kOrder] *= -1;
  }
  return result;
}

State build_state(
    const std::array<std::uint16_t, kMaximumFlips>& entries, int size,
    const Matrix& matrix, const Inverse& inverse) {
  State state;
  for (int index = 0; index < size; ++index) {
    Candidate candidate;
    candidate.entry = entries[index];
    candidate.rank = static_cast<std::uint16_t>(index);
    const long double schur =
        extension_schur(state, candidate.entry, matrix, inverse);
    if (std::fabs(schur) < kSingularTolerance) {
      throw std::runtime_error("selected flip set is numerically singular");
    }
    candidate.log_abs_ratio =
        state.log_abs_ratio + std::log(std::fabs(schur));
    state = extend_state(state, candidate, matrix, inverse);
  }
  return state;
}

bool contains_entry(const State& state, std::uint16_t entry) {
  return std::find(
             state.entries.begin(),
             state.entries.begin() + state.size,
             entry) != state.entries.begin() + state.size;
}

State refine_by_swaps(
    State state, const Matrix& baseline, const Inverse& inverse,
    Wide baseline_score, const Arguments& arguments, Statistics& statistics) {
  for (int pass = 0; pass < arguments.swap_passes; ++pass) {
    long double best_log_abs_ratio = state.log_abs_ratio;
    int best_removed = -1;
    std::uint16_t best_added = 0;

    for (int removed = 0; removed < state.size; ++removed) {
      std::array<std::uint16_t, kMaximumFlips> partial_entries{};
      int partial_size = 0;
      for (int index = 0; index < state.size; ++index) {
        if (index != removed) {
          partial_entries[partial_size++] = state.entries[index];
        }
      }
      const State partial =
          build_state(partial_entries, partial_size, baseline, inverse);
      for (int raw_entry = 0; raw_entry < kEntries; ++raw_entry) {
        const auto entry = static_cast<std::uint16_t>(raw_entry);
        if (contains_entry(state, entry)) continue;
        const int row = raw_entry / kOrder;
        const int column = raw_entry % kOrder;
        if (partial.row_counts[row] >= arguments.maximum_per_row ||
            partial.column_counts[column] >=
                arguments.maximum_per_column) {
          continue;
        }
        const long double schur =
            extension_schur(partial, entry, baseline, inverse);
        ++statistics.generated;
        ++statistics.swap_neighbors;
        if (std::fabs(schur) < kSingularTolerance) continue;
        const long double candidate_log_abs_ratio =
            partial.log_abs_ratio + std::log(std::fabs(schur));
        if (candidate_log_abs_ratio > best_log_abs_ratio) {
          best_log_abs_ratio = candidate_log_abs_ratio;
          best_removed = removed;
          best_added = entry;
        }
      }
    }

    if (best_removed < 0) break;
    std::array<std::uint16_t, kMaximumFlips> candidate_entries{};
    int candidate_size = 0;
    for (int index = 0; index < state.size; ++index) {
      if (index != best_removed) {
        candidate_entries[candidate_size++] = state.entries[index];
      }
    }
    candidate_entries[candidate_size++] = best_added;
    std::sort(
        candidate_entries.begin(),
        candidate_entries.begin() + candidate_size);
    State candidate = build_state(
        candidate_entries, candidate_size, baseline, inverse);
    const Matrix matrix = materialize(baseline, candidate);
    candidate.exact_score = absolute(exact_determinant(matrix));
    ++statistics.exact_checks;
    const long double predicted_score =
        static_cast<long double>(baseline_score) *
        std::exp(candidate.log_abs_ratio);
    const long double relative_error =
        std::fabs(
            predicted_score -
            static_cast<long double>(candidate.exact_score)) /
        std::max(1.0L, static_cast<long double>(candidate.exact_score));
    if (relative_error > 1.0e-9L) {
      throw std::runtime_error(
          "swap ranking disagrees with exact determinant");
    }
    if (candidate.exact_score <= state.exact_score) break;
    state = std::move(candidate);
    ++statistics.swap_improvements;
  }
  return state;
}

State random_state(
    int size, const Matrix& baseline, const Inverse& inverse,
    const Arguments& arguments, std::mt19937_64& randomizer) {
  std::array<std::uint16_t, kEntries> entries{};
  for (int entry = 0; entry < kEntries; ++entry) {
    entries[entry] = static_cast<std::uint16_t>(entry);
  }
  std::shuffle(entries.begin(), entries.end(), randomizer);

  std::array<std::uint16_t, kMaximumFlips> selected{};
  std::array<std::uint8_t, kOrder> row_counts{};
  std::array<std::uint8_t, kOrder> column_counts{};
  int selected_size = 0;
  for (const auto entry : entries) {
    const int row = entry / kOrder;
    const int column = entry % kOrder;
    if (row_counts[row] >= arguments.maximum_per_row ||
        column_counts[column] >= arguments.maximum_per_column) {
      continue;
    }
    selected[selected_size++] = entry;
    ++row_counts[row];
    ++column_counts[column];
    if (selected_size == size) break;
  }
  if (selected_size != size) {
    throw std::runtime_error("could not generate constrained random state");
  }
  std::sort(selected.begin(), selected.begin() + selected_size);
  return build_state(selected, selected_size, baseline, inverse);
}

void log_event(std::ofstream& log, std::string_view event,
               std::uint64_t run, int depth, std::uint64_t generated,
               std::uint64_t selected, std::uint64_t exact_checks,
               Wide best_score, double elapsed) {
  log << "{\"event\":\"" << event << "\",\"run\":" << run
      << ",\"depth\":" << depth << ",\"generated\":" << generated
      << ",\"selected\":" << selected
      << ",\"exact_checks\":" << exact_checks
      << ",\"best_score\":\"" << wide_to_string(best_score)
      << "\",\"elapsed_seconds\":" << std::fixed << std::setprecision(6)
      << elapsed << "}\n";
  log.flush();
}

std::uint64_t parse_unsigned(const std::string& text,
                             std::string_view option) {
  std::size_t consumed = 0;
  const unsigned long long value = std::stoull(text, &consumed);
  if (consumed != text.size()) {
    throw std::runtime_error("invalid value for " + std::string(option));
  }
  return static_cast<std::uint64_t>(value);
}

int parse_int(const std::string& text, std::string_view option) {
  std::size_t consumed = 0;
  const long value = std::stol(text, &consumed);
  if (consumed != text.size() ||
      value < std::numeric_limits<int>::min() ||
      value > std::numeric_limits<int>::max()) {
    throw std::runtime_error("invalid value for " + std::string(option));
  }
  return static_cast<int>(value);
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
      if (index + 1 >= argc) {
        throw std::runtime_error("missing value for " + option);
      }
      return argv[++index];
    };
    if (option == "--start") {
      arguments.start = require_value();
    } else if (option == "--output") {
      arguments.output = require_value();
    } else if (option == "--tie-output") {
      arguments.tie_output = require_value();
    } else if (option == "--log") {
      arguments.log = require_value();
    } else if (option == "--seed") {
      arguments.seed = parse_unsigned(require_value(), option);
    } else if (option == "--runs") {
      arguments.runs = parse_unsigned(require_value(), option);
    } else if (option == "--beam-width") {
      arguments.beam_width =
          static_cast<std::size_t>(parse_unsigned(require_value(), option));
    } else if (option == "--min-flips") {
      arguments.minimum_flips = parse_int(require_value(), option);
    } else if (option == "--max-flips") {
      arguments.maximum_flips = parse_int(require_value(), option);
    } else if (option == "--max-per-row") {
      arguments.maximum_per_row = parse_int(require_value(), option);
    } else if (option == "--max-per-column") {
      arguments.maximum_per_column = parse_int(require_value(), option);
    } else if (option == "--refine-states") {
      arguments.refine_states =
          static_cast<std::size_t>(parse_unsigned(require_value(), option));
    } else if (option == "--random-refinements") {
      arguments.random_refinements =
          static_cast<std::size_t>(parse_unsigned(require_value(), option));
    } else if (option == "--swap-passes") {
      arguments.swap_passes = parse_int(require_value(), option);
    } else if (option == "--seconds") {
      arguments.seconds = parse_double(require_value(), option);
    } else if (option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds = parse_double(require_value(), option);
    } else if (option == "--help") {
      std::cout
          << "usage: multiflip_beam --start MATRIX --output MATRIX --log JSONL"
          << " [--tie-output MATRIX]"
          << " [--seed N] [--runs N] [--beam-width N]"
          << " [--min-flips 4] [--max-flips 8; limit 24]"
          << " [--max-per-row 3] [--max-per-column 3]"
          << " [--refine-states 0] [--random-refinements 0]"
          << " [--swap-passes 0]"
          << " [--seconds 3600] [--heartbeat-seconds 30]\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }

  if (arguments.start.empty() || arguments.output.empty() ||
      arguments.log.empty()) {
    throw std::runtime_error("--start, --output, and --log are required");
  }
  if (arguments.runs == 0 || arguments.beam_width == 0) {
    throw std::runtime_error("--runs and --beam-width must be positive");
  }
  if (arguments.minimum_flips < 1 ||
      arguments.maximum_flips < arguments.minimum_flips ||
      arguments.maximum_flips > kMaximumFlips) {
    throw std::runtime_error(
        "flip bounds must satisfy 1 <= min <= max <= 24");
  }
  if (arguments.maximum_per_row < 1 ||
      arguments.maximum_per_row > arguments.maximum_flips ||
      arguments.maximum_per_column < 1 ||
      arguments.maximum_per_column > arguments.maximum_flips) {
    throw std::runtime_error("per-line limits must be in [1,max-flips]");
  }
  if (arguments.swap_passes < 0) {
    throw std::runtime_error("--swap-passes must be nonnegative");
  }
  if (((arguments.refine_states == 0 &&
        arguments.random_refinements == 0)) !=
      (arguments.swap_passes == 0)) {
    throw std::runtime_error(
        "refinement counts and --swap-passes must both be zero or positive");
  }
  if (!(arguments.seconds > 0.0) ||
      !(arguments.heartbeat_seconds > 0.0)) {
    throw std::runtime_error("time limits must be positive");
  }
  return arguments;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);

    const Matrix baseline = read_matrix(arguments.start);
    const Inverse inverse = invert_matrix(baseline);
    const Wide baseline_score = absolute(exact_determinant(baseline));
    Wide best_score = baseline_score;
    Matrix best_matrix = baseline;
    atomic_write_matrix(arguments.output, best_matrix, 0);
    bool tie_written = false;

    if (!arguments.log.parent_path().empty()) {
      std::filesystem::create_directories(arguments.log.parent_path());
    }
    std::ofstream log(arguments.log, std::ios::app);
    if (!log) {
      throw std::runtime_error("cannot open log: " + arguments.log.string());
    }

    const auto started = Clock::now();
    const auto deadline =
        started + std::chrono::duration<double>(arguments.seconds);
    auto next_heartbeat =
        started + std::chrono::duration<double>(arguments.heartbeat_seconds);
    Statistics statistics;
    std::mt19937_64 randomizer(arguments.seed);
    std::array<std::uint16_t, kEntries> permutation{};
    for (int entry = 0; entry < kEntries; ++entry) {
      permutation[entry] = static_cast<std::uint16_t>(entry);
    }

    auto checkpoint_frontier_tie =
        [&](const Matrix& matrix, Wide score, std::uint64_t run,
            int depth, std::uint64_t generated, std::uint64_t selected) {
          if (tie_written || arguments.tie_output.empty() ||
              score != baseline_score || matrix == baseline) {
            return;
          }
          atomic_write_matrix(
              arguments.tie_output, matrix, statistics.exact_checks);
          tie_written = true;
          const double elapsed = std::chrono::duration<double>(
              Clock::now() - started).count();
          log_event(
              log, "frontier_tie", run, depth, generated, selected,
              statistics.exact_checks, best_score, elapsed);
          std::cout << "retained non-baseline frontier tie run=" << run
                    << " flips=" << depth << '\n'
                    << std::flush;
        };

    log_event(log, "start", 0, 0, 0, 0, 0, best_score, 0.0);
    std::cout << "start |det|=" << wide_to_string(baseline_score)
              << " beam_width=" << arguments.beam_width
              << " flips=" << arguments.minimum_flips << '-'
              << arguments.maximum_flips << " runs=" << arguments.runs
              << '\n'
              << std::flush;

    bool complete = true;
    for (std::uint64_t run = 0;
         run < arguments.runs && !stop_requested; ++run) {
      std::shuffle(permutation.begin(), permutation.end(), randomizer);

      State empty;
      empty.last_rank = 0;
      empty.exact_score = baseline_score;
      std::vector<State> beam{empty};

      for (int depth = 1;
           depth <= arguments.maximum_flips && !stop_requested; ++depth) {
        std::priority_queue<
            Candidate, std::vector<Candidate>, CandidateMinHeap>
            retained;
        std::uint64_t depth_generated = 0;
        bool depth_complete = true;

        for (std::size_t parent_index = 0;
             parent_index < beam.size() && !stop_requested; ++parent_index) {
          const State& parent = beam[parent_index];
          const int first_rank = parent.size == 0
                                     ? 0
                                     : static_cast<int>(parent.last_rank) + 1;
          for (int rank = first_rank; rank < kEntries; ++rank) {
            const std::uint16_t entry = permutation[rank];
            const int row = entry / kOrder;
            const int column = entry % kOrder;
            if (parent.row_counts[row] >= arguments.maximum_per_row ||
                parent.column_counts[column] >=
                    arguments.maximum_per_column) {
              continue;
            }

            const long double schur =
                extension_schur(parent, entry, baseline, inverse);
            ++statistics.generated;
            ++depth_generated;
            if (std::fabs(schur) < kSingularTolerance) continue;
            Candidate candidate;
            candidate.parent = parent_index;
            candidate.entry = entry;
            candidate.rank = static_cast<std::uint16_t>(rank);
            candidate.log_abs_ratio =
                parent.log_abs_ratio + std::log(std::fabs(schur));
            if (retained.size() < arguments.beam_width) {
              retained.push(candidate);
            } else if (candidate.log_abs_ratio >
                       retained.top().log_abs_ratio) {
              retained.pop();
              retained.push(candidate);
            }

            const auto now = Clock::now();
            if (now >= deadline) {
              stop_requested = 1;
              depth_complete = false;
              break;
            }
            if (now >= next_heartbeat) {
              const double elapsed =
                  std::chrono::duration<double>(now - started).count();
              log_event(
                  log, "heartbeat", run, depth, depth_generated,
                  retained.size(), statistics.exact_checks, best_score,
                  elapsed);
              next_heartbeat =
                  now +
                  std::chrono::duration<double>(
                      arguments.heartbeat_seconds);
            }
          }
        }

        if (!depth_complete || stop_requested) {
          complete = false;
          break;
        }
        if (retained.empty()) {
          throw std::runtime_error("beam became empty");
        }

        std::vector<Candidate> candidates;
        candidates.reserve(retained.size());
        while (!retained.empty()) {
          candidates.push_back(retained.top());
          retained.pop();
        }
        std::reverse(candidates.begin(), candidates.end());

        std::vector<State> next_beam;
        next_beam.reserve(candidates.size());
        Wide depth_best = 0;
        for (const Candidate& candidate : candidates) {
          State state =
              extend_state(beam[candidate.parent], candidate, baseline, inverse);
          const Matrix matrix = materialize(baseline, state);
          state.exact_score = absolute(exact_determinant(matrix));
          ++statistics.exact_checks;
          const long double predicted_score =
              static_cast<long double>(baseline_score) *
              std::exp(state.log_abs_ratio);
          const long double relative_error =
              std::fabs(
                  predicted_score -
                  static_cast<long double>(state.exact_score)) /
              std::max(
                  1.0L, static_cast<long double>(state.exact_score));
          if (relative_error > 1.0e-9L) {
            throw std::runtime_error(
                "determinant-lemma ranking disagrees with exact determinant");
          }
          depth_best = std::max(depth_best, state.exact_score);
          if (depth >= arguments.minimum_flips) {
            checkpoint_frontier_tie(
                matrix, state.exact_score, run, depth, depth_generated,
                candidates.size());
          }

          if (depth >= arguments.minimum_flips &&
              state.exact_score > best_score) {
            best_score = state.exact_score;
            best_matrix = matrix;
            atomic_write_matrix(
                arguments.output, best_matrix, statistics.exact_checks);
            const double elapsed = std::chrono::duration<double>(
                Clock::now() - started).count();
            log_event(
                log, "new_best", run, depth, depth_generated,
                candidates.size(), statistics.exact_checks, best_score,
                elapsed);
            std::cout << "new best |det|=" << wide_to_string(best_score)
                      << " run=" << run << " flips=" << depth << " entries=";
            for (int index = 0; index < state.size; ++index) {
              if (index != 0) std::cout << ',';
              std::cout << state.entries[index];
            }
            std::cout << '\n' << std::flush;
          }
          next_beam.push_back(std::move(state));
        }

        std::sort(
            next_beam.begin(), next_beam.end(),
            [](const State& left, const State& right) {
              if (left.exact_score != right.exact_score) {
                return left.exact_score > right.exact_score;
              }
              return left.entries < right.entries;
            });

        const std::uint64_t swaps_before = statistics.swap_neighbors;
        std::size_t refined_count = 0;
        Wide refined_depth_best = 0;
        if (depth >= arguments.minimum_flips &&
            arguments.refine_states != 0) {
          const std::size_t beam_refined_count =
              std::min(arguments.refine_states, next_beam.size());
          const std::size_t refinement_pool = std::min(
              next_beam.size(),
              std::max(
                  beam_refined_count,
                  beam_refined_count >
                          std::numeric_limits<std::size_t>::max() / 16
                      ? next_beam.size()
                      : beam_refined_count * 16));
          for (std::size_t sample = 0;
               sample < beam_refined_count && !stop_requested; ++sample) {
            const std::size_t state_index =
                sample * refinement_pool / beam_refined_count;
            State refined = refine_by_swaps(
                next_beam[state_index], baseline, inverse, baseline_score,
                arguments, statistics);
            ++refined_count;
            refined_depth_best =
                std::max(refined_depth_best, refined.exact_score);
            checkpoint_frontier_tie(
                materialize(baseline, refined), refined.exact_score, run,
                depth, statistics.swap_neighbors - swaps_before,
                refined_count);
            if (refined.exact_score > best_score) {
              best_score = refined.exact_score;
              best_matrix = materialize(baseline, refined);
              atomic_write_matrix(
                  arguments.output, best_matrix, statistics.exact_checks);
              const double improvement_elapsed =
                  std::chrono::duration<double>(
                      Clock::now() - started).count();
              log_event(
                  log, "new_best_swap", run, depth,
                  statistics.swap_neighbors - swaps_before, refined_count,
                  statistics.exact_checks, best_score,
                  improvement_elapsed);
              std::cout
                  << "new best from fixed-cardinality swap |det|="
                  << wide_to_string(best_score) << " run=" << run
                  << " flips=" << depth << '\n'
                  << std::flush;
            }
            if (Clock::now() >= deadline) {
              stop_requested = 1;
            }
          }
        }
        if (depth >= arguments.minimum_flips &&
            arguments.random_refinements != 0) {
          for (std::size_t sample = 0;
               sample < arguments.random_refinements &&
               !stop_requested; ++sample) {
            State refined = random_state(
                depth, baseline, inverse, arguments, randomizer);
            Matrix random_matrix = materialize(baseline, refined);
            refined.exact_score =
                absolute(exact_determinant(random_matrix));
            ++statistics.exact_checks;
            refined = refine_by_swaps(
                refined, baseline, inverse, baseline_score,
                arguments, statistics);
            ++refined_count;
            refined_depth_best =
                std::max(refined_depth_best, refined.exact_score);
            checkpoint_frontier_tie(
                materialize(baseline, refined), refined.exact_score, run,
                depth, statistics.swap_neighbors - swaps_before,
                refined_count);
            if (refined.exact_score > best_score) {
              best_score = refined.exact_score;
              best_matrix = materialize(baseline, refined);
              atomic_write_matrix(
                  arguments.output, best_matrix, statistics.exact_checks);
              const double improvement_elapsed =
                  std::chrono::duration<double>(
                      Clock::now() - started).count();
              log_event(
                  log, "new_best_random_swap", run, depth,
                  statistics.swap_neighbors - swaps_before, refined_count,
                  statistics.exact_checks, best_score,
                  improvement_elapsed);
              std::cout
                  << "new best from random fixed-cardinality swap |det|="
                  << wide_to_string(best_score) << " run=" << run
                  << " flips=" << depth << '\n'
                  << std::flush;
            }
            if (Clock::now() >= deadline) {
              stop_requested = 1;
            }
          }
        }

        beam = std::move(next_beam);
        ++statistics.completed_depths;
        const double elapsed = std::chrono::duration<double>(
            Clock::now() - started).count();
        log_event(
            log, "depth_finished", run, depth, depth_generated, beam.size(),
            statistics.exact_checks, best_score, elapsed);
        std::cout << "depth finished run=" << run << " flips=" << depth
                  << " generated=" << depth_generated
                  << " retained_exact=" << beam.size()
                  << " depth_best=" << wide_to_string(depth_best)
                  << " global_best=" << wide_to_string(best_score)
                  << " refined=" << refined_count
                  << " refined_best="
                  << wide_to_string(refined_depth_best)
                  << " swap_neighbors="
                  << statistics.swap_neighbors - swaps_before
                  << " elapsed=" << std::fixed << std::setprecision(3)
                  << elapsed << "s\n"
                  << std::flush;
        if (stop_requested) {
          complete = false;
          break;
        }
      }
      if (!stop_requested) ++statistics.completed_runs;
    }

    const double elapsed = std::chrono::duration<double>(
        Clock::now() - started).count();
    const bool finished =
        complete && !stop_requested &&
        statistics.completed_runs == arguments.runs;
    log_event(
        log, finished ? "finished" : "stopped",
        statistics.completed_runs, 0, statistics.generated, 0,
        statistics.exact_checks, best_score, elapsed);
    std::cout << (finished ? "finished" : "stopped")
              << " best |det|=" << wide_to_string(best_score)
              << " generated=" << statistics.generated
              << " exact_checks=" << statistics.exact_checks
              << " swap_neighbors=" << statistics.swap_neighbors
              << " swap_improvements=" << statistics.swap_improvements
              << " completed_depths=" << statistics.completed_depths
              << " completed_runs=" << statistics.completed_runs
              << " elapsed=" << std::fixed << std::setprecision(3)
              << elapsed << "s\n";
    return 0;
  } catch (const std::exception& exception) {
    std::cerr << "error: " << exception.what() << '\n';
    return 1;
  }
}
