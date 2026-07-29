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
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

#include <fcntl.h>
#include <unistd.h>

namespace {

constexpr int kOrder = 23;
constexpr int kEntries = kOrder * kOrder;
constexpr long double kInverseZeroTolerance = 1.0e-12L;

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
  std::filesystem::path research_output;
  std::filesystem::path log;
  std::uint64_t seed = 23;
  double seconds = 3600.0;
  double heartbeat_seconds = 60.0;
  int min_flips = 1;
  int max_flips = 4;
  int kick_min = 6;
  int kick_max = 16;
  int block_rows_min = 0;
  int block_rows_max = 0;
  int block_columns_min = 0;
  int block_columns_max = 0;
  std::uint64_t restart_interval = 4096;
  std::uint64_t cooling_period = 2048;
  long double temperature_start = 0.075L;
  long double temperature_end = 0.006L;
  int floor_percent = 58;
};

struct Statistics {
  std::uint64_t proposals = 0;
  std::uint64_t row_proposals = 0;
  std::uint64_t column_proposals = 0;
  std::uint64_t accepted = 0;
  std::uint64_t accepted_uphill = 0;
  std::uint64_t accepted_downhill = 0;
  std::uint64_t accepted_equal = 0;
  std::uint64_t singular_rejections = 0;
  std::uint64_t floor_rejections = 0;
  std::uint64_t restarts = 0;
  std::uint64_t exact_checks = 0;
  std::uint64_t exact_cofactor_fallbacks = 0;
  std::uint64_t perturbed_entries = 0;
  Wide research_score = 0;
};

struct ProposalDescription {
  int first = -1;
  int second = -1;
  int optimized = -1;
  int first_flips = 0;
  int second_flips = 0;
  int block_rows = 0;
  int block_columns = 0;
  bool columns = false;
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

Matrix transpose(const Matrix& matrix) {
  Matrix result{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      result[row][column] = matrix[column][row];
    }
  }
  return result;
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
        ("." + path.filename().string() + ".triple-line-" +
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

bool rebuild_inverse(const Matrix& matrix, Inverse& inverse) {
  std::array<std::array<long double, 2 * kOrder>, kOrder> augmented{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      augmented[row][column] =
          static_cast<long double>(matrix[row][column]);
      augmented[row][column + kOrder] =
          row == column ? 1.0L : 0.0L;
    }
  }

  for (int column = 0; column < kOrder; ++column) {
    int pivot_row = column;
    for (int row = column + 1; row < kOrder; ++row) {
      if (std::fabs(augmented[row][column]) >
          std::fabs(augmented[pivot_row][column])) {
        pivot_row = row;
      }
    }
    if (std::fabs(augmented[pivot_row][column]) <
        std::numeric_limits<long double>::epsilon()) {
      return false;
    }
    if (pivot_row != column) {
      std::swap(augmented[pivot_row], augmented[column]);
    }
    const long double pivot = augmented[column][column];
    for (int inner = 0; inner < 2 * kOrder; ++inner) {
      augmented[column][inner] /= pivot;
    }
    for (int row = 0; row < kOrder; ++row) {
      if (row == column) continue;
      const long double factor = augmented[row][column];
      if (factor == 0.0L) continue;
      for (int inner = 0; inner < 2 * kOrder; ++inner) {
        augmented[row][inner] -=
            factor * augmented[column][inner];
      }
    }
  }

  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      inverse[row][column] =
          augmented[row][column + kOrder];
    }
  }
  return true;
}

Wide exact_row_cofactor(const Matrix& matrix, int row, int column) {
  Matrix basis = matrix;
  basis[row].fill(0);
  basis[row][column] = 1;
  return exact_determinant(basis);
}

void optimize_row_from_cofactors(Matrix& matrix, const Inverse& inverse,
                                 int row, Statistics& statistics) {
  std::array<long double, kOrder> ratios{};
  bool needs_exact_fallback = false;
  for (int column = 0; column < kOrder; ++column) {
    ratios[column] = inverse[column][row];
    needs_exact_fallback =
        needs_exact_fallback ||
        std::fabs(ratios[column]) <= kInverseZeroTolerance;
  }

  int determinant_sign = 0;
  if (needs_exact_fallback) {
    const Wide determinant = exact_determinant(matrix);
    ++statistics.exact_checks;
    if (determinant == 0) {
      throw std::runtime_error(
          "cofactor fallback received a singular matrix");
    }
    determinant_sign = determinant > 0 ? 1 : -1;
  }

  for (int column = 0; column < kOrder; ++column) {
    const long double ratio = ratios[column];
    if (std::fabs(ratio) > kInverseZeroTolerance) {
      // inverse[column][row] = cofactor[row][column] / det(matrix).
      // Multiplying every chosen sign by sign(det) changes only the whole-row
      // sign, so sign(inverse) is also a global |det|-maximizing pattern.
      matrix[row][column] = ratio > 0.0L ? 1 : -1;
      continue;
    }
    const Wide cofactor = exact_row_cofactor(matrix, row, column);
    ++statistics.exact_checks;
    ++statistics.exact_cofactor_fallbacks;
    if (cofactor != 0) {
      const bool inverse_sign_positive =
          (cofactor > 0) == (determinant_sign > 0);
      matrix[row][column] = inverse_sign_positive ? 1 : -1;
    }
  }
}

void perturb_row(Matrix& matrix, int row, int flips,
                 std::mt19937_64& randomizer) {
  std::array<int, kOrder> positions{};
  for (int index = 0; index < kOrder; ++index) positions[index] = index;
  for (int index = 0; index < flips; ++index) {
    std::uniform_int_distribution<int> choose(index, kOrder - 1);
    const int selected = choose(randomizer);
    std::swap(positions[index], positions[selected]);
    matrix[row][positions[index]] *= -1;
  }
}

Matrix propose_triple_move(const Matrix& state,
                           ProposalDescription& description,
                           const Arguments& arguments,
                           std::mt19937_64& randomizer,
                           Statistics& statistics) {
  std::uniform_int_distribution<int> choose_line(0, kOrder - 1);
  description.columns = (randomizer() & 1U) != 0U;

  Matrix oriented = description.columns ? transpose(state) : state;
  if (arguments.block_rows_min != 0) {
    std::uniform_int_distribution<int> choose_rows(
        arguments.block_rows_min, arguments.block_rows_max);
    std::uniform_int_distribution<int> choose_columns(
        arguments.block_columns_min, arguments.block_columns_max);
    description.block_rows = choose_rows(randomizer);
    description.block_columns = choose_columns(randomizer);
    description.first_flips = description.block_columns;
    description.second_flips = description.block_columns;

    std::array<int, kOrder> rows{};
    std::array<int, kOrder> columns{};
    for (int index = 0; index < kOrder; ++index) {
      rows[index] = index;
      columns[index] = index;
    }
    for (int index = 0; index < description.block_rows; ++index) {
      std::uniform_int_distribution<int> choose(index, kOrder - 1);
      std::swap(rows[index], rows[choose(randomizer)]);
    }
    for (int index = 0; index < description.block_columns; ++index) {
      std::uniform_int_distribution<int> choose(index, kOrder - 1);
      std::swap(columns[index], columns[choose(randomizer)]);
    }
    description.first = rows[0];
    description.second = rows[1];
    do {
      description.optimized = choose_line(randomizer);
    } while (std::find(
                 rows.begin(),
                 rows.begin() + description.block_rows,
                 description.optimized) !=
             rows.begin() + description.block_rows);
    for (int row_index = 0;
         row_index < description.block_rows; ++row_index) {
      for (int column_index = 0;
           column_index < description.block_columns; ++column_index) {
        oriented[rows[row_index]][columns[column_index]] *= -1;
      }
    }
    statistics.perturbed_entries +=
        static_cast<std::uint64_t>(description.block_rows) *
        static_cast<std::uint64_t>(description.block_columns);
  } else {
    description.first = choose_line(randomizer);
    do {
      description.second = choose_line(randomizer);
    } while (description.second == description.first);
    do {
      description.optimized = choose_line(randomizer);
    } while (description.optimized == description.first ||
             description.optimized == description.second);
    std::uniform_int_distribution<int> choose_flips(
        arguments.min_flips, arguments.max_flips);
    description.first_flips = choose_flips(randomizer);
    description.second_flips = choose_flips(randomizer);
    perturb_row(
        oriented, description.first, description.first_flips, randomizer);
    perturb_row(
        oriented, description.second, description.second_flips, randomizer);
    statistics.perturbed_entries +=
        static_cast<std::uint64_t>(
            description.first_flips + description.second_flips);
  }

  Inverse inverse{};
  if (!rebuild_inverse(oriented, inverse)) {
    ++statistics.singular_rejections;
    return state;
  }
  optimize_row_from_cofactors(
      oriented, inverse, description.optimized, statistics);
  return description.columns ? transpose(oriented) : oriented;
}

std::uint64_t parse_unsigned(std::string_view text,
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

double parse_double(std::string_view text, std::string_view option,
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
    } else if (option == "--research-output") {
      arguments.research_output = value();
    } else if (option == "--log") {
      arguments.log = value();
    } else if (option == "--seed") {
      arguments.seed = parse_unsigned(value(), option);
    } else if (option == "--seconds") {
      arguments.seconds = parse_double(value(), option, false);
    } else if (option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds =
          parse_double(value(), option, true);
    } else if (option == "--min-flips") {
      arguments.min_flips =
          static_cast<int>(parse_unsigned(value(), option));
    } else if (option == "--max-flips") {
      arguments.max_flips =
          static_cast<int>(parse_unsigned(value(), option));
    } else if (option == "--kick-min") {
      arguments.kick_min =
          static_cast<int>(parse_unsigned(value(), option));
    } else if (option == "--kick-max") {
      arguments.kick_max =
          static_cast<int>(parse_unsigned(value(), option));
    } else if (option == "--block-rows-min") {
      arguments.block_rows_min =
          static_cast<int>(parse_unsigned(value(), option));
    } else if (option == "--block-rows-max") {
      arguments.block_rows_max =
          static_cast<int>(parse_unsigned(value(), option));
    } else if (option == "--block-columns-min") {
      arguments.block_columns_min =
          static_cast<int>(parse_unsigned(value(), option));
    } else if (option == "--block-columns-max") {
      arguments.block_columns_max =
          static_cast<int>(parse_unsigned(value(), option));
    } else if (option == "--restart-interval") {
      arguments.restart_interval = parse_unsigned(value(), option);
    } else if (option == "--cooling-period") {
      arguments.cooling_period = parse_unsigned(value(), option);
    } else if (option == "--temperature-start") {
      arguments.temperature_start =
          static_cast<long double>(
              parse_double(value(), option, true));
    } else if (option == "--temperature-end") {
      arguments.temperature_end =
          static_cast<long double>(
              parse_double(value(), option, true));
    } else if (option == "--floor-percent") {
      arguments.floor_percent =
          static_cast<int>(parse_unsigned(value(), option));
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
  const auto normalized = [](const std::filesystem::path& path) {
    return std::filesystem::absolute(path).lexically_normal();
  };
  const auto start = normalized(arguments.start);
  const auto output = normalized(arguments.output);
  const auto log = normalized(arguments.log);
  if (start == output || start == log || output == log) {
    throw std::runtime_error(
        "--start, --output, and --log must be distinct");
  }
  if (!arguments.research_output.empty()) {
    const auto research_output =
        normalized(arguments.research_output);
    if (research_output == start ||
        research_output == output ||
        research_output == log) {
      throw std::runtime_error(
          "--research-output must not alias start, output, or log");
    }
  }
  if (arguments.min_flips < 1 ||
      arguments.max_flips < arguments.min_flips ||
      arguments.max_flips > kOrder) {
    throw std::runtime_error(
        "--min-flips and --max-flips must satisfy "
        "1 <= min <= max <= 23");
  }
  if (arguments.kick_min < 1 ||
      arguments.kick_max < arguments.kick_min ||
      arguments.kick_max > kEntries) {
    throw std::runtime_error(
        "--kick-min and --kick-max must satisfy "
        "1 <= min <= max <= 529");
  }
  const bool block_mode =
      arguments.block_rows_min != 0 ||
      arguments.block_rows_max != 0 ||
      arguments.block_columns_min != 0 ||
      arguments.block_columns_max != 0;
  if (block_mode &&
      (arguments.block_rows_min < 2 ||
       arguments.block_rows_max < arguments.block_rows_min ||
       arguments.block_rows_max >= kOrder ||
       arguments.block_columns_min < 1 ||
       arguments.block_columns_max <
           arguments.block_columns_min ||
       arguments.block_columns_max > kOrder)) {
    throw std::runtime_error(
        "block ranges must satisfy 2 <= row min <= row max <= 22 "
        "and 1 <= column min <= column max <= 23");
  }
  if (arguments.restart_interval == 0 ||
      arguments.cooling_period == 0) {
    throw std::runtime_error(
        "--restart-interval and --cooling-period must be positive");
  }
  if (arguments.temperature_start < arguments.temperature_end) {
    throw std::runtime_error(
        "--temperature-start must be at least --temperature-end");
  }
  if (arguments.floor_percent < 0 ||
      arguments.floor_percent > 100) {
    throw std::runtime_error(
        "--floor-percent must be between 0 and 100");
  }
  return arguments;
}

void log_record(std::ofstream& log, const Arguments& arguments,
                const Statistics& statistics,
                const ProposalDescription& proposal, const char* event,
                double elapsed, Wide best_score, Wide state_score) {
  const double throughput =
      elapsed > 0.0
          ? static_cast<double>(statistics.proposals) / elapsed
          : 0.0;
  log << "{\"absolute_determinant\":\""
      << wide_to_string(best_score)
      << "\",\"accepted\":" << statistics.accepted
      << ",\"accepted_downhill\":"
      << statistics.accepted_downhill
      << ",\"accepted_equal\":" << statistics.accepted_equal
      << ",\"accepted_uphill\":" << statistics.accepted_uphill
      << ",\"block_columns\":" << proposal.block_columns
      << ",\"block_columns_max\":" << arguments.block_columns_max
      << ",\"block_columns_min\":" << arguments.block_columns_min
      << ",\"block_rows\":" << proposal.block_rows
      << ",\"block_rows_max\":" << arguments.block_rows_max
      << ",\"block_rows_min\":" << arguments.block_rows_min
      << ",\"column_proposals\":" << statistics.column_proposals
      << ",\"cooling_period\":" << arguments.cooling_period
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(3) << elapsed
      << ",\"event\":\"" << event
      << "\",\"exact_checks\":" << statistics.exact_checks
      << ",\"exact_cofactor_fallbacks\":"
      << statistics.exact_cofactor_fallbacks
      << ",\"first_line\":" << proposal.first
      << ",\"first_line_flips\":" << proposal.first_flips
      << ",\"floor_percent\":" << arguments.floor_percent
      << ",\"floor_rejections\":" << statistics.floor_rejections
      << ",\"kick_max\":" << arguments.kick_max
      << ",\"kick_min\":" << arguments.kick_min
      << ",\"max_flips\":" << arguments.max_flips
      << ",\"min_flips\":" << arguments.min_flips
      << ",\"optimized_line\":" << proposal.optimized
      << ",\"orientation\":\""
      << (proposal.columns ? "columns" : "rows")
      << "\",\"perturbed_entries\":"
      << statistics.perturbed_entries
      << ",\"proposals\":" << statistics.proposals
      << ",\"proposals_per_second\":" << std::setprecision(3)
      << throughput
      << ",\"restarts\":" << statistics.restarts
      << ",\"restart_interval\":" << arguments.restart_interval
      << ",\"research_determinant\":\""
      << wide_to_string(statistics.research_score) << "\""
      << ",\"row_proposals\":" << statistics.row_proposals
      << ",\"second_line\":" << proposal.second
      << ",\"second_line_flips\":" << proposal.second_flips
      << ",\"seed\":" << arguments.seed
      << ",\"singular_rejections\":"
      << statistics.singular_rejections
      << ",\"temperature_end\":"
      << static_cast<double>(arguments.temperature_end)
      << ",\"temperature_start\":"
      << static_cast<double>(arguments.temperature_start)
      << ",\"state_determinant\":\""
      << wide_to_string(state_score) << "\"}\n";
  log.flush();
  if (!log) throw std::runtime_error("cannot append research log");
}

void kick_from_best(Matrix& state, Wide& state_score,
                    const Matrix& best_matrix,
                    const Arguments& arguments,
                    std::mt19937_64& randomizer,
                    Statistics& statistics) {
  std::array<int, kEntries> coordinates{};
  for (int index = 0; index < kEntries; ++index) {
    coordinates[index] = index;
  }
  std::uniform_int_distribution<int> kick_size(
      arguments.kick_min, arguments.kick_max);
  for (int attempt = 0; attempt < 32; ++attempt) {
    state = best_matrix;
    const int flips = kick_size(randomizer);
    for (int index = 0; index < flips; ++index) {
      std::uniform_int_distribution<int> choose(index, kEntries - 1);
      const int selected = choose(randomizer);
      std::swap(coordinates[index], coordinates[selected]);
      const int encoded = coordinates[index];
      state[encoded / kOrder][encoded % kOrder] *= -1;
    }
    state_score = absolute(exact_determinant(state));
    ++statistics.exact_checks;
    if (state_score != 0) return;
  }
  throw std::runtime_error("could not produce a nonsingular basin kick");
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

    Matrix best_matrix = read_matrix(arguments.start);
    Wide best_score = absolute(exact_determinant(best_matrix));
    if (best_score == 0) {
      throw std::runtime_error("start matrix must be nonsingular");
    }
    Statistics statistics;
    statistics.exact_checks = 1;
    Matrix state_matrix = best_matrix;
    Wide state_score = best_score;
    std::uint64_t checkpoint_nonce = 0;
    atomic_write_matrix(
        arguments.output, best_matrix, checkpoint_nonce++);

    std::mt19937_64 randomizer(arguments.seed);
    std::uniform_real_distribution<long double> unit(0.0L, 1.0L);
    const auto started = Clock::now();
    const auto deadline =
        started + std::chrono::duration<double>(arguments.seconds);
    auto next_heartbeat =
        started +
        std::chrono::duration<double>(arguments.heartbeat_seconds);
    std::uint64_t last_best_proposal = 0;
    std::uint64_t last_restart_proposal = 0;
    ProposalDescription last_proposal;

    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);
    log_record(log, arguments, statistics, last_proposal, "start", 0.0,
               best_score, state_score);

    while (!stop_requested && Clock::now() < deadline) {
      ++statistics.proposals;
      ProposalDescription proposal;
      Matrix candidate = propose_triple_move(
          state_matrix, proposal, arguments, randomizer, statistics);
      last_proposal = proposal;
      if (proposal.columns) {
        ++statistics.column_proposals;
      } else {
        ++statistics.row_proposals;
      }

      Wide candidate_score = 0;
      if (candidate != state_matrix) {
        candidate_score = absolute(exact_determinant(candidate));
        ++statistics.exact_checks;
      } else {
        candidate_score = state_score;
      }
      if (!arguments.research_output.empty() &&
          candidate_score < best_score &&
          candidate_score > statistics.research_score) {
        statistics.research_score = candidate_score;
        atomic_write_matrix(
            arguments.research_output, candidate, checkpoint_nonce++);
      }

      bool accept = false;
      if (candidate_score > state_score) {
        accept = true;
        ++statistics.accepted_uphill;
      } else if (candidate_score == state_score) {
        accept = true;
        ++statistics.accepted_equal;
      } else if (candidate_score != 0) {
        // Keep basin walks from collapsing too far below the best exact state.
        const bool above_floor =
            candidate_score * static_cast<Wide>(100) >=
            best_score * static_cast<Wide>(
                arguments.floor_percent);
        if (!above_floor) {
          ++statistics.floor_rejections;
        } else {
          const long double phase =
              static_cast<long double>(
                  statistics.proposals %
                  arguments.cooling_period) /
              static_cast<long double>(
                  arguments.cooling_period);
          const long double temperature =
              (arguments.temperature_start -
               arguments.temperature_end) *
                  (1.0L - phase) +
              arguments.temperature_end;
          const long double log_ratio =
              std::log(static_cast<long double>(candidate_score)) -
              std::log(static_cast<long double>(state_score));
          accept =
              std::log(std::max(unit(randomizer), 1.0e-30L)) <
              log_ratio / temperature;
          if (accept) ++statistics.accepted_downhill;
        }
      }

      if (accept) {
        state_matrix = candidate;
        state_score = candidate_score;
        ++statistics.accepted;
        if (state_score > best_score) {
          // candidate_score is the exact Bareiss authority for promotion.
          best_score = state_score;
          best_matrix = state_matrix;
          last_best_proposal = statistics.proposals;
          atomic_write_matrix(
              arguments.output, best_matrix, checkpoint_nonce++);
          const double elapsed =
              std::chrono::duration<double>(Clock::now() - started)
                  .count();
          log_record(
              log, arguments, statistics, proposal, "new_best", elapsed,
              best_score, state_score);
          std::cout << "new best |det|="
                    << wide_to_string(best_score)
                    << " proposal=" << statistics.proposals << '\n'
                    << std::flush;
        }
      }

      if (statistics.proposals - last_best_proposal >=
              arguments.restart_interval &&
          statistics.proposals - last_restart_proposal >=
              arguments.restart_interval) {
        kick_from_best(
            state_matrix, state_score, best_matrix, arguments, randomizer,
            statistics);
        last_restart_proposal = statistics.proposals;
        ++statistics.restarts;
        const double elapsed =
            std::chrono::duration<double>(Clock::now() - started)
                .count();
        log_record(log, arguments, statistics, proposal, "restart_kick",
                   elapsed, best_score, state_score);
      }

      const auto now = Clock::now();
      if (arguments.heartbeat_seconds > 0.0 &&
          now >= next_heartbeat) {
        const double elapsed =
            std::chrono::duration<double>(now - started).count();
        log_record(log, arguments, statistics, proposal, "heartbeat",
                   elapsed, best_score, state_score);
        next_heartbeat =
            now +
            std::chrono::duration<double>(
                arguments.heartbeat_seconds);
      }
    }

    const double elapsed =
        std::chrono::duration<double>(Clock::now() - started).count();
    log_record(
        log, arguments, statistics, last_proposal,
        stop_requested ? "stopped" : "finished", elapsed, best_score,
        state_score);
    const double throughput =
        elapsed > 0.0
            ? static_cast<double>(statistics.proposals) / elapsed
            : 0.0;
    std::cout << "finished |det|=" << wide_to_string(best_score)
              << " proposals=" << statistics.proposals
              << " accepted=" << statistics.accepted
              << " proposals/s=" << std::fixed << std::setprecision(1)
              << throughput << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "triple_line_search: " << error.what() << '\n';
    return 2;
  }
}
