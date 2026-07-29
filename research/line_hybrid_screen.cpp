#include <algorithm>
#include <array>
#include <bit>
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
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <unordered_set>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

namespace {

constexpr int kOrder = 23;
using Clock = std::chrono::steady_clock;
using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using Wide = __int128_t;

volatile std::sig_atomic_t stop_requested = 0;

void request_stop(int) { stop_requested = 1; }

struct Arguments {
  std::filesystem::path first;
  std::filesystem::path second;
  std::filesystem::path output;
  std::filesystem::path research_output;
  std::filesystem::path log;
  std::filesystem::path tie_dir;
  std::string orientation = "both";
  double heartbeat_seconds = 30.0;
};

struct Statistics {
  std::uint64_t assignments = 0;
  std::uint64_t exact_checks = 0;
  std::uint64_t ties = 0;
  std::uint64_t unique_ties = 0;
  std::uint64_t promotions = 0;
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
    throw std::runtime_error("cannot open matrix: " + path.string());
  }
  Matrix matrix{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      if (!(input >> matrix[row][column]) ||
          (matrix[row][column] != -1 &&
           matrix[row][column] != 1)) {
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

std::string matrix_sign_bits_hex(const Matrix& matrix) {
  constexpr char digits[] = "0123456789abcdef";
  std::string result;
  int nibble = 0;
  int bits = 0;
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      nibble = (nibble << 1) | (matrix[row][column] == 1 ? 1 : 0);
      ++bits;
      if (bits == 4) {
        result.push_back(digits[nibble]);
        nibble = 0;
        bits = 0;
      }
    }
  }
  if (bits != 0) {
    nibble <<= 4 - bits;
    result.push_back(digits[nibble]);
  }
  return result;
}

std::string json_escape(std::string_view text) {
  constexpr char digits[] = "0123456789abcdef";
  std::string result;
  result.reserve(text.size());
  for (const unsigned char character : text) {
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
          result += "\\u00";
          result.push_back(digits[character >> 4U]);
          result.push_back(digits[character & 0x0fU]);
        } else {
          result.push_back(static_cast<char>(character));
        }
        break;
    }
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
                         const Matrix& matrix,
                         std::uint64_t nonce) {
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
        ("." + path.filename().string() + ".line-hybrid-" +
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
          "cannot create checkpoint temporary file: " +
          std::string(std::strerror(errno)));
    }
  }
  if (descriptor < 0) {
    throw std::runtime_error(
        "cannot allocate a unique checkpoint temporary file");
  }
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

double strict_double(std::string_view text, std::string_view option) {
  std::size_t consumed = 0;
  const double value = std::stod(std::string(text), &consumed);
  if (consumed != text.size() || !std::isfinite(value) ||
      value <= 0.0) {
    throw std::runtime_error(
        std::string(option) + " must be finite and positive");
  }
  return value;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    auto value = [&]() -> std::string {
      if (++index >= argc) {
        throw std::runtime_error("missing value for " + option);
      }
      return argv[index];
    };
    if (option == "--first") {
      arguments.first = value();
    } else if (option == "--second") {
      arguments.second = value();
    } else if (option == "--output") {
      arguments.output = value();
    } else if (option == "--research-output") {
      arguments.research_output = value();
    } else if (option == "--log") {
      arguments.log = value();
    } else if (option == "--tie-dir") {
      arguments.tie_dir = value();
    } else if (option == "--orientation") {
      arguments.orientation = value();
    } else if (option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds = strict_double(value(), option);
    } else if (option == "--help" || option == "-h") {
      std::cout
          << "usage: line_hybrid_screen --first MATRIX --second MATRIX "
             "--output MATRIX --research-output MATRIX --log JSONL "
             "[--orientation rows|columns|both] "
             "[--tie-dir FRESH_DIRECTORY] "
             "[--heartbeat-seconds 30]\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (arguments.first.empty() || arguments.second.empty() ||
      arguments.output.empty() || arguments.research_output.empty() ||
      arguments.log.empty()) {
    throw std::runtime_error(
        "--first, --second, --output, --research-output, and --log "
        "are required");
  }
  if (arguments.orientation != "rows" &&
      arguments.orientation != "columns" &&
      arguments.orientation != "both") {
    throw std::runtime_error(
        "--orientation must be rows, columns, or both");
  }

  const auto normalized = [](const std::filesystem::path& path) {
    std::error_code error;
    const auto canonical = std::filesystem::weakly_canonical(
        std::filesystem::absolute(path), error);
    return error ? std::filesystem::absolute(path).lexically_normal()
                 : canonical;
  };
  std::vector<std::filesystem::path> paths = {
      arguments.first,
      arguments.second,
      arguments.output,
      arguments.research_output,
      arguments.log,
  };
  if (!arguments.tie_dir.empty()) {
    paths.push_back(arguments.tie_dir);
  }
  for (auto& path : paths) path = normalized(path);
  for (std::size_t first = 0; first < paths.size(); ++first) {
    for (std::size_t second = first + 1;
         second < paths.size();
         ++second) {
      if (paths[first] == paths[second]) {
        throw std::runtime_error(
            "input, output, research output, and log paths must be distinct");
      }
    }
  }
  for (const auto& path :
       {arguments.output, arguments.research_output, arguments.log}) {
    std::error_code error;
    const auto status = std::filesystem::symlink_status(path, error);
    const bool missing =
        error == std::errc::no_such_file_or_directory;
    if (missing) error.clear();
    if (error) {
      throw std::runtime_error(
          "cannot inspect fresh output path: " + path.string());
    }
    if (!missing &&
        status.type() != std::filesystem::file_type::not_found) {
      throw std::runtime_error(
          "refusing to overwrite existing output: " + path.string());
    }
  }
  if (!arguments.tie_dir.empty()) {
    std::error_code error;
    const auto status =
        std::filesystem::symlink_status(arguments.tie_dir, error);
    const bool missing =
        error == std::errc::no_such_file_or_directory;
    if (missing) error.clear();
    if (error) {
      throw std::runtime_error(
          "cannot inspect fresh tie directory: " +
          arguments.tie_dir.string());
    }
    if (!missing &&
        status.type() != std::filesystem::file_type::not_found) {
      throw std::runtime_error(
          "refusing to reuse existing tie directory: " +
          arguments.tie_dir.string());
    }
  }
  return arguments;
}

void log_record(std::ofstream& log,
                const char* event,
                std::string_view orientation,
                int differing_lines,
                const Statistics& statistics,
                Wide best_score,
                Wide research_score,
                double elapsed_seconds,
                bool complete) {
  log << "{\"absolute_determinant\":\""
      << wide_to_string(best_score)
      << "\",\"assignments\":" << statistics.assignments
      << ",\"complete\":" << (complete ? "true" : "false")
      << ",\"differing_lines\":" << differing_lines
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(6) << elapsed_seconds
      << ",\"event\":\"" << event
      << "\",\"exact_checks\":" << statistics.exact_checks
      << ",\"orientation\":\"" << orientation
      << "\",\"promotions\":" << statistics.promotions
      << ",\"research_determinant\":\""
      << wide_to_string(research_score)
      << "\",\"ties\":" << statistics.ties
      << ",\"unique_ties\":" << statistics.unique_ties << "}\n";
  log.flush();
  if (!log) throw std::runtime_error("cannot append research log");
}

std::vector<int> differing_lines(const Matrix& first,
                                 const Matrix& second,
                                 bool rows) {
  std::vector<int> lines;
  for (int line = 0; line < kOrder; ++line) {
    bool differs = false;
    for (int entry = 0; entry < kOrder; ++entry) {
      const int left = rows ? first[line][entry] : first[entry][line];
      const int right =
          rows ? second[line][entry] : second[entry][line];
      differs = differs || left != right;
    }
    if (differs) lines.push_back(line);
  }
  return lines;
}

void set_line(Matrix& target,
              const Matrix& source,
              int line,
              bool rows) {
  for (int entry = 0; entry < kOrder; ++entry) {
    if (rows) {
      target[line][entry] = source[line][entry];
    } else {
      target[entry][line] = source[entry][line];
    }
  }
}

struct Campaign {
  const Arguments& arguments;
  const Matrix& first;
  const Matrix& second;
  std::ofstream& log;
  Clock::time_point started;
  Clock::time_point next_heartbeat;
  Matrix best_matrix;
  Matrix research_matrix;
  Wide best_score = 0;
  Wide research_score = 0;
  Statistics statistics;
  std::uint64_t checkpoint_nonce = 1;
  std::ofstream* tie_manifest = nullptr;
  std::unordered_set<std::string> tie_bits;

  void archive_tie(const Matrix& matrix,
                   std::string_view orientation,
                   std::uint64_t orientation_code,
                   std::uint64_t subset_mask) {
    if (tie_manifest == nullptr) return;
    const std::string bits = matrix_sign_bits_hex(matrix);
    if (!tie_bits.insert(bits).second) return;

    const std::uint64_t unique_index = statistics.unique_ties + 1;
    const std::string filename =
        "tie-" + std::to_string(unique_index) + "-" +
        std::string(orientation) + "-assignment-" +
        std::to_string(orientation_code + 1) + ".matrix.txt";
    atomic_write_matrix(
        arguments.tie_dir / filename, matrix, checkpoint_nonce++);
    ++statistics.unique_ties;
    *tie_manifest
        << "{\"absolute_determinant\":\""
        << wide_to_string(best_score)
        << "\",\"event\":\"unique_tie\""
        << ",\"global_assignment\":" << statistics.assignments
        << ",\"matrix\":\"" << json_escape(filename)
        << "\",\"orientation\":\"" << orientation
        << "\",\"orientation_assignment\":"
        << orientation_code + 1
        << ",\"row_major_sign_bits_hex\":\"" << bits
        << "\",\"subset_mask\":" << subset_mask << "}\n";
    tie_manifest->flush();
    if (!*tie_manifest) {
      throw std::runtime_error("cannot append tie manifest");
    }
  }

  bool screen(bool rows) {
    const std::string orientation = rows ? "rows" : "columns";
    const std::vector<int> lines =
        differing_lines(first, second, rows);
    if (lines.size() >= 64) {
      throw std::runtime_error("too many differing lines for subset mask");
    }
    const std::uint64_t total =
        std::uint64_t{1} << lines.size();
    Matrix current = first;
    std::uint64_t previous_gray = 0;
    log_record(
        log, "orientation_start", orientation,
        static_cast<int>(lines.size()), statistics, best_score,
        research_score,
        std::chrono::duration<double>(Clock::now() - started).count(),
        false);

    for (std::uint64_t code = 0;
         code < total && !stop_requested;
         ++code) {
      const std::uint64_t gray = code ^ (code >> 1U);
      if (code != 0) {
        const std::uint64_t changed = gray ^ previous_gray;
        const int bit =
            static_cast<int>(std::countr_zero(changed));
        const bool use_second =
            ((gray >> bit) & std::uint64_t{1}) != 0;
        set_line(
            current,
            use_second ? second : first,
            lines[static_cast<std::size_t>(bit)],
            rows);
      }
      previous_gray = gray;

      const Wide score = absolute(exact_determinant(current));
      ++statistics.assignments;
      ++statistics.exact_checks;
      if (score > best_score) {
        if (best_score > research_score) {
          research_score = best_score;
          research_matrix = best_matrix;
          atomic_write_matrix(
              arguments.research_output,
              research_matrix,
              checkpoint_nonce++);
        }
        best_score = score;
        best_matrix = current;
        ++statistics.promotions;
        atomic_write_matrix(
            arguments.output, best_matrix, checkpoint_nonce++);
        log_record(
            log, "new_best", orientation,
            static_cast<int>(lines.size()), statistics, best_score,
            research_score,
            std::chrono::duration<double>(
                Clock::now() - started).count(),
            false);
        std::cout << "new best |det|=" << wide_to_string(best_score)
                  << " orientation=" << orientation
                  << " assignment=" << statistics.assignments << '\n'
                  << std::flush;
      } else if (score == best_score) {
        ++statistics.ties;
        archive_tie(current, orientation, code, gray);
      } else if (score > research_score) {
        research_score = score;
        research_matrix = current;
        atomic_write_matrix(
            arguments.research_output,
            research_matrix,
            checkpoint_nonce++);
      }

      const auto now = Clock::now();
      if (now >= next_heartbeat) {
        log_record(
            log, "heartbeat", orientation,
            static_cast<int>(lines.size()), statistics, best_score,
            research_score,
            std::chrono::duration<double>(now - started).count(),
            false);
        next_heartbeat =
            now + std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(
                          arguments.heartbeat_seconds));
      }
    }
    const bool complete = !stop_requested;
    log_record(
        log, "orientation_finished", orientation,
        static_cast<int>(lines.size()), statistics, best_score,
        research_score,
        std::chrono::duration<double>(Clock::now() - started).count(),
        complete);
    return complete;
  }
};

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);

    const Matrix first = read_matrix(arguments.first);
    const Matrix second = read_matrix(arguments.second);
    const Wide first_score = absolute(exact_determinant(first));
    const Wide second_score = absolute(exact_determinant(second));

    if (!arguments.log.parent_path().empty()) {
      std::filesystem::create_directories(arguments.log.parent_path());
    }
    std::ofstream log(arguments.log, std::ios::out | std::ios::trunc);
    if (!log) {
      throw std::runtime_error(
          "cannot create research log: " + arguments.log.string());
    }
    const auto started = Clock::now();
    const Matrix& initial_best =
        first_score >= second_score ? first : second;
    Campaign campaign{
        arguments,
        first,
        second,
        log,
        started,
        started + std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(
                          arguments.heartbeat_seconds)),
        initial_best,
        Matrix{},
        std::max(first_score, second_score),
        0,
        Statistics{},
        1,
        nullptr,
        {},
    };
    atomic_write_matrix(arguments.output, campaign.best_matrix, 0);

    std::ofstream tie_manifest;
    if (!arguments.tie_dir.empty()) {
      std::error_code error;
      const bool created =
          std::filesystem::create_directories(arguments.tie_dir, error);
      if (error || !created) {
        throw std::runtime_error(
            "cannot create fresh tie directory: " +
            arguments.tie_dir.string());
      }
      tie_manifest.open(
          arguments.tie_dir / "manifest.jsonl",
          std::ios::out | std::ios::trunc);
      if (!tie_manifest) {
        throw std::runtime_error("cannot create tie manifest");
      }
      tie_manifest
          << "{\"dedupe\":\"row_major_sign_bits_hex\""
          << ",\"event\":\"start\",\"first\":\""
          << json_escape(arguments.first.string())
          << "\",\"frontier_absolute_determinant\":\""
          << wide_to_string(campaign.best_score)
          << "\",\"second\":\""
          << json_escape(arguments.second.string()) << "\"}\n";
      tie_manifest.flush();
      if (!tie_manifest) {
        throw std::runtime_error("cannot initialize tie manifest");
      }
      campaign.tie_manifest = &tie_manifest;
    }

    log << "{\"event\":\"start\",\"first\":\""
        << json_escape(arguments.first.string())
        << "\",\"first_absolute_determinant\":\""
        << wide_to_string(first_score)
        << "\",\"first_row_major_sign_bits_hex\":\""
        << matrix_sign_bits_hex(first)
        << "\",\"orientation\":\""
        << json_escape(arguments.orientation)
        << "\",\"second\":\""
        << json_escape(arguments.second.string())
        << "\",\"second_absolute_determinant\":\""
        << wide_to_string(second_score)
        << "\",\"second_row_major_sign_bits_hex\":\""
        << matrix_sign_bits_hex(second) << "\",\"tie_dir\":";
    if (arguments.tie_dir.empty()) {
      log << "null";
    } else {
      log << "\"" << json_escape(arguments.tie_dir.string()) << "\"";
    }
    log << "}\n";
    log.flush();
    if (!log) throw std::runtime_error("cannot append start record");

    bool complete = true;
    if (arguments.orientation == "rows" ||
        arguments.orientation == "both") {
      complete = campaign.screen(true);
    }
    if (complete &&
        (arguments.orientation == "columns" ||
         arguments.orientation == "both")) {
      complete = campaign.screen(false);
    }

    const double elapsed =
        std::chrono::duration<double>(Clock::now() - started).count();
    log_record(
        log, complete ? "finished" : "stopped", arguments.orientation,
        -1, campaign.statistics, campaign.best_score,
        campaign.research_score, elapsed, complete);
    if (campaign.tie_manifest != nullptr) {
      tie_manifest
          << "{\"assignments\":" << campaign.statistics.assignments
          << ",\"complete\":" << (complete ? "true" : "false")
          << ",\"event\":\"" << (complete ? "finished" : "stopped")
          << "\",\"ties\":" << campaign.statistics.ties
          << ",\"unique_ties\":"
          << campaign.statistics.unique_ties << "}\n";
      tie_manifest.flush();
      if (!tie_manifest) {
        throw std::runtime_error("cannot finalize tie manifest");
      }
    }
    std::cout << (complete ? "finished" : "stopped")
              << " |det|=" << wide_to_string(campaign.best_score)
              << " research="
              << wide_to_string(campaign.research_score)
              << " assignments=" << campaign.statistics.assignments
              << " exact_checks=" << campaign.statistics.exact_checks
              << " ties=" << campaign.statistics.ties
              << " unique_ties="
              << campaign.statistics.unique_ties << '\n';
    return complete ? 0 : 3;
  } catch (const std::exception& error) {
    std::cerr << "line_hybrid_screen: " << error.what() << '\n';
    return 2;
  }
}
