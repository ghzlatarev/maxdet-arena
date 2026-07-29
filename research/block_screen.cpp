#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <csignal>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

#include <fcntl.h>
#include <unistd.h>

namespace {

constexpr int kOrder = 23;
using Wide = __int128_t;
using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using Clock = std::chrono::steady_clock;

volatile std::sig_atomic_t stop_requested = 0;
void request_stop(int) { stop_requested = 1; }

struct Arguments {
  std::filesystem::path start;
  std::filesystem::path output;
  std::filesystem::path research_output;
  std::filesystem::path log;
  int min_size = 2;
  int max_size = 5;
  int line_count = 2;
  std::string orientation = "both";
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

Wide exact_determinant(const Matrix& input) {
  std::array<std::array<Wide, kOrder>, kOrder> matrix{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      matrix[row][column] = input[row][column];
    }
  }
  Wide previous_pivot = 1;
  int sign = 1;
  for (int column = 0; column < kOrder - 1; ++column) {
    int pivot_row = column;
    while (pivot_row < kOrder && matrix[pivot_row][column] == 0) {
      ++pivot_row;
    }
    if (pivot_row == kOrder) return 0;
    if (pivot_row != column) {
      std::swap(matrix[pivot_row], matrix[column]);
      sign = -sign;
    }
    const Wide pivot = matrix[column][column];
    for (int row = column + 1; row < kOrder; ++row) {
      for (int inner = column + 1; inner < kOrder; ++inner) {
        const Wide numerator =
            matrix[row][inner] * pivot -
            matrix[row][column] * matrix[column][inner];
        if (column != 0 && numerator % previous_pivot != 0) {
          throw std::runtime_error("exact Bareiss division failed");
        }
        matrix[row][inner] =
            column == 0 ? numerator : numerator / previous_pivot;
      }
      matrix[row][column] = 0;
    }
    previous_pivot = pivot;
  }
  return static_cast<Wide>(sign) *
         matrix[kOrder - 1][kOrder - 1];
}

Matrix read_matrix(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("cannot open start matrix: " + path.string());
  }
  Matrix matrix{};
  for (auto& row : matrix) {
    for (int& value : row) {
      if (!(input >> value) || (value != -1 && value != 1)) {
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
  const std::filesystem::path temporary =
      directory /
      ("." + path.filename().string() + ".block2-" +
       std::to_string(static_cast<long long>(::getpid())) + "-" +
       std::to_string(nonce) + ".tmp");
  int flags = O_WRONLY | O_CREAT | O_EXCL;
#ifdef O_CLOEXEC
  flags |= O_CLOEXEC;
#endif
  int descriptor = ::open(temporary.c_str(), flags, 0600);
  if (descriptor < 0) {
    throw std::runtime_error(
        "cannot create checkpoint: " + std::string(std::strerror(errno)));
  }
  bool renamed = false;
  try {
    write_all(descriptor, matrix_bytes(matrix));
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

int parse_integer(std::string_view text, std::string_view option) {
  std::size_t consumed = 0;
  const int result = std::stoi(std::string(text), &consumed);
  if (consumed != text.size()) {
    throw std::runtime_error(std::string(option) + " must be an integer");
  }
  return result;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    auto value = [&]() -> std::string_view {
      if (++index >= argc) {
        throw std::runtime_error("missing value for " + option);
      }
      return argv[index];
    };
    if (option == "--start") arguments.start = value();
    else if (option == "--output") arguments.output = value();
    else if (option == "--research-output") {
      arguments.research_output = value();
    } else if (option == "--log") arguments.log = value();
    else if (option == "--min-size") {
      arguments.min_size = parse_integer(value(), option);
    } else if (option == "--max-size") {
      arguments.max_size = parse_integer(value(), option);
    } else if (option == "--line-count") {
      arguments.line_count = parse_integer(value(), option);
    } else if (option == "--orientation") {
      arguments.orientation = value();
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (arguments.start.empty() || arguments.output.empty() ||
      arguments.log.empty()) {
    throw std::runtime_error("--start, --output, and --log are required");
  }
  if (arguments.min_size < 1 ||
      arguments.max_size < arguments.min_size ||
      arguments.max_size > kOrder) {
    throw std::runtime_error(
        "block sizes must satisfy 1 <= min <= max <= 23");
  }
  if (arguments.line_count != 2 && arguments.line_count != 3) {
    throw std::runtime_error("--line-count must be 2 or 3");
  }
  if (arguments.orientation != "rows" &&
      arguments.orientation != "columns" &&
      arguments.orientation != "both") {
    throw std::runtime_error(
        "--orientation must be rows, columns, or both");
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
    const auto research = normalized(arguments.research_output);
    if (research == start || research == output || research == log) {
      throw std::runtime_error(
          "--research-output must not alias another path");
    }
  }
  return arguments;
}

void append_log(std::ofstream& log, const char* event,
                const std::string& orientation, int first, int second,
                int third,
                std::uint64_t assignments, std::uint64_t exact_checks,
                Wide best_score, Wide research_score, double elapsed) {
  log << "{\"absolute_determinant\":\""
      << wide_to_string(best_score)
      << "\",\"assignments\":" << assignments
      << ",\"elapsed_seconds\":" << std::fixed << std::setprecision(3)
      << elapsed
      << ",\"event\":\"" << event
      << "\",\"exact_checks\":" << exact_checks
      << ",\"first\":" << first
      << ",\"orientation\":\"" << orientation
      << "\",\"research_determinant\":\""
      << wide_to_string(research_score)
      << "\",\"second\":" << second
      << ",\"third\":" << third << "}\n";
  log.flush();
  if (!log) throw std::runtime_error("cannot append research log");
}

struct Campaign {
  const Arguments& arguments;
  const Matrix& original;
  std::ofstream& log;
  Wide best_score;
  Wide research_score = 0;
  Matrix best_matrix{};
  std::uint64_t assignments = 0;
  std::uint64_t exact_checks = 0;
  std::uint64_t checkpoint_nonce = 0;
  Clock::time_point started;

  void screen_orientation(bool columns) {
    const std::string label = columns ? "columns" : "rows";
    const Matrix base = columns ? transpose(original) : original;
    const Wide determinant = exact_determinant(base);
    ++exact_checks;
    for (int first = 0; first < kOrder - 1 && !stop_requested; ++first) {
      for (int second = first + 1;
           second < kOrder && !stop_requested; ++second) {
        std::array<Wide, kOrder> first_cofactors{};
        std::array<Wide, kOrder> second_cofactors{};
        for (int column = 0; column < kOrder; ++column) {
          Matrix basis = base;
          basis[first].fill(0);
          basis[first][column] = 1;
          first_cofactors[column] = exact_determinant(basis);
          basis = base;
          basis[second].fill(0);
          basis[second][column] = 1;
          second_cofactors[column] = exact_determinant(basis);
          exact_checks += 2;
        }

        std::array<std::array<Wide, kOrder>, kOrder> weight{};
        Wide reconstructed = 0;
        for (int left = 0; left < kOrder; ++left) {
          for (int right = left + 1; right < kOrder; ++right) {
            const Wide forward_numerator =
                first_cofactors[left] * second_cofactors[right] -
                second_cofactors[left] * first_cofactors[right];
            const Wide reverse_numerator =
                first_cofactors[right] * second_cofactors[left] -
                second_cofactors[right] * first_cofactors[left];
            if (forward_numerator % determinant != 0 ||
                reverse_numerator % determinant != 0) {
              throw std::runtime_error(
                  "second-cofactor division failed");
            }
            const Wide forward =
                forward_numerator / determinant;
            const Wide reverse =
                reverse_numerator / determinant;
            weight[left][right] =
                static_cast<Wide>(base[first][left]) *
                    static_cast<Wide>(base[second][right]) * forward +
                static_cast<Wide>(base[first][right]) *
                    static_cast<Wide>(base[second][left]) * reverse;
            reconstructed += weight[left][right];
          }
        }
        if (reconstructed != determinant) {
          throw std::runtime_error(
              "pair coefficient reconstruction disagrees with determinant");
        }

        std::array<int, kOrder> selection{};
        std::array<bool, kOrder> selected{};
        const auto enumerate =
            [&](auto&& self, int size, int depth, int next) -> void {
          if (stop_requested) return;
          if (depth != size) {
            for (int value = next;
                 value <= kOrder - (size - depth); ++value) {
              selection[depth] = value;
              selected[value] = true;
              self(self, size, depth + 1, value + 1);
              selected[value] = false;
              if (stop_requested) return;
            }
            return;
          }
          Wide cut = 0;
          for (int index = 0; index < size; ++index) {
            const int inside = selection[index];
            for (int outside = 0; outside < kOrder; ++outside) {
              if (selected[outside]) continue;
              cut += inside < outside
                         ? weight[inside][outside]
                         : weight[outside][inside];
            }
          }
          const Wide candidate_determinant = determinant - 2 * cut;
          const Wide candidate_score = absolute(candidate_determinant);
          ++assignments;
          if (candidate_score <= research_score ||
              candidate_score == best_score) {
            return;
          }
          Matrix candidate = base;
          for (int index = 0; index < size; ++index) {
            const int column = selection[index];
            candidate[first][column] *= -1;
            candidate[second][column] *= -1;
          }
          const Wide checked = absolute(exact_determinant(candidate));
          ++exact_checks;
          if (checked != candidate_score) {
            throw std::runtime_error(
                "block formula disagrees with exact Bareiss determinant");
          }
          const Matrix unoriented =
              columns ? transpose(candidate) : candidate;
          if (candidate_score > best_score) {
            best_score = candidate_score;
            best_matrix = unoriented;
            atomic_write_matrix(
                arguments.output, best_matrix, checkpoint_nonce++);
            append_log(
                log, "new_best", label, first, second, -1, assignments,
                exact_checks, best_score, research_score,
                std::chrono::duration<double>(Clock::now() - started)
                    .count());
            std::cout << "new best |det|="
                      << wide_to_string(best_score)
                      << " orientation=" << label << " pair=" << first
                      << ',' << second << " size=" << size << '\n'
                      << std::flush;
          } else if (candidate_score > research_score) {
            research_score = candidate_score;
            if (!arguments.research_output.empty()) {
              atomic_write_matrix(
                  arguments.research_output, unoriented,
                  checkpoint_nonce++);
            }
          }
        };
        for (int size = arguments.min_size;
             size <= arguments.max_size && !stop_requested; ++size) {
          enumerate(enumerate, size, 0, 0);
        }
        append_log(
            log, stop_requested ? "pair_interrupted" : "pair_finished",
            label, first, second, -1, assignments,
            exact_checks, best_score, research_score,
            std::chrono::duration<double>(Clock::now() - started).count());
      }
    }
  }

  void screen_orientation_three(bool columns) {
    const std::string label = columns ? "columns" : "rows";
    const Matrix base = columns ? transpose(original) : original;
    for (int first = 0; first < kOrder - 2 && !stop_requested; ++first) {
      for (int second = first + 1;
           second < kOrder - 1 && !stop_requested; ++second) {
        for (int third = second + 1;
             third < kOrder && !stop_requested; ++third) {
          std::array<int, kOrder> selection{};
          const auto enumerate =
              [&](auto&& self, int size, int depth, int next) -> void {
            if (stop_requested) return;
            if (depth != size) {
              for (int value = next;
                   value <= kOrder - (size - depth); ++value) {
                selection[depth] = value;
                self(self, size, depth + 1, value + 1);
                if (stop_requested) return;
              }
              return;
            }
            Matrix candidate = base;
            for (int index = 0; index < size; ++index) {
              const int column = selection[index];
              candidate[first][column] *= -1;
              candidate[second][column] *= -1;
              candidate[third][column] *= -1;
            }
            const Wide candidate_score =
                absolute(exact_determinant(candidate));
            ++assignments;
            ++exact_checks;
            if (candidate_score <= research_score ||
                candidate_score == best_score) {
              return;
            }
            const Matrix unoriented =
                columns ? transpose(candidate) : candidate;
            if (candidate_score > best_score) {
              best_score = candidate_score;
              best_matrix = unoriented;
              atomic_write_matrix(
                  arguments.output, best_matrix, checkpoint_nonce++);
              append_log(
                  log, "new_best", label, first, second, third,
                  assignments, exact_checks, best_score, research_score,
                  std::chrono::duration<double>(
                      Clock::now() - started).count());
              std::cout << "new best |det|="
                        << wide_to_string(best_score)
                        << " orientation=" << label << " triple="
                        << first << ',' << second << ',' << third
                        << " size=" << size << '\n' << std::flush;
            } else if (candidate_score > research_score) {
              research_score = candidate_score;
              if (!arguments.research_output.empty()) {
                atomic_write_matrix(
                    arguments.research_output, unoriented,
                    checkpoint_nonce++);
              }
            }
          };
          for (int size = arguments.min_size;
               size <= arguments.max_size && !stop_requested; ++size) {
            enumerate(enumerate, size, 0, 0);
          }
          append_log(
              log,
              stop_requested ? "triple_interrupted" : "triple_finished",
              label, first, second, third,
              assignments, exact_checks, best_score, research_score,
              std::chrono::duration<double>(
                  Clock::now() - started).count());
        }
      }
    }
  }
};

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    if (!arguments.log.parent_path().empty()) {
      std::filesystem::create_directories(arguments.log.parent_path());
    }
    std::ofstream log(arguments.log, std::ios::app);
    if (!log) throw std::runtime_error("cannot open research log");
    const Matrix start = read_matrix(arguments.start);
    Campaign campaign{
        arguments,
        start,
        log,
        absolute(exact_determinant(start)),
        0,
        start,
        0,
        1,
        0,
        Clock::now(),
    };
    atomic_write_matrix(arguments.output, start, campaign.checkpoint_nonce++);
    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);
    append_log(
        log, "start", arguments.orientation, -1, -1, -1, 0,
        campaign.exact_checks, campaign.best_score, 0, 0.0);
    if (arguments.orientation == "rows" ||
        arguments.orientation == "both") {
      if (arguments.line_count == 2) {
        campaign.screen_orientation(false);
      } else {
        campaign.screen_orientation_three(false);
      }
    }
    if (!stop_requested &&
        (arguments.orientation == "columns" ||
         arguments.orientation == "both")) {
      if (arguments.line_count == 2) {
        campaign.screen_orientation(true);
      } else {
        campaign.screen_orientation_three(true);
      }
    }
    const double elapsed =
        std::chrono::duration<double>(Clock::now() - campaign.started)
            .count();
    append_log(
        log, stop_requested ? "stopped" : "finished",
        arguments.orientation, -1, -1, -1, campaign.assignments,
        campaign.exact_checks, campaign.best_score,
        campaign.research_score, elapsed);
    std::cout << (stop_requested ? "stopped" : "finished")
              << " |det|=" << wide_to_string(campaign.best_score)
              << " research=" << wide_to_string(campaign.research_score)
              << " assignments=" << campaign.assignments
              << " exact_checks=" << campaign.exact_checks
              << " elapsed=" << std::fixed << std::setprecision(3)
              << elapsed << "s\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "block_screen: " << error.what() << '\n';
    return 2;
  }
}
