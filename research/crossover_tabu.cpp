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

namespace {

constexpr int kOrder = 23;
constexpr int kEntries = kOrder * kOrder;
constexpr int kMinimumTenure = 7;
constexpr int kMaximumTenure = 128;
constexpr std::size_t kVisitTableSize = 1U << 18U;
constexpr long double kSingularRatio = 1.0e-18L;
constexpr long double kTieTolerance = 1.0e-18L;
constexpr long double kAspirationSlack = 1.0e-9L;
constexpr std::string_view kMethodVersion = "entrywise-crossover-tabu-v1";
constexpr std::string_view kExactMethodVersion =
    "entrywise-recombinant-exhaustive-v1";

using Clock = std::chrono::steady_clock;
using Wide = __int128_t;
using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using Inverse = std::array<std::array<long double, kOrder>, kOrder>;

volatile std::sig_atomic_t stop_requested = 0;

void request_stop(int) { stop_requested = 1; }

struct Coordinate {
  int row = 0;
  int column = 0;
};

struct State {
  Matrix matrix{};
  Inverse inverse{};
  std::vector<unsigned char> endpoint_bits;
  long double log_abs_determinant =
      -std::numeric_limits<long double>::infinity();
  bool nonsingular = false;
  std::uint64_t moves_since_rebuild = 0;
};

struct Arguments {
  std::filesystem::path first;
  std::filesystem::path second;
  std::filesystem::path output;
  std::filesystem::path log;
  std::string mode = "tabu";
  std::uint64_t seed = 23;
  double seconds = 3600.0;
  double heartbeat_seconds = 30.0;
  std::uint64_t maximum_iterations = 0;
  std::uint64_t restart_interval = 4096;
  std::uint64_t rebuild_interval = 64;
  std::uint64_t exact_interval = 64;
  std::uint64_t explore_period = 17;
  std::uint64_t kick_size = 12;
  double temperature = 0.04;
  double maximum_log_drop = 0.20;
  double exact_log_window = 1.0e-9;
};

struct Move {
  std::size_t difference_index = 0;
  long double projected_log_abs_determinant =
      -std::numeric_limits<long double>::infinity();
  bool aspiration = false;
  bool exploratory = false;
};

struct Candidate {
  std::size_t difference_index = 0;
  long double projected_log_abs_determinant =
      -std::numeric_limits<long double>::infinity();
  bool aspiration = false;
};

struct Visit {
  std::uint64_t hash = 0;
  std::uint64_t iteration = 0;
  std::uint64_t epoch = 0;
};

struct Statistics {
  std::uint64_t iterations = 0;
  std::uint64_t moves = 0;
  std::uint64_t inverse_rebuilds = 0;
  std::uint64_t exact_checks = 0;
  std::uint64_t promotions = 0;
  std::uint64_t cycles = 0;
  std::uint64_t singular_moves_rejected = 0;
  std::uint64_t exploratory_moves = 0;
  std::uint64_t downhill_moves = 0;
  std::uint64_t aspiration_probes = 0;
  std::uint64_t aspiration_moves = 0;
  std::uint64_t restarts = 0;
  std::uint64_t first_endpoint_restarts = 0;
  std::uint64_t second_endpoint_restarts = 0;
  std::uint64_t random_subset_restarts = 0;
  std::uint64_t kicks = 0;
  std::uint64_t restart_singular_retries = 0;
  std::uint64_t exhaustive_assignments = 0;
  std::uint64_t exhaustive_total = 0;
  std::uint64_t exhaustive_ties = 0;
};

std::string_view method_name(const Arguments& arguments) {
  return arguments.mode == "exact" ? kExactMethodVersion
                                    : kMethodVersion;
}

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
  constexpr char digits[] = "0123456789abcdef";
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
          result += "\\u00";
          result.push_back(digits[character >> 4U]);
          result.push_back(digits[character & 0x0fU]);
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

std::string read_file_bytes(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("cannot open input file: " + path.string());
  }
  std::string bytes{
      std::istreambuf_iterator<char>(input),
      std::istreambuf_iterator<char>()};
  if (input.bad()) {
    throw std::runtime_error("cannot read input file: " + path.string());
  }
  return bytes;
}

Matrix parse_matrix(std::string_view bytes,
                    const std::filesystem::path& path) {
  std::istringstream input{std::string(bytes)};
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

class Sha256 {
 public:
  void update(std::string_view bytes) {
    for (const unsigned char byte : bytes) {
      buffer_[buffer_size_++] = byte;
      bit_count_ += 8;
      if (buffer_size_ == buffer_.size()) {
        transform();
        buffer_size_ = 0;
      }
    }
  }

  std::string finish() {
    buffer_[buffer_size_++] = 0x80U;
    if (buffer_size_ > 56) {
      while (buffer_size_ < buffer_.size()) buffer_[buffer_size_++] = 0;
      transform();
      buffer_size_ = 0;
    }
    while (buffer_size_ < 56) buffer_[buffer_size_++] = 0;
    for (int shift = 56; shift >= 0; shift -= 8) {
      buffer_[buffer_size_++] =
          static_cast<unsigned char>(bit_count_ >> shift);
    }
    transform();
    std::ostringstream output;
    output << std::hex << std::setfill('0');
    for (const std::uint32_t word : state_) {
      output << std::setw(8) << word;
    }
    return output.str();
  }

 private:
  static std::uint32_t rotate_right(std::uint32_t value, int shift) {
    return (value >> shift) | (value << (32 - shift));
  }

  void transform() {
    static constexpr std::array<std::uint32_t, 64> constants = {
        0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U,
        0x3956c25bU, 0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U,
        0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
        0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U,
        0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU,
        0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
        0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U,
        0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U,
        0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
        0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
        0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U,
        0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
        0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U,
        0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
        0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
        0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U,
    };
    std::array<std::uint32_t, 64> words{};
    for (int index = 0; index < 16; ++index) {
      words[index] =
          (static_cast<std::uint32_t>(buffer_[4 * index]) << 24U) |
          (static_cast<std::uint32_t>(buffer_[4 * index + 1]) << 16U) |
          (static_cast<std::uint32_t>(buffer_[4 * index + 2]) << 8U) |
          static_cast<std::uint32_t>(buffer_[4 * index + 3]);
    }
    for (int index = 16; index < 64; ++index) {
      const std::uint32_t first =
          rotate_right(words[index - 15], 7) ^
          rotate_right(words[index - 15], 18) ^
          (words[index - 15] >> 3U);
      const std::uint32_t second =
          rotate_right(words[index - 2], 17) ^
          rotate_right(words[index - 2], 19) ^
          (words[index - 2] >> 10U);
      words[index] =
          words[index - 16] + first + words[index - 7] + second;
    }
    std::uint32_t a = state_[0];
    std::uint32_t b = state_[1];
    std::uint32_t c = state_[2];
    std::uint32_t d = state_[3];
    std::uint32_t e = state_[4];
    std::uint32_t f = state_[5];
    std::uint32_t g = state_[6];
    std::uint32_t h = state_[7];
    for (int index = 0; index < 64; ++index) {
      const std::uint32_t sum_one =
          rotate_right(e, 6) ^ rotate_right(e, 11) ^
          rotate_right(e, 25);
      const std::uint32_t choice = (e & f) ^ ((~e) & g);
      const std::uint32_t temporary_one =
          h + sum_one + choice + constants[index] + words[index];
      const std::uint32_t sum_zero =
          rotate_right(a, 2) ^ rotate_right(a, 13) ^
          rotate_right(a, 22);
      const std::uint32_t majority =
          (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t temporary_two = sum_zero + majority;
      h = g;
      g = f;
      f = e;
      e = d + temporary_one;
      d = c;
      c = b;
      b = a;
      a = temporary_one + temporary_two;
    }
    state_[0] += a;
    state_[1] += b;
    state_[2] += c;
    state_[3] += d;
    state_[4] += e;
    state_[5] += f;
    state_[6] += g;
    state_[7] += h;
  }

  std::array<std::uint32_t, 8> state_ = {
      0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
      0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U,
  };
  std::array<unsigned char, 64> buffer_{};
  std::size_t buffer_size_ = 0;
  std::uint64_t bit_count_ = 0;
};

std::string sha256(std::string_view bytes) {
  Sha256 hash;
  hash.update(bytes);
  return hash.finish();
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
        ("." + path.filename().string() + ".crossover-tabu-" +
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

std::filesystem::path resolved_target(
    const std::filesystem::path& path) {
  const std::filesystem::path parent =
      path.parent_path().empty() ? std::filesystem::path(".")
                                 : path.parent_path();
  return std::filesystem::weakly_canonical(parent) / path.filename();
}

void ensure_fresh_outputs(const Arguments& arguments) {
  if (std::filesystem::exists(arguments.output)) {
    throw std::runtime_error(
        "refusing to overwrite existing output: " +
        arguments.output.string());
  }
  if (std::filesystem::exists(arguments.log)) {
    throw std::runtime_error(
        "refusing to overwrite existing log: " + arguments.log.string());
  }
  const std::filesystem::path first =
      std::filesystem::canonical(arguments.first);
  const std::filesystem::path second =
      std::filesystem::canonical(arguments.second);
  const std::filesystem::path output =
      resolved_target(arguments.output);
  const std::filesystem::path log = resolved_target(arguments.log);
  if (output == first || output == second ||
      log == first || log == second) {
    throw std::runtime_error(
        "output and log paths must not alias either input");
  }
  if (output == log) {
    throw std::runtime_error("output and log paths must be different");
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
        augmented[row][inner] -=
            factor * augmented[column][inner];
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

long double flip_ratio(const State& state,
                       const Coordinate& coordinate) {
  const long double delta =
      -2.0L *
      static_cast<long double>(
          state.matrix[coordinate.row][coordinate.column]);
  return 1.0L +
         delta * state.inverse[coordinate.column][coordinate.row];
}

bool apply_flip(State& state, const Coordinate& coordinate,
                bool& rebuilt_inverse) {
  rebuilt_inverse = false;
  const long double delta =
      -2.0L *
      static_cast<long double>(
          state.matrix[coordinate.row][coordinate.column]);
  const long double ratio =
      1.0L +
      delta * state.inverse[coordinate.column][coordinate.row];
  if (std::fabs(ratio) < kSingularRatio) {
    state.matrix[coordinate.row][coordinate.column] *= -1;
    rebuilt_inverse = true;
    return rebuild_inverse(state);
  }

  std::array<long double, kOrder> inverse_column{};
  std::array<long double, kOrder> inverse_row{};
  for (int index = 0; index < kOrder; ++index) {
    inverse_column[index] =
        state.inverse[index][coordinate.row];
    inverse_row[index] =
        state.inverse[coordinate.column][index];
  }
  const long double factor = delta / ratio;
  for (int inner_row = 0; inner_row < kOrder; ++inner_row) {
    for (int inner_column = 0; inner_column < kOrder; ++inner_column) {
      state.inverse[inner_row][inner_column] -=
          factor * inverse_column[inner_row] *
          inverse_row[inner_column];
    }
  }
  state.matrix[coordinate.row][coordinate.column] *= -1;
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

std::vector<std::uint64_t> make_zobrist(
    std::size_t dimensions, std::uint64_t seed) {
  std::vector<std::uint64_t> values(dimensions);
  std::uint64_t state = seed ^ 0x43726f7373546162ULL;
  for (std::uint64_t& value : values) {
    state = splitmix64(state);
    value = state;
  }
  return values;
}

std::uint64_t subset_hash(
    const std::vector<unsigned char>& bits,
    const std::vector<std::uint64_t>& zobrist) {
  std::uint64_t hash = 0;
  for (std::size_t index = 0; index < bits.size(); ++index) {
    if (bits[index] != 0) hash ^= zobrist[index];
  }
  return hash;
}

std::uint64_t strict_unsigned(std::string_view text,
                              std::string_view option) {
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
    if (option == "--first") {
      arguments.first = value();
    } else if (option == "--second") {
      arguments.second = value();
    } else if (option == "--output") {
      arguments.output = value();
    } else if (option == "--log") {
      arguments.log = value();
    } else if (option == "--mode") {
      arguments.mode = value();
    } else if (option == "--seed") {
      arguments.seed = strict_unsigned(value(), option);
    } else if (option == "--seconds") {
      arguments.seconds = strict_double(value(), option, false);
    } else if (option == "--heartbeat" ||
               option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds =
          strict_double(value(), option, true);
    } else if (option == "--iterations") {
      arguments.maximum_iterations =
          strict_unsigned(value(), option);
    } else if (option == "--restart-interval") {
      arguments.restart_interval =
          strict_unsigned(value(), option);
    } else if (option == "--rebuild-interval") {
      arguments.rebuild_interval =
          strict_unsigned(value(), option);
    } else if (option == "--exact-interval") {
      arguments.exact_interval =
          strict_unsigned(value(), option);
    } else if (option == "--explore-period") {
      arguments.explore_period =
          strict_unsigned(value(), option);
    } else if (option == "--kick-size") {
      arguments.kick_size =
          strict_unsigned(value(), option);
    } else if (option == "--temperature") {
      arguments.temperature =
          strict_double(value(), option, false);
    } else if (option == "--max-log-drop") {
      arguments.maximum_log_drop =
          strict_double(value(), option, true);
    } else if (option == "--exact-log-window") {
      arguments.exact_log_window =
          strict_double(value(), option, true);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (arguments.first.empty()) {
    throw std::runtime_error("--first is required");
  }
  if (arguments.second.empty()) {
    throw std::runtime_error("--second is required");
  }
  if (arguments.output.empty()) {
    throw std::runtime_error("--output is required");
  }
  if (arguments.log.empty()) {
    throw std::runtime_error("--log is required");
  }
  if (arguments.mode != "tabu" && arguments.mode != "exact") {
    throw std::runtime_error("--mode must be tabu or exact");
  }
  if (arguments.rebuild_interval == 0) {
    throw std::runtime_error("--rebuild-interval must be positive");
  }
  if (arguments.exact_interval == 0) {
    throw std::runtime_error("--exact-interval must be positive");
  }
  return arguments;
}

std::size_t subset_weight(const State& state) {
  return static_cast<std::size_t>(
      std::count(
          state.endpoint_bits.begin(), state.endpoint_bits.end(),
          static_cast<unsigned char>(1)));
}

void log_record(
    std::ofstream& log, const Arguments& arguments,
    const Statistics& statistics, const State& state,
    std::size_t dimensions, const char* event, const char* reason,
    double elapsed_seconds, Wide best_score, Wide last_exact_score,
    int tenure) {
  const std::size_t weight = subset_weight(state);
  log << "{\"allowed_hamming_dimension\":" << dimensions
      << ",\"aspiration_moves\":" << statistics.aspiration_moves
      << ",\"aspiration_probes\":" << statistics.aspiration_probes
      << ",\"best_absolute_determinant\":\""
      << wide_to_string(best_score)
      << "\",\"cycles\":" << statistics.cycles
      << ",\"downhill_moves\":" << statistics.downhill_moves
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(6) << elapsed_seconds
      << ",\"event\":\"" << event
      << "\",\"exact_checks\":" << statistics.exact_checks
      << ",\"exhaustive_assignments\":"
      << statistics.exhaustive_assignments
      << ",\"exhaustive_ties\":" << statistics.exhaustive_ties
      << ",\"exhaustive_total\":" << statistics.exhaustive_total
      << ",\"exploratory_moves\":" << statistics.exploratory_moves
      << ",\"first_endpoint_restarts\":"
      << statistics.first_endpoint_restarts
      << ",\"hamming_from_first\":" << weight
      << ",\"hamming_from_second\":" << (dimensions - weight)
      << ",\"inverse_rebuilds\":"
      << statistics.inverse_rebuilds
      << ",\"iterations\":" << statistics.iterations
      << ",\"kicks\":" << statistics.kicks
      << ",\"last_exact_determinant\":\""
      << wide_to_string(last_exact_score)
      << "\",\"method\":\"" << method_name(arguments)
      << "\",\"moves\":" << statistics.moves
      << ",\"promotions\":" << statistics.promotions
      << ",\"random_subset_restarts\":"
      << statistics.random_subset_restarts
      << ",\"reason\":\"" << reason
      << "\",\"restart_singular_retries\":"
      << statistics.restart_singular_retries
      << ",\"restarts\":" << statistics.restarts
      << ",\"second_endpoint_restarts\":"
      << statistics.second_endpoint_restarts
      << ",\"seed\":" << arguments.seed
      << ",\"singular_moves_rejected\":"
      << statistics.singular_moves_rejected
      << ",\"tenure\":" << tenure << "}\n";
  log.flush();
  if (!log) {
    throw std::runtime_error("cannot append research log");
  }
}

bool exact_check_and_promote(
    const State& state, Matrix& best_matrix,
    std::vector<unsigned char>& best_bits, Wide& best_score,
    Wide& last_exact_score, const Arguments& arguments,
    Statistics& statistics, std::ofstream& log,
    const Clock::time_point& started, std::size_t dimensions,
    int tenure, std::uint64_t& checkpoint_nonce) {
  const Wide exact_score = absolute(exact_determinant(state.matrix));
  ++statistics.exact_checks;
  last_exact_score = exact_score;
  if (exact_score <= best_score) return false;

  // Floating inverse scores are proposal-only. This exact Bareiss comparison
  // is the sole promotion and checkpoint gate.
  best_score = exact_score;
  best_matrix = state.matrix;
  best_bits = state.endpoint_bits;
  ++statistics.promotions;
  atomic_write_matrix(
      arguments.output, best_matrix, checkpoint_nonce++);
  const double elapsed =
      std::chrono::duration<double>(Clock::now() - started).count();
  log_record(
      log, arguments, statistics, state, dimensions, "new_best",
      "exact_bareiss_promotion", elapsed, best_score,
      last_exact_score, tenure);
  std::cout << "new best |det|=" << wide_to_string(best_score)
            << " iteration=" << statistics.iterations << '\n'
            << std::flush;
  return true;
}

long double random_unit(std::mt19937_64& randomizer) {
  constexpr long double denominator =
      static_cast<long double>(std::uint64_t{1} << 53U);
  return static_cast<long double>(randomizer() >> 11U) / denominator;
}

Move choose_move(
    const State& state,
    const std::vector<Coordinate>& differences,
    const std::vector<std::uint64_t>& tabu_until,
    std::uint64_t iteration, Wide best_score,
    const Arguments& arguments, std::mt19937_64& randomizer,
    Statistics& statistics) {
  const long double best_exact_log =
      std::log(static_cast<long double>(best_score));
  std::vector<Candidate> candidates;
  candidates.reserve(differences.size());
  std::vector<Candidate> forbidden;
  forbidden.reserve(differences.size());
  long double best_projected =
      -std::numeric_limits<long double>::infinity();

  for (std::size_t index = 0; index < differences.size(); ++index) {
    const long double ratio =
        flip_ratio(state, differences[index]);
    const long double magnitude = std::fabs(ratio);
    const long double projected =
        magnitude < kSingularRatio
            ? -std::numeric_limits<long double>::infinity()
            : state.log_abs_determinant + std::log(magnitude);
    const bool tabu = iteration < tabu_until[index];
    bool aspiration = false;
    if (tabu &&
        projected + kAspirationSlack >= best_exact_log) {
      Matrix candidate_matrix = state.matrix;
      const Coordinate coordinate = differences[index];
      candidate_matrix[coordinate.row][coordinate.column] *= -1;
      ++statistics.aspiration_probes;
      ++statistics.exact_checks;
      aspiration =
          absolute(exact_determinant(candidate_matrix)) > best_score;
      if (aspiration) {
        return Move{index, projected, true, false};
      }
    }
    const Candidate candidate{index, projected, aspiration};
    if (!tabu || aspiration) {
      candidates.push_back(candidate);
      if (projected > best_projected) {
        best_projected = projected;
      }
    } else {
      forbidden.push_back(candidate);
    }
  }

  if (candidates.empty()) {
    if (forbidden.empty()) {
      throw std::runtime_error("recombinant hypercube has no moves");
    }
    const auto earliest = std::min_element(
        forbidden.begin(), forbidden.end(),
        [&](const Candidate& left, const Candidate& right) {
          const std::uint64_t left_until =
              tabu_until[left.difference_index];
          const std::uint64_t right_until =
              tabu_until[right.difference_index];
          if (left_until != right_until) return left_until < right_until;
          return left.projected_log_abs_determinant >
                 right.projected_log_abs_determinant;
        });
    return Move{
        earliest->difference_index,
        earliest->projected_log_abs_determinant,
        false,
        false};
  }

  const bool explore =
      arguments.explore_period != 0 &&
      iteration % arguments.explore_period == 0 &&
      candidates.size() > 1;
  if (explore && std::isfinite(best_projected)) {
    std::vector<std::size_t> pool;
    std::vector<long double> weights;
    long double total = 0.0L;
    for (std::size_t index = 0; index < candidates.size(); ++index) {
      const long double loss =
          best_projected -
          candidates[index].projected_log_abs_determinant;
      if (!std::isfinite(loss) ||
          loss > static_cast<long double>(
                     arguments.maximum_log_drop)) {
        continue;
      }
      const long double weight =
          std::exp(
              -loss /
              static_cast<long double>(arguments.temperature));
      if (!(weight > 0.0L) || !std::isfinite(weight)) continue;
      pool.push_back(index);
      weights.push_back(weight);
      total += weight;
    }
    if (!pool.empty() && total > 0.0L) {
      long double choice = random_unit(randomizer) * total;
      std::size_t selected = pool.back();
      for (std::size_t offset = 0; offset < pool.size(); ++offset) {
        if (choice <= weights[offset]) {
          selected = pool[offset];
          break;
        }
        choice -= weights[offset];
      }
      return Move{
          candidates[selected].difference_index,
          candidates[selected].projected_log_abs_determinant,
          candidates[selected].aspiration,
          true};
    }
  }

  Candidate selected = candidates.front();
  for (std::size_t index = 1; index < candidates.size(); ++index) {
    const Candidate& candidate = candidates[index];
    const bool better =
        candidate.projected_log_abs_determinant >
        selected.projected_log_abs_determinant;
    const bool tied =
        std::fabs(
            candidate.projected_log_abs_determinant -
            selected.projected_log_abs_determinant) <=
        kTieTolerance;
    if (better || (tied && (randomizer() & 1U) != 0U)) {
      selected = candidate;
    }
  }
  return Move{
      selected.difference_index,
      selected.projected_log_abs_determinant,
      selected.aspiration,
      false};
}

void install_subset(
    State& state, const Matrix& first, const Matrix& second,
    const std::vector<Coordinate>& differences,
    const std::vector<unsigned char>& bits) {
  state.matrix = first;
  state.endpoint_bits = bits;
  for (std::size_t index = 0; index < differences.size(); ++index) {
    if (bits[index] == 0) continue;
    const Coordinate coordinate = differences[index];
    state.matrix[coordinate.row][coordinate.column] =
        second[coordinate.row][coordinate.column];
  }
}

std::string restart_state(
    State& state, const Matrix& first, const Matrix& second,
    const std::vector<Coordinate>& differences,
    const std::vector<unsigned char>& best_bits,
    const Arguments& arguments, Statistics& statistics,
    std::mt19937_64& randomizer) {
  const std::size_t dimensions = differences.size();
  const std::uint64_t restart_number = statistics.restarts++;
  std::vector<unsigned char> bits(dimensions, 0);
  std::string mode;
  switch (restart_number % 4U) {
    case 0:
      mode = "first_endpoint";
      ++statistics.first_endpoint_restarts;
      break;
    case 1:
      mode = "second_endpoint";
      std::fill(bits.begin(), bits.end(), 1);
      ++statistics.second_endpoint_restarts;
      break;
    case 2:
      mode = "random_subset";
      for (unsigned char& bit : bits) {
        bit = static_cast<unsigned char>(randomizer() & 1U);
      }
      ++statistics.random_subset_restarts;
      break;
    default: {
      mode = "best_subset_kick";
      bits = best_bits;
      std::vector<std::size_t> order(dimensions);
      std::iota(order.begin(), order.end(), 0);
      const std::size_t count =
          std::min<std::size_t>(
              dimensions,
              static_cast<std::size_t>(arguments.kick_size));
      for (std::size_t index = 0; index < count; ++index) {
        std::uniform_int_distribution<std::size_t> distribution(
            index, dimensions - 1);
        const std::size_t selected = distribution(randomizer);
        std::swap(order[index], order[selected]);
        bits[order[index]] ^= 1U;
      }
      ++statistics.kicks;
      break;
    }
  }

  install_subset(state, first, second, differences, bits);
  ++statistics.inverse_rebuilds;
  if (rebuild_inverse(state)) return mode;

  // A random recombinant or kick can be singular. Retry random subsets before
  // falling back to a guaranteed nonsingular exact endpoint.
  for (int attempt = 0; attempt < 32; ++attempt) {
    ++statistics.restart_singular_retries;
    for (unsigned char& bit : bits) {
      bit = static_cast<unsigned char>(randomizer() & 1U);
    }
    install_subset(state, first, second, differences, bits);
    ++statistics.inverse_rebuilds;
    if (rebuild_inverse(state)) {
      return mode + "_singular_retry";
    }
  }

  bits.assign(dimensions, 0);
  install_subset(state, first, second, differences, bits);
  ++statistics.inverse_rebuilds;
  if (!rebuild_inverse(state)) {
    throw std::runtime_error(
        "cannot recover a nonsingular endpoint after restart");
  }
  return mode + "_fallback_first";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    ensure_fresh_outputs(arguments);

    const std::string first_raw = read_file_bytes(arguments.first);
    const std::string second_raw = read_file_bytes(arguments.second);
    const Matrix first = parse_matrix(first_raw, arguments.first);
    const Matrix second = parse_matrix(second_raw, arguments.second);
    const Wide first_score = absolute(exact_determinant(first));
    const Wide second_score = absolute(exact_determinant(second));
    if (first_score == 0 || second_score == 0) {
      throw std::runtime_error("both endpoints must be nonsingular");
    }

    std::vector<Coordinate> differences;
    differences.reserve(kEntries);
    for (int row = 0; row < kOrder; ++row) {
      for (int column = 0; column < kOrder; ++column) {
        if (first[row][column] != second[row][column]) {
          differences.push_back(Coordinate{row, column});
        }
      }
    }
    if (differences.empty()) {
      throw std::runtime_error(
          "endpoints are identical; recombinant dimension is zero");
    }
    if (arguments.mode == "exact" && differences.size() >= 63) {
      throw std::runtime_error(
          "exact mode supports at most 62 differing entries");
    }

    if (!arguments.log.parent_path().empty()) {
      std::filesystem::create_directories(
          arguments.log.parent_path());
    }
    std::ofstream log(
        arguments.log, std::ios::out | std::ios::trunc);
    if (!log) {
      throw std::runtime_error(
          "cannot create research log: " + arguments.log.string());
    }

    State state;
    state.matrix = first;
    state.endpoint_bits.assign(differences.size(), 0);
    if (!rebuild_inverse(state)) {
      throw std::runtime_error("first endpoint inverse rebuild failed");
    }

    Statistics statistics;
    statistics.inverse_rebuilds = 1;
    statistics.exact_checks = 2;
    if (arguments.mode == "exact") {
      statistics.exhaustive_total =
          std::uint64_t{1} << differences.size();
    }
    Matrix best_matrix =
        first_score >= second_score ? first : second;
    std::vector<unsigned char> best_bits(
        differences.size(),
        static_cast<unsigned char>(first_score >= second_score ? 0 : 1));
    Wide best_score = std::max(first_score, second_score);
    Wide last_exact_score = first_score;
    std::uint64_t checkpoint_nonce = 0;
    atomic_write_matrix(
        arguments.output, best_matrix, checkpoint_nonce++);

    std::mt19937_64 randomizer(arguments.seed);
    const std::vector<std::uint64_t> zobrist =
        make_zobrist(differences.size(), arguments.seed);
    std::uint64_t current_hash =
        subset_hash(state.endpoint_bits, zobrist);
    std::vector<Visit> visits(kVisitTableSize);
    std::uint64_t visit_epoch = 1;
    visits[current_hash & (kVisitTableSize - 1U)] =
        Visit{current_hash, 0, visit_epoch};
    std::vector<std::uint64_t> tabu_until(
        differences.size(), 0);
    const int baseline_tenure =
        kMinimumTenure +
        static_cast<int>(arguments.seed % 5U);
    int tenure = baseline_tenure;
    std::uint64_t last_cycle_iteration = 0;

    const auto started = Clock::now();
    const auto deadline =
        started +
        std::chrono::duration_cast<Clock::duration>(
            std::chrono::duration<double>(arguments.seconds));
    auto next_heartbeat =
        started +
        std::chrono::duration_cast<Clock::duration>(
            std::chrono::duration<double>(
                arguments.heartbeat_seconds));

    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);

    log << "{\"allowed_hamming_dimension\":"
        << differences.size()
        << ",\"event\":\"start\""
        << ",\"exact_interval\":" << arguments.exact_interval
        << ",\"exact_log_window\":" << arguments.exact_log_window
        << ",\"exhaustive_total\":"
        << statistics.exhaustive_total
        << ",\"explore_period\":" << arguments.explore_period
        << ",\"first\":\""
        << json_escape(arguments.first.string())
        << "\",\"first_absolute_determinant\":\""
        << wide_to_string(first_score)
        << "\",\"first_parsed_matrix_sha256\":\""
        << sha256(matrix_bytes(first))
        << "\",\"first_raw_sha256\":\"" << sha256(first_raw)
        << "\",\"first_row_major_sign_bits_hex\":\""
        << matrix_sign_bits_hex(first)
        << "\",\"kick_size\":" << arguments.kick_size
        << ",\"maximum_iterations\":"
        << arguments.maximum_iterations
        << ",\"maximum_log_drop\":"
        << arguments.maximum_log_drop
        << ",\"method\":\"" << method_name(arguments)
        << "\",\"mode\":\"" << arguments.mode
        << "\",\"output\":\""
        << json_escape(arguments.output.string())
        << "\",\"rebuild_interval\":"
        << arguments.rebuild_interval
        << ",\"restart_interval\":"
        << arguments.restart_interval
        << ",\"second\":\""
        << json_escape(arguments.second.string())
        << "\",\"second_absolute_determinant\":\""
        << wide_to_string(second_score)
        << "\",\"second_parsed_matrix_sha256\":\""
        << sha256(matrix_bytes(second))
        << "\",\"second_raw_sha256\":\"" << sha256(second_raw)
        << "\",\"second_row_major_sign_bits_hex\":\""
        << matrix_sign_bits_hex(second)
        << "\",\"seconds\":" << arguments.seconds
        << ",\"seed\":" << arguments.seed
        << ",\"temperature\":" << arguments.temperature
        << "}\n";
    log.flush();
    if (!log) {
      throw std::runtime_error("cannot append start record");
    }
    std::cout << "start dimension=" << differences.size()
              << " first=" << wide_to_string(first_score)
              << " second=" << wide_to_string(second_score)
              << " seed=" << arguments.seed << '\n'
              << std::flush;

    if (arguments.mode == "exact") {
      std::uint64_t previous_gray = 0;
      while (!stop_requested &&
             Clock::now() < deadline &&
             statistics.exhaustive_assignments <
                 statistics.exhaustive_total &&
             (arguments.maximum_iterations == 0 ||
              statistics.exhaustive_assignments <
                  arguments.maximum_iterations)) {
        const std::uint64_t ordinal =
            statistics.exhaustive_assignments;
        const std::uint64_t gray = ordinal ^ (ordinal >> 1U);
        if (ordinal != 0) {
          std::uint64_t changed = gray ^ previous_gray;
          std::size_t difference_index = 0;
          while ((changed & 1U) == 0U) {
            changed >>= 1U;
            ++difference_index;
          }
          const Coordinate coordinate =
              differences[difference_index];
          state.matrix[coordinate.row][coordinate.column] *= -1;
          state.endpoint_bits[difference_index] ^= 1U;
          ++statistics.moves;
        }
        previous_gray = gray;

        const Wide score =
            absolute(exact_determinant(state.matrix));
        ++statistics.exact_checks;
        ++statistics.exhaustive_assignments;
        statistics.iterations =
            statistics.exhaustive_assignments;
        last_exact_score = score;
        if (score > best_score) {
          best_score = score;
          best_matrix = state.matrix;
          best_bits = state.endpoint_bits;
          statistics.exhaustive_ties = 1;
          ++statistics.promotions;
          atomic_write_matrix(
              arguments.output, best_matrix, checkpoint_nonce++);
          log_record(
              log, arguments, statistics, state,
              differences.size(), "new_best",
              "exact_bareiss_promotion",
              std::chrono::duration<double>(
                  Clock::now() - started).count(),
              best_score, last_exact_score, 0);
          std::cout << "new best |det|="
                    << wide_to_string(best_score)
                    << " assignment="
                    << statistics.exhaustive_assignments << '\n'
                    << std::flush;
        } else if (score == best_score) {
          ++statistics.exhaustive_ties;
        }

        const auto now = Clock::now();
        if (arguments.heartbeat_seconds > 0.0 &&
            now >= next_heartbeat) {
          log_record(
              log, arguments, statistics, state,
              differences.size(), "heartbeat",
              "exact_enumeration",
              std::chrono::duration<double>(
                  now - started).count(),
              best_score, last_exact_score, 0);
          next_heartbeat =
              now +
              std::chrono::duration_cast<Clock::duration>(
                  std::chrono::duration<double>(
                      arguments.heartbeat_seconds));
        }
      }

      const bool complete =
          statistics.exhaustive_assignments ==
          statistics.exhaustive_total;
      const char* reason =
          complete
              ? "complete"
              : (stop_requested
                     ? "signal"
                     : ((arguments.maximum_iterations != 0 &&
                         statistics.exhaustive_assignments >=
                             arguments.maximum_iterations)
                            ? "assignment_limit"
                            : "time_limit"));
      log_record(
          log, arguments, statistics, state, differences.size(),
          complete ? "finished" : "stopped", reason,
          std::chrono::duration<double>(
              Clock::now() - started).count(),
          best_score, last_exact_score, 0);
      std::cout
          << (complete ? "finished" : "stopped")
          << " |det|=" << wide_to_string(best_score)
          << " dimension=" << differences.size()
          << " assignments=" << statistics.exhaustive_assignments
          << "/" << statistics.exhaustive_total
          << " exact_checks=" << statistics.exact_checks
          << " ties=" << statistics.exhaustive_ties << '\n';
      return 0;
    }

    const auto iteration_limit_reached = [&]() {
      return arguments.maximum_iterations != 0 &&
             statistics.iterations >=
                 arguments.maximum_iterations;
    };

    while (!stop_requested &&
           Clock::now() < deadline &&
           !iteration_limit_reached()) {
      if (arguments.restart_interval != 0 &&
          statistics.iterations != 0 &&
          statistics.iterations % arguments.restart_interval == 0) {
        const std::string reason = restart_state(
            state, first, second, differences, best_bits,
            arguments, statistics, randomizer);
        exact_check_and_promote(
            state, best_matrix, best_bits, best_score,
            last_exact_score, arguments, statistics, log, started,
            differences.size(), tenure, checkpoint_nonce);
        current_hash =
            subset_hash(state.endpoint_bits, zobrist);
        std::fill(tabu_until.begin(), tabu_until.end(), 0);
        tenure = baseline_tenure;
        last_cycle_iteration = statistics.iterations;
        ++visit_epoch;
        visits[current_hash & (kVisitTableSize - 1U)] =
            Visit{
                current_hash, statistics.iterations, visit_epoch};
        log_record(
            log, arguments, statistics, state, differences.size(),
            "restart", reason.c_str(),
            std::chrono::duration<double>(
                Clock::now() - started).count(),
            best_score, last_exact_score, tenure);
      }

      ++statistics.iterations;
      const Move move = choose_move(
          state, differences, tabu_until,
          statistics.iterations, best_score, arguments,
          randomizer, statistics);
      const Coordinate coordinate =
          differences[move.difference_index];
      const Matrix before = state.matrix;
      const long double before_log = state.log_abs_determinant;
      bool rebuilt_during_flip = false;

      if (!apply_flip(state, coordinate, rebuilt_during_flip)) {
        if (rebuilt_during_flip) ++statistics.inverse_rebuilds;
        state.matrix = before;
        if (!rebuild_inverse(state)) {
          throw std::runtime_error(
              "cannot recover inverse after singular move");
        }
        ++statistics.inverse_rebuilds;
        ++statistics.singular_moves_rejected;
        tabu_until[move.difference_index] =
            statistics.iterations +
            static_cast<std::uint64_t>(kMaximumTenure);
        continue;
      }
      if (rebuilt_during_flip) ++statistics.inverse_rebuilds;
      state.endpoint_bits[move.difference_index] ^= 1U;
      ++statistics.moves;
      if (move.exploratory) ++statistics.exploratory_moves;
      if (move.aspiration) ++statistics.aspiration_moves;
      if (move.projected_log_abs_determinant + kTieTolerance <
          before_log) {
        ++statistics.downhill_moves;
      }
      current_hash ^= zobrist[move.difference_index];
      const int jitter = static_cast<int>(randomizer() % 5U);
      tabu_until[move.difference_index] =
          statistics.iterations +
          static_cast<std::uint64_t>(tenure + jitter + 1);

      Visit& visit =
          visits[current_hash & (kVisitTableSize - 1U)];
      if (visit.epoch == visit_epoch &&
          visit.hash == current_hash &&
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
                      std::min<std::uint64_t>(
                          cycle_length / 8U, 8U)));
        }
      }
      visit =
          Visit{current_hash, statistics.iterations, visit_epoch};
      if (statistics.iterations - last_cycle_iteration >= 512U &&
          (statistics.iterations & 127U) == 0U &&
          tenure > baseline_tenure) {
        --tenure;
      }

      const long double best_exact_log =
          std::log(static_cast<long double>(best_score));
      const bool periodic_exact =
          statistics.iterations % arguments.exact_interval == 0;
      const bool near_best =
          move.aspiration ||
          state.log_abs_determinant +
                  static_cast<long double>(
                      arguments.exact_log_window) >=
              best_exact_log;
      if (periodic_exact || near_best) {
        exact_check_and_promote(
            state, best_matrix, best_bits, best_score,
            last_exact_score, arguments, statistics, log, started,
            differences.size(), tenure, checkpoint_nonce);
      }

      if (state.moves_since_rebuild >=
          arguments.rebuild_interval) {
        ++statistics.inverse_rebuilds;
        if (!rebuild_inverse(state)) {
          ++statistics.singular_moves_rejected;
          const std::string reason = restart_state(
              state, first, second, differences, best_bits,
              arguments, statistics, randomizer);
          current_hash =
              subset_hash(state.endpoint_bits, zobrist);
          std::fill(tabu_until.begin(), tabu_until.end(), 0);
          tenure = baseline_tenure;
          last_cycle_iteration = statistics.iterations;
          ++visit_epoch;
          visits[current_hash & (kVisitTableSize - 1U)] =
              Visit{
                  current_hash, statistics.iterations, visit_epoch};
          log_record(
              log, arguments, statistics, state,
              differences.size(), "restart", reason.c_str(),
              std::chrono::duration<double>(
                  Clock::now() - started).count(),
              best_score, last_exact_score, tenure);
        }
        exact_check_and_promote(
            state, best_matrix, best_bits, best_score,
            last_exact_score, arguments, statistics, log, started,
            differences.size(), tenure, checkpoint_nonce);
      }

      const auto now = Clock::now();
      if (arguments.heartbeat_seconds > 0.0 &&
          now >= next_heartbeat) {
        log_record(
            log, arguments, statistics, state,
            differences.size(), "heartbeat", "running",
            std::chrono::duration<double>(now - started).count(),
            best_score, last_exact_score, tenure);
        next_heartbeat =
            now +
            std::chrono::duration_cast<Clock::duration>(
                std::chrono::duration<double>(
                    arguments.heartbeat_seconds));
      }
    }

    exact_check_and_promote(
        state, best_matrix, best_bits, best_score,
        last_exact_score, arguments, statistics, log, started,
        differences.size(), tenure, checkpoint_nonce);
    const bool hit_iteration_limit = iteration_limit_reached();
    const char* reason =
        stop_requested
            ? "signal"
            : (hit_iteration_limit ? "iteration_limit"
                                   : "time_limit");
    const char* event = stop_requested ? "stopped" : "finished";
    const double elapsed =
        std::chrono::duration<double>(
            Clock::now() - started).count();
    log_record(
        log, arguments, statistics, state, differences.size(),
        event, reason, elapsed, best_score, last_exact_score,
        tenure);
    std::cout << event
              << " |det|=" << wide_to_string(best_score)
              << " dimension=" << differences.size()
              << " iterations=" << statistics.iterations
              << " exact_checks=" << statistics.exact_checks
              << " promotions=" << statistics.promotions
              << " restarts=" << statistics.restarts << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "crossover_tabu: " << error.what() << '\n';
    return 2;
  }
}
