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
constexpr int kMinimumRadius = 9;
constexpr int kMaximumRadius = 48;
constexpr std::size_t kVisitTableSize = 1U << 18U;
constexpr long double kSingularThreshold = 1.0e-20L;
constexpr long double kAspirationSlack = 1.0e-10L;
constexpr long double kProposalWarningTolerance = 1.0e-9L;
constexpr std::string_view kEngine = "hamming-sphere-tabu-v1";

using Clock = std::chrono::steady_clock;
using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using Inverse =
    std::array<std::array<long double, kOrder>, kOrder>;
using Wide = __int128_t;
using UnsignedWide = __uint128_t;

volatile std::sig_atomic_t stop_requested = 0;

void request_stop(int) { stop_requested = 1; }

struct Arguments {
  std::filesystem::path start;
  std::filesystem::path output;
  std::filesystem::path tie_output;
  std::filesystem::path checkpoint;
  std::filesystem::path log;
  std::vector<int> radii{12, 16, 24, 32};
  std::uint64_t seed = 23;
  std::uint64_t maximum_iterations = 0;
  std::uint64_t restart_iterations = 2048;
  std::size_t swap_samples = 0;
  int baseline_tenure = 13;
  int maximum_tenure = 192;
  double seconds = 3600.0;
  double heartbeat_seconds = 30.0;
  Wide score_floor = 0;
  bool score_floor_was_set = false;
};

struct State {
  Matrix matrix{};
  Inverse inverse{};
  std::array<unsigned char, kEntries> selected{};
  Wide exact_score = 0;
  long double log_abs_determinant =
      -std::numeric_limits<long double>::infinity();
  std::uint64_t subset_hash = 0;
  int radius = 0;
};

struct Move {
  int removed = -1;
  int added = -1;
  long double determinant_ratio = 0.0L;
  long double projected_log_abs_determinant =
      -std::numeric_limits<long double>::infinity();
  bool aspiration = false;
};

struct Visit {
  std::uint64_t hash = 0;
  std::uint64_t iteration = 0;
  bool occupied = false;
};

struct Statistics {
  std::uint64_t iterations = 0;
  std::uint64_t accepted_moves = 0;
  std::uint64_t downhill_moves = 0;
  std::uint64_t equal_moves = 0;
  std::uint64_t candidates_available = 0;
  std::uint64_t candidates_evaluated = 0;
  std::uint64_t exact_checks = 0;
  std::uint64_t inverse_rebuilds = 0;
  std::uint64_t restarts = 0;
  std::uint64_t restart_retries = 0;
  std::uint64_t singular_proposals = 0;
  std::uint64_t numerical_rejections = 0;
  std::uint64_t proposal_warnings = 0;
  std::uint64_t cycles = 0;
  std::uint64_t tabu_resets = 0;
  std::uint64_t aspiration_moves = 0;
  std::uint64_t strict_promotions = 0;
  std::uint64_t frontier_ties = 0;
  long double maximum_proposal_relative_error = 0.0L;
};

Wide absolute(Wide value) { return value < 0 ? -value : value; }

std::string wide_to_string(Wide value) {
  if (value == 0) return "0";
  const bool negative = value < 0;
  UnsignedWide magnitude =
      negative ? static_cast<UnsignedWide>(-(value + 1)) + 1U
               : static_cast<UnsignedWide>(value);
  std::string result;
  while (magnitude != 0) {
    result.push_back(
        static_cast<char>('0' + static_cast<int>(magnitude % 10U)));
    magnitude /= 10U;
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

template <typename Predicate>
std::string bits_hex(Predicate bit) {
  constexpr char digits[] = "0123456789abcdef";
  std::string result;
  result.reserve((kEntries + 3) / 4);
  unsigned nibble = 0;
  int used = 0;
  for (int entry = 0; entry < kEntries; ++entry) {
    nibble =
        (nibble << 1U) | static_cast<unsigned>(bit(entry) ? 1U : 0U);
    ++used;
    if (used == 4) {
      result.push_back(digits[nibble]);
      nibble = 0;
      used = 0;
    }
  }
  if (used != 0) {
    nibble <<= static_cast<unsigned>(4 - used);
    result.push_back(digits[nibble]);
  }
  return result;
}

std::string matrix_sign_bits_hex(const Matrix& matrix) {
  return bits_hex([&](int entry) {
    return matrix[entry / kOrder][entry % kOrder] == 1;
  });
}

std::string flip_bits_hex(
    const std::array<unsigned char, kEntries>& selected) {
  return bits_hex(
      [&](int entry) { return selected[entry] != 0; });
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

std::string hexadecimal_u64(std::uint64_t value) {
  std::ostringstream output;
  output << std::hex << std::setfill('0') << std::setw(16) << value;
  return output.str();
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

void atomic_write(const std::filesystem::path& path,
                  std::string_view bytes, std::string_view tag,
                  std::uint64_t nonce) {
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
    throw std::runtime_error("cannot allocate atomic output temporary");
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

  for (int column = 0; column < kOrder; ++column) {
    int pivot_row = column;
    for (int row = column + 1; row < kOrder; ++row) {
      if (std::fabs(augmented[row][column]) >
          std::fabs(augmented[pivot_row][column])) {
        pivot_row = row;
      }
    }
    if (std::fabs(augmented[pivot_row][column]) <
        kSingularThreshold) {
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
        augmented[row][inner] -= factor * augmented[column][inner];
      }
    }
  }

  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      state.inverse[row][column] =
          augmented[row][column + kOrder];
    }
  }
  return true;
}

std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

std::array<std::uint64_t, kEntries> make_zobrist(
    std::uint64_t seed) {
  std::array<std::uint64_t, kEntries> values{};
  std::uint64_t state = seed ^ 0x48616d5370686572ULL;
  for (std::uint64_t& value : values) {
    state = splitmix64(state);
    value = state;
  }
  return values;
}

std::uint64_t calculate_subset_hash(
    const std::array<unsigned char, kEntries>& selected,
    const std::array<std::uint64_t, kEntries>& zobrist) {
  std::uint64_t hash = 0;
  for (int entry = 0; entry < kEntries; ++entry) {
    if (selected[entry] != 0) hash ^= zobrist[entry];
  }
  return hash;
}

void verify_state(
    const State& state, const Matrix& baseline,
    const std::array<std::uint64_t, kEntries>& zobrist) {
  int distance = 0;
  for (int entry = 0; entry < kEntries; ++entry) {
    const int row = entry / kOrder;
    const int column = entry % kOrder;
    const bool differs =
        state.matrix[row][column] != baseline[row][column];
    if (differs != (state.selected[entry] != 0)) {
      throw std::runtime_error(
          "internal flip bitset disagrees with current matrix");
    }
    if (differs) ++distance;
  }
  if (distance != state.radius) {
    throw std::runtime_error(
        "fixed-radius invariant was violated");
  }
  if (calculate_subset_hash(state.selected, zobrist) !=
      state.subset_hash) {
    throw std::runtime_error(
        "incremental subset hash disagrees with flip bitset");
  }
}

long double rank_two_ratio(const State& state, int removed, int added) {
  const int first_row = removed / kOrder;
  const int first_column = removed % kOrder;
  const int second_row = added / kOrder;
  const int second_column = added % kOrder;
  const long double first_delta =
      -2.0L *
      static_cast<long double>(
          state.matrix[first_row][first_column]);
  const long double second_delta =
      -2.0L *
      static_cast<long double>(
          state.matrix[second_row][second_column]);

  // A' = A + U V^T with columns
  // U = [d1 e_r1, d2 e_r2] and V = [e_c1, e_c2].
  // det(A') / det(A) = det(I_2 + V^T A^{-1} U).
  const long double first_first =
      1.0L +
      first_delta *
          state.inverse[first_column][first_row];
  const long double first_second =
      second_delta *
      state.inverse[first_column][second_row];
  const long double second_first =
      first_delta *
      state.inverse[second_column][first_row];
  const long double second_second =
      1.0L +
      second_delta *
          state.inverse[second_column][second_row];
  return first_first * second_second -
         first_second * second_first;
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

int strict_integer(std::string_view text, std::string_view option) {
  const std::uint64_t parsed = strict_unsigned(text, option);
  if (parsed >
      static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error(
        std::string(option) + " is too large");
  }
  return static_cast<int>(parsed);
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

Wide strict_wide(std::string_view text, std::string_view option) {
  if (text.empty() ||
      !std::all_of(
          text.begin(), text.end(),
          [](unsigned char character) {
            return character >= '0' && character <= '9';
          })) {
    throw std::runtime_error(
        std::string(option) + " must be a non-negative integer");
  }
  const UnsignedWide maximum =
      ~static_cast<UnsignedWide>(0) >> 1U;
  UnsignedWide value = 0;
  for (const char character : text) {
    const unsigned digit =
        static_cast<unsigned>(character - '0');
    if (value > (maximum - digit) / 10U) {
      throw std::runtime_error(
          std::string(option) + " is too large");
    }
    value = value * 10U + digit;
  }
  return static_cast<Wide>(value);
}

std::vector<int> parse_radii(std::string_view text) {
  if (text.empty()) {
    throw std::runtime_error("--radii must not be empty");
  }
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
          "--radii entries must be between 9 and 48");
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
    } else if (option == "--tie-output") {
      arguments.tie_output = value();
    } else if (option == "--checkpoint") {
      arguments.checkpoint = value();
    } else if (option == "--log") {
      arguments.log = value();
    } else if (option == "--radii") {
      arguments.radii = parse_radii(value());
    } else if (option == "--radius") {
      arguments.radii =
          parse_radii(value());
    } else if (option == "--seed") {
      arguments.seed = strict_unsigned(value(), option);
    } else if (option == "--max-iterations") {
      arguments.maximum_iterations =
          strict_unsigned(value(), option);
    } else if (option == "--restart-iterations") {
      arguments.restart_iterations =
          strict_unsigned(value(), option);
    } else if (option == "--swap-samples") {
      arguments.swap_samples =
          static_cast<std::size_t>(
              strict_unsigned(value(), option));
    } else if (option == "--tabu-tenure") {
      arguments.baseline_tenure =
          strict_integer(value(), option);
    } else if (option == "--max-tabu-tenure") {
      arguments.maximum_tenure =
          strict_integer(value(), option);
    } else if (option == "--seconds") {
      arguments.seconds =
          strict_double(value(), option, false);
    } else if (option == "--heartbeat" ||
               option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds =
          strict_double(value(), option, true);
    } else if (option == "--score-floor") {
      arguments.score_floor = strict_wide(value(), option);
      arguments.score_floor_was_set = true;
    } else if (option == "--help") {
      std::cout
          << "usage: hamming_sphere_tabu --start MATRIX --output MATRIX "
             "--checkpoint JSON --log JSONL [options]\n\n"
          << "  --tie-output MATRIX       first exact non-baseline floor tie\n"
          << "  --score-floor INTEGER     exact promotion floor "
             "(default start score)\n"
          << "  --radii 12,16,24,32       fixed Hamming spheres; each in [9,48]\n"
          << "  --seed N                  deterministic random seed\n"
          << "  --seconds S               wall-clock bound (default 3600)\n"
          << "  --max-iterations N        deterministic bound; 0 disables\n"
          << "  --restart-iterations N    accepted moves per random restart; "
             "0 disables\n"
          << "  --swap-samples N          sampled swaps per move; 0 scores all\n"
          << "  --tabu-tenure N           baseline reactive tenure (default 13)\n"
          << "  --max-tabu-tenure N       maximum reactive tenure (default 192)\n"
          << "  --heartbeat-seconds S     0 disables periodic heartbeats\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }

  if (arguments.start.empty() || arguments.output.empty() ||
      arguments.checkpoint.empty() || arguments.log.empty()) {
    throw std::runtime_error(
        "--start, --output, --checkpoint, and --log are required");
  }
  if (arguments.baseline_tenure <= 0 ||
      arguments.maximum_tenure < arguments.baseline_tenure) {
    throw std::runtime_error(
        "tabu tenures must satisfy 1 <= baseline <= maximum");
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

void ensure_fresh_paths(const Arguments& arguments) {
  const std::filesystem::path start =
      std::filesystem::canonical(arguments.start);
  std::vector<std::pair<std::string, std::filesystem::path>> targets = {
      {"output", arguments.output},
      {"checkpoint", arguments.checkpoint},
      {"log", arguments.log},
  };
  if (!arguments.tie_output.empty()) {
    targets.emplace_back("tie output", arguments.tie_output);
  }

  std::vector<std::filesystem::path> resolved;
  resolved.reserve(targets.size());
  for (const auto& [label, target] : targets) {
    if (std::filesystem::exists(target)) {
      throw std::runtime_error(
          "refusing to overwrite existing " + label + ": " +
          target.string());
    }
    const std::filesystem::path normalized = resolved_target(target);
    if (normalized == start) {
      throw std::runtime_error(
          label + " path aliases the start matrix");
    }
    if (std::find(resolved.begin(), resolved.end(), normalized) !=
        resolved.end()) {
      throw std::runtime_error("output paths must be distinct");
    }
    resolved.push_back(normalized);
  }
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

std::string selected_entries_json(const State& state) {
  std::ostringstream output;
  output << '[';
  bool first = true;
  for (int entry = 0; entry < kEntries; ++entry) {
    if (state.selected[entry] == 0) continue;
    if (!first) output << ',';
    output << entry;
    first = false;
  }
  output << ']';
  return output.str();
}

std::string make_record(
    std::string_view event, const Arguments& arguments,
    const Statistics& statistics, const State& state,
    const Matrix& baseline, std::string_view input_raw_sha256,
    const Matrix& best_matrix, Wide best_score, int tenure,
    double elapsed_seconds) {
  const std::string current_bytes = matrix_bytes(state.matrix);
  const std::string best_bytes = matrix_bytes(best_matrix);
  std::ostringstream output;
  output << "{\"accepted_moves\":" << statistics.accepted_moves
         << ",\"aspiration_moves\":" << statistics.aspiration_moves
         << ",\"baseline_matrix_sha256\":\""
         << sha256(matrix_bytes(baseline))
         << "\",\"best_matrix_sha256\":\"" << sha256(best_bytes)
         << "\",\"best_score\":\"" << wide_to_string(best_score)
         << "\",\"candidates_available\":"
         << statistics.candidates_available
         << ",\"candidates_evaluated\":"
         << statistics.candidates_evaluated
         << ",\"checkpoint_path\":\""
         << json_escape(arguments.checkpoint.string())
         << "\",\"current_flip_bits_hex\":\""
         << flip_bits_hex(state.selected)
         << "\",\"current_matrix_sha256\":\""
         << sha256(current_bytes)
         << "\",\"current_score\":\""
         << wide_to_string(state.exact_score)
         << "\",\"current_sign_bits_hex\":\""
         << matrix_sign_bits_hex(state.matrix)
         << "\",\"current_subset_hash64\":\""
         << hexadecimal_u64(state.subset_hash)
         << "\",\"cycles\":" << statistics.cycles
         << ",\"downhill_moves\":" << statistics.downhill_moves
         << ",\"elapsed_seconds\":" << std::fixed
         << std::setprecision(6) << elapsed_seconds
         << ",\"engine\":\"" << kEngine
         << "\",\"equal_moves\":" << statistics.equal_moves
         << ",\"event\":\"" << json_escape(event)
         << "\",\"exact_checks\":" << statistics.exact_checks
         << ",\"frontier_ties\":" << statistics.frontier_ties
         << ",\"input_raw_sha256\":\"" << input_raw_sha256
         << "\",\"inverse_rebuilds\":"
         << statistics.inverse_rebuilds
         << ",\"iterations\":" << statistics.iterations
         << ",\"max_iterations\":" << arguments.maximum_iterations
         << ",\"max_proposal_relative_error\":"
         << std::scientific << std::setprecision(8)
         << static_cast<double>(
                statistics.maximum_proposal_relative_error)
         << ",\"numerical_rejections\":"
         << statistics.numerical_rejections
         << ",\"output_path\":\""
         << json_escape(arguments.output.string())
         << "\",\"proposal_warnings\":"
         << statistics.proposal_warnings
         << ",\"radii\":" << radii_json(arguments.radii)
         << ",\"radius\":" << state.radius
         << ",\"restart_iterations\":"
         << arguments.restart_iterations
         << ",\"restart_retries\":" << statistics.restart_retries
         << ",\"restarts\":" << statistics.restarts
         << ",\"score_floor\":\""
         << wide_to_string(arguments.score_floor)
         << "\",\"seconds\":" << std::fixed << std::setprecision(6)
         << arguments.seconds
         << ",\"seed\":" << arguments.seed
         << ",\"selected_entries\":"
         << selected_entries_json(state)
         << ",\"singular_proposals\":"
         << statistics.singular_proposals
         << ",\"start_path\":\""
         << json_escape(arguments.start.string())
         << "\",\"strict_promotions\":"
         << statistics.strict_promotions
         << ",\"swap_samples\":" << arguments.swap_samples
         << ",\"tabu_resets\":" << statistics.tabu_resets
         << ",\"tenure\":" << tenure << '}';
  return output.str();
}

void emit_record(
    std::ofstream& log, const Arguments& arguments,
    const Statistics& statistics, std::string_view event,
    const State& state, const Matrix& baseline,
    std::string_view input_raw_sha256, const Matrix& best_matrix,
    Wide best_score, int tenure, const Clock::time_point& started,
    std::uint64_t& checkpoint_nonce) {
  const double elapsed =
      std::chrono::duration<double>(Clock::now() - started).count();
  const std::string record =
      make_record(
          event, arguments, statistics, state, baseline,
          input_raw_sha256, best_matrix, best_score, tenure, elapsed);
  log << record << '\n';
  log.flush();
  if (!log) {
    throw std::runtime_error("cannot append research log");
  }
  atomic_write(
      arguments.checkpoint, record + "\n", "hamming-sphere-checkpoint",
      checkpoint_nonce++);
}

State random_state(
    const Matrix& baseline, int radius,
    const std::array<std::uint64_t, kEntries>& zobrist,
    std::mt19937_64& randomizer, Statistics& statistics) {
  std::array<int, kEntries> entries{};
  std::iota(entries.begin(), entries.end(), 0);

  for (std::uint64_t attempt = 0; attempt < 4096; ++attempt) {
    State state;
    state.matrix = baseline;
    state.radius = radius;
    std::shuffle(entries.begin(), entries.end(), randomizer);
    for (int index = 0; index < radius; ++index) {
      const int entry = entries[index];
      state.selected[entry] = 1;
      state.matrix[entry / kOrder][entry % kOrder] *= -1;
      state.subset_hash ^= zobrist[entry];
    }
    state.exact_score = absolute(exact_determinant(state.matrix));
    ++statistics.exact_checks;
    if (state.exact_score == 0 || !rebuild_inverse(state)) {
      ++statistics.restart_retries;
      continue;
    }
    ++statistics.inverse_rebuilds;
    state.log_abs_determinant =
        std::log(static_cast<long double>(state.exact_score));
    verify_state(state, baseline, zobrist);
    return state;
  }
  throw std::runtime_error(
      "could not construct a nonsingular random fixed-radius state");
}

Move choose_move(
    const State& state,
    const std::array<std::uint64_t, kEntries>& add_tabu_until,
    std::uint64_t iteration, Wide best_score,
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

  const std::size_t total =
      selected.size() * unselected.size();
  statistics.candidates_available +=
      static_cast<std::uint64_t>(total);
  const std::size_t evaluated =
      swap_samples == 0 ? total : std::min(swap_samples, total);
  std::vector<std::size_t> sampled_indices;
  if (evaluated != total) {
    sampled_indices.resize(total);
    std::iota(
        sampled_indices.begin(), sampled_indices.end(),
        static_cast<std::size_t>(0));
    for (std::size_t index = 0; index < evaluated; ++index) {
      const std::size_t remaining = total - index;
      const std::size_t chosen =
          index +
          static_cast<std::size_t>(
              randomizer() % static_cast<std::uint64_t>(remaining));
      std::swap(sampled_indices[index], sampled_indices[chosen]);
    }
    sampled_indices.resize(evaluated);
  }

  Move best;
  std::uint64_t equal_best_count = 0;
  const long double best_log =
      std::log(static_cast<long double>(best_score));
  auto score_index = [&](std::size_t flat_index) {
    const int removed =
        selected[flat_index / unselected.size()];
    const int added =
        unselected[flat_index % unselected.size()];
    const long double ratio =
        rank_two_ratio(state, removed, added);
    ++statistics.candidates_evaluated;
    if (std::fabs(ratio) < kSingularThreshold) return;
    const long double projected =
        state.log_abs_determinant +
        std::log(std::fabs(ratio));
    const bool tabu =
        !ignore_tabu && iteration < add_tabu_until[added];
    const bool aspiration =
        tabu && projected > best_log + kAspirationSlack;
    if (tabu && !aspiration) return;

    const bool better =
        best.removed < 0 ||
        projected >
            best.projected_log_abs_determinant + 1.0e-18L;
    const bool equal =
        best.removed >= 0 &&
        std::fabs(
            projected -
            best.projected_log_abs_determinant) <= 1.0e-18L;
    if (better) {
      best =
          Move{removed, added, ratio, projected, aspiration};
      equal_best_count = 1;
    } else if (equal) {
      ++equal_best_count;
      if (randomizer() % equal_best_count == 0) {
        best =
            Move{removed, added, ratio, projected, aspiration};
      }
    }
  };

  if (evaluated == total) {
    for (std::size_t index = 0; index < total; ++index) {
      score_index(index);
      if ((index & 1023U) == 0U &&
          (stop_requested || Clock::now() >= deadline)) {
        halted = true;
        break;
      }
    }
  } else {
    for (std::size_t index = 0; index < sampled_indices.size(); ++index) {
      score_index(sampled_indices[index]);
      if ((index & 1023U) == 0U &&
          (stop_requested || Clock::now() >= deadline)) {
        halted = true;
        break;
      }
    }
  }
  return best;
}

bool apply_exact_move(
    State& state, const Move& move, const Matrix& baseline,
    const std::array<std::uint64_t, kEntries>& zobrist,
    Statistics& statistics) {
  State candidate = state;
  candidate.matrix[move.removed / kOrder][move.removed % kOrder] *= -1;
  candidate.matrix[move.added / kOrder][move.added % kOrder] *= -1;
  candidate.selected[move.removed] = 0;
  candidate.selected[move.added] = 1;
  candidate.subset_hash ^=
      zobrist[move.removed] ^ zobrist[move.added];

  candidate.exact_score =
      absolute(exact_determinant(candidate.matrix));
  ++statistics.exact_checks;
  const long double predicted =
      static_cast<long double>(state.exact_score) *
      std::fabs(move.determinant_ratio);
  const long double denominator =
      std::max(
          1.0L, static_cast<long double>(candidate.exact_score));
  const long double relative_error =
      std::fabs(
          predicted -
          static_cast<long double>(candidate.exact_score)) /
      denominator;
  statistics.maximum_proposal_relative_error =
      std::max(
          statistics.maximum_proposal_relative_error,
          relative_error);
  if (relative_error > kProposalWarningTolerance) {
    ++statistics.proposal_warnings;
  }
  if (candidate.exact_score == 0) {
    ++statistics.singular_proposals;
    return false;
  }
  if (!rebuild_inverse(candidate)) {
    ++statistics.numerical_rejections;
    return false;
  }
  ++statistics.inverse_rebuilds;
  candidate.log_abs_determinant =
      std::log(static_cast<long double>(candidate.exact_score));
  verify_state(candidate, baseline, zobrist);
  state = std::move(candidate);
  return true;
}

void consider_promotion(
    const State& state, const Matrix& baseline,
    const Arguments& arguments, Statistics& statistics,
    Matrix& best_matrix, Wide& best_score, bool& tie_written,
    std::ofstream& log, std::string_view input_raw_sha256,
    int tenure, const Clock::time_point& started,
    std::uint64_t& checkpoint_nonce,
    std::uint64_t& matrix_nonce) {
  if (state.exact_score > best_score) {
    best_score = state.exact_score;
    best_matrix = state.matrix;
    ++statistics.strict_promotions;
    atomic_write(
        arguments.output, matrix_bytes(best_matrix),
        "hamming-sphere-best", matrix_nonce++);
    emit_record(
        log, arguments, statistics, "new_best", state, baseline,
        input_raw_sha256, best_matrix, best_score, tenure, started,
        checkpoint_nonce);
    std::cout << "new best |det|=" << wide_to_string(best_score)
              << " radius=" << state.radius
              << " iteration=" << statistics.iterations << '\n'
              << std::flush;
  }

  if (!tie_written && !arguments.tie_output.empty() &&
      state.exact_score == arguments.score_floor &&
      state.matrix != baseline) {
    atomic_write(
        arguments.tie_output, matrix_bytes(state.matrix),
        "hamming-sphere-tie", matrix_nonce++);
    tie_written = true;
    ++statistics.frontier_ties;
    emit_record(
        log, arguments, statistics, "frontier_tie", state, baseline,
        input_raw_sha256, best_matrix, best_score, tenure, started,
        checkpoint_nonce);
    std::cout << "retained exact floor tie |det|="
              << wide_to_string(state.exact_score)
              << " radius=" << state.radius
              << " iteration=" << statistics.iterations << '\n'
              << std::flush;
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments parsed_arguments = parse_arguments(argc, argv);
    Arguments arguments = parsed_arguments;
    ensure_fresh_paths(arguments);
    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);

    const std::string input_bytes =
        read_file_bytes(arguments.start);
    const std::string input_raw_sha256 = sha256(input_bytes);
    const Matrix baseline =
        parse_matrix(input_bytes, arguments.start);
    const Wide baseline_score =
        absolute(exact_determinant(baseline));
    if (baseline_score == 0) {
      throw std::runtime_error(
          "start matrix must have nonzero determinant");
    }
    if (!arguments.score_floor_was_set) {
      arguments.score_floor = baseline_score;
    }
    if (baseline_score < arguments.score_floor) {
      throw std::runtime_error(
          "exact start score is below --score-floor");
    }

    for (const std::filesystem::path& path :
         {arguments.output, arguments.checkpoint, arguments.log}) {
      if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
      }
    }
    if (!arguments.tie_output.empty() &&
        !arguments.tie_output.parent_path().empty()) {
      std::filesystem::create_directories(
          arguments.tie_output.parent_path());
    }
    std::ofstream log(arguments.log, std::ios::out);
    if (!log) {
      throw std::runtime_error(
          "cannot create research log: " +
          arguments.log.string());
    }

    Statistics statistics;
    Matrix best_matrix = baseline;
    Wide best_score = baseline_score;
    std::uint64_t matrix_nonce = 0;
    std::uint64_t checkpoint_nonce = 0;
    atomic_write(
        arguments.output, matrix_bytes(best_matrix),
        "hamming-sphere-best", matrix_nonce++);

    std::mt19937_64 randomizer(arguments.seed);
    const auto zobrist = make_zobrist(arguments.seed);
    std::size_t radius_index = 0;
    State state =
        random_state(
            baseline, arguments.radii[radius_index], zobrist,
            randomizer, statistics);
    ++statistics.restarts;

    // An accepted exchange removes one baseline-relative flip.  Tabuing that
    // entry only in the future "add" role blocks the exact reverse exchange
    // without exhausting the much smaller selected side of the sphere.
    std::array<std::uint64_t, kEntries> add_tabu_until{};
    std::vector<Visit> visits(kVisitTableSize);
    visits[state.subset_hash & (kVisitTableSize - 1U)] =
        Visit{state.subset_hash, 0, true};
    int tenure = arguments.baseline_tenure;
    std::uint64_t last_cycle_iteration = 0;
    std::uint64_t restart_start_accepted_moves = 0;
    bool tie_written = false;

    const auto started = Clock::now();
    const Clock::time_point deadline =
        started +
        std::chrono::duration_cast<Clock::duration>(
            std::chrono::duration<double>(arguments.seconds));
    Clock::time_point next_heartbeat =
        started +
        std::chrono::duration_cast<Clock::duration>(
            std::chrono::duration<double>(
                arguments.heartbeat_seconds));

    emit_record(
        log, arguments, statistics, "start", state, baseline,
        input_raw_sha256, best_matrix, best_score, tenure, started,
        checkpoint_nonce);
    consider_promotion(
        state, baseline, arguments, statistics, best_matrix, best_score,
        tie_written, log, input_raw_sha256, tenure, started,
        checkpoint_nonce, matrix_nonce);
    std::cout << "start |det|=" << wide_to_string(baseline_score)
              << " score_floor=" << wide_to_string(arguments.score_floor)
              << " radii=" << radii_json(arguments.radii)
              << " swap_samples=" << arguments.swap_samples << '\n'
              << std::flush;

    while (!stop_requested) {
      if (Clock::now() >= deadline) break;
      if (arguments.maximum_iterations != 0 &&
          statistics.iterations >= arguments.maximum_iterations) {
        break;
      }

      if (arguments.restart_iterations != 0 &&
          statistics.accepted_moves - restart_start_accepted_moves >=
              arguments.restart_iterations) {
        radius_index =
            (radius_index + 1U) % arguments.radii.size();
        state =
            random_state(
                baseline, arguments.radii[radius_index], zobrist,
                randomizer, statistics);
        ++statistics.restarts;
        restart_start_accepted_moves =
            statistics.accepted_moves;
        add_tabu_until.fill(0);
        std::fill(visits.begin(), visits.end(), Visit{});
        visits[state.subset_hash & (kVisitTableSize - 1U)] =
            Visit{state.subset_hash, statistics.iterations, true};
        tenure = arguments.baseline_tenure;
        last_cycle_iteration = statistics.iterations;
        consider_promotion(
            state, baseline, arguments, statistics, best_matrix,
            best_score, tie_written, log, input_raw_sha256, tenure,
            started, checkpoint_nonce, matrix_nonce);
        emit_record(
            log, arguments, statistics, "restart", state, baseline,
            input_raw_sha256, best_matrix, best_score, tenure, started,
            checkpoint_nonce);
      }

      ++statistics.iterations;
      bool halted = false;
      Move move =
          choose_move(
              state, add_tabu_until, statistics.iterations, best_score,
              arguments.swap_samples, randomizer, statistics, false,
              deadline, halted);
      if (halted) {
        break;
      }
      if (move.removed < 0) {
        add_tabu_until.fill(0);
        ++statistics.tabu_resets;
        move =
            choose_move(
                state, add_tabu_until, statistics.iterations, best_score,
                arguments.swap_samples, randomizer, statistics, true,
                deadline, halted);
      }
      if (halted) {
        break;
      }
      if (move.removed < 0) {
        throw std::runtime_error(
            "no nonsingular exchange proposal was available");
      }

      const Wide previous_score = state.exact_score;
      if (!apply_exact_move(
              state, move, baseline, zobrist, statistics)) {
        const std::uint64_t expires =
            statistics.iterations +
            static_cast<std::uint64_t>(
                arguments.maximum_tenure);
        add_tabu_until[move.added] = expires;
        continue;
      }
      ++statistics.accepted_moves;
      if (state.exact_score < previous_score) {
        ++statistics.downhill_moves;
      } else if (state.exact_score == previous_score) {
        ++statistics.equal_moves;
      }
      if (move.aspiration) ++statistics.aspiration_moves;

      const int jitter =
          static_cast<int>(randomizer() % 7U);
      const std::uint64_t expires =
          statistics.iterations +
          static_cast<std::uint64_t>(tenure + jitter + 1);
      add_tabu_until[move.removed] = expires;

      Visit& visit =
          visits[state.subset_hash & (kVisitTableSize - 1U)];
      if (visit.occupied && visit.hash == state.subset_hash &&
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
      visit =
          Visit{state.subset_hash, statistics.iterations, true};
      if (statistics.iterations - last_cycle_iteration >= 2048U &&
          (statistics.iterations & 127U) == 0U &&
          tenure > arguments.baseline_tenure) {
        --tenure;
      }

      consider_promotion(
          state, baseline, arguments, statistics, best_matrix,
          best_score, tie_written, log, input_raw_sha256, tenure,
          started, checkpoint_nonce, matrix_nonce);

      const auto now = Clock::now();
      if (arguments.heartbeat_seconds > 0.0 &&
          now >= next_heartbeat) {
        emit_record(
            log, arguments, statistics, "heartbeat", state, baseline,
            input_raw_sha256, best_matrix, best_score, tenure, started,
            checkpoint_nonce);
        std::cout
            << "heartbeat radius=" << state.radius
            << " iterations=" << statistics.iterations
            << " accepted=" << statistics.accepted_moves
            << " candidates=" << statistics.candidates_evaluated
            << " exact_checks=" << statistics.exact_checks
            << " current=" << wide_to_string(state.exact_score)
            << " best=" << wide_to_string(best_score)
            << " tenure=" << tenure << '\n'
            << std::flush;
        next_heartbeat =
            now +
            std::chrono::duration_cast<Clock::duration>(
                std::chrono::duration<double>(
                    arguments.heartbeat_seconds));
      }
    }

    const std::string_view final_event =
        stop_requested ? "stopped" : "finished";
    emit_record(
        log, arguments, statistics, final_event, state, baseline,
        input_raw_sha256, best_matrix, best_score, tenure, started,
        checkpoint_nonce);
    const double elapsed =
        std::chrono::duration<double>(Clock::now() - started).count();
    std::cout << final_event
              << " best |det|=" << wide_to_string(best_score)
              << " iterations=" << statistics.iterations
              << " accepted=" << statistics.accepted_moves
              << " candidates=" << statistics.candidates_evaluated
              << " exact_checks=" << statistics.exact_checks
              << " restarts=" << statistics.restarts
              << " elapsed=" << std::fixed << std::setprecision(3)
              << elapsed << "s\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "hamming_sphere_tabu: " << error.what() << '\n';
    return 2;
  }
}
