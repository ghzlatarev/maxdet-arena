// Exhaustive exact search of an entry-flip cube by fast principal minors.
//
// For an invertible base matrix A and distinct flipped coordinates
// (r_i, c_i), write delta_i = -2 A[r_i,c_i].  The matrix determinant lemma
// gives, for every selector set S,
//
//   det(A with S flipped) / det(A) = det(M[S,S]),
//   M[i,j] = 1[i=j] + delta_i (A^-1)[c_i,r_j].
//
// Thus all 2^m candidates are obtained by computing all principal minors of
// one m-by-m matrix over p = 2^32 - 5.  The elimination below is a fresh
// implementation of the fast-principal-minor recursion.  Exact zero pivots
// are temporarily shifted by one and corrected in reverse dependency order.
//
// Every order-23 sign determinant is divisible by 2^22, while
//
//   floor(sqrt(23^23))/2^22 = 1,089,457,290 < p/2.
//
// Consequently one centered residue recovers the signed integer determinant.
// Any promoted matrix is independently checked with integer Bareiss
// elimination before its checkpoint is atomically installed.

#include <algorithm>
#include <array>
#include <bit>
#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <optional>
#include <queue>
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
constexpr int kMaximumCubeDimension = 27;
constexpr std::size_t kFrontierTieMaskLimit = 32;
constexpr std::size_t kMaximumTopK = 256;
constexpr std::size_t kMaximumTopKPerWeight = 16;
constexpr std::uint64_t kPrime64 = UINT64_C(4294967291);
constexpr std::uint64_t kTwo22 = UINT64_C(1) << 22U;
constexpr std::uint64_t kHadamardQuotientBound = UINT64_C(1089457290);
constexpr std::uint64_t kFrontierFloor = UINT64_C(2779447296000000);
constexpr std::string_view kEngine =
    "fast-principal-minor-entry-cube-v1";

static_assert(2U * kHadamardQuotientBound < kPrime64);

using Clock = std::chrono::steady_clock;
using Mod = std::uint32_t;
using Wide = __int128_t;
using Matrix = std::array<std::array<int, kOrder>, kOrder>;

struct Coordinate {
  int row = 0;
  int column = 0;
};

struct Arguments {
  bool self_test = false;
  std::filesystem::path start;
  std::optional<std::filesystem::path> second;
  std::optional<std::filesystem::path> coordinates;
  std::filesystem::path output;
  std::filesystem::path log;
  std::filesystem::path report;
  std::optional<std::filesystem::path> tie_output;
  std::size_t top_k = 0;
  std::size_t top_k_per_weight = 0;
  std::optional<std::filesystem::path> top_k_output_directory;
};

struct ModFactorization {
  Mod determinant = 0;
  std::vector<Mod> inverse;
  bool nonsingular = false;
};

struct Correction {
  std::size_t mask = 0;
};

struct PrincipalMinorResult {
  std::vector<Mod> minors;
  std::uint64_t zero_pivot_corrections = 0;
};

struct RecoveredDeterminant {
  std::int64_t signed_value = 0;
  std::uint64_t absolute_value = 0;
};

struct TopCandidate {
  std::uint64_t absolute_value = 0;
  std::int64_t signed_value = 0;
  std::size_t mask = 0;
};

bool top_candidate_better(
    const TopCandidate& first, const TopCandidate& second) {
  if (first.absolute_value != second.absolute_value) {
    return first.absolute_value > second.absolute_value;
  }
  return first.mask < second.mask;
}

struct BetterTopCandidate {
  bool operator()(
      const TopCandidate& first,
      const TopCandidate& second) const {
    // priority_queue puts the element considered last at the top.  Treating
    // a better candidate as "less" makes top() the deterministic worst
    // retained candidate, ready for bounded replacement.
    return top_candidate_better(first, second);
  }
};

using TopCandidateHeap = std::priority_queue<
    TopCandidate, std::vector<TopCandidate>, BetterTopCandidate>;

void retain_top_candidate(
    TopCandidateHeap& heap, std::size_t limit,
    const TopCandidate& candidate) {
  if (limit == 0) return;
  if (heap.size() < limit) {
    heap.push(candidate);
  } else if (top_candidate_better(candidate, heap.top())) {
    heap.pop();
    heap.push(candidate);
  }
}

std::vector<TopCandidate> sorted_top_candidates(
    TopCandidateHeap heap) {
  std::vector<TopCandidate> candidates;
  candidates.reserve(heap.size());
  while (!heap.empty()) {
    candidates.push_back(heap.top());
    heap.pop();
  }
  std::sort(
      candidates.begin(), candidates.end(),
      top_candidate_better);
  return candidates;
}

Wide wide_absolute(Wide value) { return value < 0 ? -value : value; }

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

Mod add_mod(Mod first, Mod second) {
  const std::uint64_t sum =
      static_cast<std::uint64_t>(first) + second;
  return static_cast<Mod>(sum >= kPrime64 ? sum - kPrime64 : sum);
}

Mod sub_mod(Mod first, Mod second) {
  return first >= second
             ? static_cast<Mod>(first - second)
             : static_cast<Mod>(
                   static_cast<std::uint64_t>(first) + kPrime64 -
                   second);
}

Mod reduce_mod(std::uint64_t value) {
  std::uint64_t reduced =
      static_cast<std::uint32_t>(value) + 5U * (value >> 32U);
  reduced =
      static_cast<std::uint32_t>(reduced) + 5U * (reduced >> 32U);
  if (reduced >= kPrime64) reduced -= kPrime64;
  return static_cast<Mod>(reduced);
}

Mod mul_mod(Mod first, Mod second) {
  return reduce_mod(
      static_cast<std::uint64_t>(first) *
      static_cast<std::uint64_t>(second));
}

Mod neg_mod(Mod value) {
  return value == 0 ? 0 : static_cast<Mod>(kPrime64 - value);
}

Mod pow_mod(Mod base, std::uint64_t exponent) {
  Mod result = 1;
  while (exponent != 0) {
    if ((exponent & 1U) != 0) result = mul_mod(result, base);
    base = mul_mod(base, base);
    exponent >>= 1U;
  }
  return result;
}

Mod inverse_mod(Mod value) {
  if (value == 0) {
    throw std::runtime_error("attempted to invert zero modulo p");
  }
  return pow_mod(value, kPrime64 - 2U);
}

Mod signed_mod(Wide value) {
  Wide residue = value % static_cast<Wide>(kPrime64);
  if (residue < 0) residue += static_cast<Wide>(kPrime64);
  return static_cast<Mod>(residue);
}

std::size_t checked_product(
    std::size_t first, std::size_t second,
    std::string_view description) {
  if (first != 0 &&
      second > std::numeric_limits<std::size_t>::max() / first) {
    throw std::runtime_error(
        "size overflow while allocating " + std::string(description));
  }
  return first * second;
}

std::string read_file_bytes(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("cannot open input: " + path.string());
  }
  std::string bytes{
      std::istreambuf_iterator<char>(input),
      std::istreambuf_iterator<char>()};
  if (input.bad()) {
    throw std::runtime_error("cannot read input: " + path.string());
  }
  return bytes;
}

Matrix parse_matrix(
    std::string_view bytes, const std::filesystem::path& path) {
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
      nibble =
          (nibble << 1U) | static_cast<unsigned>(value == 1);
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
    for (int word_index = 0; word_index < 16; ++word_index) {
      words[word_index] =
          (static_cast<std::uint32_t>(buffer_[4 * word_index]) << 24U) |
          (static_cast<std::uint32_t>(
               buffer_[4 * word_index + 1])
           << 16U) |
          (static_cast<std::uint32_t>(
               buffer_[4 * word_index + 2])
           << 8U) |
          static_cast<std::uint32_t>(
              buffer_[4 * word_index + 3]);
    }
    for (int word_index = 16; word_index < 64; ++word_index) {
      const std::uint32_t first =
          rotate_right(words[word_index - 15], 7) ^
          rotate_right(words[word_index - 15], 18) ^
          (words[word_index - 15] >> 3U);
      const std::uint32_t second =
          rotate_right(words[word_index - 2], 17) ^
          rotate_right(words[word_index - 2], 19) ^
          (words[word_index - 2] >> 10U);
      words[word_index] =
          words[word_index - 16] + first +
          words[word_index - 7] + second;
    }
    std::uint32_t a = state_[0];
    std::uint32_t b = state_[1];
    std::uint32_t c = state_[2];
    std::uint32_t d = state_[3];
    std::uint32_t e = state_[4];
    std::uint32_t f = state_[5];
    std::uint32_t g = state_[6];
    std::uint32_t h = state_[7];
    for (int round = 0; round < 64; ++round) {
      const std::uint32_t sum_one =
          rotate_right(e, 6) ^ rotate_right(e, 11) ^
          rotate_right(e, 25);
      const std::uint32_t choice = (e & f) ^ ((~e) & g);
      const std::uint32_t temporary_one =
          h + sum_one + choice + constants[round] + words[round];
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

class AtomicJournal {
 public:
  explicit AtomicJournal(std::filesystem::path path)
      : path_(std::move(path)) {}

  void append(std::string record) {
    if (record.empty() || record.back() != '\n') record.push_back('\n');
    contents_ += record;
    atomic_write(path_, contents_, "journal", nonce_++);
  }

 private:
  std::filesystem::path path_;
  std::string contents_;
  std::uint64_t nonce_ = 0;
};

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

Wide exact_determinant_dynamic(
    const std::vector<int>& matrix, int order) {
  if (order == 0) return 1;
  std::vector<Wide> work(
      checked_product(
          static_cast<std::size_t>(order),
          static_cast<std::size_t>(order), "dynamic Bareiss matrix"));
  for (std::size_t index = 0; index < work.size(); ++index) {
    work[index] = matrix[index];
  }

  Wide previous_pivot = 1;
  int sign = 1;
  for (int column = 0; column < order - 1; ++column) {
    int pivot_row = column;
    while (pivot_row < order &&
           work[static_cast<std::size_t>(pivot_row) * order + column] ==
               0) {
      ++pivot_row;
    }
    if (pivot_row == order) return 0;
    if (pivot_row != column) {
      for (int inner = 0; inner < order; ++inner) {
        std::swap(
            work[static_cast<std::size_t>(pivot_row) * order + inner],
            work[static_cast<std::size_t>(column) * order + inner]);
      }
      sign = -sign;
    }
    const Wide pivot =
        work[static_cast<std::size_t>(column) * order + column];
    for (int row = column + 1; row < order; ++row) {
      for (int inner = column + 1; inner < order; ++inner) {
        const Wide numerator =
            work[static_cast<std::size_t>(row) * order + inner] * pivot -
            work[static_cast<std::size_t>(row) * order + column] *
                work[static_cast<std::size_t>(column) * order + inner];
        if (column != 0 && numerator % previous_pivot != 0) {
          throw std::runtime_error(
              "dynamic exact Bareiss division failed");
        }
        work[static_cast<std::size_t>(row) * order + inner] =
            column == 0 ? numerator : numerator / previous_pivot;
      }
      work[static_cast<std::size_t>(row) * order + column] = 0;
    }
    previous_pivot = pivot;
  }
  return static_cast<Wide>(sign) *
         work[static_cast<std::size_t>(order - 1) * order + order - 1];
}

ModFactorization factorize_mod(
    const std::vector<Mod>& matrix, int order) {
  if (order <= 0 ||
      matrix.size() !=
          static_cast<std::size_t>(order) *
              static_cast<std::size_t>(order)) {
    throw std::runtime_error("invalid modular matrix dimensions");
  }
  const int width = 2 * order;
  std::vector<Mod> work(
      static_cast<std::size_t>(order) *
      static_cast<std::size_t>(width));
  for (int row = 0; row < order; ++row) {
    for (int column = 0; column < order; ++column) {
      work[static_cast<std::size_t>(row) * width + column] =
          matrix[static_cast<std::size_t>(row) * order + column];
      work[static_cast<std::size_t>(row) * width + order + column] =
          row == column ? 1U : 0U;
    }
  }

  Mod determinant = 1;
  bool negative = false;
  for (int column = 0; column < order; ++column) {
    int pivot_row = column;
    while (pivot_row < order &&
           work[static_cast<std::size_t>(pivot_row) * width + column] ==
               0) {
      ++pivot_row;
    }
    if (pivot_row == order) return ModFactorization{};
    if (pivot_row != column) {
      for (int inner = 0; inner < width; ++inner) {
        std::swap(
            work[static_cast<std::size_t>(pivot_row) * width + inner],
            work[static_cast<std::size_t>(column) * width + inner]);
      }
      negative = !negative;
    }
    const Mod pivot =
        work[static_cast<std::size_t>(column) * width + column];
    determinant = mul_mod(determinant, pivot);
    const Mod pivot_inverse = inverse_mod(pivot);
    for (int inner = 0; inner < width; ++inner) {
      Mod& entry =
          work[static_cast<std::size_t>(column) * width + inner];
      entry = mul_mod(entry, pivot_inverse);
    }
    for (int row = 0; row < order; ++row) {
      if (row == column) continue;
      const Mod multiplier =
          work[static_cast<std::size_t>(row) * width + column];
      if (multiplier == 0) continue;
      for (int inner = 0; inner < width; ++inner) {
        Mod& entry =
            work[static_cast<std::size_t>(row) * width + inner];
        entry = sub_mod(
            entry,
            mul_mod(
                multiplier,
                work[static_cast<std::size_t>(column) * width +
                     inner]));
      }
    }
  }
  if (negative) determinant = neg_mod(determinant);

  ModFactorization result;
  result.determinant = determinant;
  result.nonsingular = true;
  result.inverse.resize(
      static_cast<std::size_t>(order) *
      static_cast<std::size_t>(order));
  for (int row = 0; row < order; ++row) {
    for (int column = 0; column < order; ++column) {
      result.inverse[
          static_cast<std::size_t>(row) * order + column] =
          work[static_cast<std::size_t>(row) * width + order + column];
    }
  }
  return result;
}

Mod determinant_mod(std::vector<Mod> matrix, int order) {
  if (order == 0) return 1;
  Mod determinant = 1;
  bool negative = false;
  for (int column = 0; column < order; ++column) {
    int pivot_row = column;
    while (pivot_row < order &&
           matrix[static_cast<std::size_t>(pivot_row) * order + column] ==
               0) {
      ++pivot_row;
    }
    if (pivot_row == order) return 0;
    if (pivot_row != column) {
      for (int inner = column; inner < order; ++inner) {
        std::swap(
            matrix[static_cast<std::size_t>(pivot_row) * order + inner],
            matrix[static_cast<std::size_t>(column) * order + inner]);
      }
      negative = !negative;
    }
    const Mod pivot =
        matrix[static_cast<std::size_t>(column) * order + column];
    determinant = mul_mod(determinant, pivot);
    const Mod pivot_inverse = inverse_mod(pivot);
    for (int row = column + 1; row < order; ++row) {
      const Mod entry =
          matrix[static_cast<std::size_t>(row) * order + column];
      if (entry == 0) continue;
      const Mod multiplier = mul_mod(entry, pivot_inverse);
      for (int inner = column + 1; inner < order; ++inner) {
        Mod& target =
            matrix[static_cast<std::size_t>(row) * order + inner];
        target = sub_mod(
            target,
            mul_mod(
                multiplier,
                matrix[static_cast<std::size_t>(column) * order +
                       inner]));
      }
    }
  }
  return negative ? neg_mod(determinant) : determinant;
}

std::vector<Mod> build_ratio_matrix(
    const std::vector<int>& base, int order,
    const ModFactorization& factorization,
    const std::vector<Coordinate>& support) {
  const std::size_t dimension = support.size();
  std::vector<Mod> ratio_matrix(
      checked_product(dimension, dimension, "ratio matrix"));
  for (std::size_t first = 0; first < dimension; ++first) {
    const Coordinate left = support[first];
    const int base_entry =
        base[static_cast<std::size_t>(left.row) * order + left.column];
    const Mod delta =
        base_entry == 1 ? static_cast<Mod>(kPrime64 - 2U) : 2U;
    for (std::size_t second = 0; second < dimension; ++second) {
      const Coordinate right = support[second];
      Mod value = mul_mod(
          delta,
          factorization.inverse[
              static_cast<std::size_t>(left.column) * order +
              right.row]);
      if (first == second) value = add_mod(value, 1U);
      ratio_matrix[first * dimension + second] = value;
    }
  }
  return ratio_matrix;
}

using LevelCallback = std::function<void(
    int, std::size_t, std::size_t, std::uint64_t, double)>;

PrincipalMinorResult all_principal_minors(
    const std::vector<Mod>& matrix, int dimension,
    const LevelCallback& callback = LevelCallback{}) {
  if (dimension <= 0 || dimension > kMaximumCubeDimension) {
    throw std::runtime_error("principal-minor dimension must be 1..27");
  }
  const std::size_t initial_cells =
      static_cast<std::size_t>(dimension) *
      static_cast<std::size_t>(dimension);
  if (matrix.size() != initial_cells) {
    throw std::runtime_error("principal-minor matrix has wrong size");
  }
  const std::size_t total =
      std::size_t{1} << static_cast<unsigned>(dimension);
  PrincipalMinorResult result;
  result.minors.assign(total, 0);
  result.minors[0] = 1;
  std::vector<Mod> current = matrix;
  std::vector<Correction> corrections;
  const auto started = Clock::now();

  for (int level = 0; level < dimension; ++level) {
    const int current_dimension = dimension - level;
    const std::size_t nodes =
        std::size_t{1} << static_cast<unsigned>(level);
    const std::size_t node_cells =
        static_cast<std::size_t>(current_dimension) *
        static_cast<std::size_t>(current_dimension);
    if (current.size() != checked_product(
                              nodes, node_cells,
                              "current principal-minor level")) {
      throw std::runtime_error(
          "internal principal-minor level size mismatch");
    }
    const std::size_t selector_bit =
        std::size_t{1} << static_cast<unsigned>(level);
    if (current_dimension == 1) {
      // No Schur complement follows the scalar level, so a zero is a valid
      // final multiplier.  Consuming it directly avoids a 2^(m-1)-element
      // batch-inversion buffer at the largest level.
      for (std::size_t node = 0; node < nodes; ++node) {
        result.minors[node + selector_bit] =
            mul_mod(result.minors[node], current[node]);
      }
      if (callback) {
        callback(
            level + 1, nodes, 0,
            static_cast<std::uint64_t>(corrections.size()),
            std::chrono::duration<double>(
                Clock::now() - started).count());
      }
      current.clear();
      break;
    }
    std::vector<Mod> prefixes(nodes + 1U);
    prefixes[0] = 1;
    for (std::size_t node = 0; node < nodes; ++node) {
      const std::size_t offset = node * node_cells;
      Mod pivot = current[offset];
      if (pivot == 0) {
        // Adding one to this transformed diagonal entry makes elimination
        // legal.  Reverse correction below removes its linear contribution
        // from exactly the affected principal minors.
        pivot = 1;
        current[offset] = pivot;
        corrections.push_back(Correction{node + selector_bit});
      }
      prefixes[node + 1U] = mul_mod(prefixes[node], pivot);
      result.minors[node + selector_bit] =
          mul_mod(result.minors[node], pivot);
    }

    const int next_dimension = current_dimension - 1;
    const std::size_t next_node_cells =
        static_cast<std::size_t>(next_dimension) *
        static_cast<std::size_t>(next_dimension);
    std::vector<Mod> next(
        checked_product(
            2U * nodes, next_node_cells,
            "next principal-minor level"));
    Mod inverse_suffix = inverse_mod(prefixes[nodes]);
    for (std::size_t reverse_node = nodes;
         reverse_node-- > 0;) {
      const std::size_t offset = reverse_node * node_cells;
      const Mod pivot = current[offset];
      const Mod pivot_inverse =
          mul_mod(inverse_suffix, prefixes[reverse_node]);
      inverse_suffix = mul_mod(inverse_suffix, pivot);
      if (next_dimension == 0) continue;
      const std::size_t left_offset =
          reverse_node * next_node_cells;
      const std::size_t right_offset =
          (nodes + reverse_node) * next_node_cells;
      for (int row = 0; row < next_dimension; ++row) {
        const Mod column_factor = mul_mod(
            current[
                offset +
                static_cast<std::size_t>(row + 1) *
                    current_dimension],
            pivot_inverse);
        for (int column = 0; column < next_dimension; ++column) {
          const Mod trailing =
              current[
                  offset +
                  static_cast<std::size_t>(row + 1) *
                      current_dimension +
                  column + 1];
          next[
              left_offset +
              static_cast<std::size_t>(row) * next_dimension +
              column] = trailing;
          next[
              right_offset +
              static_cast<std::size_t>(row) * next_dimension +
              column] =
              sub_mod(
                  trailing,
                  mul_mod(
                      column_factor,
                      current[offset + column + 1]));
        }
      }
    }
    current.swap(next);
    if (callback) {
      callback(
          level + 1, nodes, current.size(),
          static_cast<std::uint64_t>(corrections.size()),
          std::chrono::duration<double>(Clock::now() - started).count());
    }
  }

  // A correction recorded at bit d fixes a lower-bit prefix.  Adding one at
  // that pivot adds minor(T without d) to every affected minor T.  Descendant
  // corrections are undone first, so the dependency on the (possibly also
  // corrected) parent remains valid until its own turn.
  for (auto iterator = corrections.rbegin();
       iterator != corrections.rend(); ++iterator) {
    const std::size_t mask = iterator->mask;
    std::size_t selector_bit = 1;
    while ((selector_bit << 1U) <= mask) selector_bit <<= 1U;
    const std::size_t parent = mask - selector_bit;
    const Mod parent_minor =
        parent == 0 ? 1U : result.minors[parent];
    result.minors[mask] =
        sub_mod(result.minors[mask], parent_minor);
    const std::size_t stride = selector_bit << 1U;
    for (std::size_t affected = mask + stride;
         affected < total; affected += stride) {
      result.minors[affected] = sub_mod(
          result.minors[affected],
          result.minors[affected - selector_bit]);
    }
  }
  result.zero_pivot_corrections =
      static_cast<std::uint64_t>(corrections.size());
  return result;
}

RecoveredDeterminant recover_determinant(Mod determinant_residue) {
  static const Mod inverse_two22 =
      inverse_mod(static_cast<Mod>(kTwo22 % kPrime64));
  const Mod quotient_residue =
      mul_mod(determinant_residue, inverse_two22);
  const std::int64_t quotient =
      quotient_residue <= kPrime64 / 2U
          ? static_cast<std::int64_t>(quotient_residue)
          : static_cast<std::int64_t>(quotient_residue) -
                static_cast<std::int64_t>(kPrime64);
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
  return RecoveredDeterminant{signed_value, absolute_value};
}

void validate_support(
    const std::vector<Coordinate>& support, int order,
    int maximum_dimension) {
  if (support.empty()) {
    throw std::runtime_error("cube support must not be empty");
  }
  if (support.size() >
      static_cast<std::size_t>(maximum_dimension)) {
    throw std::runtime_error(
        "cube dimension exceeds the maximum of " +
        std::to_string(maximum_dimension));
  }
  std::vector<unsigned char> seen(
      static_cast<std::size_t>(order) *
          static_cast<std::size_t>(order),
      0);
  for (const Coordinate coordinate : support) {
    if (coordinate.row < 0 || coordinate.row >= order ||
        coordinate.column < 0 || coordinate.column >= order) {
      throw std::runtime_error("cube coordinate is out of range");
    }
    const std::size_t index =
        static_cast<std::size_t>(coordinate.row) * order +
        coordinate.column;
    if (seen[index] != 0) {
      throw std::runtime_error(
          "duplicate cube coordinate: " +
          std::to_string(coordinate.row + 1) + " " +
          std::to_string(coordinate.column + 1));
    }
    seen[index] = 1;
  }
}

std::vector<Coordinate> parse_coordinates(
    std::string_view bytes, const std::filesystem::path& path) {
  std::vector<Coordinate> support;
  std::istringstream input{std::string(bytes)};
  std::string line;
  int line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    const std::size_t comment = line.find('#');
    if (comment != std::string::npos) line.resize(comment);
    std::istringstream line_input(line);
    int row = 0;
    int column = 0;
    if (!(line_input >> row)) continue;
    if (!(line_input >> column)) {
      throw std::runtime_error(
          "coordinate line must contain row and column: " +
          path.string() + ":" + std::to_string(line_number));
    }
    std::string extra;
    if (line_input >> extra) {
      throw std::runtime_error(
          "coordinate line contains extra data: " +
          path.string() + ":" + std::to_string(line_number));
    }
    support.push_back(Coordinate{row - 1, column - 1});
  }
  validate_support(support, kOrder, kMaximumCubeDimension);
  return support;
}

std::vector<Coordinate> endpoint_differences(
    const Matrix& first, const Matrix& second) {
  std::vector<Coordinate> support;
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      if (first[row][column] != second[row][column]) {
        support.push_back(Coordinate{row, column});
      }
    }
  }
  validate_support(support, kOrder, kMaximumCubeDimension);
  return support;
}

Matrix apply_mask(
    const Matrix& base, const std::vector<Coordinate>& support,
    std::size_t mask) {
  Matrix candidate = base;
  for (std::size_t index = 0; index < support.size(); ++index) {
    if ((mask & (std::size_t{1} << index)) == 0) continue;
    const Coordinate coordinate = support[index];
    candidate[coordinate.row][coordinate.column] *= -1;
  }
  return candidate;
}

std::string support_json(
    const std::vector<Coordinate>& support) {
  std::ostringstream output;
  output << '[';
  for (std::size_t index = 0; index < support.size(); ++index) {
    if (index != 0) output << ',';
    output << '[' << support[index].row + 1 << ','
           << support[index].column + 1 << ']';
  }
  output << ']';
  return output.str();
}

std::filesystem::path resolved_target(
    const std::filesystem::path& path) {
  const std::filesystem::path parent =
      path.parent_path().empty() ? std::filesystem::path(".")
                                 : path.parent_path();
  return std::filesystem::weakly_canonical(parent) / path.filename();
}

void ensure_safe_paths(const Arguments& arguments) {
  std::vector<std::filesystem::path> inputs = {
      std::filesystem::canonical(arguments.start)};
  if (arguments.second.has_value()) {
    inputs.push_back(std::filesystem::canonical(*arguments.second));
  }
  if (arguments.coordinates.has_value()) {
    inputs.push_back(std::filesystem::canonical(*arguments.coordinates));
  }
  std::vector<std::filesystem::path> outputs = {
      resolved_target(arguments.output),
      resolved_target(arguments.log),
      resolved_target(arguments.report)};
  if (arguments.tie_output.has_value()) {
    outputs.push_back(resolved_target(*arguments.tie_output));
    if (std::filesystem::exists(*arguments.tie_output)) {
      throw std::runtime_error(
          "refusing to overwrite existing tie output: " +
          arguments.tie_output->string());
    }
  }
  if (arguments.top_k_output_directory.has_value()) {
    outputs.push_back(
        resolved_target(*arguments.top_k_output_directory));
    if (std::filesystem::exists(
            *arguments.top_k_output_directory)) {
      throw std::runtime_error(
          "refusing to reuse existing top-K output directory: " +
          arguments.top_k_output_directory->string());
    }
  }
  for (const auto& output : outputs) {
    for (const auto& input : inputs) {
      if (output == input) {
        throw std::runtime_error(
            "output, log, and report must not alias an input");
      }
    }
  }
  for (std::size_t first = 0; first < outputs.size(); ++first) {
    for (std::size_t second = first + 1;
         second < outputs.size(); ++second) {
      if (outputs[first] == outputs[second]) {
        throw std::runtime_error(
            "output, log, report, tie output, and top-K directory "
            "paths must be distinct");
      }
    }
  }
}

std::size_t theoretical_peak_bytes(int dimension) {
  const std::size_t total =
      std::size_t{1} << static_cast<unsigned>(dimension);
  std::size_t peak_elements = total;
  for (int level = 0; level < dimension; ++level) {
    const std::size_t nodes =
        std::size_t{1} << static_cast<unsigned>(level);
    const std::size_t current_dimension =
        static_cast<std::size_t>(dimension - level);
    const std::size_t next_dimension = current_dimension - 1U;
    const std::size_t level_elements =
        current_dimension == 1U
            ? total + nodes
            : total +
                  nodes * current_dimension * current_dimension +
                  (nodes + 1U) +
                  2U * nodes * next_dimension * next_dimension;
    peak_elements = std::max(peak_elements, level_elements);
  }
  return checked_product(
      peak_elements, sizeof(Mod), "theoretical memory estimate");
}

void expect(bool condition, std::string_view message) {
  if (!condition) {
    throw std::runtime_error(
        "self-test failure: " + std::string(message));
  }
}

Mod direct_principal_minor(
    const std::vector<Mod>& matrix, int dimension,
    std::size_t mask) {
  std::vector<int> selected;
  for (int index = 0; index < dimension; ++index) {
    if ((mask & (std::size_t{1} << index)) != 0) {
      selected.push_back(index);
    }
  }
  const int order = static_cast<int>(selected.size());
  std::vector<Mod> principal(
      static_cast<std::size_t>(order) *
      static_cast<std::size_t>(order));
  for (int row = 0; row < order; ++row) {
    for (int column = 0; column < order; ++column) {
      principal[static_cast<std::size_t>(row) * order + column] =
          matrix[
              static_cast<std::size_t>(selected[row]) * dimension +
              selected[column]];
    }
  }
  return determinant_mod(std::move(principal), order);
}

void test_principal_minors() {
  std::vector<std::vector<Mod>> hand_cases = {
      {0},
      {0, 1, 1, 0},
      {0, 1, 0, 1, 0, 1, 0, 1, 0},
      {1, 1, 1, 1, 1, 1, 1, 1, 1},
  };
  const std::array<int, 4> hand_dimensions = {1, 2, 3, 3};
  std::uint64_t correction_count = 0;
  for (std::size_t case_index = 0;
       case_index < hand_cases.size(); ++case_index) {
    const int dimension = hand_dimensions[case_index];
    const PrincipalMinorResult fast =
        all_principal_minors(hand_cases[case_index], dimension);
    correction_count += fast.zero_pivot_corrections;
    for (std::size_t mask = 0; mask < fast.minors.size(); ++mask) {
      expect(
          fast.minors[mask] ==
              direct_principal_minor(
                  hand_cases[case_index], dimension, mask),
          "hand principal minor mismatch");
    }
  }
  expect(
      correction_count > 0,
      "zero-pivot correction was not exercised");

  std::mt19937_64 randomizer(UINT64_C(0x2a5f9d31c4e87701));
  constexpr std::array<Mod, 7> values = {
      0U, 0U, 1U, 2U,
      static_cast<Mod>(kPrime64 - 1U),
      static_cast<Mod>(kPrime64 - 2U), 7U};
  for (int dimension = 1; dimension <= 8; ++dimension) {
    for (int trial = 0; trial < 32; ++trial) {
      std::vector<Mod> matrix(
          static_cast<std::size_t>(dimension) *
          static_cast<std::size_t>(dimension));
      for (Mod& value : matrix) {
        value = values[
            static_cast<std::size_t>(randomizer() % values.size())];
      }
      if ((trial & 1) == 0) {
        for (int index = 0; index < dimension; ++index) {
          matrix[
              static_cast<std::size_t>(index) * dimension + index] = 0;
        }
      }
      const PrincipalMinorResult fast =
          all_principal_minors(matrix, dimension);
      for (std::size_t mask = 0; mask < fast.minors.size(); ++mask) {
        expect(
            fast.minors[mask] ==
                direct_principal_minor(matrix, dimension, mask),
            "random principal minor mismatch");
      }
    }
  }
}

void test_entry_flip_cube() {
  int singular_candidates = 0;
  const auto check_cube =
      [&](const std::vector<int>& base, int order,
          const std::vector<Coordinate>& support) {
        std::vector<Mod> base_mod(base.size());
        for (std::size_t index = 0; index < base.size(); ++index) {
          base_mod[index] =
              base[index] == 1
                  ? 1U
                  : static_cast<Mod>(kPrime64 - 1U);
        }
        const ModFactorization factorization =
            factorize_mod(base_mod, order);
        expect(
            factorization.nonsingular,
            "entry-flip test base must be invertible");
        const std::vector<Mod> ratio_matrix =
            build_ratio_matrix(base, order, factorization, support);
        const PrincipalMinorResult ratios =
            all_principal_minors(
                ratio_matrix, static_cast<int>(support.size()));
        for (std::size_t mask = 0; mask < ratios.minors.size(); ++mask) {
          std::vector<int> candidate = base;
          for (std::size_t index = 0; index < support.size(); ++index) {
            if ((mask & (std::size_t{1} << index)) == 0) continue;
            const Coordinate coordinate = support[index];
            candidate[
                static_cast<std::size_t>(coordinate.row) * order +
                coordinate.column] *= -1;
          }
          const Wide exact = exact_determinant_dynamic(candidate, order);
          if (exact == 0) ++singular_candidates;
          const Mod predicted =
              mul_mod(
                  factorization.determinant, ratios.minors[mask]);
          expect(
              predicted == signed_mod(exact),
              "entry-flip Bareiss differential mismatch");
        }
      };

  check_cube(
      {1, 1, 1, -1}, 2,
      {{0, 0}, {0, 1}, {1, 0}, {1, 1}});

  std::mt19937_64 randomizer(UINT64_C(0x89ac1419f735e20d));
  for (int order = 2; order <= 7; ++order) {
    for (int trial = 0; trial < 12; ++trial) {
      std::vector<int> base(
          static_cast<std::size_t>(order) *
          static_cast<std::size_t>(order));
      ModFactorization factorization;
      do {
        std::vector<Mod> modular(base.size());
        for (std::size_t index = 0; index < base.size(); ++index) {
          base[index] = (randomizer() & 1U) == 0 ? -1 : 1;
          modular[index] =
              base[index] == 1
                  ? 1U
                  : static_cast<Mod>(kPrime64 - 1U);
        }
        factorization = factorize_mod(modular, order);
      } while (!factorization.nonsingular);

      std::vector<int> coordinate_indices(base.size());
      for (std::size_t index = 0;
           index < coordinate_indices.size(); ++index) {
        coordinate_indices[index] = static_cast<int>(index);
      }
      std::shuffle(
          coordinate_indices.begin(), coordinate_indices.end(),
          randomizer);
      const int dimension =
          std::min(8, static_cast<int>(coordinate_indices.size()));
      std::vector<Coordinate> support;
      support.reserve(static_cast<std::size_t>(dimension));
      for (int index = 0; index < dimension; ++index) {
        const int coordinate = coordinate_indices[index];
        support.push_back(
            Coordinate{coordinate / order, coordinate % order});
      }
      check_cube(base, order, support);
    }
  }
  expect(
      singular_candidates > 0,
      "singular entry-flip candidates were not exercised");
}

void test_support_validation() {
  bool duplicate_rejected = false;
  try {
    validate_support({{0, 0}, {0, 0}}, 2, 4);
  } catch (const std::runtime_error&) {
    duplicate_rejected = true;
  }
  expect(duplicate_rejected, "duplicate coordinate was accepted");
}

void test_quotient_recovery() {
  const std::array<std::int64_t, 7> quotients = {
      -static_cast<std::int64_t>(kHadamardQuotientBound),
      -662671875, -1, 0, 1, 662671875,
      static_cast<std::int64_t>(kHadamardQuotientBound)};
  for (const std::int64_t quotient : quotients) {
    const Wide determinant =
        static_cast<Wide>(quotient) * static_cast<Wide>(kTwo22);
    const RecoveredDeterminant recovered =
        recover_determinant(signed_mod(determinant));
    expect(
        recovered.signed_value ==
            static_cast<std::int64_t>(determinant),
        "centered quotient recovery mismatch");
  }
}

void test_top_k_retention() {
  std::mt19937_64 randomizer(UINT64_C(0x4e70c31a91d2b608));
  std::vector<TopCandidate> candidates;
  for (std::size_t index = 0; index < 257; ++index) {
    const std::uint64_t score =
        UINT64_C(1000000) + randomizer() % 37U;
    candidates.push_back(
        TopCandidate{
            score,
            (index & 1U) == 0
                ? static_cast<std::int64_t>(score)
                : -static_cast<std::int64_t>(score),
            index + 1U});
  }
  std::vector<TopCandidate> expected = candidates;
  std::sort(
      expected.begin(), expected.end(), top_candidate_better);
  constexpr std::size_t limit = 17;
  expected.resize(limit);

  std::shuffle(candidates.begin(), candidates.end(), randomizer);
  TopCandidateHeap heap;
  for (const TopCandidate& candidate : candidates) {
    retain_top_candidate(heap, limit, candidate);
  }
  const std::vector<TopCandidate> retained =
      sorted_top_candidates(std::move(heap));
  expect(retained.size() == expected.size(), "top-K size mismatch");
  for (std::size_t index = 0; index < expected.size(); ++index) {
    expect(
        retained[index].absolute_value ==
                expected[index].absolute_value &&
            retained[index].signed_value ==
                expected[index].signed_value &&
            retained[index].mask == expected[index].mask,
        "top-K deterministic ordering mismatch");
  }

  TopCandidateHeap empty;
  retain_top_candidate(
      empty, 0, TopCandidate{1, 1, 1});
  expect(empty.empty(), "zero-sized top-K retained a candidate");
}

void run_self_test() {
  const auto started = Clock::now();
  expect(
      sha256("abc") ==
          "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb"
          "410ff61f20015ad",
      "SHA-256 known-answer mismatch");
  test_principal_minors();
  test_entry_flip_cube();
  test_support_validation();
  test_quotient_recovery();
  test_top_k_retention();
  const double elapsed =
      std::chrono::duration<double>(Clock::now() - started).count();
  std::cout << "{\"engine\":\"" << kEngine
            << "\",\"event\":\"self_test\",\"passed\":true"
            << ",\"random_principal_minor_cases\":256"
            << ",\"maximum_random_dimension\":8"
            << ",\"singular_and_zero_pivot_cases\":true"
            << ",\"duplicate_coordinate_rejection\":true"
            << ",\"entry_flip_bareiss_differential\":true"
            << ",\"random_entry_flip_bareiss_cubes\":72"
            << ",\"deterministic_top_k_heap\":true"
            << ",\"seconds\":" << std::fixed << std::setprecision(6)
            << elapsed << "}\n";
}

std::size_t strict_top_k(std::string_view text) {
  if (text.empty() ||
      !std::all_of(
          text.begin(), text.end(),
          [](unsigned char character) {
            return character >= '0' && character <= '9';
          })) {
    throw std::runtime_error(
        "--top-k must be a non-negative integer");
  }
  std::size_t consumed = 0;
  const unsigned long long parsed =
      std::stoull(std::string(text), &consumed);
  if (consumed != text.size() || parsed > kMaximumTopK) {
    throw std::runtime_error(
        "--top-k must be between 0 and " +
        std::to_string(kMaximumTopK));
  }
  return static_cast<std::size_t>(parsed);
}

std::size_t strict_top_k_per_weight(std::string_view text) {
  if (text.empty() ||
      !std::all_of(
          text.begin(), text.end(),
          [](unsigned char character) {
            return character >= '0' && character <= '9';
          })) {
    throw std::runtime_error(
        "--top-k-per-weight must be a non-negative integer");
  }
  std::size_t consumed = 0;
  const unsigned long long parsed =
      std::stoull(std::string(text), &consumed);
  if (consumed != text.size() ||
      parsed > kMaximumTopKPerWeight) {
    throw std::runtime_error(
        "--top-k-per-weight must be between 0 and " +
        std::to_string(kMaximumTopKPerWeight));
  }
  return static_cast<std::size_t>(parsed);
}

void print_usage(std::ostream& output) {
  output
      << "usage:\n"
      << "  fast_principal_cube --self-test\n"
      << "  fast_principal_cube --start MATRIX "
         "(--second MATRIX | --coordinates FILE)\\\n"
      << "    --output MATRIX --log JSONL --report JSON "
         "[--tie-output MATRIX]\\\n"
      << "    [--top-k N [--top-k-output-dir DIRECTORY]] "
         "[--top-k-per-weight N]\n"
      << "\nCoordinates are one-based row-column pairs; blank lines and # "
         "comments are allowed.\n";
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string_view option = argv[index];
    if (option == "--help") {
      print_usage(std::cout);
      std::exit(0);
    }
    if (option == "--self-test") {
      arguments.self_test = true;
      continue;
    }
    if (index + 1 >= argc) {
      throw std::runtime_error(
          std::string(option) + " requires a value");
    }
    const std::string value = argv[++index];
    if (option == "--start") {
      arguments.start = value;
    } else if (option == "--second") {
      arguments.second = value;
    } else if (option == "--coordinates") {
      arguments.coordinates = value;
    } else if (option == "--output") {
      arguments.output = value;
    } else if (option == "--log") {
      arguments.log = value;
    } else if (option == "--report") {
      arguments.report = value;
    } else if (option == "--tie-output") {
      arguments.tie_output = value;
    } else if (option == "--top-k") {
      arguments.top_k = strict_top_k(value);
    } else if (option == "--top-k-per-weight") {
      arguments.top_k_per_weight =
          strict_top_k_per_weight(value);
    } else if (option == "--top-k-output-dir") {
      arguments.top_k_output_directory =
          std::filesystem::path(value);
    } else {
      throw std::runtime_error("unknown option: " + std::string(option));
    }
  }
  if (arguments.self_test) {
    if (argc != 2) {
      throw std::runtime_error("--self-test does not accept other options");
    }
    return arguments;
  }
  if (arguments.start.empty() || arguments.output.empty() ||
      arguments.log.empty() || arguments.report.empty()) {
    throw std::runtime_error(
        "--start, --output, --log, and --report are required");
  }
  if (arguments.second.has_value() ==
      arguments.coordinates.has_value()) {
    throw std::runtime_error(
        "provide exactly one of --second or --coordinates");
  }
  if (arguments.top_k_output_directory.has_value() &&
      arguments.top_k == 0) {
    throw std::runtime_error(
        "--top-k-output-dir requires a positive --top-k");
  }
  return arguments;
}

std::vector<int> flatten_matrix(const Matrix& matrix) {
  std::vector<int> flat(kEntries);
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      flat[static_cast<std::size_t>(row) * kOrder + column] =
          matrix[row][column];
    }
  }
  return flat;
}

std::vector<Mod> modular_matrix(const Matrix& matrix) {
  std::vector<Mod> modular(kEntries);
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      modular[static_cast<std::size_t>(row) * kOrder + column] =
          matrix[row][column] == 1
              ? 1U
              : static_cast<Mod>(kPrime64 - 1U);
    }
  }
  return modular;
}

void verify_exact_prediction(
    const Matrix& candidate,
    const RecoveredDeterminant& predicted) {
  const Wide exact = exact_determinant(candidate);
  if (exact != static_cast<Wide>(predicted.signed_value)) {
    throw std::runtime_error(
        "independent Bareiss determinant disagrees with modular recovery");
  }
  if (exact % static_cast<Wide>(kTwo22) != 0) {
    throw std::runtime_error(
        "candidate determinant violates 2^22 divisibility");
  }
}

int run_search(const Arguments& arguments) {
  ensure_safe_paths(arguments);
  const auto started = Clock::now();
  const std::string start_raw = read_file_bytes(arguments.start);
  const Matrix start = parse_matrix(start_raw, arguments.start);
  const Wide start_determinant_wide = exact_determinant(start);
  if (start_determinant_wide == 0) {
    throw std::runtime_error("start matrix must be nonsingular");
  }
  if (start_determinant_wide % static_cast<Wide>(kTwo22) != 0) {
    throw std::runtime_error(
        "start determinant violates order-23 divisibility");
  }
  if (wide_absolute(start_determinant_wide) >
      static_cast<Wide>(kHadamardQuotientBound) *
          static_cast<Wide>(kTwo22)) {
    throw std::runtime_error(
        "start determinant violates the Hadamard bound");
  }
  const std::int64_t start_determinant =
      static_cast<std::int64_t>(start_determinant_wide);
  const std::uint64_t start_score =
      static_cast<std::uint64_t>(wide_absolute(start_determinant_wide));

  std::optional<std::string> second_raw;
  std::optional<Matrix> second_matrix;
  std::optional<Wide> second_determinant;
  std::optional<std::string> coordinate_raw;
  std::vector<Coordinate> support;
  std::string support_source;
  if (arguments.second.has_value()) {
    second_raw = read_file_bytes(*arguments.second);
    second_matrix = parse_matrix(*second_raw, *arguments.second);
    second_determinant = exact_determinant(*second_matrix);
    support = endpoint_differences(start, *second_matrix);
    support_source = "endpoint_differences";
  } else {
    coordinate_raw = read_file_bytes(*arguments.coordinates);
    support =
        parse_coordinates(*coordinate_raw, *arguments.coordinates);
    support_source = "coordinate_file";
  }
  const int dimension = static_cast<int>(support.size());
  const std::size_t assignments =
      std::size_t{1} << static_cast<unsigned>(dimension);

  const ModFactorization factorization =
      factorize_mod(modular_matrix(start), kOrder);
  if (!factorization.nonsingular) {
    throw std::runtime_error(
        "start matrix is unexpectedly singular modulo p");
  }
  if (factorization.determinant !=
      signed_mod(start_determinant_wide)) {
    throw std::runtime_error(
        "modular factorization determinant disagrees with Bareiss");
  }
  const std::vector<Mod> ratio_matrix =
      build_ratio_matrix(
          flatten_matrix(start), kOrder, factorization, support);

  std::uint64_t output_nonce = 0;
  std::uint64_t report_nonce = 0;
  std::uint64_t tie_output_nonce = 0;
  atomic_write(
      arguments.output, matrix_bytes(start), "cube-matrix",
      output_nonce++);
  AtomicJournal journal(arguments.log);
  const std::string support_as_json = support_json(support);
  {
    std::ostringstream record;
    record
        << "{\"assignments\":" << assignments
        << ",\"coordinate_indexing\":\"one_based\""
        << ",\"dimension\":" << dimension
        << ",\"engine\":\"" << kEngine << "\""
        << ",\"event\":\"start\""
        << ",\"modulus\":" << kPrime64
        << ",\"output\":\""
        << json_escape(arguments.output.string()) << "\""
        << ",\"report\":\""
        << json_escape(arguments.report.string()) << "\""
        << ",\"start\":\""
        << json_escape(arguments.start.string()) << "\""
        << ",\"start_absolute_determinant\":\"" << start_score << "\""
        << ",\"start_parsed_matrix_sha256\":\""
        << sha256(matrix_bytes(start)) << "\""
        << ",\"start_raw_sha256\":\"" << sha256(start_raw) << "\""
        << ",\"support\":" << support_as_json
        << ",\"support_source\":\"" << support_source << "\""
        << ",\"frontier_tie_mask_limit\":" << kFrontierTieMaskLimit;
    if (arguments.tie_output.has_value()) {
      record
          << ",\"tie_output\":\""
          << json_escape(arguments.tie_output->string()) << "\"";
    }
    if (arguments.top_k != 0) {
      record << ",\"top_k\":" << arguments.top_k;
      if (arguments.top_k_output_directory.has_value()) {
        record
            << ",\"top_k_output_directory\":\""
            << json_escape(
                   arguments.top_k_output_directory->string())
            << "\"";
      }
    }
    if (arguments.top_k_per_weight != 0) {
      record
          << ",\"top_k_per_weight\":"
          << arguments.top_k_per_weight;
    }
    record
        << ",\"theoretical_peak_bytes\":"
        << theoretical_peak_bytes(dimension) << '}';
    journal.append(record.str());
  }
  {
    std::ostringstream in_progress;
    in_progress
        << "{\"complete\":false,\"dimension\":" << dimension
        << ",\"engine\":\"" << kEngine
        << "\",\"event\":\"in_progress\"}\n";
    atomic_write(
        arguments.report, in_progress.str(), "cube-report",
        report_nonce++);
  }

  const auto principal_started = Clock::now();
  const PrincipalMinorResult principal = all_principal_minors(
      ratio_matrix, dimension,
      [&](int completed_levels, std::size_t nodes,
          std::size_t retained_elements, std::uint64_t corrections,
          double elapsed) {
        std::ostringstream record;
        record
            << "{\"completed_levels\":" << completed_levels
            << ",\"dimension\":" << dimension
            << ",\"elapsed_seconds\":" << std::fixed
            << std::setprecision(6) << elapsed
            << ",\"engine\":\"" << kEngine << "\""
            << ",\"event\":\"principal_minor_level\""
            << ",\"nodes\":" << nodes
            << ",\"retained_elements\":" << retained_elements
            << ",\"zero_pivot_corrections\":" << corrections << '}';
        journal.append(record.str());
      });
  const double principal_seconds =
      std::chrono::duration<double>(
          Clock::now() - principal_started).count();

  const auto scan_started = Clock::now();
  std::uint64_t best_score = start_score;
  std::int64_t best_signed = start_determinant;
  std::size_t best_mask = 0;
  std::uint64_t best_ties = 0;
  std::uint64_t promotions = 0;
  std::uint64_t bareiss_verifications = 0;
  std::uint64_t nonzero_frontier_ties = 0;
  std::vector<std::size_t> frontier_tie_masks;
  frontier_tie_masks.reserve(kFrontierTieMaskLimit);
  std::optional<std::string> first_tie_bytes;
  TopCandidateHeap top_candidates_heap;
  std::vector<TopCandidateHeap> top_candidates_by_weight(
      static_cast<std::size_t>(dimension) + 1U);
  std::uint64_t top_k_eligible_candidates = 0;
  for (std::size_t mask = 0; mask < assignments; ++mask) {
    const Mod determinant_residue =
        mul_mod(factorization.determinant, principal.minors[mask]);
    const RecoveredDeterminant recovered =
        recover_determinant(determinant_residue);
    if (mask != 0 &&
        recovered.absolute_value == kFrontierFloor) {
      ++nonzero_frontier_ties;
      if (frontier_tie_masks.size() < kFrontierTieMaskLimit) {
        frontier_tie_masks.push_back(mask);
      }
      if (arguments.tie_output.has_value() &&
          !first_tie_bytes.has_value()) {
        const Matrix tie_matrix = apply_mask(start, support, mask);
        verify_exact_prediction(tie_matrix, recovered);
        ++bareiss_verifications;
        first_tie_bytes = matrix_bytes(tie_matrix);
        atomic_write(
            *arguments.tie_output, *first_tie_bytes,
            "cube-frontier-tie", tie_output_nonce++);
        std::ostringstream tie_record;
        tie_record
            << "{\"absolute_determinant\":\""
            << recovered.absolute_value << "\""
            << ",\"bareiss_verified\":true"
            << ",\"engine\":\"" << kEngine << "\""
            << ",\"event\":\"frontier_tie_artifact\""
            << ",\"mask_decimal\":\"" << mask << "\""
            << ",\"output\":\""
            << json_escape(arguments.tie_output->string()) << "\"}";
        journal.append(tie_record.str());
      }
    }
    if ((arguments.top_k != 0 ||
         arguments.top_k_per_weight != 0) &&
        mask != 0 &&
        recovered.absolute_value < kFrontierFloor) {
      ++top_k_eligible_candidates;
      const TopCandidate candidate{
          recovered.absolute_value, recovered.signed_value, mask};
      retain_top_candidate(
          top_candidates_heap, arguments.top_k, candidate);
      if (arguments.top_k_per_weight != 0) {
        const std::size_t weight =
            static_cast<std::size_t>(std::popcount(mask));
        retain_top_candidate(
            top_candidates_by_weight[weight],
            arguments.top_k_per_weight, candidate);
      }
    }
    if (recovered.absolute_value > best_score) {
      const Matrix candidate = apply_mask(start, support, mask);
      verify_exact_prediction(candidate, recovered);
      ++bareiss_verifications;
      best_score = recovered.absolute_value;
      best_signed = recovered.signed_value;
      best_mask = mask;
      best_ties = 1;
      ++promotions;
      atomic_write(
          arguments.output, matrix_bytes(candidate), "cube-matrix",
          output_nonce++);
      std::ostringstream record;
      record
          << "{\"absolute_determinant\":\"" << best_score << "\""
          << ",\"bareiss_verified\":true"
          << ",\"engine\":\"" << kEngine << "\""
          << ",\"event\":\"promotion\""
          << ",\"mask_decimal\":\"" << mask << "\""
          << ",\"signed_determinant\":\"" << best_signed << "\"}";
      journal.append(record.str());
    } else if (recovered.absolute_value == best_score) {
      ++best_ties;
    }
  }
  const double scan_seconds =
      std::chrono::duration<double>(Clock::now() - scan_started).count();

  std::vector<TopCandidate> top_candidates =
      sorted_top_candidates(std::move(top_candidates_heap));
  std::vector<std::vector<TopCandidate>>
      sorted_top_candidates_by_weight(
          static_cast<std::size_t>(dimension) + 1U);
  if (arguments.top_k_per_weight != 0) {
    for (int weight = 1; weight <= dimension; ++weight) {
      sorted_top_candidates_by_weight[
          static_cast<std::size_t>(weight)] =
          sorted_top_candidates(
              std::move(
                  top_candidates_by_weight[
                      static_cast<std::size_t>(weight)]));
    }
  }
  std::vector<std::string> top_candidate_artifact_paths(
      top_candidates.size());
  std::vector<std::string> top_candidate_artifact_hashes(
      top_candidates.size());
  if (arguments.top_k_output_directory.has_value()) {
    std::filesystem::create_directories(
        *arguments.top_k_output_directory);
    for (std::size_t index = 0;
         index < top_candidates.size(); ++index) {
      const TopCandidate candidate = top_candidates[index];
      const Matrix candidate_matrix =
          apply_mask(start, support, candidate.mask);
      verify_exact_prediction(
          candidate_matrix,
          RecoveredDeterminant{
              candidate.signed_value, candidate.absolute_value});
      ++bareiss_verifications;
      const std::string candidate_bytes =
          matrix_bytes(candidate_matrix);
      std::ostringstream filename;
      filename << "rank-" << std::setfill('0') << std::setw(3)
               << index + 1U << "-mask-" << candidate.mask
               << ".matrix.txt";
      const std::filesystem::path artifact_path =
          *arguments.top_k_output_directory / filename.str();
      atomic_write(
          artifact_path, candidate_bytes, "cube-top-k",
          static_cast<std::uint64_t>(index));
      top_candidate_artifact_paths[index] = artifact_path.string();
      top_candidate_artifact_hashes[index] =
          sha256(candidate_bytes);
    }
  }

  const Matrix best_matrix = apply_mask(start, support, best_mask);
  const RecoveredDeterminant final_prediction{
      best_signed, best_score};
  verify_exact_prediction(best_matrix, final_prediction);
  ++bareiss_verifications;
  atomic_write(
      arguments.output, matrix_bytes(best_matrix), "cube-matrix",
      output_nonce++);
  const std::string best_bytes = matrix_bytes(best_matrix);
  const double elapsed_seconds =
      std::chrono::duration<double>(Clock::now() - started).count();
  const double assignments_per_second =
      static_cast<double>(assignments) /
      std::max(elapsed_seconds, std::numeric_limits<double>::min());
  const Wide frontier_gain =
      static_cast<Wide>(best_score) - static_cast<Wide>(kFrontierFloor);

  std::ostringstream report;
  report
      << "{\"all_assignments_bound_checked\":true"
      << ",\"assignments\":" << assignments
      << ",\"assignments_per_second\":" << std::fixed
      << std::setprecision(3) << assignments_per_second
      << ",\"bareiss_verifications\":" << bareiss_verifications
      << ",\"best_absolute_determinant\":\"" << best_score << "\""
      << ",\"best_mask_decimal\":\"" << best_mask << "\""
      << ",\"best_matrix_sha256\":\"" << sha256(best_bytes) << "\""
      << ",\"best_row_major_sign_bits_hex\":\""
      << matrix_sign_bits_hex(best_matrix) << "\""
      << ",\"best_signed_determinant\":\"" << best_signed << "\""
      << ",\"best_ties\":" << best_ties
      << ",\"complete\":true"
      << ",\"coordinate_indexing\":\"one_based\""
      << ",\"dimension\":" << dimension
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(6) << elapsed_seconds
      << ",\"engine\":\"" << kEngine << "\""
      << ",\"frontier_floor\":\"" << kFrontierFloor << "\""
      << ",\"frontier_gain\":\"" << wide_to_string(frontier_gain) << "\""
      << ",\"frontier_nonzero_ties\":" << nonzero_frontier_ties
      << ",\"frontier_tie_mask_limit\":" << kFrontierTieMaskLimit
      << ",\"frontier_tie_masks_decimal\":[";
  for (std::size_t index = 0;
       index < frontier_tie_masks.size(); ++index) {
    if (index != 0) report << ',';
    report << frontier_tie_masks[index];
  }
  report
      << "]"
      << ",\"frontier_tie_masks_truncated\":"
      << (nonzero_frontier_ties > frontier_tie_masks.size()
              ? "true"
              : "false")
      << ",\"modulus\":" << kPrime64
      << ",\"output\":\"" << json_escape(arguments.output.string())
      << "\""
      << ",\"output_raw_sha256\":\"" << sha256(best_bytes) << "\""
      << ",\"principal_minor_seconds\":" << std::fixed
      << std::setprecision(6) << principal_seconds
      << ",\"promotions\":" << promotions
      << ",\"quotient_bound\":" << kHadamardQuotientBound
      << ",\"quotient_divisor\":" << kTwo22
      << ",\"scan_seconds\":" << std::fixed << std::setprecision(6)
      << scan_seconds
      << ",\"start\":\"" << json_escape(arguments.start.string())
      << "\""
      << ",\"start_absolute_determinant\":\"" << start_score << "\""
      << ",\"start_parsed_matrix_sha256\":\""
      << sha256(matrix_bytes(start)) << "\""
      << ",\"start_raw_sha256\":\"" << sha256(start_raw) << "\""
      << ",\"support\":" << support_as_json
      << ",\"support_source\":\"" << support_source << "\""
      << ",\"theoretical_peak_bytes\":"
      << theoretical_peak_bytes(dimension)
      << ",\"zero_pivot_corrections\":"
      << principal.zero_pivot_corrections;
  if (arguments.tie_output.has_value()) {
    report
        << ",\"tie_output\":\""
        << json_escape(arguments.tie_output->string()) << "\""
        << ",\"tie_output_written\":"
        << (first_tie_bytes.has_value() ? "true" : "false");
    if (first_tie_bytes.has_value()) {
      report
          << ",\"tie_output_raw_sha256\":\""
          << sha256(*first_tie_bytes) << "\"";
    }
  }
  if (arguments.top_k != 0) {
    report
        << ",\"top_k_captured\":" << top_candidates.size()
        << ",\"top_k_eligible_candidates\":"
        << top_k_eligible_candidates
        << ",\"top_k_nonfrontier_only\":true"
        << ",\"top_k_nonzero_masks_only\":true"
        << ",\"top_k_requested\":" << arguments.top_k
        << ",\"top_k_candidates\":[";
    for (std::size_t index = 0;
         index < top_candidates.size(); ++index) {
      if (index != 0) report << ',';
      const TopCandidate candidate = top_candidates[index];
      report
          << "{\"absolute_determinant\":\""
          << candidate.absolute_value << "\""
          << ",\"hamming_weight\":"
          << std::popcount(candidate.mask)
          << ",\"mask_decimal\":\"" << candidate.mask << "\""
          << ",\"rank\":" << index + 1U
          << ",\"signed_determinant\":\""
          << candidate.signed_value << "\"";
      if (arguments.top_k_output_directory.has_value()) {
        report
            << ",\"artifact\":\""
            << json_escape(top_candidate_artifact_paths[index])
            << "\""
            << ",\"artifact_bareiss_verified\":true"
            << ",\"artifact_raw_sha256\":\""
            << top_candidate_artifact_hashes[index] << "\"";
      }
      report << '}';
    }
    report << ']';
    if (arguments.top_k_output_directory.has_value()) {
      report
          << ",\"top_k_output_directory\":\""
          << json_escape(
                 arguments.top_k_output_directory->string())
          << "\"";
    }
  }
  if (arguments.top_k_per_weight != 0) {
    report
        << ",\"top_k_per_weight_eligible_candidates\":"
        << top_k_eligible_candidates
        << ",\"top_k_per_weight_nonfrontier_only\":true"
        << ",\"top_k_per_weight_nonzero_masks_only\":true"
        << ",\"top_k_per_weight_requested\":"
        << arguments.top_k_per_weight
        << ",\"top_k_per_weight_candidates\":[";
    bool first_weight = true;
    for (int weight = 1; weight <= dimension; ++weight) {
      const auto& candidates =
          sorted_top_candidates_by_weight[
              static_cast<std::size_t>(weight)];
      if (candidates.empty()) continue;
      if (!first_weight) report << ',';
      first_weight = false;
      report
          << "{\"hamming_weight\":" << weight
          << ",\"candidates\":[";
      for (std::size_t index = 0;
           index < candidates.size(); ++index) {
        if (index != 0) report << ',';
        const TopCandidate candidate = candidates[index];
        report
            << "{\"absolute_determinant\":\""
            << candidate.absolute_value << "\""
            << ",\"mask_decimal\":\"" << candidate.mask << "\""
            << ",\"rank_within_weight\":" << index + 1U
            << ",\"signed_determinant\":\""
            << candidate.signed_value << "\"}";
      }
      report << "]}";
    }
    report << ']';
  }
  if (arguments.coordinates.has_value()) {
    report
        << ",\"coordinate_file\":\""
        << json_escape(arguments.coordinates->string()) << "\""
        << ",\"coordinate_file_raw_sha256\":\""
        << sha256(*coordinate_raw) << "\"";
  }
  if (arguments.second.has_value()) {
    report
        << ",\"second\":\""
        << json_escape(arguments.second->string()) << "\""
        << ",\"second_parsed_matrix_sha256\":\""
        << sha256(matrix_bytes(*second_matrix)) << "\""
        << ",\"second_raw_sha256\":\"" << sha256(*second_raw) << "\""
        << ",\"second_signed_determinant\":\""
        << wide_to_string(*second_determinant) << "\"";
  }
  report << "}\n";
  atomic_write(
      arguments.report, report.str(), "cube-report", report_nonce++);

  {
    std::ostringstream record;
    record
        << "{\"assignments\":" << assignments
        << ",\"bareiss_verifications\":" << bareiss_verifications
        << ",\"best_absolute_determinant\":\"" << best_score << "\""
        << ",\"best_mask_decimal\":\"" << best_mask << "\""
        << ",\"best_ties\":" << best_ties
        << ",\"complete\":true"
        << ",\"elapsed_seconds\":" << std::fixed
        << std::setprecision(6) << elapsed_seconds
        << ",\"engine\":\"" << kEngine << "\""
        << ",\"event\":\"finish\""
        << ",\"frontier_gain\":\"" << wide_to_string(frontier_gain)
        << "\",\"frontier_nonzero_ties\":"
        << nonzero_frontier_ties
        << ",\"promotions\":" << promotions
        << ",\"zero_pivot_corrections\":"
        << principal.zero_pivot_corrections;
    if (arguments.top_k != 0) {
      record
          << ",\"top_k_captured\":" << top_candidates.size()
          << ",\"top_k_eligible_candidates\":"
          << top_k_eligible_candidates;
    }
    if (arguments.top_k_per_weight != 0) {
      record
          << ",\"top_k_per_weight\":"
          << arguments.top_k_per_weight
          << ",\"top_k_per_weight_eligible_candidates\":"
          << top_k_eligible_candidates;
    }
    record << '}';
    journal.append(record.str());
  }

  std::cout << "complete dimension=" << dimension
            << " assignments=" << assignments
            << " best=" << best_score
            << " ties=" << best_ties
            << " corrections=" << principal.zero_pivot_corrections
            << " seconds=" << std::fixed << std::setprecision(6)
            << elapsed_seconds << '\n';
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    if (arguments.self_test) {
      run_self_test();
      return 0;
    }
    return run_search(arguments);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
