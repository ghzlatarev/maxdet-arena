// Exact Hamming-radius <= 3 audit for frontier-score order-23 sign matrices.
//
// For an invertible base A, flipping entries (r_i,c_i) changes A by U V^T,
// where U_i = -2 A[r_i,c_i] e_{r_i} and V_i = e_{c_i}.  The matrix
// determinant lemma gives
//
//   det(A') = det(A) det(I + V^T A^-1 U).
//
// We evaluate this identity modulo p = 2^32 - 5.  Every determinant of an
// order-23 sign matrix is divisible by 2^22, and Hadamard's bound gives
//
//   |det(A)| / 2^22 <= 1,089,457,290 < p/2.
//
// Therefore the centered residue of det(A') / 2^22 uniquely recovers the
// exact signed integer determinant.  This is exact modular arithmetic, not a
// floating screen.  The implementation differentially checks deterministic
// and pseudorandom radius-1/2/3 samples with fraction-free Bareiss, and every
// retained frontier tie or strict improvement is independently checked with
// Bareiss before being written atomically.

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

namespace {

constexpr int kOrder = 23;
constexpr int kEntries = kOrder * kOrder;
constexpr std::uint64_t kPrime64 = UINT64_C(4294967291);
constexpr std::uint64_t kTwo22 = UINT64_C(1) << 22U;
constexpr std::uint64_t kHadamardQuotientBound = UINT64_C(1089457290);
constexpr std::uint64_t kDefaultFrontier = UINT64_C(2779447296000000);
constexpr std::uint64_t kRadius1Count = 529;
constexpr std::uint64_t kRadius2Count = 139656;
constexpr std::uint64_t kRadius3Count = 24532904;
constexpr std::uint64_t kPerRepresentativeCount =
    kRadius1Count + kRadius2Count + kRadius3Count;
constexpr std::string_view kEngine = "frontier-portal-radius3-v1";

static_assert(2U * kHadamardQuotientBound < kPrime64);
static_assert(kPerRepresentativeCount == UINT64_C(24673089));

using Clock = std::chrono::steady_clock;
using Mod = std::uint32_t;
using Wide = __int128_t;
using Matrix = std::array<std::array<int, kOrder>, kOrder>;

struct Arguments {
  std::filesystem::path manifest;
  std::filesystem::path authority_audit;
  std::filesystem::path source;
  std::filesystem::path binary;
  std::filesystem::path output;
  std::filesystem::path ties_directory;
  std::string build_command;
  std::uint64_t frontier = kDefaultFrontier;
  unsigned threads = 1;
  bool differential_only = false;
};

struct ManifestEntry {
  std::string ht_certificate_sha256;
  std::string expected_raw_sha256;
  std::filesystem::path path;
};

struct ModFactorization {
  Mod determinant = 0;
  std::array<Mod, kEntries> inverse{};
  bool nonsingular = false;
};

struct CoordinateSet {
  int weight = 0;
  std::array<int, 3> entries = {-1, -1, -1};
};

struct RadiusStatistics {
  std::uint64_t candidates = 0;
  std::uint64_t frontier_ties = 0;
  std::uint64_t strict_improvements = 0;
  std::uint64_t singular = 0;
  std::uint64_t best_absolute_determinant = 0;
  std::int64_t best_signed_determinant = 0;
  CoordinateSet best_coordinates;
};

struct RetainedArtifact {
  std::string ht_certificate_sha256;
  CoordinateSet coordinates;
  std::int64_t signed_determinant = 0;
  std::filesystem::path path;
  std::string raw_sha256;
};

struct RepresentativeResult {
  ManifestEntry manifest;
  std::string canonical_matrix_sha256;
  std::int64_t base_signed_determinant = 0;
  std::uint64_t base_absolute_determinant = 0;
  std::uint64_t bareiss_differential_samples = 0;
  std::array<RadiusStatistics, 3> radius;
  std::vector<RetainedArtifact> artifacts;
  double differential_seconds = 0.0;
  double enumeration_seconds = 0.0;
  bool complete = false;
};

std::mutex output_mutex;

std::uint64_t absolute_i64(std::int64_t value) {
  return value < 0 ? static_cast<std::uint64_t>(-value)
                   : static_cast<std::uint64_t>(value);
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
          (static_cast<std::uint32_t>(buffer_[4 * word_index + 1])
           << 16U) |
          (static_cast<std::uint32_t>(buffer_[4 * word_index + 2])
           << 8U) |
          static_cast<std::uint32_t>(buffer_[4 * word_index + 3]);
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
    std::string_view tag) {
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
         std::to_string(attempt) + ".tmp");
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

ModFactorization factorize_mod(const Matrix& matrix) {
  constexpr int width = 2 * kOrder;
  std::array<std::array<Mod, width>, kOrder> work{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      work[row][column] =
          matrix[row][column] == 1 ? 1U
                                   : static_cast<Mod>(kPrime64 - 1U);
      work[row][kOrder + column] = row == column ? 1U : 0U;
    }
  }

  Mod determinant = 1;
  bool negative = false;
  for (int column = 0; column < kOrder; ++column) {
    int pivot_row = column;
    while (pivot_row < kOrder && work[pivot_row][column] == 0) {
      ++pivot_row;
    }
    if (pivot_row == kOrder) return ModFactorization{};
    if (pivot_row != column) {
      std::swap(work[pivot_row], work[column]);
      negative = !negative;
    }
    const Mod pivot = work[column][column];
    determinant = mul_mod(determinant, pivot);
    const Mod pivot_inverse = inverse_mod(pivot);
    for (int inner = 0; inner < width; ++inner) {
      work[column][inner] =
          mul_mod(work[column][inner], pivot_inverse);
    }
    for (int row = 0; row < kOrder; ++row) {
      if (row == column) continue;
      const Mod multiplier = work[row][column];
      if (multiplier == 0) continue;
      for (int inner = 0; inner < width; ++inner) {
        work[row][inner] =
            sub_mod(
                work[row][inner],
                mul_mod(multiplier, work[column][inner]));
      }
    }
  }
  if (negative) determinant = neg_mod(determinant);

  ModFactorization result;
  result.determinant = determinant;
  result.nonsingular = true;
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      result.inverse[row * kOrder + column] =
          work[row][kOrder + column];
    }
  }
  return result;
}

std::int64_t recover_determinant_from_quotient(Mod quotient) {
  const std::int64_t centered =
      quotient <= kPrime64 / 2U
          ? static_cast<std::int64_t>(quotient)
          : static_cast<std::int64_t>(quotient) -
                static_cast<std::int64_t>(kPrime64);
  if (absolute_i64(centered) > kHadamardQuotientBound) {
    throw std::runtime_error(
        "centered determinant quotient exceeds Hadamard bound");
  }
  return centered * static_cast<std::int64_t>(kTwo22);
}

std::vector<ManifestEntry> parse_manifest(std::string_view bytes) {
  std::vector<ManifestEntry> entries;
  std::istringstream input{std::string(bytes)};
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    if (line.empty() || line[0] == '#') continue;
    std::array<std::string, 3> fields;
    std::size_t start = 0;
    for (int field = 0; field < 3; ++field) {
      const std::size_t tab =
          field == 2 ? std::string::npos : line.find('\t', start);
      if (field != 2 && tab == std::string::npos) {
        throw std::runtime_error(
            "manifest line " + std::to_string(line_number) +
            " must contain three tab-separated fields");
      }
      fields[field] =
          line.substr(
              start,
              tab == std::string::npos ? tab : tab - start);
      start = tab == std::string::npos ? line.size() : tab + 1;
    }
    if (fields[0].size() != 64 || fields[1].size() != 64 ||
        fields[2].empty()) {
      throw std::runtime_error(
          "invalid manifest field on line " +
          std::to_string(line_number));
    }
    entries.push_back(
        ManifestEntry{fields[0], fields[1], fields[2]});
  }
  if (entries.empty()) {
    throw std::runtime_error(
        "manifest must contain at least one representative");
  }
  std::sort(
      entries.begin(), entries.end(),
      [](const ManifestEntry& first, const ManifestEntry& second) {
        return first.ht_certificate_sha256 <
               second.ht_certificate_sha256;
      });
  for (std::size_t index = 1; index < entries.size(); ++index) {
    if (entries[index - 1].ht_certificate_sha256 ==
        entries[index].ht_certificate_sha256) {
      throw std::runtime_error("duplicate HT certificate in manifest");
    }
  }
  return entries;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  const unsigned hardware = std::thread::hardware_concurrency();
  arguments.threads = hardware == 0 ? 1U : std::min(6U, hardware);
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    const auto value = [&]() -> std::string {
      if (index + 1 >= argc) {
        throw std::runtime_error("missing value for " + option);
      }
      return argv[++index];
    };
    if (option == "--manifest") {
      arguments.manifest = value();
    } else if (option == "--authority-audit") {
      arguments.authority_audit = value();
    } else if (option == "--source") {
      arguments.source = value();
    } else if (option == "--binary") {
      arguments.binary = value();
    } else if (option == "--output") {
      arguments.output = value();
    } else if (option == "--ties-dir") {
      arguments.ties_directory = value();
    } else if (option == "--build-command") {
      arguments.build_command = value();
    } else if (option == "--frontier") {
      arguments.frontier = std::stoull(value());
    } else if (option == "--threads") {
      arguments.threads =
          static_cast<unsigned>(std::stoul(value()));
    } else if (option == "--differential-only") {
      arguments.differential_only = true;
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (arguments.manifest.empty() ||
      arguments.authority_audit.empty() ||
      arguments.source.empty() ||
      arguments.binary.empty() ||
      arguments.output.empty() ||
      arguments.ties_directory.empty() ||
      arguments.build_command.empty()) {
    throw std::runtime_error(
        "usage: frontier_portal_radius3 --manifest FILE "
        "--authority-audit FILE --source FILE --binary FILE "
        "--output FILE --ties-dir DIR --build-command COMMAND "
        "[--frontier N] [--threads N] "
        "[--differential-only]");
  }
  if (arguments.frontier == 0 || arguments.threads == 0) {
    throw std::runtime_error("frontier and thread count must be positive");
  }
  return arguments;
}

void validate_authority(
    std::string_view authority,
    const std::vector<ManifestEntry>& entries) {
  for (const ManifestEntry& entry : entries) {
    if (authority.find(entry.ht_certificate_sha256) ==
            std::string_view::npos ||
        authority.find(entry.path.string()) == std::string_view::npos) {
      throw std::runtime_error(
          "manifest representative absent from authority audit: " +
          entry.ht_certificate_sha256);
    }
  }
}

std::uint64_t parse_unsigned_json_field(
    std::string_view json, std::string_view field) {
  const std::string key = "\"" + std::string(field) + "\"";
  std::size_t position = json.find(key);
  if (position == std::string_view::npos) {
    throw std::runtime_error(
        "authority audit lacks required field: " +
        std::string(field));
  }
  position = json.find(':', position + key.size());
  if (position == std::string_view::npos) {
    throw std::runtime_error(
        "authority audit has malformed field: " +
        std::string(field));
  }
  ++position;
  while (position < json.size() &&
         (json[position] == ' ' || json[position] == '\t' ||
          json[position] == '\r' || json[position] == '\n')) {
    ++position;
  }
  if (position == json.size() ||
      json[position] < '0' || json[position] > '9') {
    throw std::runtime_error(
        "authority audit field is not unsigned: " +
        std::string(field));
  }
  std::uint64_t value = 0;
  while (position < json.size() &&
         json[position] >= '0' && json[position] <= '9') {
    const unsigned digit =
        static_cast<unsigned>(json[position] - '0');
    if (value >
        (std::numeric_limits<std::uint64_t>::max() - digit) / 10U) {
      throw std::runtime_error(
          "authority audit field overflows uint64: " +
          std::string(field));
    }
    value = value * 10U + digit;
    ++position;
  }
  return value;
}

bool is_original_six_class_manifest(
    const std::vector<RepresentativeResult>& results) {
  constexpr std::array<std::string_view, 6> certificates = {
      "9035bdf2a85b8a2a600a76c6d55af36f627327694e904e7b483e88993716c91b",
      "b584c923ea12af6634bb5681f1eb903196e5a7cedadf5d1ea321bc63613c986e",
      "b64c33090c9aac4dc7386213917cec66ac28d2da9685319b585067e95d0d63f6",
      "df0b940533f84c9d61ec8df73000b4fe3a646f073a6fd32f3826d72331ebedc0",
      "eb138a06ec638735c34bdacf77bd1cdd869c5d2fbc3450be25d63af0cde1a134",
      "ff1b5d3735bd5735fd4648ba24627fa997045ce004d6ea08a9a5ff7567919490",
  };
  if (results.size() != certificates.size()) return false;
  for (std::size_t index = 0; index < certificates.size(); ++index) {
    if (results[index].manifest.ht_certificate_sha256 !=
        certificates[index]) {
      return false;
    }
  }
  return true;
}

struct ExpansionTables {
  Mod base_quotient = 0;
  std::vector<Mod> b;
  std::array<Mod, kEntries> linear{};
  std::vector<Mod> pair;
};

ExpansionTables build_expansion_tables(
    const Matrix& matrix, const ModFactorization& factorization) {
  ExpansionTables tables;
  tables.base_quotient =
      mul_mod(
          factorization.determinant,
          inverse_mod(static_cast<Mod>(kTwo22)));
  tables.b.resize(
      static_cast<std::size_t>(kEntries) * kEntries);
  for (int source = 0; source < kEntries; ++source) {
    const int source_column = source % kOrder;
    for (int target = 0; target < kEntries; ++target) {
      const int target_row = target / kOrder;
      const int target_column = target % kOrder;
      const int delta = -2 * matrix[target_row][target_column];
      const Mod delta_mod =
          delta >= 0 ? static_cast<Mod>(delta)
                     : static_cast<Mod>(kPrime64 + delta);
      tables.b[
          static_cast<std::size_t>(source) * kEntries + target] =
          mul_mod(
              delta_mod,
              factorization.inverse[
                  source_column * kOrder + target_row]);
    }
  }
  tables.pair.resize(
      static_cast<std::size_t>(kEntries) * kEntries);
  for (int first = 0; first < kEntries; ++first) {
    const Mod diagonal =
        tables.b[static_cast<std::size_t>(first) * kEntries + first];
    tables.linear[first] =
        mul_mod(tables.base_quotient, diagonal);
    for (int second = first + 1; second < kEntries; ++second) {
      const Mod other_diagonal =
          tables.b[
              static_cast<std::size_t>(second) * kEntries + second];
      const Mod determinant =
          sub_mod(
              mul_mod(diagonal, other_diagonal),
              mul_mod(
                  tables.b[
                      static_cast<std::size_t>(first) * kEntries +
                      second],
                  tables.b[
                      static_cast<std::size_t>(second) * kEntries +
                      first]));
      const Mod coefficient =
          mul_mod(tables.base_quotient, determinant);
      tables.pair[
          static_cast<std::size_t>(first) * kEntries + second] =
          coefficient;
      tables.pair[
          static_cast<std::size_t>(second) * kEntries + first] =
          coefficient;
    }
  }
  return tables;
}

Mod determinant3_coefficient(
    const ExpansionTables& tables, int first, int second, int third) {
  const auto b = [&](int row, int column) -> Mod {
    return tables.b[
        static_cast<std::size_t>(row) * kEntries + column];
  };
  const Mod minor0 =
      sub_mod(
          mul_mod(b(second, second), b(third, third)),
          mul_mod(b(second, third), b(third, second)));
  const Mod minor1 =
      sub_mod(
          mul_mod(b(second, first), b(third, third)),
          mul_mod(b(second, third), b(third, first)));
  const Mod minor2 =
      sub_mod(
          mul_mod(b(second, first), b(third, second)),
          mul_mod(b(second, second), b(third, first)));
  const Mod determinant =
      add_mod(
          sub_mod(
              mul_mod(b(first, first), minor0),
              mul_mod(b(first, second), minor1)),
          mul_mod(b(first, third), minor2));
  return mul_mod(tables.base_quotient, determinant);
}

Mod expansion_quotient(
    const ExpansionTables& tables, const CoordinateSet& coordinates) {
  Mod result = tables.base_quotient;
  for (int index = 0; index < coordinates.weight; ++index) {
    result = add_mod(result, tables.linear[coordinates.entries[index]]);
  }
  for (int first = 0; first < coordinates.weight; ++first) {
    for (int second = first + 1;
         second < coordinates.weight; ++second) {
      result =
          add_mod(
              result,
              tables.pair[
                  static_cast<std::size_t>(
                      coordinates.entries[first]) *
                      kEntries +
                  coordinates.entries[second]]);
    }
  }
  if (coordinates.weight == 3) {
    const int first = coordinates.entries[0];
    const int second = coordinates.entries[1];
    const int third = coordinates.entries[2];
    const bool distinct_rows =
        first / kOrder != second / kOrder &&
        first / kOrder != third / kOrder &&
        second / kOrder != third / kOrder;
    const bool distinct_columns =
        first % kOrder != second % kOrder &&
        first % kOrder != third % kOrder &&
        second % kOrder != third % kOrder;
    if (distinct_rows && distinct_columns) {
      result =
          add_mod(
              result,
              determinant3_coefficient(tables, first, second, third));
    }
  }
  return result;
}

Mod direct_quotient(
    const ExpansionTables& tables, const CoordinateSet& coordinates) {
  std::array<std::array<Mod, 3>, 3> kernel{};
  for (int row = 0; row < coordinates.weight; ++row) {
    for (int column = 0; column < coordinates.weight; ++column) {
      kernel[row][column] =
          tables.b[
              static_cast<std::size_t>(coordinates.entries[row]) *
                  kEntries +
              coordinates.entries[column]];
      if (row == column) {
        kernel[row][column] = add_mod(kernel[row][column], 1U);
      }
    }
  }
  Mod determinant = 0;
  if (coordinates.weight == 1) {
    determinant = kernel[0][0];
  } else if (coordinates.weight == 2) {
    determinant =
        sub_mod(
            mul_mod(kernel[0][0], kernel[1][1]),
            mul_mod(kernel[0][1], kernel[1][0]));
  } else if (coordinates.weight == 3) {
    const Mod minor0 =
        sub_mod(
            mul_mod(kernel[1][1], kernel[2][2]),
            mul_mod(kernel[1][2], kernel[2][1]));
    const Mod minor1 =
        sub_mod(
            mul_mod(kernel[1][0], kernel[2][2]),
            mul_mod(kernel[1][2], kernel[2][0]));
    const Mod minor2 =
        sub_mod(
            mul_mod(kernel[1][0], kernel[2][1]),
            mul_mod(kernel[1][1], kernel[2][0]));
    determinant =
        add_mod(
            sub_mod(
                mul_mod(kernel[0][0], minor0),
                mul_mod(kernel[0][1], minor1)),
            mul_mod(kernel[0][2], minor2));
  } else {
    throw std::runtime_error("direct quotient requires weight 1, 2, or 3");
  }
  return mul_mod(tables.base_quotient, determinant);
}

Matrix flipped_matrix(
    Matrix matrix, const CoordinateSet& coordinates) {
  for (int index = 0; index < coordinates.weight; ++index) {
    const int entry = coordinates.entries[index];
    matrix[entry / kOrder][entry % kOrder] *= -1;
  }
  return matrix;
}

std::string coordinate_string(const CoordinateSet& coordinates) {
  std::ostringstream output;
  for (int index = 0; index < coordinates.weight; ++index) {
    if (index != 0) output << '-';
    output << 'e' << std::setw(3) << std::setfill('0')
           << coordinates.entries[index];
  }
  return output.str();
}

std::uint64_t run_differential_checks(
    const Matrix& matrix, const ExpansionTables& tables,
    std::uint64_t seed) {
  std::vector<CoordinateSet> samples;
  samples.reserve(1297);
  for (int entry = 0; entry < kEntries; ++entry) {
    samples.push_back(CoordinateSet{1, {entry, -1, -1}});
  }
  for (int index = 0; index < 256; ++index) {
    int first = (index * 149 + 17) % kEntries;
    int second = (index * 263 + 311) % kEntries;
    if (second == first) second = (second + 1) % kEntries;
    if (first > second) std::swap(first, second);
    samples.push_back(CoordinateSet{2, {first, second, -1}});
  }
  std::mt19937_64 randomizer(seed);
  for (int index = 0; index < 512; ++index) {
    std::array<int, 3> entries{};
    do {
      entries[0] = static_cast<int>(randomizer() % kEntries);
      entries[1] = static_cast<int>(randomizer() % kEntries);
      entries[2] = static_cast<int>(randomizer() % kEntries);
    } while (
        entries[0] == entries[1] ||
        entries[0] == entries[2] ||
        entries[1] == entries[2]);
    std::sort(entries.begin(), entries.end());
    samples.push_back(CoordinateSet{3, entries});
  }

  for (const CoordinateSet& coordinates : samples) {
    const Mod expanded = expansion_quotient(tables, coordinates);
    const Mod direct = direct_quotient(tables, coordinates);
    if (expanded != direct) {
      throw std::runtime_error(
          "principal-minor expansion/direct determinant mismatch");
    }
    const std::int64_t predicted =
        recover_determinant_from_quotient(expanded);
    const Wide exact =
        exact_determinant(flipped_matrix(matrix, coordinates));
    if (exact != static_cast<Wide>(predicted)) {
      throw std::runtime_error(
          "modular rank-k prediction/Bareiss mismatch");
    }
  }
  return samples.size();
}

bool coordinate_less(
    const CoordinateSet& first, const CoordinateSet& second) {
  if (first.weight != second.weight) return first.weight < second.weight;
  return first.entries < second.entries;
}

void update_best(
    RadiusStatistics& statistics, std::int64_t signed_determinant,
    const CoordinateSet& coordinates) {
  const std::uint64_t score = absolute_i64(signed_determinant);
  if (score > statistics.best_absolute_determinant ||
      (score == statistics.best_absolute_determinant &&
       coordinate_less(coordinates, statistics.best_coordinates))) {
    statistics.best_absolute_determinant = score;
    statistics.best_signed_determinant = signed_determinant;
    statistics.best_coordinates = coordinates;
  }
}

RetainedArtifact retain_artifact(
    const Arguments& arguments, const ManifestEntry& manifest,
    const Matrix& base, const CoordinateSet& coordinates,
    std::int64_t signed_determinant) {
  const Matrix candidate = flipped_matrix(base, coordinates);
  const Wide independent = exact_determinant(candidate);
  if (independent != static_cast<Wide>(signed_determinant)) {
    throw std::runtime_error(
        "retained artifact failed independent Bareiss validation");
  }
  const std::string bytes = matrix_bytes(candidate);
  const std::string filename =
      "ht-" + manifest.ht_certificate_sha256.substr(0, 12) +
      "-r" + std::to_string(coordinates.weight) + "-" +
      coordinate_string(coordinates) + ".matrix.txt";
  const std::filesystem::path path =
      arguments.ties_directory / filename;
  atomic_write(path, bytes, "matrix");
  return RetainedArtifact{
      manifest.ht_certificate_sha256,
      coordinates,
      signed_determinant,
      path,
      sha256(bytes)};
}

void evaluate_candidate(
    const Arguments& arguments, const ManifestEntry& manifest,
    const Matrix& matrix, RadiusStatistics& statistics,
    std::vector<RetainedArtifact>& artifacts,
    const CoordinateSet& coordinates, Mod quotient) {
  const std::int64_t signed_determinant =
      recover_determinant_from_quotient(quotient);
  const std::uint64_t score = absolute_i64(signed_determinant);
  ++statistics.candidates;
  if (score == 0) ++statistics.singular;
  update_best(statistics, signed_determinant, coordinates);
  if (score == arguments.frontier) {
    ++statistics.frontier_ties;
    artifacts.push_back(
        retain_artifact(
            arguments, manifest, matrix, coordinates,
            signed_determinant));
  } else if (score > arguments.frontier) {
    ++statistics.strict_improvements;
    artifacts.push_back(
        retain_artifact(
            arguments, manifest, matrix, coordinates,
            signed_determinant));
  }
}

RepresentativeResult audit_representative(
    const Arguments& arguments, const ManifestEntry& manifest) {
  RepresentativeResult result;
  result.manifest = manifest;
  const std::string raw = read_file_bytes(manifest.path);
  const std::string raw_hash = sha256(raw);
  if (raw_hash != manifest.expected_raw_sha256) {
    throw std::runtime_error(
        "input raw SHA-256 mismatch for " + manifest.path.string());
  }
  const Matrix matrix = parse_matrix(raw, manifest.path);
  result.canonical_matrix_sha256 = sha256(matrix_bytes(matrix));
  const Wide base_exact = exact_determinant(matrix);
  if (base_exact < std::numeric_limits<std::int64_t>::min() ||
      base_exact > std::numeric_limits<std::int64_t>::max()) {
    throw std::runtime_error("base determinant does not fit int64");
  }
  result.base_signed_determinant =
      static_cast<std::int64_t>(base_exact);
  result.base_absolute_determinant =
      absolute_i64(result.base_signed_determinant);
  if (result.base_absolute_determinant != arguments.frontier) {
    throw std::runtime_error(
        "manifest input is not at the requested frontier: " +
        manifest.path.string());
  }
  const ModFactorization factorization = factorize_mod(matrix);
  if (!factorization.nonsingular ||
      factorization.determinant != signed_mod(base_exact)) {
    throw std::runtime_error(
        "modular factorization/Bareiss base mismatch");
  }
  const ExpansionTables tables =
      build_expansion_tables(matrix, factorization);
  if (recover_determinant_from_quotient(tables.base_quotient) !=
      result.base_signed_determinant) {
    throw std::runtime_error("base determinant quotient recovery mismatch");
  }

  const auto differential_started = Clock::now();
  std::uint64_t seed = UINT64_C(0xc7e31f07a524bd19);
  for (const unsigned char character :
       manifest.ht_certificate_sha256) {
    seed = seed * UINT64_C(0x100000001b3) ^ character;
  }
  result.bareiss_differential_samples =
      run_differential_checks(matrix, tables, seed);
  result.differential_seconds =
      std::chrono::duration<double>(
          Clock::now() - differential_started).count();
  if (arguments.differential_only) {
    result.complete = true;
    return result;
  }

  const auto enumeration_started = Clock::now();
  for (int first = 0; first < kEntries; ++first) {
    const CoordinateSet coordinates{1, {first, -1, -1}};
    evaluate_candidate(
        arguments, manifest, matrix, result.radius[0],
        result.artifacts, coordinates,
        add_mod(tables.base_quotient, tables.linear[first]));
  }
  for (int first = 0; first < kEntries; ++first) {
    for (int second = first + 1; second < kEntries; ++second) {
      Mod quotient =
          add_mod(tables.base_quotient, tables.linear[first]);
      quotient = add_mod(quotient, tables.linear[second]);
      quotient =
          add_mod(
              quotient,
              tables.pair[
                  static_cast<std::size_t>(first) * kEntries + second]);
      const CoordinateSet coordinates{
          2, {first, second, -1}};
      evaluate_candidate(
          arguments, manifest, matrix, result.radius[1],
          result.artifacts, coordinates, quotient);
    }
  }

  std::uint64_t triples_since_heartbeat = 0;
  for (int first = 0; first < kEntries; ++first) {
    for (int second = first + 1; second < kEntries; ++second) {
      Mod prefix =
          add_mod(tables.base_quotient, tables.linear[first]);
      prefix = add_mod(prefix, tables.linear[second]);
      prefix =
          add_mod(
              prefix,
              tables.pair[
                  static_cast<std::size_t>(first) * kEntries + second]);
      for (int third = second + 1; third < kEntries; ++third) {
        Mod quotient = add_mod(prefix, tables.linear[third]);
        quotient =
            add_mod(
                quotient,
                tables.pair[
                    static_cast<std::size_t>(first) * kEntries + third]);
        quotient =
            add_mod(
                quotient,
                tables.pair[
                    static_cast<std::size_t>(second) * kEntries + third]);
        const bool distinct_rows =
            first / kOrder != second / kOrder &&
            first / kOrder != third / kOrder &&
            second / kOrder != third / kOrder;
        const bool distinct_columns =
            first % kOrder != second % kOrder &&
            first % kOrder != third % kOrder &&
            second % kOrder != third % kOrder;
        if (distinct_rows && distinct_columns) {
          quotient =
              add_mod(
                  quotient,
                  determinant3_coefficient(
                      tables, first, second, third));
        }
        const CoordinateSet coordinates{
            3, {first, second, third}};
        evaluate_candidate(
            arguments, manifest, matrix, result.radius[2],
            result.artifacts, coordinates, quotient);
        ++triples_since_heartbeat;
        if (triples_since_heartbeat == UINT64_C(2097152)) {
          triples_since_heartbeat = 0;
          const double elapsed =
              std::chrono::duration<double>(
                  Clock::now() - enumeration_started).count();
          std::lock_guard<std::mutex> lock(output_mutex);
          std::cout
              << "{\"engine\":\"" << kEngine
              << "\",\"event\":\"heartbeat\",\"ht_certificate\":\""
              << manifest.ht_certificate_sha256
              << "\",\"radius3_completed\":"
              << result.radius[2].candidates
              << ",\"radius3_total\":" << kRadius3Count
              << ",\"frontier_ties\":"
              << result.radius[2].frontier_ties
              << ",\"strict_improvements\":"
              << result.radius[2].strict_improvements
              << ",\"elapsed_seconds\":" << std::fixed
              << std::setprecision(6) << elapsed << "}\n";
          std::cout.flush();
        }
      }
    }
  }
  result.enumeration_seconds =
      std::chrono::duration<double>(
          Clock::now() - enumeration_started).count();

  if (result.radius[0].candidates != kRadius1Count ||
      result.radius[1].candidates != kRadius2Count ||
      result.radius[2].candidates != kRadius3Count) {
    throw std::runtime_error("enumeration completion count mismatch");
  }
  result.complete = true;
  {
    std::lock_guard<std::mutex> lock(output_mutex);
    std::cout
        << "{\"engine\":\"" << kEngine
        << "\",\"event\":\"representative_finished\""
        << ",\"ht_certificate\":\""
        << manifest.ht_certificate_sha256
        << "\",\"candidates\":" << kPerRepresentativeCount
        << ",\"frontier_ties\":"
        << result.radius[0].frontier_ties +
               result.radius[1].frontier_ties +
               result.radius[2].frontier_ties
        << ",\"strict_improvements\":"
        << result.radius[0].strict_improvements +
               result.radius[1].strict_improvements +
               result.radius[2].strict_improvements
        << ",\"enumeration_seconds\":" << std::fixed
        << std::setprecision(6) << result.enumeration_seconds << "}\n";
    std::cout.flush();
  }
  return result;
}

void append_coordinates_json(
    std::ostringstream& output, const CoordinateSet& coordinates) {
  output << "{\"weight\":" << coordinates.weight
         << ",\"entry_indices\":[";
  for (int index = 0; index < coordinates.weight; ++index) {
    if (index != 0) output << ',';
    output << coordinates.entries[index];
  }
  output << "],\"row_column_zero_based\":[";
  for (int index = 0; index < coordinates.weight; ++index) {
    if (index != 0) output << ',';
    const int entry = coordinates.entries[index];
    output << '[' << entry / kOrder << ',' << entry % kOrder << ']';
  }
  output << "]}";
}

std::string build_report(
    const Arguments& arguments,
    const std::vector<RepresentativeResult>& results,
    std::string_view manifest_raw,
    std::string_view authority_raw,
    double elapsed_seconds) {
  std::uint64_t total_candidates = 0;
  std::uint64_t total_ties = 0;
  std::uint64_t total_improvements = 0;
  std::uint64_t total_singular = 0;
  std::uint64_t total_bareiss_samples = 0;
  std::vector<RetainedArtifact> all_artifacts;
  for (const RepresentativeResult& result : results) {
    total_bareiss_samples += result.bareiss_differential_samples;
    for (const RadiusStatistics& radius : result.radius) {
      total_candidates += radius.candidates;
      total_ties += radius.frontier_ties;
      total_improvements += radius.strict_improvements;
      total_singular += radius.singular;
    }
    all_artifacts.insert(
        all_artifacts.end(),
        result.artifacts.begin(), result.artifacts.end());
  }
  std::sort(
      all_artifacts.begin(), all_artifacts.end(),
      [](const RetainedArtifact& first,
         const RetainedArtifact& second) {
        if (first.ht_certificate_sha256 !=
            second.ht_certificate_sha256) {
          return first.ht_certificate_sha256 <
                 second.ht_certificate_sha256;
        }
        return coordinate_less(
            first.coordinates, second.coordinates);
      });
  const bool original_six_class_selection =
      is_original_six_class_manifest(results);
  const std::uint64_t authority_ht_class_count =
      parse_unsigned_json_field(authority_raw, "ht_class_count");

  std::ostringstream output;
  output << "{\"schema_version\":1"
         << ",\"challenge_id\":\"maxdet-23-v1\""
         << ",\"engine\":\"" << kEngine << "\""
         << ",\"complete\":true"
         << ",\"mode\":\""
         << (arguments.differential_only
                 ? "differential-only"
                 : original_six_class_selection
                       ? "six-uncovered-ht-radius-1-2-3"
                       : "manifest-radius-1-2-3")
         << "\""
         << ",\"frontier\":\"" << arguments.frontier << "\""
         << ",\"proof\":{\"modulus\":\"4294967291\""
         << ",\"modulus_description\":\"prime 2^32-5\""
         << ",\"determinant_divisor\":\"4194304\""
         << ",\"hadamard_quotient_bound\":\""
         << kHadamardQuotientBound << "\""
         << ",\"centered_recovery_unique\":true"
         << ",\"identity\":\"det(A+UV^T)=det(A)det(I+V^TA^-1U)\""
         << ",\"floating_point_ranking\":false"
         << ",\"bareiss_differential_samples\":"
         << total_bareiss_samples
         << ",\"retained_artifacts_bareiss_verified\":true}"
         << ",\"authority\":{\"path\":\""
         << json_escape(arguments.authority_audit.string()) << "\""
         << ",\"sha256\":\"" << sha256(authority_raw) << "\""
         << ",\"ht_class_count\":" << authority_ht_class_count
         << ",\"previously_completed_ht_certificates\":["
         << "\"4072ead420f0688da998d030955fecc419fc96e0494c6771fe6a92f74615305f\","
         << "\"db2cddf4b8f12da99a32b1563689ad693b89d0a9ae343a4aedde56361fdc5a81\""
         << "]}"
         << ",\"inputs\":{\"manifest_path\":\""
         << json_escape(arguments.manifest.string()) << "\""
         << ",\"manifest_sha256\":\"" << sha256(manifest_raw) << "\"}"
         << ",\"build\":{\"source_path\":\""
         << json_escape(arguments.source.string()) << "\""
         << ",\"source_sha256\":\""
         << sha256(read_file_bytes(arguments.source)) << "\""
         << ",\"binary_path\":\""
         << json_escape(arguments.binary.string()) << "\""
         << ",\"binary_sha256\":\""
         << sha256(read_file_bytes(arguments.binary)) << "\""
         << ",\"compile_command\":\""
         << json_escape(arguments.build_command) << "\"}"
         << ",\"coverage\":{\"representatives\":" << results.size()
         << ",\"original_six_class_selection\":"
         << (original_six_class_selection ? "true" : "false")
         << ",\"radius1_per_representative\":" << kRadius1Count
         << ",\"radius2_per_representative\":" << kRadius2Count
         << ",\"radius3_per_representative\":" << kRadius3Count
         << ",\"per_representative\":" << kPerRepresentativeCount
         << ",\"total_candidates\":" << total_candidates
         << ",\"expected_total_candidates\":"
         << (arguments.differential_only
                 ? 0
                 : kPerRepresentativeCount * results.size())
         << ",\"exact_completion_counts_match\":"
         << (arguments.differential_only ||
                     total_candidates ==
                         kPerRepresentativeCount * results.size()
                 ? "true"
                 : "false")
         << "}"
         << ",\"statistics\":{\"frontier_ties\":" << total_ties
         << ",\"strict_improvements\":" << total_improvements
         << ",\"singular_candidates\":" << total_singular
         << ",\"retained_artifact_count\":" << all_artifacts.size()
         << "}"
         << ",\"representatives\":[";
  for (std::size_t result_index = 0;
       result_index < results.size(); ++result_index) {
    if (result_index != 0) output << ',';
    const RepresentativeResult& result = results[result_index];
    output << "{\"ht_certificate_sha256\":\""
           << result.manifest.ht_certificate_sha256 << "\""
           << ",\"path\":\""
           << json_escape(result.manifest.path.string()) << "\""
           << ",\"raw_sha256\":\""
           << result.manifest.expected_raw_sha256 << "\""
           << ",\"canonical_matrix_sha256\":\""
           << result.canonical_matrix_sha256 << "\""
           << ",\"base_signed_determinant\":\""
           << result.base_signed_determinant << "\""
           << ",\"base_absolute_determinant\":\""
           << result.base_absolute_determinant << "\""
           << ",\"bareiss_differential_samples\":"
           << result.bareiss_differential_samples
           << ",\"differential_seconds\":" << std::fixed
           << std::setprecision(6) << result.differential_seconds
           << ",\"enumeration_seconds\":" << std::fixed
           << std::setprecision(6) << result.enumeration_seconds
           << ",\"complete\":"
           << (result.complete ? "true" : "false")
           << ",\"radii\":[";
    for (std::size_t radius_index = 0;
         radius_index < result.radius.size(); ++radius_index) {
      if (radius_index != 0) output << ',';
      const RadiusStatistics& radius = result.radius[radius_index];
      output << "{\"radius\":" << radius_index + 1
             << ",\"candidates\":" << radius.candidates
             << ",\"frontier_ties\":" << radius.frontier_ties
             << ",\"strict_improvements\":"
             << radius.strict_improvements
             << ",\"singular\":" << radius.singular
             << ",\"best_absolute_determinant\":\""
             << radius.best_absolute_determinant << "\""
             << ",\"best_signed_determinant\":\""
             << radius.best_signed_determinant << "\""
             << ",\"best_coordinates\":";
      append_coordinates_json(output, radius.best_coordinates);
      output << '}';
    }
    output << "]}";
  }
  output << "],\"retained_artifacts\":[";
  for (std::size_t index = 0; index < all_artifacts.size(); ++index) {
    if (index != 0) output << ',';
    const RetainedArtifact& artifact = all_artifacts[index];
    output << "{\"ht_certificate_sha256\":\""
           << artifact.ht_certificate_sha256 << "\""
           << ",\"coordinates\":";
    append_coordinates_json(output, artifact.coordinates);
    output << ",\"signed_determinant\":\""
           << artifact.signed_determinant << "\""
           << ",\"absolute_determinant\":\""
           << absolute_i64(artifact.signed_determinant) << "\""
           << ",\"path\":\"" << json_escape(artifact.path.string())
           << "\",\"raw_sha256\":\"" << artifact.raw_sha256
           << "\"}";
  }
  output << "]"
         << ",\"elapsed_seconds\":" << std::fixed
         << std::setprecision(6) << elapsed_seconds
         << ",\"claim_boundary\":\"Exact only for all entry-flip "
            "subsets of weights 1, 2, or 3 around the manifest-pinned "
            "representatives. H/HT certificates are inherited from the "
            "frozen authority audit; retained matrices require a new "
            "equivalence audit before any novelty claim.\""
         << "}\n";
  return output.str();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto started = Clock::now();
    const Arguments arguments = parse_arguments(argc, argv);
    if (sha256("abc") !=
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb"
        "410ff61f20015ad") {
      throw std::runtime_error("SHA-256 known-answer test failed");
    }
    const std::string manifest_raw =
        read_file_bytes(arguments.manifest);
    const std::vector<ManifestEntry> manifest =
        parse_manifest(manifest_raw);
    const std::string authority_raw =
        read_file_bytes(arguments.authority_audit);
    validate_authority(authority_raw, manifest);
    if (arguments.binary != std::filesystem::path(argv[0])) {
      throw std::runtime_error(
          "--binary must equal the invoked binary path");
    }
    if (std::filesystem::exists(arguments.output)) {
      throw std::runtime_error(
          "refusing to overwrite existing report: " +
          arguments.output.string());
    }
    if (std::filesystem::exists(arguments.ties_directory)) {
      if (!std::filesystem::is_directory(arguments.ties_directory)) {
        throw std::runtime_error(
            "ties output exists and is not a directory: " +
            arguments.ties_directory.string());
      }
      if (std::filesystem::directory_iterator(
              arguments.ties_directory) !=
          std::filesystem::directory_iterator()) {
        throw std::runtime_error(
            "refusing to use non-empty ties directory: " +
            arguments.ties_directory.string());
      }
    }
    std::filesystem::create_directories(arguments.ties_directory);

    std::vector<RepresentativeResult> results(manifest.size());
    std::atomic<std::size_t> next_index{0};
    std::atomic<bool> failed{false};
    std::mutex error_mutex;
    std::string error_message;
    const unsigned worker_count =
        std::min<unsigned>(
            arguments.threads,
            static_cast<unsigned>(manifest.size()));
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    for (unsigned worker = 0; worker < worker_count; ++worker) {
      workers.emplace_back([&]() {
        while (!failed.load(std::memory_order_relaxed)) {
          const std::size_t index =
              next_index.fetch_add(1, std::memory_order_relaxed);
          if (index >= manifest.size()) return;
          try {
            results[index] =
                audit_representative(arguments, manifest[index]);
          } catch (const std::exception& error) {
            failed.store(true, std::memory_order_relaxed);
            std::lock_guard<std::mutex> lock(error_mutex);
            if (error_message.empty()) error_message = error.what();
            return;
          }
        }
      });
    }
    for (std::thread& worker : workers) worker.join();
    if (failed.load(std::memory_order_relaxed)) {
      throw std::runtime_error(error_message);
    }

    const double elapsed_seconds =
        std::chrono::duration<double>(Clock::now() - started).count();
    const std::string report =
        build_report(
            arguments, results, manifest_raw, authority_raw,
            elapsed_seconds);
    atomic_write(arguments.output, report, "report");
    std::cout
        << "{\"engine\":\"" << kEngine
        << "\",\"event\":\"finished\",\"complete\":true"
        << ",\"report\":\""
        << json_escape(arguments.output.string()) << "\""
        << ",\"report_sha256\":\"" << sha256(report) << "\""
        << ",\"elapsed_seconds\":" << std::fixed
        << std::setprecision(6) << elapsed_seconds << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "frontier_portal_radius3: " << error.what() << '\n';
    return 1;
  }
}
