#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cctype>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <random>
#include <set>
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
using Wide = __int128_t;
using UnsignedWide = __uint128_t;
using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using Gram = std::array<std::array<int, kOrder>, kOrder>;
using Clock = std::chrono::steady_clock;

volatile std::sig_atomic_t stop_requested = 0;

void request_stop(int) { stop_requested = 1; }

Wide absolute(Wide value) { return value < 0 ? -value : value; }

Wide wide_gcd(Wide left, Wide right) {
  left = absolute(left);
  right = absolute(right);
  while (right != 0) {
    const Wide remainder = left % right;
    left = right;
    right = remainder;
  }
  return left;
}

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

class Rational {
 public:
  Rational() = default;
  explicit Rational(Wide numerator) : numerator_(numerator) {}
  Rational(Wide numerator, Wide denominator)
      : numerator_(numerator), denominator_(denominator) {
    normalize();
  }

  bool is_zero() const { return numerator_ == 0; }
  bool is_integer() const { return denominator_ == 1; }
  Wide numerator() const { return numerator_; }
  Wide denominator() const { return denominator_; }

  Rational operator-() const {
    return Rational(-numerator_, denominator_);
  }

  friend Rational operator+(const Rational& left,
                            const Rational& right) {
    const Wide divisor =
        wide_gcd(left.denominator_, right.denominator_);
    const Wide left_scale = right.denominator_ / divisor;
    const Wide right_scale = left.denominator_ / divisor;
    return Rational(
        left.numerator_ * left_scale +
            right.numerator_ * right_scale,
        left.denominator_ * left_scale);
  }

  friend Rational operator-(const Rational& left,
                            const Rational& right) {
    return left + (-right);
  }

  friend Rational operator*(const Rational& left,
                            const Rational& right) {
    const Wide first =
        wide_gcd(left.numerator_, right.denominator_);
    const Wide second =
        wide_gcd(right.numerator_, left.denominator_);
    return Rational(
        (left.numerator_ / first) *
            (right.numerator_ / second),
        (left.denominator_ / second) *
            (right.denominator_ / first));
  }

  friend Rational operator/(const Rational& left,
                            const Rational& right) {
    if (right.numerator_ == 0) {
      throw std::runtime_error("division by zero in exact elimination");
    }
    return left *
           Rational(right.denominator_, right.numerator_);
  }

  Rational& operator+=(const Rational& other) {
    *this = *this + other;
    return *this;
  }

  Rational& operator-=(const Rational& other) {
    *this = *this - other;
    return *this;
  }

 private:
  Wide numerator_ = 0;
  Wide denominator_ = 1;

  void normalize() {
    if (denominator_ == 0) {
      throw std::runtime_error("zero rational denominator");
    }
    if (numerator_ == 0) {
      denominator_ = 1;
      return;
    }
    if (denominator_ < 0) {
      numerator_ = -numerator_;
      denominator_ = -denominator_;
    }
    const Wide divisor = wide_gcd(numerator_, denominator_);
    numerator_ /= divisor;
    denominator_ /= divisor;
  }
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
      const Wide left = work[row][column];
      for (int inner = column + 1; inner < kOrder; ++inner) {
        const Wide numerator =
            work[row][inner] * pivot -
            left * work[column][inner];
        if (column != 0 && numerator % previous_pivot != 0) {
          throw std::runtime_error(
              "exact Bareiss division failed");
        }
        work[row][inner] =
            column == 0 ? numerator
                        : numerator / previous_pivot;
      }
      work[row][column] = 0;
    }
    previous_pivot = pivot;
  }
  return static_cast<Wide>(sign) *
         work[kOrder - 1][kOrder - 1];
}

std::uint64_t modular_power(
    std::uint64_t base,
    std::uint64_t exponent,
    std::uint64_t modulus) {
  std::uint64_t result = 1;
  while (exponent != 0) {
    if ((exponent & 1U) != 0U) {
      result = (result * base) % modulus;
    }
    base = (base * base) % modulus;
    exponent >>= 1U;
  }
  return result;
}

std::uint64_t gram_determinant_modulo(
    const Gram& matrix, int order, std::uint64_t prime) {
  std::array<std::array<std::uint64_t, kOrder>, kOrder> work{};
  for (int row = 0; row < order; ++row) {
    for (int column = 0; column < order; ++column) {
      const int value = matrix[row][column];
      work[row][column] =
          value >= 0
              ? static_cast<std::uint64_t>(value) % prime
              : prime -
                    static_cast<std::uint64_t>(-value) % prime;
    }
  }

  std::uint64_t determinant = 1;
  for (int column = 0; column < order; ++column) {
    int pivot_row = column;
    while (pivot_row < order &&
           work[pivot_row][column] == 0) {
      ++pivot_row;
    }
    if (pivot_row == order) return 0;
    if (pivot_row != column) {
      std::swap(work[pivot_row], work[column]);
      determinant =
          determinant == 0 ? 0 : prime - determinant;
    }
    const std::uint64_t pivot = work[column][column];
    determinant = (determinant * pivot) % prime;
    const std::uint64_t inverse =
        modular_power(pivot, prime - 2U, prime);
    for (int row = column + 1; row < order; ++row) {
      if (work[row][column] == 0) continue;
      const std::uint64_t factor =
          (work[row][column] * inverse) % prime;
      for (int inner = column + 1; inner < order; ++inner) {
        const std::uint64_t product =
            (factor * work[column][inner]) % prime;
        work[row][inner] =
            work[row][inner] >= product
                ? work[row][inner] - product
                : work[row][inner] + prime - product;
      }
      work[row][column] = 0;
    }
  }
  return determinant;
}

Wide exact_normalized_gram_determinant(
    const Gram& matrix, int order = kOrder) {
  if (order <= 0 || order > kOrder) {
    throw std::runtime_error(
        "invalid exact Gram determinant order");
  }
  for (int row = 0; row < order; ++row) {
    for (int column = 0; column < order; ++column) {
      const int expected_diagonal = row == column ? kOrder : 0;
      if ((row == column &&
           matrix[row][column] != expected_diagonal) ||
          (row != column && matrix[row][column] != -1 &&
           matrix[row][column] != 3)) {
        throw std::runtime_error(
            "CRT determinant requires G=24I-J+4A");
      }
    }
  }

  static constexpr std::array<std::uint64_t, 4> primes{
      1'000'000'007ULL,
      1'000'000'009ULL,
      1'000'000'033ULL,
      1'000'000'087ULL};
  // A normalized Gram row has squared norm at most
  // 23^2 + 22*3^2 = 727. Hadamard therefore bounds every
  // principal determinant by 727^(23/2) < 10^33, while the
  // following CRT modulus exceeds 10^36. Symmetric
  // reconstruction is consequently unique.
  UnsignedWide reconstructed = 0;
  UnsignedWide modulus = 1;
  for (const std::uint64_t prime : primes) {
    const std::uint64_t residue =
        gram_determinant_modulo(matrix, order, prime);
    const std::uint64_t current =
        static_cast<std::uint64_t>(reconstructed % prime);
    const std::uint64_t difference =
        residue >= current ? residue - current
                           : residue + prime - current;
    const std::uint64_t inverse = modular_power(
        static_cast<std::uint64_t>(modulus % prime),
        prime - 2U, prime);
    const std::uint64_t multiplier =
        (difference * inverse) % prime;
    reconstructed += modulus * multiplier;
    modulus *= prime;
  }
  if (reconstructed > modulus / 2U) {
    return static_cast<Wide>(reconstructed) -
           static_cast<Wide>(modulus);
  }
  return static_cast<Wide>(reconstructed);
}

bool exact_positive_definite(const Gram& matrix) {
  for (int order = 1; order <= kOrder; ++order) {
    if (exact_normalized_gram_determinant(matrix, order) <= 0) {
      return false;
    }
  }
  return true;
}

Matrix read_matrix(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error(
        "cannot open matrix: " + path.string());
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

Gram gram(const Matrix& matrix) {
  Gram result{};
  for (int row = 0; row < kOrder; ++row) {
    for (int other = row; other < kOrder; ++other) {
      int product = 0;
      for (int column = 0; column < kOrder; ++column) {
        product +=
            matrix[row][column] * matrix[other][column];
      }
      result[row][other] = product;
      result[other][row] = product;
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

Matrix sign_normalize(const Matrix& matrix) {
  Matrix normalized{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      normalized[row][column] =
          matrix[row][column] *
          matrix[row][0] *
          matrix[0][column] *
          matrix[0][0];
    }
  }
  return normalized;
}

Matrix column_canonicalize(const Matrix& matrix) {
  std::array<std::array<int, kOrder>, kOrder> columns{};
  for (int column = 0; column < kOrder; ++column) {
    const int sign = matrix[0][column];
    for (int row = 0; row < kOrder; ++row) {
      columns[column][row] = sign * matrix[row][column];
    }
  }
  std::sort(columns.begin(), columns.end());

  Matrix canonical{};
  for (int column = 0; column < kOrder; ++column) {
    for (int row = 0; row < kOrder; ++row) {
      canonical[row][column] = columns[column][row];
    }
  }
  return canonical;
}

class Sha256 {
 public:
  void update(std::string_view bytes) {
    for (const unsigned char byte : bytes) {
      buffer_[buffer_size_++] = byte;
      ++total_bytes_;
      if (buffer_size_ == buffer_.size()) {
        transform(buffer_.data());
        buffer_size_ = 0;
      }
    }
  }

  std::string hex_digest() {
    const std::uint64_t bit_length = total_bytes_ * 8U;
    buffer_[buffer_size_++] = 0x80U;
    if (buffer_size_ > 56U) {
      while (buffer_size_ < buffer_.size()) {
        buffer_[buffer_size_++] = 0;
      }
      transform(buffer_.data());
      buffer_size_ = 0;
    }
    while (buffer_size_ < 56U) buffer_[buffer_size_++] = 0;
    for (int shift = 56; shift >= 0; shift -= 8) {
      buffer_[buffer_size_++] =
          static_cast<std::uint8_t>(bit_length >> shift);
    }
    transform(buffer_.data());
    buffer_size_ = 0;

    std::ostringstream output;
    output << std::hex << std::setfill('0');
    for (const std::uint32_t word : state_) {
      output << std::setw(8) << word;
    }
    return output.str();
  }

 private:
  std::array<std::uint32_t, 8> state_{
      0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U,
      0xa54ff53aU, 0x510e527fU, 0x9b05688cU,
      0x1f83d9abU, 0x5be0cd19U};
  std::array<std::uint8_t, 64> buffer_{};
  std::size_t buffer_size_ = 0;
  std::uint64_t total_bytes_ = 0;

  static std::uint32_t rotate_right(std::uint32_t value,
                                    unsigned amount) {
    return (value >> amount) |
           (value << (32U - amount));
  }

  void transform(const std::uint8_t* block) {
    static constexpr std::array<std::uint32_t, 64> constants{
        0x428a2f98U, 0x71374491U, 0xb5c0fbcfU,
        0xe9b5dba5U, 0x3956c25bU, 0x59f111f1U,
        0x923f82a4U, 0xab1c5ed5U, 0xd807aa98U,
        0x12835b01U, 0x243185beU, 0x550c7dc3U,
        0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U,
        0xc19bf174U, 0xe49b69c1U, 0xefbe4786U,
        0x0fc19dc6U, 0x240ca1ccU, 0x2de92c6fU,
        0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
        0x983e5152U, 0xa831c66dU, 0xb00327c8U,
        0xbf597fc7U, 0xc6e00bf3U, 0xd5a79147U,
        0x06ca6351U, 0x14292967U, 0x27b70a85U,
        0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
        0x650a7354U, 0x766a0abbU, 0x81c2c92eU,
        0x92722c85U, 0xa2bfe8a1U, 0xa81a664bU,
        0xc24b8b70U, 0xc76c51a3U, 0xd192e819U,
        0xd6990624U, 0xf40e3585U, 0x106aa070U,
        0x19a4c116U, 0x1e376c08U, 0x2748774cU,
        0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU,
        0x5b9cca4fU, 0x682e6ff3U, 0x748f82eeU,
        0x78a5636fU, 0x84c87814U, 0x8cc70208U,
        0x90befffaU, 0xa4506cebU, 0xbef9a3f7U,
        0xc67178f2U};

    std::array<std::uint32_t, 64> words{};
    for (int index = 0; index < 16; ++index) {
      const int offset = 4 * index;
      words[index] =
          (static_cast<std::uint32_t>(block[offset]) << 24U) |
          (static_cast<std::uint32_t>(block[offset + 1]) << 16U) |
          (static_cast<std::uint32_t>(block[offset + 2]) << 8U) |
          static_cast<std::uint32_t>(block[offset + 3]);
    }
    for (int index = 16; index < 64; ++index) {
      const std::uint32_t first =
          rotate_right(words[index - 15], 7U) ^
          rotate_right(words[index - 15], 18U) ^
          (words[index - 15] >> 3U);
      const std::uint32_t second =
          rotate_right(words[index - 2], 17U) ^
          rotate_right(words[index - 2], 19U) ^
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
      const std::uint32_t sigma_one =
          rotate_right(e, 6U) ^ rotate_right(e, 11U) ^
          rotate_right(e, 25U);
      const std::uint32_t choice =
          (e & f) ^ ((~e) & g);
      const std::uint32_t temporary_one =
          h + sigma_one + choice + constants[index] + words[index];
      const std::uint32_t sigma_zero =
          rotate_right(a, 2U) ^ rotate_right(a, 13U) ^
          rotate_right(a, 22U);
      const std::uint32_t majority =
          (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t temporary_two =
          sigma_zero + majority;
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
};

std::string sha256(std::string_view bytes) {
  Sha256 hash;
  hash.update(bytes);
  return hash.hex_digest();
}

std::string read_file_bytes(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error(
        "cannot open file: " + path.string());
  }
  std::string bytes{
      std::istreambuf_iterator<char>(input),
      std::istreambuf_iterator<char>()};
  if (input.bad()) {
    throw std::runtime_error(
        "cannot read file: " + path.string());
  }
  return bytes;
}

std::size_t json_field_start(
    std::string_view text, std::string_view key) {
  const std::string token = "\"" + std::string(key) + "\"";
  const std::size_t key_position = text.find(token);
  if (key_position == std::string_view::npos) {
    throw std::runtime_error(
        "missing JSON field: " + std::string(key));
  }
  std::size_t position = key_position + token.size();
  while (position < text.size() &&
         std::isspace(
             static_cast<unsigned char>(text[position])) != 0) {
    ++position;
  }
  if (position >= text.size() || text[position] != ':') {
    throw std::runtime_error(
        "invalid JSON field separator: " + std::string(key));
  }
  ++position;
  while (position < text.size() &&
         std::isspace(
             static_cast<unsigned char>(text[position])) != 0) {
    ++position;
  }
  return position;
}

std::string json_string_field(
    std::string_view text, std::string_view key) {
  std::size_t position = json_field_start(text, key);
  if (position >= text.size() || text[position] != '"') {
    throw std::runtime_error(
        "JSON field is not a string: " + std::string(key));
  }
  ++position;
  std::string value;
  bool escaped = false;
  while (position < text.size()) {
    const char character = text[position++];
    if (escaped) {
      switch (character) {
        case '"':
        case '\\':
        case '/':
          value.push_back(character);
          break;
        case 'b':
          value.push_back('\b');
          break;
        case 'f':
          value.push_back('\f');
          break;
        case 'n':
          value.push_back('\n');
          break;
        case 'r':
          value.push_back('\r');
          break;
        case 't':
          value.push_back('\t');
          break;
        default:
          throw std::runtime_error(
              "unsupported JSON string escape in field: " +
              std::string(key));
      }
      escaped = false;
    } else if (character == '\\') {
      escaped = true;
    } else if (character == '"') {
      return value;
    } else {
      value.push_back(character);
    }
  }
  throw std::runtime_error(
      "unterminated JSON string field: " + std::string(key));
}

bool json_bool_field(
    std::string_view text, std::string_view key) {
  const std::size_t position = json_field_start(text, key);
  if (text.substr(position, 4) == "true") return true;
  if (text.substr(position, 5) == "false") return false;
  throw std::runtime_error(
      "JSON field is not boolean: " + std::string(key));
}

std::uint64_t json_unsigned_field(
    std::string_view text, std::string_view key) {
  std::size_t position = json_field_start(text, key);
  const std::size_t begin = position;
  while (position < text.size() &&
         text[position] >= '0' && text[position] <= '9') {
    ++position;
  }
  if (position == begin) {
    throw std::runtime_error(
        "JSON field is not unsigned: " + std::string(key));
  }
  std::size_t consumed = 0;
  const std::string digits{text.substr(begin, position - begin)};
  const std::uint64_t value = std::stoull(digits, &consumed);
  if (consumed != digits.size()) {
    throw std::runtime_error(
        "invalid unsigned JSON field: " + std::string(key));
  }
  return value;
}

std::size_t matching_json_delimiter(
    std::string_view text, std::size_t start,
    char opening, char closing) {
  if (start >= text.size() || text[start] != opening) {
    throw std::runtime_error("invalid JSON delimiter start");
  }
  int depth = 0;
  bool in_string = false;
  bool escaped = false;
  for (std::size_t position = start;
       position < text.size(); ++position) {
    const char character = text[position];
    if (in_string) {
      if (escaped) {
        escaped = false;
      } else if (character == '\\') {
        escaped = true;
      } else if (character == '"') {
        in_string = false;
      }
      continue;
    }
    if (character == '"') {
      in_string = true;
    } else if (character == opening) {
      ++depth;
    } else if (character == closing) {
      --depth;
      if (depth == 0) return position;
      if (depth < 0) break;
    }
  }
  throw std::runtime_error("unterminated JSON delimiter");
}

std::string_view json_array_object(
    std::string_view text, std::string_view key,
    std::uint64_t requested_index) {
  std::size_t position = json_field_start(text, key);
  if (position >= text.size() || text[position] != '[') {
    throw std::runtime_error(
        "JSON field is not an array: " + std::string(key));
  }
  ++position;
  std::uint64_t index = 0;
  for (;;) {
    while (position < text.size() &&
           (std::isspace(
                static_cast<unsigned char>(text[position])) != 0 ||
            text[position] == ',')) {
      ++position;
    }
    if (position >= text.size() || text[position] == ']') {
      throw std::runtime_error(
          "Gram snapshot hit index is out of range");
    }
    if (text[position] != '{') {
      throw std::runtime_error(
          "Gram snapshot hits must contain objects");
    }
    const std::size_t end =
        matching_json_delimiter(text, position, '{', '}');
    if (index == requested_index) {
      return text.substr(position, end - position + 1);
    }
    ++index;
    position = end + 1;
  }
}

Wide positive_wide_decimal(
    std::string_view text, std::string_view label) {
  if (text.empty()) {
    throw std::runtime_error(
        std::string(label) + " must be a positive decimal integer");
  }
  constexpr UnsignedWide maximum =
      (static_cast<UnsignedWide>(~UnsignedWide{0})) >> 1U;
  UnsignedWide value = 0;
  for (const char character : text) {
    if (character < '0' || character > '9') {
      throw std::runtime_error(
          std::string(label) +
          " must be a positive decimal integer");
    }
    const unsigned digit =
        static_cast<unsigned>(character - '0');
    if (value > (maximum - digit) / 10U) {
      throw std::runtime_error(
          std::string(label) + " exceeds signed 128-bit range");
    }
    value = value * 10U + digit;
  }
  if (value == 0) {
    throw std::runtime_error(
        std::string(label) + " must be positive");
  }
  return static_cast<Wide>(value);
}

std::string gram_bytes(const Gram& matrix) {
  std::ostringstream output;
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      if (column != 0) output << ' ';
      output << matrix[row][column];
    }
    output << '\n';
  }
  return output.str();
}

struct SnapshotTarget {
  Gram gram{};
  Wide root = 0;
  Wide determinant = 0;
  std::size_t edge_count = 0;
  std::string snapshot_sha256;
};

SnapshotTarget read_snapshot_target(
    const std::filesystem::path& path,
    std::uint64_t hit_index) {
  const std::string bytes = read_file_bytes(path);
  if (json_string_field(bytes, "engine") != "gram-tabu") {
    throw std::runtime_error(
        "Gram snapshot engine must be gram-tabu");
  }
  if (json_string_field(bytes, "normalization") !=
      "G=24I-J+4A") {
    throw std::runtime_error(
        "unsupported Gram snapshot normalization");
  }
  const std::string_view hit =
      json_array_object(bytes, "hits", hit_index);
  if (!json_bool_field(hit, "qualified") ||
      !json_bool_field(hit, "divisible_by_2_22") ||
      !json_bool_field(hit, "positive_definite")) {
    throw std::runtime_error(
        "selected Gram snapshot hit is not qualified");
  }

  SnapshotTarget target;
  target.snapshot_sha256 = sha256(bytes);
  target.root = positive_wide_decimal(
      json_string_field(hit, "square_root"), "square_root");
  target.determinant = positive_wide_decimal(
      json_string_field(hit, "determinant"), "determinant");
  if (target.root >
      static_cast<Wide>(
          (static_cast<UnsignedWide>(~UnsignedWide{0}) >> 1U)) /
          target.root ||
      target.root * target.root != target.determinant) {
    throw std::runtime_error(
        "snapshot determinant is not the square of square_root");
  }
  if (target.root % (Wide{1} << 22) != 0) {
    throw std::runtime_error(
        "snapshot square_root is not divisible by 2^22");
  }

  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      target.gram[row][column] =
          row == column ? kOrder : -1;
    }
  }

  const std::size_t edges_start =
      json_field_start(hit, "edges");
  if (edges_start >= hit.size() || hit[edges_start] != '[') {
    throw std::runtime_error("snapshot edges field is not an array");
  }
  const std::size_t edges_end =
      matching_json_delimiter(hit, edges_start, '[', ']');
  const std::string_view edges =
      hit.substr(edges_start, edges_end - edges_start + 1);
  std::vector<unsigned> numbers;
  for (std::size_t position = 0; position < edges.size();) {
    const unsigned char character =
        static_cast<unsigned char>(edges[position]);
    if (std::isdigit(character) != 0) {
      unsigned value = 0;
      while (position < edges.size() &&
             std::isdigit(static_cast<unsigned char>(
                 edges[position])) != 0) {
        value = value * 10U +
                static_cast<unsigned>(edges[position] - '0');
        ++position;
      }
      numbers.push_back(value);
    } else if (std::isspace(character) != 0 ||
               edges[position] == '[' ||
               edges[position] == ']' ||
               edges[position] == ',') {
      ++position;
    } else {
      throw std::runtime_error(
          "invalid character in snapshot edge list");
    }
  }
  if ((numbers.size() & 1U) != 0U) {
    throw std::runtime_error(
        "snapshot edge list has an odd number of vertices");
  }
  std::set<std::pair<unsigned, unsigned>> unique_edges;
  for (std::size_t index = 0; index < numbers.size(); index += 2) {
    unsigned first = numbers[index];
    unsigned second = numbers[index + 1];
    if (first == 0 || first > kOrder ||
        second == 0 || second > kOrder || first == second) {
      throw std::runtime_error(
          "snapshot edge vertex is out of range or repeated");
    }
    if (first > second) std::swap(first, second);
    if (!unique_edges.emplace(first, second).second) {
      throw std::runtime_error("snapshot contains a duplicate edge");
    }
    target.gram[first - 1][second - 1] = 3;
    target.gram[second - 1][first - 1] = 3;
  }
  target.edge_count = unique_edges.size();
  if (json_unsigned_field(hit, "edge_count") !=
      target.edge_count) {
    throw std::runtime_error(
        "snapshot edge_count disagrees with edges");
  }
  if (exact_normalized_gram_determinant(target.gram) !=
      target.determinant) {
    throw std::runtime_error(
        "snapshot Gram determinant failed exact recomputation");
  }
  if (!exact_positive_definite(target.gram)) {
    throw std::runtime_error(
        "snapshot Gram failed exact positive-definiteness");
  }
  return target;
}

std::string json_escape(std::string_view text) {
  std::string result;
  result.reserve(text.size() + 8);
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
          static constexpr char digits[] = "0123456789abcdef";
          result += "\\u00";
          result.push_back(digits[(character >> 4U) & 0xfU]);
          result.push_back(digits[character & 0xfU]);
        } else {
          result.push_back(static_cast<char>(character));
        }
    }
  }
  return result;
}

void write_all(int descriptor, std::string_view bytes) {
  std::size_t offset = 0;
  while (offset < bytes.size()) {
    const ssize_t written =
        ::write(descriptor, bytes.data() + offset,
                bytes.size() - offset);
    if (written < 0) {
      if (errno == EINTR) continue;
      throw std::runtime_error(
          "cannot write factor: " +
          std::string(std::strerror(errno)));
    }
    if (written == 0) {
      throw std::runtime_error(
          "short write while writing factor");
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
  const int result = ::fsync(descriptor);
  const int saved_errno = errno;
  ::close(descriptor);
  if (result != 0) {
    throw std::runtime_error(
        "cannot sync output directory: " +
        std::string(std::strerror(saved_errno)));
  }
}

void atomic_write(const std::filesystem::path& path,
                  std::string_view bytes,
                  std::uint64_t nonce) {
  const std::filesystem::path directory =
      path.parent_path().empty() ? std::filesystem::path(".")
                                 : path.parent_path();
  std::filesystem::create_directories(directory);

  if (std::filesystem::exists(path)) {
    std::ifstream existing(path, std::ios::binary);
    const std::string existing_bytes{
        std::istreambuf_iterator<char>(existing),
        std::istreambuf_iterator<char>()};
    if (existing_bytes == bytes) return;
    throw std::runtime_error(
        "refusing to replace nonmatching existing factor: " +
        path.string());
  }

  std::filesystem::path temporary;
  int descriptor = -1;
  for (std::uint64_t attempt = 0; attempt < 64; ++attempt) {
    temporary =
        directory /
        ("." + path.filename().string() + ".gram-decompose-" +
         std::to_string(static_cast<long long>(::getpid())) + "-" +
         std::to_string(nonce) + "-" + std::to_string(attempt) +
         ".tmp");
    int flags = O_WRONLY | O_CREAT | O_EXCL;
#ifdef O_CLOEXEC
    flags |= O_CLOEXEC;
#endif
    descriptor = ::open(temporary.c_str(), flags, 0600);
    if (descriptor >= 0) break;
    if (errno != EEXIST) {
      throw std::runtime_error(
          "cannot create factor temporary file: " +
          std::string(std::strerror(errno)));
    }
  }
  if (descriptor < 0) {
    throw std::runtime_error(
        "cannot allocate unique factor temporary file");
  }

  bool renamed = false;
  try {
    write_all(descriptor, bytes);
    if (::fsync(descriptor) != 0) {
      throw std::runtime_error(
          "cannot sync factor: " +
          std::string(std::strerror(errno)));
    }
    if (::close(descriptor) != 0) {
      descriptor = -1;
      throw std::runtime_error(
          "cannot close factor: " +
          std::string(std::strerror(errno)));
    }
    descriptor = -1;
    std::filesystem::rename(temporary, path);
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

struct Arguments {
  std::filesystem::path input;
  std::filesystem::path gram_snapshot;
  std::filesystem::path output_directory;
  std::filesystem::path log;
  std::vector<std::filesystem::path> check_factors;
  std::vector<std::filesystem::path> anchor_factors;
  std::uint64_t seed = 23;
  std::uint64_t hit_index = 0;
  double seconds = 60.0;
  double heartbeat_seconds = 15.0;
  double branch_probability = 0.3;
  std::uint64_t restarts = 0;
  std::uint64_t max_solutions = 10;
  bool shuffle_rows = true;
  bool anchor_input = true;
};

std::uint64_t strict_unsigned(std::string_view text,
                              std::string_view option) {
  if (text.empty() ||
      !std::all_of(
          text.begin(), text.end(),
          [](unsigned char character) {
            return character >= '0' && character <= '9';
          })) {
    throw std::runtime_error(
        std::string(option) +
        " must be a non-negative integer");
  }
  std::size_t consumed = 0;
  const std::uint64_t result =
      std::stoull(std::string(text), &consumed);
  if (consumed != text.size()) {
    throw std::runtime_error(
        std::string(option) +
        " must be a non-negative integer");
  }
  return result;
}

double strict_double(std::string_view text,
                     std::string_view option) {
  std::size_t consumed = 0;
  const double result =
      std::stod(std::string(text), &consumed);
  if (consumed != text.size() || !std::isfinite(result) ||
      result < 0.0) {
    throw std::runtime_error(
        std::string(option) +
        " must be finite and non-negative");
  }
  return result;
}

void print_usage(std::ostream& output) {
  output
      << "usage: gram_decompose (--input MATRIX | "
         "--gram-snapshot FILE) --output-dir DIR "
         "--log FILE [options]\n"
      << "\n"
      << "Randomly decompose an exact order-23 Gram target using "
         "framed row backtracking.\n"
      << "\n"
      << "options:\n"
      << "  --input MATRIX              derive G from a known sign "
         "factor\n"
      << "  --gram-snapshot FILE        read a qualified gram_tabu "
         "checkpoint hit\n"
      << "  --hit-index N               zero-based hit in a gram_tabu "
         "snapshot (default 0)\n"
      << "  --seed N                    PRNG seed (default 23)\n"
      << "  --seconds S                 wall limit; 0 needs finite "
         "--restarts (default 60)\n"
      << "  --restarts N                restart limit; 0 means "
         "until time limit (default 0)\n"
      << "  --max-solutions N           distinct factors to emit "
         "(default 10)\n"
      << "  --branch-probability P      probability of sampling a "
         "second child (default 0.3)\n"
      << "  --heartbeat-seconds S       JSONL heartbeat interval "
         "(default 15)\n"
      << "  --check-factor MATRIX       verify another exact factor; "
         "repeatable\n"
      << "  --anchor-factor MATRIX      add a known factor as a "
         "guaranteed search branch; repeatable\n"
      << "  --no-anchor                 disable the default input-factor "
         "anchor\n"
      << "  --no-shuffle-rows           retain target row order on "
         "every restart\n"
      << "  --help                      show this text\n";
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    auto value = [&]() -> std::string_view {
      ++index;
      if (index >= argc) {
        throw std::runtime_error(
            "missing value for " + option);
      }
      return argv[index];
    };
    if (option == "--input" || option == "--start") {
      arguments.input = value();
    } else if (option == "--gram-snapshot") {
      arguments.gram_snapshot = value();
    } else if (option == "--output-dir") {
      arguments.output_directory = value();
    } else if (option == "--log") {
      arguments.log = value();
    } else if (option == "--check-factor") {
      arguments.check_factors.emplace_back(value());
    } else if (option == "--anchor-factor") {
      arguments.anchor_factors.emplace_back(value());
    } else if (option == "--seed") {
      arguments.seed = strict_unsigned(value(), option);
    } else if (option == "--hit-index") {
      arguments.hit_index = strict_unsigned(value(), option);
    } else if (option == "--seconds") {
      arguments.seconds = strict_double(value(), option);
    } else if (option == "--heartbeat" ||
               option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds =
          strict_double(value(), option);
    } else if (option == "--branch-probability") {
      arguments.branch_probability =
          strict_double(value(), option);
    } else if (option == "--restarts") {
      arguments.restarts = strict_unsigned(value(), option);
    } else if (option == "--max-solutions") {
      arguments.max_solutions =
          strict_unsigned(value(), option);
    } else if (option == "--shuffle-rows") {
      arguments.shuffle_rows = true;
    } else if (option == "--no-shuffle-rows") {
      arguments.shuffle_rows = false;
    } else if (option == "--no-anchor") {
      arguments.anchor_input = false;
    } else if (option == "--help" || option == "-h") {
      print_usage(std::cout);
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }

  if (arguments.input.empty() == arguments.gram_snapshot.empty()) {
    throw std::runtime_error(
        "choose exactly one of --input or --gram-snapshot");
  }
  if (arguments.output_directory.empty()) {
    throw std::runtime_error("--output-dir is required");
  }
  if (arguments.log.empty()) {
    throw std::runtime_error("--log is required");
  }
  if (arguments.max_solutions == 0) {
    throw std::runtime_error("--max-solutions must be positive");
  }
  if (arguments.branch_probability > 1.0) {
    throw std::runtime_error(
        "--branch-probability must be at most 1");
  }
  if (arguments.seconds == 0.0 && arguments.restarts == 0) {
    throw std::runtime_error(
        "zero --seconds requires a finite --restarts limit");
  }

  const auto normalized = [](const std::filesystem::path& path) {
    return std::filesystem::absolute(path).lexically_normal();
  };
  const std::filesystem::path target_path =
      arguments.input.empty() ? arguments.gram_snapshot
                              : arguments.input;
  if (normalized(target_path) == normalized(arguments.log)) {
    throw std::runtime_error(
        "target input and --log must be distinct");
  }
  return arguments;
}

struct Frame {
  std::uint32_t pattern = 0;
  int width = 0;
};

int frame_sign(const Frame& frame, int row) {
  return ((frame.pattern >> row) & 1U) != 0U ? 1 : -1;
}

std::vector<std::vector<Rational>> rational_matrix(
    const std::vector<Frame>& frames,
    int rows,
    const std::vector<int>& columns) {
  std::vector<std::vector<Rational>> result(
      static_cast<std::size_t>(rows),
      std::vector<Rational>(columns.size()));
  for (int row = 0; row < rows; ++row) {
    for (std::size_t column = 0; column < columns.size();
         ++column) {
      result[static_cast<std::size_t>(row)][column] =
          Rational(frame_sign(
              frames[static_cast<std::size_t>(
                  columns[column])],
              row));
    }
  }
  return result;
}

int exact_rank(std::vector<std::vector<Rational>> matrix) {
  if (matrix.empty() || matrix.front().empty()) return 0;
  const int rows = static_cast<int>(matrix.size());
  const int columns =
      static_cast<int>(matrix.front().size());
  int pivot_row = 0;
  for (int column = 0;
       column < columns && pivot_row < rows; ++column) {
    int selected = pivot_row;
    while (selected < rows &&
           matrix[static_cast<std::size_t>(selected)]
                 [static_cast<std::size_t>(column)]
                     .is_zero()) {
      ++selected;
    }
    if (selected == rows) continue;
    if (selected != pivot_row) {
      std::swap(
          matrix[static_cast<std::size_t>(selected)],
          matrix[static_cast<std::size_t>(pivot_row)]);
    }
    const Rational pivot =
        matrix[static_cast<std::size_t>(pivot_row)]
              [static_cast<std::size_t>(column)];
    for (int inner = column; inner < columns; ++inner) {
      matrix[static_cast<std::size_t>(pivot_row)]
            [static_cast<std::size_t>(inner)] =
          matrix[static_cast<std::size_t>(pivot_row)]
                [static_cast<std::size_t>(inner)] /
          pivot;
    }
    for (int row = pivot_row + 1; row < rows; ++row) {
      const Rational factor =
          matrix[static_cast<std::size_t>(row)]
                [static_cast<std::size_t>(column)];
      if (factor.is_zero()) continue;
      for (int inner = column; inner < columns; ++inner) {
        matrix[static_cast<std::size_t>(row)]
              [static_cast<std::size_t>(inner)] -=
            factor *
            matrix[static_cast<std::size_t>(pivot_row)]
                  [static_cast<std::size_t>(inner)];
      }
    }
    ++pivot_row;
  }
  return pivot_row;
}

std::vector<int> select_basis(
    const std::vector<Frame>& frames,
    int rows) {
  std::vector<int> order(frames.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(
      order.begin(), order.end(),
      [&](int left, int right) {
        if (frames[static_cast<std::size_t>(left)].width !=
            frames[static_cast<std::size_t>(right)].width) {
          return frames[static_cast<std::size_t>(left)].width >
                 frames[static_cast<std::size_t>(right)].width;
        }
        return frames[static_cast<std::size_t>(left)].pattern <
               frames[static_cast<std::size_t>(right)].pattern;
      });

  std::vector<int> basis;
  int rank = 0;
  for (const int column : order) {
    std::vector<int> trial = basis;
    trial.push_back(column);
    const int trial_rank =
        exact_rank(rational_matrix(frames, rows, trial));
    if (trial_rank > rank) {
      basis.push_back(column);
      rank = trial_rank;
      if (rank == rows) break;
    }
  }
  if (static_cast<int>(basis.size()) != rows) {
    throw std::runtime_error(
        "partial factor lost full row rank");
  }
  return basis;
}

std::vector<std::vector<Rational>> exact_inverse(
    const std::vector<std::vector<Rational>>& matrix) {
  const int order = static_cast<int>(matrix.size());
  if (order == 0) return {};
  for (const auto& row : matrix) {
    if (static_cast<int>(row.size()) != order) {
      throw std::runtime_error(
          "exact inverse requires a square matrix");
    }
  }

  std::vector<std::vector<Rational>> augmented(
      static_cast<std::size_t>(order),
      std::vector<Rational>(
          static_cast<std::size_t>(2 * order)));
  for (int row = 0; row < order; ++row) {
    for (int column = 0; column < order; ++column) {
      augmented[static_cast<std::size_t>(row)]
               [static_cast<std::size_t>(column)] =
          matrix[static_cast<std::size_t>(row)]
                [static_cast<std::size_t>(column)];
      augmented[static_cast<std::size_t>(row)]
               [static_cast<std::size_t>(column + order)] =
          Rational(row == column ? 1 : 0);
    }
  }

  for (int column = 0; column < order; ++column) {
    int pivot_row = column;
    while (pivot_row < order &&
           augmented[static_cast<std::size_t>(pivot_row)]
                    [static_cast<std::size_t>(column)]
                        .is_zero()) {
      ++pivot_row;
    }
    if (pivot_row == order) {
      throw std::runtime_error(
          "selected frame basis is singular");
    }
    if (pivot_row != column) {
      std::swap(
          augmented[static_cast<std::size_t>(pivot_row)],
          augmented[static_cast<std::size_t>(column)]);
    }

    const Rational pivot =
        augmented[static_cast<std::size_t>(column)]
                 [static_cast<std::size_t>(column)];
    for (int inner = 0; inner < 2 * order; ++inner) {
      augmented[static_cast<std::size_t>(column)]
               [static_cast<std::size_t>(inner)] =
          augmented[static_cast<std::size_t>(column)]
                   [static_cast<std::size_t>(inner)] /
          pivot;
    }
    for (int row = 0; row < order; ++row) {
      if (row == column) continue;
      const Rational factor =
          augmented[static_cast<std::size_t>(row)]
                   [static_cast<std::size_t>(column)];
      if (factor.is_zero()) continue;
      for (int inner = 0; inner < 2 * order; ++inner) {
        augmented[static_cast<std::size_t>(row)]
                 [static_cast<std::size_t>(inner)] -=
            factor *
            augmented[static_cast<std::size_t>(column)]
                     [static_cast<std::size_t>(inner)];
      }
    }
  }

  std::vector<std::vector<Rational>> inverse(
      static_cast<std::size_t>(order),
      std::vector<Rational>(static_cast<std::size_t>(order)));
  for (int row = 0; row < order; ++row) {
    for (int column = 0; column < order; ++column) {
      inverse[static_cast<std::size_t>(row)]
             [static_cast<std::size_t>(column)] =
          augmented[static_cast<std::size_t>(row)]
                   [static_cast<std::size_t>(column + order)];
    }
  }
  return inverse;
}

struct IntegerLinearForm {
  std::vector<Wide> coefficients;
  Wide denominator = 1;
};

std::vector<IntegerLinearForm> integer_forms(
    const std::vector<std::vector<Rational>>& inverse) {
  std::vector<IntegerLinearForm> forms;
  forms.reserve(inverse.size());
  for (const auto& row : inverse) {
    Wide common_denominator = 1;
    for (const Rational& value : row) {
      const Wide divisor =
          wide_gcd(common_denominator, value.denominator());
      common_denominator =
          (common_denominator / divisor) * value.denominator();
    }
    IntegerLinearForm form;
    form.denominator = common_denominator;
    form.coefficients.reserve(row.size());
    for (const Rational& value : row) {
      form.coefficients.push_back(
          value.numerator() *
          (common_denominator / value.denominator()));
    }
    forms.push_back(std::move(form));
  }
  return forms;
}

struct Statistics {
  std::uint64_t restarts_started = 0;
  std::uint64_t restarts_completed = 0;
  std::uint64_t nodes = 0;
  std::uint64_t node_assignments = 0;
  std::uint64_t feasible_children = 0;
  std::uint64_t sampled_children = 0;
  std::uint64_t dead_ends = 0;
  std::uint64_t leaves = 0;
  std::uint64_t exact_factor_checks = 0;
  std::uint64_t duplicate_factors = 0;
  std::uint64_t factors_written = 0;
  int maximum_depth = 1;
};

class Search;

struct SampledChild {
  std::vector<int> values;
  bool follows_anchor = false;
};

struct SampledChildren {
  std::vector<SampledChild> children;
  std::uint64_t feasible = 0;
  bool complete = true;
};

class Search {
 public:
  Search(const Arguments& arguments,
         const Gram& target_gram,
         Wide target_score,
         std::vector<Matrix> anchors,
         std::set<std::string> known_factor_hashes,
         std::ofstream& log)
      : arguments_(arguments),
        target_gram_(target_gram),
        target_score_(target_score),
        anchors_(std::move(anchors)),
        log_(log),
        randomizer_(arguments.seed),
        started_(Clock::now()),
        deadline_(
            arguments.seconds == 0.0
                ? Clock::time_point::max()
                : started_ +
                      std::chrono::duration_cast<Clock::duration>(
                          std::chrono::duration<double>(
                              arguments.seconds))),
        next_heartbeat_(
            arguments.heartbeat_seconds == 0.0
                ? Clock::time_point::max()
                : started_ +
                      std::chrono::duration_cast<Clock::duration>(
                          std::chrono::duration<double>(
                              arguments.heartbeat_seconds))),
        factor_hashes_(std::move(known_factor_hashes)) {}

  Statistics& statistics() { return statistics_; }

  bool stopped() const {
    return stop_requested != 0 ||
           Clock::now() >= deadline_ ||
           statistics_.factors_written >=
               arguments_.max_solutions;
  }

  bool time_or_signal_stopped() const {
    return stop_requested != 0 || Clock::now() >= deadline_;
  }

  double elapsed() const {
    return std::chrono::duration<double>(
               Clock::now() - started_)
        .count();
  }

  void write_event(std::string_view event,
                   std::string_view extra = {}) {
    log_ << "{\"dead_ends\":" << statistics_.dead_ends
         << ",\"duplicate_factors\":"
         << statistics_.duplicate_factors
         << ",\"elapsed_seconds\":" << std::fixed
         << std::setprecision(3) << elapsed()
         << ",\"event\":\"" << json_escape(event) << "\""
         << ",\"exact_factor_checks\":"
         << statistics_.exact_factor_checks
         << ",\"factors_written\":"
         << statistics_.factors_written
         << ",\"feasible_children\":"
         << statistics_.feasible_children
         << ",\"leaves\":" << statistics_.leaves
         << ",\"maximum_depth\":" << statistics_.maximum_depth
         << ",\"node_assignments\":"
         << statistics_.node_assignments
         << ",\"nodes\":" << statistics_.nodes
         << ",\"restarts_completed\":"
         << statistics_.restarts_completed
         << ",\"restarts_started\":"
         << statistics_.restarts_started
         << ",\"sampled_children\":"
         << statistics_.sampled_children
         << ",\"seed\":" << arguments_.seed;
    if (!extra.empty()) log_ << ',' << extra;
    log_ << "}\n";
    log_.flush();
    if (!log_) {
      throw std::runtime_error(
          "cannot append Gram decomposition log");
    }
  }

  void maybe_heartbeat(bool force = false) {
    const auto now = Clock::now();
    if (!force && now < next_heartbeat_) return;
    write_event("heartbeat");
    if (arguments_.heartbeat_seconds == 0.0) {
      next_heartbeat_ = Clock::time_point::max();
    } else {
      next_heartbeat_ =
          now +
          std::chrono::duration_cast<Clock::duration>(
              std::chrono::duration<double>(
                  arguments_.heartbeat_seconds));
    }
  }

  void run() {
    while (!stopped() &&
           (arguments_.restarts == 0 ||
            statistics_.restarts_started <
                arguments_.restarts)) {
      std::array<int, kOrder> permutation{};
      std::iota(permutation.begin(), permutation.end(), 0);
      if (arguments_.shuffle_rows &&
          statistics_.restarts_started != 0) {
        std::shuffle(
            permutation.begin(), permutation.end(), randomizer_);
      }
      Gram permuted{};
      for (int row = 0; row < kOrder; ++row) {
        for (int column = 0; column < kOrder; ++column) {
          permuted[row][column] =
              target_gram_[permutation[row]]
                          [permutation[column]];
        }
      }

      ++statistics_.restarts_started;
      const std::uint64_t restart_number =
          statistics_.restarts_started;
      Matrix internal_anchor{};
      const Matrix* anchor = nullptr;
      if (!anchors_.empty()) {
        const Matrix& selected =
            anchors_[static_cast<std::size_t>(
                (restart_number - 1U) % anchors_.size())];
        for (int row = 0; row < kOrder; ++row) {
          for (int column = 0; column < kOrder; ++column) {
            internal_anchor[row][column] =
                selected[permutation[row]][column] *
                selected[permutation[0]][column];
          }
        }
        anchor = &internal_anchor;
      }
      const std::uint64_t before_nodes = statistics_.nodes;
      const std::uint64_t before_assignments =
          statistics_.node_assignments;
      std::vector<Frame> frames{{1U, kOrder}};
      explore(1, frames, permuted, permutation, anchor);
      const bool solution_limit =
          statistics_.factors_written >=
          arguments_.max_solutions;
      if (!time_or_signal_stopped() && !solution_limit) {
        ++statistics_.restarts_completed;
      }

      std::ostringstream extra;
      extra << "\"restart\":" << restart_number
            << ",\"restart_assignments\":"
            << statistics_.node_assignments - before_assignments
            << ",\"restart_nodes\":"
            << statistics_.nodes - before_nodes;
      write_event(
          time_or_signal_stopped()
              ? "restart_interrupted"
              : solution_limit ? "restart_solution_limit"
                               : "restart_complete",
          extra.str());
    }
  }

 private:
  const Arguments& arguments_;
  const Gram& target_gram_;
  Wide target_score_;
  std::vector<Matrix> anchors_;
  std::ofstream& log_;
  std::mt19937_64 randomizer_;
  Clock::time_point started_;
  Clock::time_point deadline_;
  Clock::time_point next_heartbeat_;
  Statistics statistics_;
  std::set<std::string> factor_hashes_;
  std::uint64_t output_nonce_ = 0;

  SampledChildren sample_children(
      int depth,
      const std::vector<Frame>& frames,
      const Gram& target,
      const Matrix* anchor) {
    SampledChildren result;
    const int variable_count =
        static_cast<int>(frames.size());
    const std::vector<int> basis =
        select_basis(frames, depth);
    std::vector<bool> is_basic(
        frames.size(), false);
    for (const int column : basis) {
      is_basic[static_cast<std::size_t>(column)] = true;
    }
    std::vector<int> nonbasic;
    for (int column = 0; column < variable_count; ++column) {
      if (!is_basic[static_cast<std::size_t>(column)]) {
        nonbasic.push_back(column);
      }
    }

    const auto basis_matrix =
        rational_matrix(frames, depth, basis);
    const auto forms =
        integer_forms(exact_inverse(basis_matrix));
    std::vector<Wide> constant(
        static_cast<std::size_t>(depth));
    for (int row = 0; row < depth; ++row) {
      int weighted_sum = 0;
      for (const Frame& frame : frames) {
        weighted_sum +=
            frame_sign(frame, row) * frame.width;
      }
      const int numerator =
          target[row][depth] + weighted_sum;
      if ((numerator & 1) != 0) {
        throw std::runtime_error(
            "Gram parity is incompatible with a sign factor");
      }
      constant[static_cast<std::size_t>(row)] =
          numerator / 2;
    }

    std::vector<int> anchor_values;
    if (anchor != nullptr) {
      anchor_values.assign(frames.size(), 0);
      std::vector<int> anchor_widths(frames.size(), 0);
      const std::uint32_t mask =
          (std::uint32_t{1} << depth) - 1U;
      for (int column = 0; column < kOrder; ++column) {
        std::uint32_t pattern = 0;
        for (int row = 0; row < depth; ++row) {
          if ((*anchor)[row][column] == 1) {
            pattern |= std::uint32_t{1} << row;
          }
        }
        pattern &= mask;
        const auto iterator = std::lower_bound(
            frames.begin(), frames.end(), pattern,
            [](const Frame& frame, std::uint32_t value) {
              return frame.pattern < value;
            });
        if (iterator == frames.end() ||
            iterator->pattern != pattern) {
          throw std::runtime_error(
              "anchor prefix is absent from the current framing");
        }
        const std::size_t index =
            static_cast<std::size_t>(iterator - frames.begin());
        ++anchor_widths[index];
        if ((*anchor)[depth][column] == 1) {
          ++anchor_values[index];
        }
      }
      for (std::size_t index = 0; index < frames.size(); ++index) {
        if (anchor_widths[index] != frames[index].width) {
          throw std::runtime_error(
              "anchor prefix width disagrees with the framing");
        }
      }
    }

    const bool sample_second =
        std::bernoulli_distribution(
            arguments_.branch_probability)(randomizer_);
    const int random_slots =
        anchor == nullptr ? 1 + (sample_second ? 1 : 0)
                          : (sample_second ? 1 : 0);
    std::vector<std::vector<int>> random_children;
    std::uint64_t random_feasible = 0;
    bool anchor_found = false;
    std::vector<int> nonbasic_values(nonbasic.size(), 0);
    bool more = true;
    std::uint64_t assignments_since_check = 0;
    while (more) {
      ++statistics_.node_assignments;
      ++assignments_since_check;

      std::vector<Wide> right = constant;
      for (std::size_t position = 0;
           position < nonbasic.size(); ++position) {
        const int column = nonbasic[position];
        const int value = nonbasic_values[position];
        if (value == 0) continue;
        for (int row = 0; row < depth; ++row) {
          right[static_cast<std::size_t>(row)] -=
              static_cast<Wide>(
                  frame_sign(
                      frames[static_cast<std::size_t>(column)],
                      row)) *
              value;
        }
      }

      std::vector<int> values(
          static_cast<std::size_t>(variable_count), 0);
      for (std::size_t position = 0;
           position < nonbasic.size(); ++position) {
        values[static_cast<std::size_t>(
            nonbasic[position])] =
            nonbasic_values[position];
      }
      bool feasible = true;
      for (int basic_position = 0;
           basic_position < depth; ++basic_position) {
        const IntegerLinearForm& form =
            forms[static_cast<std::size_t>(basic_position)];
        Wide numerator = 0;
        for (int row = 0; row < depth; ++row) {
          numerator +=
              form.coefficients[static_cast<std::size_t>(row)] *
              right[static_cast<std::size_t>(row)];
        }
        if (numerator % form.denominator != 0) {
          feasible = false;
          break;
        }
        const Wide value = numerator / form.denominator;
        const int column = basis[static_cast<std::size_t>(
            basic_position)];
        if (value < 0 ||
            value >
                frames[static_cast<std::size_t>(column)].width) {
          feasible = false;
          break;
        }
        values[static_cast<std::size_t>(column)] =
            static_cast<int>(value);
      }

      if (feasible) {
        ++result.feasible;
        if (anchor != nullptr && values == anchor_values) {
          anchor_found = true;
        } else if (random_slots != 0) {
          ++random_feasible;
          if (static_cast<int>(random_children.size()) <
              random_slots) {
            random_children.push_back(std::move(values));
          } else {
            std::uniform_int_distribution<std::uint64_t> choose(
                0, random_feasible - 1);
            const std::uint64_t selected = choose(randomizer_);
            if (selected <
                static_cast<std::uint64_t>(random_slots)) {
              random_children[static_cast<std::size_t>(selected)] =
                  std::move(values);
            }
          }
        }
      }

      if (assignments_since_check >= 1024U) {
        assignments_since_check = 0;
        maybe_heartbeat();
        if (time_or_signal_stopped()) {
          result.complete = false;
          break;
        }
      }

      if (nonbasic_values.empty()) {
        more = false;
      } else {
        std::size_t position = 0;
        while (position < nonbasic_values.size()) {
          const int column = nonbasic[position];
          if (nonbasic_values[position] <
              frames[static_cast<std::size_t>(column)].width) {
            ++nonbasic_values[position];
            break;
          }
          nonbasic_values[position] = 0;
          ++position;
        }
        if (position == nonbasic_values.size()) more = false;
      }
    }

    statistics_.feasible_children += result.feasible;
    if (!result.complete) return result;
    if (anchor != nullptr) {
      if (!anchor_found) {
        throw std::runtime_error(
            "known anchor was not found among feasible children");
      }
      result.children.push_back(
          SampledChild{std::move(anchor_values), true});
    }
    for (auto& values : random_children) {
      result.children.push_back(
          SampledChild{std::move(values), false});
    }
    statistics_.sampled_children += result.children.size();
    return result;
  }

  std::vector<Frame> refine(
      int new_row,
      const std::vector<Frame>& frames,
      const std::vector<int>& positive_counts) {
    std::vector<Frame> refined;
    refined.reserve(frames.size() * 2U);
    for (std::size_t index = 0; index < frames.size(); ++index) {
      const Frame frame = frames[index];
      const int positive = positive_counts[index];
      const int negative = frame.width - positive;
      if (negative > 0) {
        refined.push_back(Frame{frame.pattern, negative});
      }
      if (positive > 0) {
        refined.push_back(
            Frame{
                frame.pattern |
                    (std::uint32_t{1} << new_row),
                positive});
      }
    }
    std::sort(
        refined.begin(), refined.end(),
        [](const Frame& left, const Frame& right) {
          return left.pattern < right.pattern;
        });
    return refined;
  }

  void explore(
      int depth,
      const std::vector<Frame>& frames,
      const Gram& target,
      const std::array<int, kOrder>& permutation,
      const Matrix* anchor) {
    if (stopped()) return;
    ++statistics_.nodes;
    statistics_.maximum_depth =
        std::max(statistics_.maximum_depth, depth);
    maybe_heartbeat();

    if (depth == kOrder) {
      ++statistics_.leaves;
      process_leaf(frames, permutation);
      return;
    }

    SampledChildren children =
        sample_children(depth, frames, target, anchor);
    if (!children.complete || stopped()) return;
    if (children.feasible == 0) {
      ++statistics_.dead_ends;
      return;
    }
    for (const SampledChild& child : children.children) {
      if (stopped()) return;
      explore(
          depth + 1, refine(depth, frames, child.values),
          target, permutation,
          child.follows_anchor ? anchor : nullptr);
    }
  }

  void process_leaf(
      const std::vector<Frame>& frames,
      const std::array<int, kOrder>& permutation) {
    Matrix internal{};
    int column = 0;
    for (const Frame& frame : frames) {
      for (int copy = 0; copy < frame.width; ++copy) {
        if (column >= kOrder) {
          throw std::runtime_error(
              "leaf contains too many columns");
        }
        for (int row = 0; row < kOrder; ++row) {
          internal[row][column] = frame_sign(frame, row);
        }
        ++column;
      }
    }
    if (column != kOrder) {
      throw std::runtime_error(
          "leaf contains too few columns");
    }

    Matrix factor{};
    for (int internal_row = 0;
         internal_row < kOrder; ++internal_row) {
      factor[permutation[internal_row]] =
          internal[internal_row];
    }
    factor = column_canonicalize(factor);

    ++statistics_.exact_factor_checks;
    if (gram(factor) != target_gram_) {
      throw std::runtime_error(
          "generated factor failed exact Gram reconstruction");
    }
    const Wide determinant = exact_determinant(factor);
    if (absolute(determinant) != target_score_) {
      throw std::runtime_error(
          "generated factor determinant disagrees with target");
    }

    const std::string bytes = matrix_bytes(factor);
    const std::string canonical_hash = sha256(bytes);
    if (!factor_hashes_.insert(canonical_hash).second) {
      ++statistics_.duplicate_factors;
      return;
    }
    const std::string normalized_hash =
        sha256(matrix_bytes(sign_normalize(factor)));
    const std::filesystem::path output =
        arguments_.output_directory /
        ("solution-" + canonical_hash + ".matrix.txt");
    atomic_write(output, bytes, output_nonce_++);
    ++statistics_.factors_written;

    std::ostringstream extra;
    extra << "\"column_canonical_sha256\":\""
          << canonical_hash
          << "\",\"determinant\":\""
          << wide_to_string(determinant)
          << "\",\"normalized_sha256\":\""
          << normalized_hash
          << "\",\"output\":\""
          << json_escape(output.string())
          << "\",\"score\":\""
          << wide_to_string(absolute(determinant)) << "\"";
    write_event("solution", extra.str());
    std::cout << "solution "
              << statistics_.factors_written
              << " score=" << wide_to_string(absolute(determinant))
              << " sha256=" << canonical_hash
              << " output=" << output.string() << '\n'
              << std::flush;
  }
};

std::string factor_description(
    const Matrix& factor,
    const std::filesystem::path& path) {
  const Matrix canonical = column_canonicalize(factor);
  std::ostringstream output;
  output << "\"column_canonical_sha256\":\""
         << sha256(matrix_bytes(canonical))
         << "\",\"determinant\":\""
         << wide_to_string(exact_determinant(factor))
         << "\",\"factor\":\""
         << json_escape(path.string())
         << "\",\"normalized_sha256\":\""
         << sha256(matrix_bytes(sign_normalize(factor)))
         << "\",\"raw_sha256\":\""
         << sha256(read_file_bytes(path)) << "\"";
  return output.str();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    std::filesystem::create_directories(
        arguments.output_directory);
    if (!arguments.log.parent_path().empty()) {
      std::filesystem::create_directories(
          arguments.log.parent_path());
    }
    std::ofstream log(arguments.log, std::ios::app);
    if (!log) {
      throw std::runtime_error(
          "cannot open Gram decomposition log");
    }

    Matrix input{};
    bool has_input = !arguments.input.empty();
    Gram target_gram{};
    Wide target_score = 0;
    std::string snapshot_sha256;
    std::size_t snapshot_edge_count = 0;
    if (has_input) {
      input = read_matrix(arguments.input);
      target_score = absolute(exact_determinant(input));
      if (target_score == 0) {
        throw std::runtime_error(
            "input matrix must be nonsingular");
      }
      target_gram = gram(input);
    } else {
      const SnapshotTarget snapshot =
          read_snapshot_target(
              arguments.gram_snapshot, arguments.hit_index);
      target_gram = snapshot.gram;
      target_score = snapshot.root;
      snapshot_sha256 = snapshot.snapshot_sha256;
      snapshot_edge_count = snapshot.edge_count;
    }
    for (int index = 0; index < kOrder; ++index) {
      if (target_gram[index][index] != kOrder) {
        throw std::runtime_error(
            "input Gram diagonal is invalid");
      }
    }

    std::set<std::string> known_factor_hashes;
    if (has_input) {
      known_factor_hashes.insert(
          sha256(matrix_bytes(column_canonicalize(input))));
    }
    std::vector<std::pair<std::filesystem::path, Matrix>>
        checked_factors;
    for (const auto& path : arguments.check_factors) {
      const Matrix factor = read_matrix(path);
      if (gram(factor) != target_gram) {
        throw std::runtime_error(
            "check factor has a different Gram matrix: " +
            path.string());
      }
      if (absolute(exact_determinant(factor)) != target_score) {
        throw std::runtime_error(
            "check factor has a different determinant: " +
            path.string());
      }
      known_factor_hashes.insert(
          sha256(matrix_bytes(column_canonicalize(factor))));
      checked_factors.emplace_back(path, factor);
    }

    std::vector<Matrix> anchors;
    if (has_input && arguments.anchor_input) {
      anchors.push_back(input);
    }
    for (const auto& path : arguments.anchor_factors) {
      const Matrix factor = read_matrix(path);
      if (gram(factor) != target_gram ||
          absolute(exact_determinant(factor)) != target_score) {
        throw std::runtime_error(
            "anchor factor does not match the target Gram and "
            "determinant: " +
            path.string());
      }
      known_factor_hashes.insert(
          sha256(matrix_bytes(column_canonicalize(factor))));
      anchors.push_back(factor);
    }

    Search search(
        arguments, target_gram, target_score,
        std::move(anchors), std::move(known_factor_hashes), log);
    std::ostringstream start_extra;
    start_extra
        << "\"branch_probability\":"
        << std::setprecision(17)
        << arguments.branch_probability
        << ",\"max_solutions\":"
        << arguments.max_solutions
        << ",\"anchor_count\":"
        << (has_input && arguments.anchor_input ? 1U : 0U) +
               arguments.anchor_factors.size()
        << ",\"restart_limit\":" << arguments.restarts
        << ",\"row_shuffle\":"
        << (arguments.shuffle_rows ? "true" : "false")
        << ",\"score\":\"" << wide_to_string(target_score)
        << "\",\"seconds\":" << arguments.seconds;
    if (has_input) {
      start_extra
          << ",\"target_mode\":\"factor\""
          << ",\"input\":\""
          << json_escape(arguments.input.string())
          << "\",\"input_column_canonical_sha256\":\""
          << sha256(matrix_bytes(column_canonicalize(input)))
          << "\",\"input_normalized_sha256\":\""
          << sha256(matrix_bytes(sign_normalize(input)))
          << "\",\"input_raw_sha256\":\""
          << sha256(read_file_bytes(arguments.input)) << "\"";
    } else {
      start_extra
          << ",\"target_mode\":\"gram_snapshot\""
          << ",\"gram_snapshot\":\""
          << json_escape(arguments.gram_snapshot.string())
          << "\",\"gram_snapshot_sha256\":\""
          << snapshot_sha256
          << "\",\"gram_sha256\":\""
          << sha256(gram_bytes(target_gram))
          << "\",\"hit_index\":" << arguments.hit_index
          << ",\"edge_count\":" << snapshot_edge_count;
    }
    search.write_event("start", start_extra.str());

    for (const auto& [path, factor] : checked_factors) {
      ++search.statistics().exact_factor_checks;
      search.write_event(
          "check_factor", factor_description(factor, path));
    }

    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);
    search.run();

    std::string termination;
    if (stop_requested != 0) {
      termination = "signal";
    } else if (search.statistics().factors_written >=
               arguments.max_solutions) {
      termination = "max_solutions";
    } else if (search.time_or_signal_stopped()) {
      termination = "time_limit";
    } else if (arguments.restarts != 0 &&
               search.statistics().restarts_started >=
                   arguments.restarts) {
      termination = "restart_limit";
    } else {
      termination = "time_limit";
    }
    search.write_event(
        "finished",
        "\"termination\":\"" + json_escape(termination) + "\"");
    std::cout << "finished factors="
              << search.statistics().factors_written
              << " restarts="
              << search.statistics().restarts_started
              << " nodes=" << search.statistics().nodes
              << " assignments="
              << search.statistics().node_assignments
              << " termination=" << termination << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "gram_decompose: " << error.what() << '\n';
    return 2;
  }
}
