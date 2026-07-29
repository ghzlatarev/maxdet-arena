#include <algorithm>
#include <array>
#include <bit>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

constexpr int kOrder = 23;
constexpr int kEdgeCount = kOrder * (kOrder - 1) / 2;
constexpr int kFeatureCount = 1 + kEdgeCount;
constexpr int kFeatureWords = (kFeatureCount + 63) / 64;
constexpr int kSpanPrime = 3;
constexpr std::uint64_t kFrontierRoot =
    2'779'447'296'000'000ULL;
using Clock = std::chrono::steady_clock;
using Wide = __int128_t;
using UnsignedWide = __uint128_t;
using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using ExactMatrix =
    std::array<std::array<Wide, kOrder>, kOrder>;

struct Edge {
  int first = 0;
  int second = 0;
};

const std::array<Edge, kEdgeCount> kEdges = [] {
  std::array<Edge, kEdgeCount> result{};
  int next = 0;
  for (int first = 0; first < kOrder; ++first) {
    for (int second = first + 1; second < kOrder; ++second) {
      result[next++] = {first, second};
    }
  }
  return result;
}();

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
          std::string(label) + " must contain decimal digits only");
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

std::uint64_t determinant_modulo(
    const Matrix& matrix, int order, std::uint64_t prime) {
  std::array<std::array<std::uint64_t, kOrder>, kOrder> work{};
  for (int row = 0; row < order; ++row) {
    for (int column = 0; column < order; ++column) {
      const int value = matrix[row][column];
      if (value >= 0) {
        work[row][column] =
            static_cast<std::uint64_t>(value) % prime;
      } else {
        const std::uint64_t magnitude =
            static_cast<std::uint64_t>(-value) % prime;
        work[row][column] =
            magnitude == 0 ? 0 : prime - magnitude;
      }
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
      determinant = prime - determinant;
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

Wide exact_determinant(const Matrix& matrix, int order) {
  if (order <= 0 || order > kOrder) {
    throw std::runtime_error("invalid determinant order");
  }
  static constexpr std::array<std::uint64_t, 4> primes{
      2'147'483'647ULL,
      2'147'483'629ULL,
      2'147'483'587ULL,
      2'147'483'579ULL};

  UnsignedWide reconstructed = 0;
  UnsignedWide modulus = 1;
  for (const std::uint64_t prime : primes) {
    const std::uint64_t residue =
        determinant_modulo(matrix, order, prime);
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

  // A supported normalized Gram row has squared norm at most
  // 23^2 + 22*7^2 = 1607. Every full, principal, or cofactor
  // determinant used here is therefore below 1607^(23/2) < 10^37.
  // The CRT modulus exceeds 2*10^37, so symmetric reconstruction is
  // exact and unique.
  if (reconstructed > modulus / 2U) {
    return static_cast<Wide>(reconstructed) -
           static_cast<Wide>(modulus);
  }
  return static_cast<Wide>(reconstructed);
}

bool exact_positive_definite(const Matrix& matrix) {
  Matrix leading{};
  for (int order = 1; order <= kOrder; ++order) {
    for (int row = 0; row < order; ++row) {
      for (int column = 0; column < order; ++column) {
        leading[row][column] = matrix[row][column];
      }
    }
    if (exact_determinant(leading, order) <= 0) return false;
  }
  return true;
}

Wide integer_square_root(Wide value) {
  if (value < 0) {
    throw std::runtime_error("square root of a negative integer");
  }
  if (value == 0) return 0;
  const UnsignedWide input = static_cast<UnsignedWide>(value);
  unsigned bits = 0;
  for (UnsignedWide copy = input; copy != 0; copy >>= 1U) ++bits;
  UnsignedWide estimate =
      UnsignedWide{1} << ((bits + 1U) / 2U);
  for (;;) {
    const UnsignedWide next =
        (estimate + input / estimate) >> 1U;
    if (next >= estimate) return static_cast<Wide>(estimate);
    estimate = next;
  }
}

std::string read_file_bytes(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("cannot open file: " + path.string());
  }
  std::string bytes{
      std::istreambuf_iterator<char>(input),
      std::istreambuf_iterator<char>()};
  if (input.bad()) {
    throw std::runtime_error("cannot read file: " + path.string());
  }
  return bytes;
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

  static std::uint32_t rotate_right(
      std::uint32_t value, unsigned amount) {
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
  const std::string digits{text.substr(begin, position - begin)};
  std::size_t consumed = 0;
  const std::uint64_t value = std::stoull(digits, &consumed);
  if (consumed != digits.size()) {
    throw std::runtime_error(
        "invalid unsigned JSON field: " + std::string(key));
  }
  return value;
}

std::size_t matching_json_delimiter(
    std::string_view text,
    std::size_t start,
    char opening,
    char closing) {
  if (start >= text.size() || text[start] != opening) {
    throw std::runtime_error("invalid JSON delimiter start");
  }
  int depth = 0;
  bool in_string = false;
  bool escaped = false;
  for (std::size_t position = start;
       position < text.size();
       ++position) {
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

std::vector<std::string_view> json_array_objects(
    std::string_view text, std::string_view key) {
  std::size_t position = json_field_start(text, key);
  if (position >= text.size() || text[position] != '[') {
    throw std::runtime_error(
        "JSON field is not an array: " + std::string(key));
  }
  ++position;
  std::vector<std::string_view> objects;
  for (;;) {
    while (position < text.size() &&
           (std::isspace(
                static_cast<unsigned char>(text[position])) != 0 ||
            text[position] == ',')) {
      ++position;
    }
    if (position >= text.size()) {
      throw std::runtime_error("unterminated JSON object array");
    }
    if (text[position] == ']') return objects;
    if (text[position] != '{') {
      throw std::runtime_error(
          "JSON array must contain objects");
    }
    const std::size_t end =
        matching_json_delimiter(text, position, '{', '}');
    objects.push_back(text.substr(position, end - position + 1));
    position = end + 1;
  }
}

std::string json_escape(std::string_view text) {
  std::string result;
  result.push_back('"');
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
        if (character < 0x20) {
          const char digits[] = "0123456789abcdef";
          result += "\\u00";
          result.push_back(digits[character >> 4]);
          result.push_back(digits[character & 15]);
        } else {
          result.push_back(static_cast<char>(character));
        }
    }
  }
  result.push_back('"');
  return result;
}

Matrix parse_sign_matrix(
    std::string_view bytes, const std::filesystem::path& path) {
  std::istringstream input{std::string(bytes)};
  Matrix matrix{};
  for (auto& row : matrix) {
    for (int& value : row) {
      if (!(input >> value) || (value != -1 && value != 1)) {
        throw std::runtime_error(
            "matrix must contain exactly 23x23 signs: " + path.string());
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

Matrix gram(const Matrix& factor) {
  Matrix result{};
  for (int row = 0; row < kOrder; ++row) {
    for (int other = row; other < kOrder; ++other) {
      int value = 0;
      for (int column = 0; column < kOrder; ++column) {
        value += factor[row][column] * factor[other][column];
      }
      result[row][other] = value;
      result[other][row] = value;
    }
  }
  return result;
}

Matrix row_switch_normalize(Matrix factor) {
  const Matrix original_gram = gram(factor);
  for (int row = 1; row < kOrder; ++row) {
    const int value = original_gram[0][row];
    const bool already_normalized =
        value == -5 || value == -1 || value == 3 || value == 7;
    const bool normalized_after_switch =
        -value == -5 || -value == -1 || -value == 3 || -value == 7;
    if (!already_normalized && !normalized_after_switch) {
      throw std::runtime_error(
          "matrix Gram cannot be row-switched into {-5,-1,3,7}");
    }
    if (normalized_after_switch) {
      for (int& entry : factor[row]) entry = -entry;
    }
  }
  return factor;
}

void validate_normalized_gram(const Matrix& matrix) {
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      if (row == column) {
        if (matrix[row][column] != kOrder) {
          throw std::runtime_error("Gram diagonal is not 23");
        }
      } else if (
          matrix[row][column] != 7 &&
          matrix[row][column] != -5 &&
          matrix[row][column] != -1 &&
          matrix[row][column] != 3) {
        throw std::runtime_error(
            "Gram off-diagonal is outside {-5,-1,3,7}");
      }
    }
  }
}

struct SnapshotHit {
  Matrix gram{};
  Wide determinant = 0;
  Wide root = 0;
  std::size_t edge_count = 0;
};

std::set<std::pair<unsigned, unsigned>> parse_edge_array(
    std::string_view hit,
    std::string_view field,
    std::string_view count_field) {
  const std::size_t start = json_field_start(hit, field);
  if (start >= hit.size() || hit[start] != '[') {
    throw std::runtime_error(
        "snapshot " + std::string(field) + " field is not an array");
  }
  const std::size_t end =
      matching_json_delimiter(hit, start, '[', ']');
  const std::string_view edges =
      hit.substr(start, end - start + 1);
  std::set<std::pair<unsigned, unsigned>> unique;
  std::size_t position = 1;
  auto skip_space = [&]() {
    while (position < edges.size() &&
           std::isspace(
               static_cast<unsigned char>(edges[position])) != 0) {
      ++position;
    }
  };
  auto parse_vertex = [&]() {
    skip_space();
    if (position >= edges.size() ||
        std::isdigit(
            static_cast<unsigned char>(edges[position])) == 0) {
      throw std::runtime_error(
          "snapshot edge vertex must be unsigned");
    }
    unsigned value = 0;
    while (position < edges.size() &&
           std::isdigit(
               static_cast<unsigned char>(edges[position])) != 0) {
      value = value * 10U +
              static_cast<unsigned>(edges[position] - '0');
      if (value > 1000) {
        throw std::runtime_error(
            "snapshot edge vertex is unreasonably large");
      }
      ++position;
    }
    return value;
  };
  skip_space();
  while (position < edges.size() && edges[position] != ']') {
    if (edges[position++] != '[') {
      throw std::runtime_error(
          "snapshot edge entry must be a two-element array");
    }
    unsigned first = parse_vertex();
    skip_space();
    if (position >= edges.size() || edges[position++] != ',') {
      throw std::runtime_error(
          "snapshot edge entry is missing its comma");
    }
    unsigned second = parse_vertex();
    skip_space();
    if (position >= edges.size() || edges[position++] != ']') {
      throw std::runtime_error(
          "snapshot edge entry must contain exactly two vertices");
    }
    if (first == 0 || first > kOrder ||
        second == 0 || second > kOrder || first == second) {
      throw std::runtime_error(
          "snapshot edge vertex is invalid");
    }
    if (first > second) std::swap(first, second);
    if (!unique.emplace(first, second).second) {
      throw std::runtime_error(
          "snapshot contains a duplicate edge in " +
          std::string(field));
    }
    skip_space();
    if (position < edges.size() && edges[position] == ',') {
      ++position;
      skip_space();
      if (position < edges.size() && edges[position] == ']') {
        throw std::runtime_error(
            "snapshot edge list has a trailing comma");
      }
    } else if (
        position >= edges.size() || edges[position] != ']') {
      throw std::runtime_error(
          "snapshot edge entries must be comma-separated");
    }
  }
  if (position >= edges.size() || edges[position] != ']') {
    throw std::runtime_error(
        "snapshot edge list is unterminated");
  }
  ++position;
  skip_space();
  if (position != edges.size()) {
    throw std::runtime_error(
        "snapshot edge list contains trailing data");
  }
  if (json_unsigned_field(hit, count_field) != unique.size()) {
    throw std::runtime_error(
        "snapshot " + std::string(count_field) +
        " disagrees with " + std::string(field));
  }
  return unique;
}

SnapshotHit parse_snapshot_hit(
    std::string_view hit, std::string_view normalization,
    bool allow_subfrontier_research) {
  if (!json_bool_field(hit, "qualified") ||
      !json_bool_field(hit, "divisible_by_2_22") ||
      !json_bool_field(hit, "positive_definite")) {
    throw std::runtime_error(
        "selected snapshot hit is not qualified");
  }

  SnapshotHit result;
  result.root = positive_wide_decimal(
      json_string_field(hit, "square_root"), "square_root");
  result.determinant = positive_wide_decimal(
      json_string_field(hit, "determinant"), "determinant");
  constexpr UnsignedWide maximum =
      (static_cast<UnsignedWide>(~UnsignedWide{0})) >> 1U;
  if (static_cast<UnsignedWide>(result.root) >
          maximum / static_cast<UnsignedWide>(result.root) ||
      result.root * result.root != result.determinant) {
    throw std::runtime_error(
        "snapshot determinant is not square_root squared");
  }
  if (integer_square_root(result.determinant) != result.root) {
    throw std::runtime_error(
        "snapshot square_root failed exact recomputation");
  }
  if (!allow_subfrontier_research &&
      result.root <= static_cast<Wide>(kFrontierRoot)) {
    throw std::runtime_error(
        "snapshot square_root is not above the frontier");
  }
  if (result.root % (Wide{1} << 22) != 0) {
    throw std::runtime_error(
        "snapshot square_root is not divisible by 2^22");
  }

  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      result.gram[row][column] =
          row == column ? kOrder : -1;
    }
  }

  std::set<std::pair<unsigned, unsigned>> minus5_edges;
  std::set<std::pair<unsigned, unsigned>> plus3_edges;
  std::set<std::pair<unsigned, unsigned>> plus7_edges;
  if (normalization == "G=24I-J+4A") {
    plus3_edges = parse_edge_array(
        hit, "edges", "edge_count");
  } else if (
      normalization == "diag=23;offdiag={-5,-1,3}" ||
      normalization == "diag=23;offdiag={-5,-1,3,7}") {
    minus5_edges = parse_edge_array(
        hit, "minus5_edges", "minus5_count");
    plus3_edges = parse_edge_array(
        hit, "plus3_edges", "plus3_count");
    if (normalization == "diag=23;offdiag={-5,-1,3,7}") {
      plus7_edges = parse_edge_array(
          hit, "plus7_edges", "plus7_count");
    }
    std::vector<std::pair<unsigned, unsigned>> overlap;
    std::set_intersection(
        minus5_edges.begin(),
        minus5_edges.end(),
        plus3_edges.begin(),
        plus3_edges.end(),
        std::back_inserter(overlap));
    if (!overlap.empty()) {
      throw std::runtime_error(
          "snapshot -5 and +3 edge arrays overlap");
    }
    overlap.clear();
    std::set_intersection(
        minus5_edges.begin(),
        minus5_edges.end(),
        plus7_edges.begin(),
        plus7_edges.end(),
        std::back_inserter(overlap));
    if (!overlap.empty()) {
      throw std::runtime_error(
          "snapshot -5 and +7 edge arrays overlap");
    }
    overlap.clear();
    std::set_intersection(
        plus3_edges.begin(),
        plus3_edges.end(),
        plus7_edges.begin(),
        plus7_edges.end(),
        std::back_inserter(overlap));
    if (!overlap.empty()) {
      throw std::runtime_error(
          "snapshot +3 and +7 edge arrays overlap");
    }
  } else {
    throw std::runtime_error(
        "unsupported snapshot normalization");
  }
  for (const auto& [first, second] : minus5_edges) {
    result.gram[first - 1][second - 1] = -5;
    result.gram[second - 1][first - 1] = -5;
  }
  for (const auto& [first, second] : plus3_edges) {
    result.gram[first - 1][second - 1] = 3;
    result.gram[second - 1][first - 1] = 3;
  }
  for (const auto& [first, second] : plus7_edges) {
    result.gram[first - 1][second - 1] = 7;
    result.gram[second - 1][first - 1] = 7;
  }
  result.edge_count =
      minus5_edges.size() + plus3_edges.size() + plus7_edges.size();
  validate_normalized_gram(result.gram);
  if (exact_determinant(result.gram, kOrder) !=
      result.determinant) {
    throw std::runtime_error(
        "snapshot determinant failed exact recomputation");
  }
  if (!exact_positive_definite(result.gram)) {
    throw std::runtime_error(
        "snapshot Gram is not exactly positive definite");
  }
  return result;
}

struct ScaledInverse {
  ExactMatrix numerator{};
  Wide denominator = 0;
};

ScaledInverse exact_scaled_inverse(
    const Matrix& matrix, Wide determinant) {
  ExactMatrix adjugate{};
  Wide divisor = absolute(determinant);
  for (int removed_row = 0;
       removed_row < kOrder;
       ++removed_row) {
    for (int removed_column = removed_row;
         removed_column < kOrder;
         ++removed_column) {
      Matrix minor{};
      int target_row = 0;
      for (int source_row = 0; source_row < kOrder; ++source_row) {
        if (source_row == removed_row) continue;
        int target_column = 0;
        for (int source_column = 0;
             source_column < kOrder;
             ++source_column) {
          if (source_column == removed_column) continue;
          minor[target_row][target_column++] =
              matrix[source_row][source_column];
        }
        ++target_row;
      }
      Wide cofactor = exact_determinant(minor, kOrder - 1);
      if (((removed_row + removed_column) & 1) != 0) {
        cofactor = -cofactor;
      }
      adjugate[removed_row][removed_column] = cofactor;
      adjugate[removed_column][removed_row] = cofactor;
      divisor = wide_gcd(divisor, cofactor);
    }
  }
  if (divisor == 0) {
    throw std::runtime_error("singular Gram has no scaled inverse");
  }

  ScaledInverse result;
  result.denominator = determinant / divisor;
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      result.numerator[row][column] =
          adjugate[row][column] / divisor;
    }
  }

  // Independently check G * (adj/gcd) = (det/gcd) I.
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      Wide value = 0;
      for (int inner = 0; inner < kOrder; ++inner) {
        value += static_cast<Wide>(matrix[row][inner]) *
                 result.numerator[inner][column];
      }
      const Wide expected =
          row == column ? result.denominator : 0;
      if (value != expected) {
        throw std::runtime_error(
            "exact scaled inverse identity check failed");
      }
    }
  }
  return result;
}

struct TernaryVector {
  std::array<std::uint64_t, kFeatureWords> one{};
  std::array<std::uint64_t, kFeatureWords> two{};
};

std::uint64_t final_word_mask() {
  constexpr int remainder = kFeatureCount % 64;
  if constexpr (remainder == 0) {
    return ~std::uint64_t{0};
  } else {
    return (std::uint64_t{1} << remainder) - 1U;
  }
}

int ternary_get(const TernaryVector& vector, int coordinate) {
  const std::uint64_t mask =
      std::uint64_t{1} << (coordinate % 64);
  if ((vector.one[coordinate / 64] & mask) != 0U) return 1;
  if ((vector.two[coordinate / 64] & mask) != 0U) return 2;
  return 0;
}

void ternary_set(
    TernaryVector& vector, int coordinate, int value) {
  const int word = coordinate / 64;
  const std::uint64_t mask =
      std::uint64_t{1} << (coordinate % 64);
  vector.one[word] &= ~mask;
  vector.two[word] &= ~mask;
  if (value == 1) vector.one[word] |= mask;
  if (value == 2) vector.two[word] |= mask;
}

bool ternary_zero(const TernaryVector& vector) {
  for (int word = 0; word < kFeatureWords; ++word) {
    if ((vector.one[word] | vector.two[word]) != 0U) return false;
  }
  return true;
}

void ternary_scale_two(TernaryVector& vector) {
  std::swap(vector.one, vector.two);
}

void ternary_add_scaled(
    TernaryVector& left,
    const TernaryVector& right,
    int scale) {
  if (scale != 1 && scale != 2) {
    throw std::runtime_error("invalid ternary scale");
  }
  for (int word = 0; word < kFeatureWords; ++word) {
    const std::uint64_t left_one = left.one[word];
    const std::uint64_t left_two = left.two[word];
    const std::uint64_t right_one =
        scale == 1 ? right.one[word] : right.two[word];
    const std::uint64_t right_two =
        scale == 1 ? right.two[word] : right.one[word];
    const std::uint64_t left_zero = ~(left_one | left_two);
    const std::uint64_t right_zero = ~(right_one | right_two);
    left.one[word] =
        (left_zero & right_one) |
        (left_one & right_zero) |
        (left_two & right_two);
    left.two[word] =
        (left_zero & right_two) |
        (left_two & right_zero) |
        (left_one & right_one);
  }
  left.one.back() &= final_word_mask();
  left.two.back() &= final_word_mask();
}

int ternary_dot(
    const TernaryVector& left,
    const std::array<unsigned char, kFeatureCount>& right) {
  int result = 0;
  for (int coordinate = 0;
       coordinate < kFeatureCount;
       ++coordinate) {
    result += ternary_get(left, coordinate) * right[coordinate];
  }
  return result % kSpanPrime;
}

TernaryVector shell_feature(
    const std::array<int, kOrder>& signs) {
  TernaryVector feature;
  ternary_set(feature, 0, 1);
  for (int edge = 0; edge < kEdgeCount; ++edge) {
    const int value =
        signs[kEdges[edge].first] == signs[kEdges[edge].second]
            ? 1
            : 2;
    ternary_set(feature, edge + 1, value);
  }
  return feature;
}

TernaryVector gram_target_feature(const Matrix& gram_matrix) {
  TernaryVector target;
  ternary_set(target, 0, kOrder % kSpanPrime);
  for (int edge = 0; edge < kEdgeCount; ++edge) {
    int value =
        gram_matrix[kEdges[edge].first][kEdges[edge].second] %
        kSpanPrime;
    if (value < 0) value += kSpanPrime;
    ternary_set(target, edge + 1, value);
  }
  return target;
}

class TernarySpan {
 public:
  void insert(TernaryVector vector) {
    for (int pivot = 0; pivot < kFeatureCount; ++pivot) {
      const int value = ternary_get(vector, pivot);
      if (value == 0) continue;
      if (!used_[pivot]) {
        if (value == 2) ternary_scale_two(vector);
        basis_[pivot] = vector;
        used_[pivot] = true;
        ++rank_;
        return;
      }
      // basis pivot is one, so add -value times the basis row.
      ternary_add_scaled(vector, basis_[pivot], value == 1 ? 2 : 1);
    }
  }

  TernaryVector reduce(TernaryVector vector) const {
    for (int pivot = 0; pivot < kFeatureCount; ++pivot) {
      const int value = ternary_get(vector, pivot);
      if (value == 0 || !used_[pivot]) continue;
      ternary_add_scaled(vector, basis_[pivot], value == 1 ? 2 : 1);
    }
    return vector;
  }

  std::array<unsigned char, kFeatureCount> separating_certificate(
      const TernaryVector& target) const {
    const TernaryVector remainder = reduce(target);
    if (ternary_zero(remainder)) {
      throw std::runtime_error(
          "cannot certify a target already in the span");
    }
    int free_coordinate = -1;
    for (int coordinate = 0;
         coordinate < kFeatureCount;
         ++coordinate) {
      if (ternary_get(remainder, coordinate) != 0) {
        if (used_[coordinate]) {
          throw std::runtime_error(
              "reduced target retains a pivot coordinate");
        }
        free_coordinate = coordinate;
        break;
      }
    }
    if (free_coordinate < 0) {
      throw std::runtime_error("missing free certificate coordinate");
    }

    std::array<unsigned char, kFeatureCount> certificate{};
    certificate[free_coordinate] = 1;
    for (int pivot = kFeatureCount - 1; pivot >= 0; --pivot) {
      if (!used_[pivot]) continue;
      int sum = 0;
      for (int coordinate = pivot + 1;
           coordinate < kFeatureCount;
           ++coordinate) {
        sum += ternary_get(basis_[pivot], coordinate) *
               certificate[coordinate];
      }
      certificate[pivot] =
          static_cast<unsigned char>((kSpanPrime - sum % kSpanPrime) %
                                     kSpanPrime);
    }

    for (int pivot = 0; pivot < kFeatureCount; ++pivot) {
      if (used_[pivot] &&
          ternary_dot(basis_[pivot], certificate) != 0) {
        throw std::runtime_error(
            "constructed certificate is not in the annihilator");
      }
    }
    if (ternary_dot(target, certificate) == 0) {
      throw std::runtime_error(
          "constructed certificate does not separate target");
    }
    return certificate;
  }

  int rank() const { return rank_; }

 private:
  std::array<TernaryVector, kFeatureCount> basis_{};
  std::array<bool, kFeatureCount> used_{};
  int rank_ = 0;
};

template <typename Callback>
std::uint64_t enumerate_shell(
    const ScaledInverse& inverse, Callback callback) {
  std::array<int, kOrder> signs{};
  signs.fill(-1);
  signs[0] = 1;
  std::array<Wide, kOrder> weighted{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      weighted[row] +=
          inverse.numerator[row][column] * signs[column];
    }
  }
  Wide quadratic = 0;
  for (int index = 0; index < kOrder; ++index) {
    quadratic += signs[index] * weighted[index];
  }

  const std::uint64_t assignments =
      std::uint64_t{1} << (kOrder - 1);
  std::uint64_t shell_size = 0;
  for (std::uint64_t code = 0; code < assignments; ++code) {
    if (quadratic == inverse.denominator) {
      ++shell_size;
      callback(signs);
    }
    if (code + 1 == assignments) break;
    const int bit =
        static_cast<int>(std::countr_zero(code + 1)) + 1;
    const int old_sign = signs[bit];
    const int delta = -2 * old_sign;
    quadratic +=
        2 * delta * weighted[bit] +
        delta * delta * inverse.numerator[bit][bit];
    for (int row = 0; row < kOrder; ++row) {
      weighted[row] += delta * inverse.numerator[row][bit];
    }
    signs[bit] = -old_sign;
  }
  return shell_size;
}

struct FilterResult {
  std::string source;
  std::string source_sha256;
  std::string normalization;
  std::uint64_t hit_index = 0;
  bool has_hit_index = false;
  Matrix gram{};
  Wide determinant = 0;
  Wide inverse_denominator = 0;
  std::uint64_t shell_size = 0;
  int span_rank = 0;
  int augmented_rank = 0;
  bool rejected = false;
  bool factor_columns_in_shell = false;
  bool has_factor = false;
  std::vector<std::uint32_t> shell_sign_masks;
  std::array<unsigned char, kFeatureCount> certificate{};
  int certificate_target_dot = 0;
  double elapsed_seconds = 0;
};

FilterResult filter_gram(
    const Matrix& gram_matrix,
    std::string source,
    std::string source_sha256,
    std::string normalization,
    bool has_hit_index,
    std::uint64_t hit_index,
    const Matrix* known_factor,
    bool include_shell_vectors) {
  const auto started = Clock::now();
  validate_normalized_gram(gram_matrix);
  const Wide determinant =
      exact_determinant(gram_matrix, kOrder);
  if (determinant <= 0 || !exact_positive_definite(gram_matrix)) {
    throw std::runtime_error(
        "shell filter requires a positive-definite Gram");
  }
  const ScaledInverse inverse =
      exact_scaled_inverse(gram_matrix, determinant);

  TernarySpan span;
  std::vector<std::uint32_t> shell_sign_masks;
  const std::uint64_t shell_size = enumerate_shell(
      inverse,
      [&](const std::array<int, kOrder>& signs) {
        span.insert(shell_feature(signs));
        if (include_shell_vectors) {
          std::uint32_t mask = 0;
          for (int index = 0; index < kOrder; ++index) {
            if (signs[index] == 1) {
              mask |= std::uint32_t{1} << index;
            }
          }
          shell_sign_masks.push_back(mask);
        }
      });
  const TernaryVector target = gram_target_feature(gram_matrix);
  const TernaryVector reduced = span.reduce(target);
  const bool rejected = !ternary_zero(reduced);

  FilterResult result;
  result.source = std::move(source);
  result.source_sha256 = std::move(source_sha256);
  result.normalization = std::move(normalization);
  result.hit_index = hit_index;
  result.has_hit_index = has_hit_index;
  result.gram = gram_matrix;
  result.determinant = determinant;
  result.inverse_denominator = inverse.denominator;
  result.shell_size = shell_size;
  result.span_rank = span.rank();
  result.augmented_rank = span.rank() + (rejected ? 1 : 0);
  result.rejected = rejected;
  result.shell_sign_masks = std::move(shell_sign_masks);

  if (rejected) {
    result.certificate = span.separating_certificate(target);
    result.certificate_target_dot =
        ternary_dot(target, result.certificate);
    std::uint64_t verified_shell_size = 0;
    verified_shell_size = enumerate_shell(
        inverse,
        [&](const std::array<int, kOrder>& signs) {
          const TernaryVector feature = shell_feature(signs);
          if (ternary_dot(feature, result.certificate) != 0) {
            throw std::runtime_error(
                "certificate does not annihilate the exact shell");
          }
        });
    if (verified_shell_size != shell_size) {
      throw std::runtime_error(
          "shell count changed during certificate verification");
    }
  }

  if (known_factor != nullptr) {
    result.has_factor = true;
    if (gram(*known_factor) != gram_matrix) {
      throw std::runtime_error(
          "known factor does not reproduce its Gram");
    }
    result.factor_columns_in_shell = true;
    for (int column = 0; column < kOrder; ++column) {
      Wide quadratic = 0;
      for (int row = 0; row < kOrder; ++row) {
        for (int other = 0; other < kOrder; ++other) {
          quadratic +=
              static_cast<Wide>((*known_factor)[row][column]) *
              inverse.numerator[row][other] *
              (*known_factor)[other][column];
        }
      }
      if (quadratic != inverse.denominator) {
        result.factor_columns_in_shell = false;
      }
    }
    if (!result.factor_columns_in_shell || rejected) {
      throw std::runtime_error(
          "known factor failed the shell-span necessary condition");
    }
  }

  result.elapsed_seconds =
      std::chrono::duration<double>(Clock::now() - started).count();
  return result;
}

struct Arguments {
  std::filesystem::path snapshot;
  std::vector<std::uint64_t> hit_indices;
  bool all_hits = false;
  bool allow_subfrontier_research = false;
  bool help = false;
  bool include_shell_vectors = false;
  std::vector<std::filesystem::path> matrices;
  std::filesystem::path output;
};

std::uint64_t parse_unsigned(
    std::string_view text, std::string_view option) {
  const std::string copy(text);
  std::size_t consumed = 0;
  const unsigned long long value = std::stoull(copy, &consumed);
  if (consumed != copy.size()) {
    throw std::runtime_error(
        "invalid integer for " + std::string(option));
  }
  return static_cast<std::uint64_t>(value);
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    auto value = [&]() -> std::string {
      if (++index >= argc) throw std::runtime_error("missing option value");
      return argv[index];
    };
    if (option == "--snapshot") {
      arguments.snapshot = value();
    } else if (option == "--hit-index") {
      arguments.hit_indices.push_back(
          parse_unsigned(value(), option));
    } else if (option == "--all-hits") {
      arguments.all_hits = true;
    } else if (option == "--allow-subfrontier-research") {
      arguments.allow_subfrontier_research = true;
    } else if (option == "--help" || option == "-h") {
      arguments.help = true;
    } else if (option == "--matrix") {
      arguments.matrices.emplace_back(value());
    } else if (option == "--include-shell-vectors") {
      arguments.include_shell_vectors = true;
    } else if (option == "--output") {
      arguments.output = value();
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }

  if (arguments.help) return arguments;
  if (arguments.snapshot.empty() &&
      arguments.matrices.empty()) {
    throw std::runtime_error(
        "provide --snapshot with hit selection and/or --matrix");
  }
  if (arguments.snapshot.empty() &&
      (!arguments.hit_indices.empty() || arguments.all_hits)) {
    throw std::runtime_error(
        "--hit-index/--all-hits requires --snapshot");
  }
  if (!arguments.snapshot.empty() &&
      arguments.hit_indices.empty() && !arguments.all_hits) {
    throw std::runtime_error(
        "--snapshot requires repeatable --hit-index or --all-hits");
  }
  if (arguments.all_hits && !arguments.hit_indices.empty()) {
    throw std::runtime_error(
        "choose --all-hits or explicit --hit-index values");
  }
  if (arguments.allow_subfrontier_research &&
      arguments.snapshot.empty()) {
    throw std::runtime_error(
        "--allow-subfrontier-research requires --snapshot");
  }
  if (!arguments.output.empty()) {
    auto normalized = [](const std::filesystem::path& path) {
      std::error_code error;
      const std::filesystem::path result =
          std::filesystem::weakly_canonical(
              std::filesystem::absolute(path), error);
      if (!error) return result;
      return std::filesystem::absolute(path).lexically_normal();
    };
    const auto output = normalized(arguments.output);
    if (!arguments.snapshot.empty() &&
        output == normalized(arguments.snapshot)) {
      throw std::runtime_error(
          "--output must not overwrite --snapshot");
    }
    for (const auto& matrix : arguments.matrices) {
      if (output == normalized(matrix)) {
        throw std::runtime_error(
            "--output must not overwrite --matrix");
      }
    }
  }
  return arguments;
}

void print_usage(std::ostream& output) {
  output
      << "Usage: gram_shell_filter "
         "[--snapshot FILE (--hit-index N ... | --all-hits)] "
         "[--matrix FILE ...] [--include-shell-vectors] [--output FILE]\n"
      << "Exact order-23 sign-column shell/span obstruction for "
         "supported normalized Gram matrices.\n"
      << "  --snapshot FILE   validated Gram-search JSON checkpoint\n"
      << "  --hit-index N     zero-based hit selector; repeatable\n"
      << "  --all-hits        filter every hit in the snapshot\n"
      << "  --allow-subfrontier-research\n"
         "                    permit qualified square roots at or below the "
         "frontier; default remains strict-above only\n"
      << "  --matrix FILE     known 23x23 sign factor check; repeatable\n"
      << "  --include-shell-vectors\n"
         "                    include normalized sign-column masks in JSON\n"
      << "  --output FILE     atomic JSON output (default: stdout)\n";
}

void append_certificate(
    std::ostream& output,
    const FilterResult& result) {
  if (!result.rejected) {
    output << "null";
    return;
  }
  output << "{\"coordinate_system\":"
         << json_escape(
                "c0=common diagonal; remaining coordinates are "
                "lexicographic one-based edges");
  output << ",\"diagonal_coefficient\":"
         << static_cast<int>(result.certificate[0]);
  output << ",\"edge_coefficients\":[";
  bool first = true;
  for (int edge = 0; edge < kEdgeCount; ++edge) {
    const int coefficient = result.certificate[edge + 1];
    if (coefficient == 0) continue;
    if (!first) output << ',';
    first = false;
    output << '[' << kEdges[edge].first + 1
           << ',' << kEdges[edge].second + 1
           << ',' << coefficient << ']';
  }
  output << ']';
  output << ",\"modulus\":" << kSpanPrime;
  output << ",\"shell_dot\":0";
  output << ",\"target_dot\":" << result.certificate_target_dot;
  output << '}';
}

std::string results_json(
    const Arguments& arguments,
    const std::vector<FilterResult>& results,
    std::string_view snapshot_sha256,
    std::string_view normalization) {
  std::ostringstream output;
  output << "{\"claim_boundary\":"
         << json_escape(
                "A rejection is an exact necessary-condition obstruction; "
                "a pass does not construct or imply a sign factor.");
  output << ",\"complete\":true";
  output << ",\"engine\":\"gram-shell-filter\"";
  output << ",\"feature_dimension\":" << kFeatureCount;
  output << ",\"modulus\":" << kSpanPrime;
  output << ",\"normalization\":" << json_escape(normalization);
  output << ",\"allow_subfrontier_research\":"
         << (arguments.allow_subfrontier_research ? "true" : "false");
  output << ",\"results\":[";
  for (std::size_t index = 0; index < results.size(); ++index) {
    if (index != 0) output << ',';
    const FilterResult& result = results[index];
    output << "{\"augmented_rank\":" << result.augmented_rank;
    output << ",\"certificate\":";
    append_certificate(output, result);
    output << ",\"determinant\":"
           << json_escape(wide_to_string(result.determinant));
    output << ",\"elapsed_seconds\":" << result.elapsed_seconds;
    if (result.has_factor) {
      output << ",\"factor_columns_in_shell\":"
             << (result.factor_columns_in_shell ? "true" : "false");
    }
    if (result.has_hit_index) {
      output << ",\"hit_index\":" << result.hit_index;
    }
    output << ",\"inverse_denominator\":"
           << json_escape(wide_to_string(result.inverse_denominator));
    output << ",\"reason\":"
           << json_escape(
                  result.rejected
                      ? result.shell_size == 0
                            ? "empty-sign-column-shell"
                            : "gram-outside-mod-3-shell-outer-product-span"
                      : "no-shell-span-obstruction");
    output << ",\"rejected\":"
           << (result.rejected ? "true" : "false");
    output << ",\"shell_size\":" << result.shell_size;
    if (arguments.include_shell_vectors) {
      output << ",\"shell_sign_masks\":[";
      for (std::size_t shell_index = 0;
           shell_index < result.shell_sign_masks.size();
           ++shell_index) {
        if (shell_index != 0) output << ',';
        output << result.shell_sign_masks[shell_index];
      }
      output << ']';
    }
    output << ",\"source\":" << json_escape(result.source);
    if (!result.source_sha256.empty()) {
      output << ",\"source_sha256\":"
             << json_escape(result.source_sha256);
    }
    output << ",\"source_normalization\":"
           << json_escape(result.normalization);
    output << ",\"span_rank\":" << result.span_rank;
    output << '}';
  }
  output << ']';
  output << ",\"schema_version\":1";
  if (!arguments.snapshot.empty()) {
    output << ",\"snapshot\":"
           << json_escape(arguments.snapshot.string());
    output << ",\"snapshot_sha256\":"
           << json_escape(snapshot_sha256);
  }
  output << "}\n";
  return output.str();
}

void write_output(
    const std::filesystem::path& path, std::string_view contents) {
  if (path.empty()) {
    std::cout << contents;
    return;
  }
  if (!path.parent_path().empty()) {
    std::filesystem::create_directories(path.parent_path());
  }
  const std::filesystem::path temporary = path.string() + ".tmp";
  {
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("cannot create output");
    output.write(
        contents.data(), static_cast<std::streamsize>(contents.size()));
    output.flush();
    if (!output) throw std::runtime_error("cannot flush output");
  }
  std::filesystem::rename(temporary, path);
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    if (arguments.help) {
      print_usage(std::cout);
      return 0;
    }
    std::vector<FilterResult> results;
    std::string snapshot_sha256;
    std::string snapshot_normalization =
        "diag=23;offdiag={-5,-1,3,7}";

    if (!arguments.snapshot.empty()) {
      const std::string snapshot_bytes =
          read_file_bytes(arguments.snapshot);
      snapshot_sha256 = sha256(snapshot_bytes);
      const std::string engine =
          json_string_field(snapshot_bytes, "engine");
      if (engine != "gram-tabu" &&
          engine != "gram-gm-switch" &&
          engine != "gram-multilevel-tabu") {
        throw std::runtime_error(
            "unsupported snapshot engine");
      }
      if (json_string_field(snapshot_bytes, "challenge_id") !=
          "maxdet-23-v1") {
        throw std::runtime_error(
            "snapshot challenge_id must be maxdet-23-v1");
      }
      snapshot_normalization =
          json_string_field(snapshot_bytes, "normalization");
      const bool supported_normalization =
          engine == "gram-tabu" || engine == "gram-gm-switch"
              ? snapshot_normalization == "G=24I-J+4A"
              : snapshot_normalization ==
                        "diag=23;offdiag={-5,-1,3}" ||
                    snapshot_normalization ==
                        "diag=23;offdiag={-5,-1,3,7}";
      if (!supported_normalization) {
        throw std::runtime_error(
            "snapshot engine and normalization disagree");
      }
      const auto hits =
          json_array_objects(snapshot_bytes, "hits");
      std::vector<std::uint64_t> selected = arguments.hit_indices;
      if (arguments.all_hits) {
        selected.resize(hits.size());
        std::iota(selected.begin(), selected.end(), 0);
      }
      std::set<std::uint64_t> unique;
      for (const std::uint64_t hit_index : selected) {
        if (!unique.insert(hit_index).second) {
          throw std::runtime_error(
              "duplicate --hit-index selection");
        }
        if (hit_index >= hits.size()) {
          throw std::runtime_error(
              "snapshot hit index is out of range");
        }
        const SnapshotHit hit =
            parse_snapshot_hit(
                hits[hit_index], snapshot_normalization,
                arguments.allow_subfrontier_research);
        results.push_back(filter_gram(
            hit.gram,
            "snapshot-hit",
            snapshot_sha256,
            snapshot_normalization,
            true,
            hit_index,
            nullptr,
            arguments.include_shell_vectors));
        if (results.back().determinant != hit.determinant) {
          throw std::runtime_error(
              "filter determinant disagrees with validated snapshot");
        }
      }
    }

    for (const auto& path : arguments.matrices) {
      const std::string matrix_bytes = read_file_bytes(path);
      const Matrix factor =
          row_switch_normalize(parse_sign_matrix(matrix_bytes, path));
      const Matrix factor_gram = gram(factor);
      bool has_plus7 = false;
      for (int row = 0; row < kOrder; ++row) {
        for (int column = row + 1; column < kOrder; ++column) {
          has_plus7 =
              has_plus7 || factor_gram[row][column] == 7;
        }
      }
      results.push_back(filter_gram(
          factor_gram,
          path.string(),
          sha256(matrix_bytes),
          has_plus7
              ? "diag=23;offdiag={-5,-1,3,7}"
              : "diag=23;offdiag={-5,-1,3}",
          false,
          0,
          &factor,
          arguments.include_shell_vectors));
    }

    write_output(
        arguments.output,
        results_json(
            arguments,
            results,
            snapshot_sha256,
            snapshot_normalization));
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "gram_shell_filter: " << error.what() << '\n';
    return 2;
  }
}
