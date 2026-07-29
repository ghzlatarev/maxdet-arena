#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
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
constexpr int kRequiredPower = 22;
constexpr std::uint64_t kFrontierRoot = 2'779'447'296'000'000ULL;
constexpr std::string_view kThreeLevelNormalization =
    "diag=23;offdiag={-5,-1,3}";
constexpr std::string_view kFourLevelNormalization =
    "diag=23;offdiag={-5,-1,3,7}";
using Clock = std::chrono::steady_clock;
using Exact = __int128_t;
using UnsignedExact = __uint128_t;
using SignMatrix = std::array<std::array<int, kOrder>, kOrder>;
using Gram = std::array<std::array<int, kOrder>, kOrder>;
using FloatingMatrix =
    std::array<std::array<long double, kOrder>, kOrder>;
using Levels = std::array<std::int8_t, kEdgeCount>;
using LevelKey = std::array<std::uint64_t, 8>;

volatile std::sig_atomic_t stop_requested = 0;

void request_stop(int) { stop_requested = 1; }

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

Exact frontier_determinant() {
  const Exact root = static_cast<Exact>(kFrontierRoot);
  return root * root;
}

std::string exact_string(Exact value) {
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

constexpr std::array<std::uint64_t, 4> kCrtPrimes = {
    2'147'483'647ULL,
    2'147'483'629ULL,
    2'147'483'587ULL,
    2'147'483'579ULL,
};

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

bool is_prime(std::uint64_t value) {
  if (value < 2) return false;
  if (value % 2 == 0) return value == 2;
  for (std::uint64_t divisor = 3;
       divisor * divisor <= value;
       divisor += 2) {
    if (value % divisor == 0) return false;
  }
  return true;
}

std::uint64_t determinant_modulo(
    const Gram& matrix, int order, std::uint64_t prime) {
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
    while (pivot_row < order && work[pivot_row][column] == 0) {
      ++pivot_row;
    }
    if (pivot_row == order) return 0;
    if (pivot_row != column) {
      std::swap(work[pivot_row], work[column]);
      determinant = determinant == 0 ? 0 : prime - determinant;
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

Exact exact_determinant(const Gram& matrix, int order = kOrder) {
  if (order <= 0 || order > kOrder) {
    throw std::runtime_error("invalid exact determinant order");
  }
  UnsignedExact reconstructed = 0;
  UnsignedExact modulus = 1;
  for (const std::uint64_t prime : kCrtPrimes) {
    const std::uint64_t residue =
        determinant_modulo(matrix, order, prime);
    const std::uint64_t current =
        static_cast<std::uint64_t>(reconstructed % prime);
    const std::uint64_t difference =
        residue >= current
            ? residue - current
            : residue + prime - current;
    const std::uint64_t inverse = modular_power(
        static_cast<std::uint64_t>(modulus % prime),
        prime - 2U,
        prime);
    const std::uint64_t multiplier =
        (difference * inverse) % prime;
    reconstructed += modulus * multiplier;
    modulus *= prime;
  }
  if (reconstructed > modulus / 2U) {
    return static_cast<Exact>(reconstructed) -
           static_cast<Exact>(modulus);
  }
  return static_cast<Exact>(reconstructed);
}

bool exact_positive_definite(const Gram& matrix) {
  for (int order = 1; order <= kOrder; ++order) {
    if (exact_determinant(matrix, order) <= 0) return false;
  }
  return true;
}

Exact integer_square_root(Exact value) {
  if (value < 0) {
    throw std::runtime_error("square root of a negative integer");
  }
  if (value == 0) return 0;
  const UnsignedExact input = static_cast<UnsignedExact>(value);
  unsigned bits = 0;
  for (UnsignedExact copy = input; copy != 0; copy >>= 1U) ++bits;
  UnsignedExact estimate =
      UnsignedExact{1} << ((bits + 1U) / 2U);
  for (;;) {
    const UnsignedExact next =
        (estimate + input / estimate) >> 1U;
    if (next >= estimate) return static_cast<Exact>(estimate);
    estimate = next;
  }
}

std::string read_file_bytes(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("cannot open input file");
  return std::string(
      std::istreambuf_iterator<char>(input),
      std::istreambuf_iterator<char>());
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
      const std::uint32_t sum1 =
          rotate_right(e, 6) ^ rotate_right(e, 11) ^
          rotate_right(e, 25);
      const std::uint32_t choice = (e & f) ^ (~e & g);
      const std::uint32_t temporary1 =
          h + sum1 + choice + constants[index] + words[index];
      const std::uint32_t sum0 =
          rotate_right(a, 2) ^ rotate_right(a, 13) ^
          rotate_right(a, 22);
      const std::uint32_t majority =
          (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t temporary2 = sum0 + majority;
      h = g;
      g = f;
      f = e;
      e = d + temporary1;
      d = c;
      c = b;
      b = a;
      a = temporary1 + temporary2;
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

SignMatrix read_sign_matrix(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open start matrix");
  SignMatrix matrix{};
  for (auto& row : matrix) {
    for (int& value : row) {
      if (!(input >> value) || (value != -1 && value != 1)) {
        throw std::runtime_error("invalid start matrix");
      }
    }
  }
  std::string extra;
  if (input >> extra) throw std::runtime_error("extra start matrix data");
  return matrix;
}

Gram gram_of(const SignMatrix& matrix) {
  Gram gram{};
  for (int row = 0; row < kOrder; ++row) {
    for (int other = 0; other < kOrder; ++other) {
      for (int column = 0; column < kOrder; ++column) {
        gram[row][other] +=
            matrix[row][column] * matrix[other][column];
      }
    }
  }
  return gram;
}

Levels levels_from_gram(const Gram& gram) {
  Levels levels{};
  for (int diagonal = 0; diagonal < kOrder; ++diagonal) {
    if (gram[diagonal][diagonal] != kOrder) {
      throw std::runtime_error("start Gram diagonal is not 23");
    }
  }
  std::array<int, kOrder> switches{};
  switches[0] = 1;
  for (int row = 1; row < kOrder; ++row) {
    const int value = gram[0][row];
    if (value == -5 || value == -1 || value == 3 || value == 7) {
      switches[row] = 1;
    } else if (
        -value == -5 || -value == -1 ||
        -value == 3 || -value == 7) {
      switches[row] = -1;
    } else {
      throw std::runtime_error(
          "start Gram cannot be row-switched into {-5,-1,3,7}");
    }
  }
  for (int index = 0; index < kEdgeCount; ++index) {
    const Edge edge = kEdges[index];
    const int value =
        switches[edge.first] * switches[edge.second] *
        gram[edge.first][edge.second];
    if (value != -5 && value != -1 && value != 3 && value != 7) {
      throw std::runtime_error(
          "start Gram is not consistently switchable into "
          "{-5,-1,3,7}");
    }
    levels[index] = static_cast<std::int8_t>(value);
  }
  return levels;
}

Gram gram_from_levels(const Levels& levels) {
  Gram gram{};
  for (int row = 0; row < kOrder; ++row) gram[row][row] = kOrder;
  for (int index = 0; index < kEdgeCount; ++index) {
    const Edge edge = kEdges[index];
    gram[edge.first][edge.second] = levels[index];
    gram[edge.second][edge.first] = levels[index];
  }
  return gram;
}

int level_count(const Levels& levels, int value) {
  return static_cast<int>(
      std::count(levels.begin(), levels.end(), value));
}

LevelKey level_key(const Levels& levels) {
  LevelKey key{};
  for (int index = 0; index < kEdgeCount; ++index) {
    const int code = (static_cast<int>(levels[index]) + 5) / 4;
    key[(2 * index) / 64] |=
        static_cast<std::uint64_t>(code) << ((2 * index) % 64);
  }
  return key;
}

struct State {
  Levels levels{};
  Gram gram{};
  FloatingMatrix inverse{};
  long double log_determinant =
      -std::numeric_limits<long double>::infinity();
  Exact exact_determinant = 0;
};

bool rebuild_numeric(State& state) {
  state.gram = gram_from_levels(state.levels);
  FloatingMatrix lower{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column <= row; ++column) {
      long double value = state.gram[row][column];
      for (int inner = 0; inner < column; ++inner) {
        value -= lower[row][inner] * lower[column][inner];
      }
      if (row == column) {
        if (!(value > 1e-16L) || !std::isfinite(value)) return false;
        lower[row][column] = std::sqrt(value);
      } else {
        lower[row][column] = value / lower[column][column];
      }
    }
  }
  state.log_determinant = 0;
  for (int index = 0; index < kOrder; ++index) {
    state.log_determinant += 2 * std::log(lower[index][index]);
  }

  FloatingMatrix inverse{};
  for (int rhs = 0; rhs < kOrder; ++rhs) {
    std::array<long double, kOrder> forward{};
    for (int row = 0; row < kOrder; ++row) {
      long double value = static_cast<long double>(row == rhs);
      for (int inner = 0; inner < row; ++inner) {
        value -= lower[row][inner] * forward[inner];
      }
      forward[row] = value / lower[row][row];
    }
    std::array<long double, kOrder> solution{};
    for (int row = kOrder - 1; row >= 0; --row) {
      long double value = forward[row];
      for (int inner = row + 1; inner < kOrder; ++inner) {
        value -= lower[inner][row] * solution[inner];
      }
      solution[row] = value / lower[row][row];
    }
    for (int row = 0; row < kOrder; ++row) {
      inverse[row][rhs] = solution[row];
    }
  }
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      state.inverse[row][column] =
          (inverse[row][column] + inverse[column][row]) / 2;
    }
  }
  return std::isfinite(state.log_determinant);
}

long double projected_log_determinant(
    const State& state, int edge_index, int new_value) {
  const Edge edge = kEdges[edge_index];
  const long double delta =
      static_cast<long double>(
          new_value - static_cast<int>(state.levels[edge_index]));
  const long double diagonal =
      1 + delta * state.inverse[edge.first][edge.second];
  const long double ratio =
      diagonal * diagonal -
      delta * delta *
          state.inverse[edge.first][edge.first] *
          state.inverse[edge.second][edge.second];
  if (!(ratio > 0) || !std::isfinite(ratio)) {
    return -std::numeric_limits<long double>::infinity();
  }
  return state.log_determinant + std::log(ratio);
}

struct Arguments {
  std::filesystem::path start;
  std::filesystem::path output;
  std::string mode = "tabu";
  std::uint64_t seed = 23;
  double seconds = 60.0;
  double heartbeat_seconds = 10.0;
  double checkpoint_seconds = 30.0;
  std::uint64_t max_iterations = 0;
  std::uint64_t restart_iterations = 5000;
  int kick_size = 8;
  int tabu_tenure = 13;
  int min_minus5 = 1;
  int max_minus5 = 8;
  bool relocate_minus5 = false;
  int min_plus7 = 0;
  int max_plus7 = 0;
  bool relocate_plus7 = false;
  std::uint64_t cooling_iterations = 10000;
  std::uint64_t max_stored_hits = 256;
  long double temperature_start = 0.05L;
  long double temperature_end = 0.001L;
};

std::uint64_t parse_unsigned(
    std::string_view text, std::string_view option) {
  const std::string copy(text);
  std::size_t parsed = 0;
  const unsigned long long value = std::stoull(copy, &parsed);
  if (parsed != copy.size()) {
    throw std::runtime_error("invalid integer for " + std::string(option));
  }
  return static_cast<std::uint64_t>(value);
}

int parse_integer(std::string_view text, std::string_view option) {
  const std::string copy(text);
  std::size_t parsed = 0;
  const long value = std::stol(copy, &parsed);
  if (parsed != copy.size() ||
      value < std::numeric_limits<int>::min() ||
      value > std::numeric_limits<int>::max()) {
    throw std::runtime_error("invalid integer for " + std::string(option));
  }
  return static_cast<int>(value);
}

double parse_double(std::string_view text, std::string_view option) {
  const std::string copy(text);
  std::size_t parsed = 0;
  const double value = std::stod(copy, &parsed);
  if (parsed != copy.size() || !std::isfinite(value)) {
    throw std::runtime_error("invalid number for " + std::string(option));
  }
  return value;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    auto value = [&]() -> std::string {
      if (++index >= argc) throw std::runtime_error("missing option value");
      return argv[index];
    };
    if (option == "--start") arguments.start = value();
    else if (option == "--output") arguments.output = value();
    else if (option == "--mode") arguments.mode = value();
    else if (option == "--seed") {
      arguments.seed = parse_unsigned(value(), option);
    } else if (option == "--seconds") {
      arguments.seconds = parse_double(value(), option);
    } else if (option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds = parse_double(value(), option);
    } else if (option == "--checkpoint-seconds") {
      arguments.checkpoint_seconds = parse_double(value(), option);
    } else if (option == "--max-iterations") {
      arguments.max_iterations = parse_unsigned(value(), option);
    } else if (option == "--restart-iterations") {
      arguments.restart_iterations = parse_unsigned(value(), option);
    } else if (option == "--kick-size") {
      arguments.kick_size = parse_integer(value(), option);
    } else if (option == "--tabu-tenure") {
      arguments.tabu_tenure = parse_integer(value(), option);
    } else if (option == "--min-minus5") {
      arguments.min_minus5 = parse_integer(value(), option);
    } else if (option == "--max-minus5") {
      arguments.max_minus5 = parse_integer(value(), option);
    } else if (option == "--relocate-minus5") {
      arguments.relocate_minus5 = true;
    } else if (option == "--min-plus7") {
      arguments.min_plus7 = parse_integer(value(), option);
    } else if (option == "--max-plus7") {
      arguments.max_plus7 = parse_integer(value(), option);
    } else if (option == "--relocate-plus7") {
      arguments.relocate_plus7 = true;
    } else if (option == "--cooling-iterations") {
      arguments.cooling_iterations = parse_unsigned(value(), option);
    } else if (option == "--max-stored-hits") {
      arguments.max_stored_hits = parse_unsigned(value(), option);
    } else if (option == "--temperature-start") {
      arguments.temperature_start = parse_double(value(), option);
    } else if (option == "--temperature-end") {
      arguments.temperature_end = parse_double(value(), option);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (arguments.start.empty()) throw std::runtime_error("--start is required");
  if (arguments.output.empty()) {
    throw std::runtime_error("--output is required");
  }
  if (arguments.mode != "hill" &&
      arguments.mode != "tabu" &&
      arguments.mode != "anneal") {
    throw std::runtime_error("--mode must be hill, tabu, or anneal");
  }
  if (!(arguments.seconds > 0)) {
    throw std::runtime_error("--seconds must be positive");
  }
  if (arguments.heartbeat_seconds < 0 ||
      arguments.checkpoint_seconds <= 0) {
    throw std::runtime_error(
        "heartbeat must be non-negative and checkpoint must be positive");
  }
  if (arguments.kick_size < 0 ||
      arguments.kick_size > kEdgeCount) {
    throw std::runtime_error("--kick-size must be between 0 and 253");
  }
  if (arguments.tabu_tenure <= 0) {
    throw std::runtime_error("--tabu-tenure must be positive");
  }
  if (arguments.min_minus5 < 0 ||
      arguments.max_minus5 < arguments.min_minus5 ||
      arguments.max_minus5 > kEdgeCount) {
    throw std::runtime_error(
        "minus5 bounds must satisfy 0 <= min <= max <= 253");
  }
  if (arguments.relocate_minus5 &&
      arguments.min_minus5 == 0) {
    throw std::runtime_error(
        "--relocate-minus5 requires a positive -5 minimum");
  }
  if (arguments.min_plus7 < 0 ||
      arguments.max_plus7 < arguments.min_plus7 ||
      arguments.max_plus7 > kEdgeCount) {
    throw std::runtime_error(
        "plus7 bounds must satisfy 0 <= min <= max <= 253");
  }
  if (arguments.relocate_plus7 &&
      arguments.min_plus7 == 0) {
    throw std::runtime_error(
        "--relocate-plus7 requires a positive +7 minimum");
  }
  if (arguments.cooling_iterations == 0) {
    throw std::runtime_error("--cooling-iterations must be positive");
  }
  if (!(arguments.temperature_start > 0) ||
      !(arguments.temperature_end > 0) ||
      arguments.temperature_end > arguments.temperature_start) {
    throw std::runtime_error(
        "temperatures must satisfy start >= end > 0");
  }
  auto normalized = [](const std::filesystem::path& path) {
    std::error_code error;
    const std::filesystem::path result =
        std::filesystem::weakly_canonical(
            std::filesystem::absolute(path), error);
    if (!error) return result;
    return std::filesystem::absolute(path).lexically_normal();
  };
  const auto start = normalized(arguments.start);
  const auto output = normalized(arguments.output);
  if (start == output) {
    throw std::runtime_error("--start and --output must be distinct");
  }
  return arguments;
}

void print_usage(std::ostream& output) {
  output
      << "usage: gram_multilevel --start MATRIX --output SNAPSHOT "
         "[options]\n\n"
      << "Search exact-SPD order-23 Gram matrices with diagonal 23 and "
         "off-diagonal levels -5, -1, 3, and optionally 7.\n\n"
      << "options:\n"
      << "  --mode MODE                 hill, tabu, or anneal "
         "(default tabu)\n"
      << "  --seed N                    PRNG seed (default 23)\n"
      << "  --seconds S                 positive wall limit\n"
      << "  --max-iterations N          0 means wall limit only\n"
      << "  --restart-iterations N      stagnant iterations before kick\n"
      << "  --kick-size N               adjacent-level changes per kick\n"
      << "  --tabu-tenure N             positive edge tabu tenure\n"
      << "  --min-minus5 N              minimum -5 edges (default 1)\n"
      << "  --max-minus5 N              maximum -5 edges (default 8)\n"
      << "  --relocate-minus5           add atomic -5-edge relocation "
         "moves (requires a positive minimum)\n"
      << "  --min-plus7 N               minimum +7 edges (default 0)\n"
      << "  --max-plus7 N               maximum +7 edges (default 0)\n"
      << "  --relocate-plus7            add atomic +7-edge relocation "
         "moves (requires a positive minimum)\n"
      << "  --cooling-iterations N      annealing cycle length\n"
      << "  --temperature-start X       annealing start temperature\n"
      << "  --temperature-end X         annealing end temperature\n"
      << "  --max-stored-hits N         bounded exact-hit array size\n"
      << "  --heartbeat-seconds S       0 disables stdout heartbeats\n"
      << "  --checkpoint-seconds S      positive snapshot interval\n"
      << "  --help                      show this text\n";
}

struct Statistics {
  std::uint64_t iterations = 0;
  std::uint64_t candidate_moves = 0;
  std::uint64_t accepted_moves = 0;
  std::uint64_t exact_neighbor_screens = 0;
  std::uint64_t above_frontier = 0;
  std::uint64_t exact_squares = 0;
  std::uint64_t exact_pd_hit_checks = 0;
  std::uint64_t exact_pd_state_checks = 0;
  std::uint64_t qualified_survivors = 0;
  std::uint64_t unrecorded_square_observations = 0;
  std::uint64_t non_pd_proposals = 0;
  std::uint64_t restarts = 0;
};

struct SquareHit {
  Levels levels{};
  Exact determinant = 0;
  Exact root = 0;
  bool divisible = false;
  bool positive_definite = false;
};

bool screen_exact_candidate(
    const Levels& levels,
    const Gram& gram,
    Exact determinant,
    Statistics& statistics,
    std::vector<SquareHit>& hits,
    std::set<LevelKey>& recorded_hits,
    std::uint64_t max_stored_hits) {
  ++statistics.exact_neighbor_screens;
  if (determinant <= frontier_determinant()) return false;
  ++statistics.above_frontier;
  const Exact root = integer_square_root(determinant);
  if (root * root != determinant) return false;
  ++statistics.exact_squares;
  const bool divisible = root % (Exact{1} << kRequiredPower) == 0;
  ++statistics.exact_pd_hit_checks;
  const bool positive_definite = exact_positive_definite(gram);
  const bool qualified = divisible && positive_definite;
  if (qualified) ++statistics.qualified_survivors;

  const LevelKey key = level_key(levels);
  const bool first_observation = recorded_hits.insert(key).second;
  if (!first_observation) return false;
  if (hits.size() >= max_stored_hits) {
    ++statistics.unrecorded_square_observations;
    return false;
  }
  hits.push_back(
      {levels, determinant, root, divisible, positive_definite});
  return true;
}

std::string json_escape(std::string_view value) {
  std::string result = "\"";
  for (const unsigned char character : value) {
    switch (character) {
      case '"':
        result += "\\\"";
        break;
      case '\\':
        result += "\\\\";
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

void append_level_edges(
    std::ostream& output, const Levels& levels, int selected_value) {
  output << '[';
  bool first = true;
  for (int index = 0; index < kEdgeCount; ++index) {
    if (levels[index] != selected_value) continue;
    if (!first) output << ',';
    first = false;
    output << '[' << kEdges[index].first + 1
           << ',' << kEdges[index].second + 1 << ']';
  }
  output << ']';
}

void append_level_state(std::ostream& output, const State& state) {
  output << "\"exact_determinant\":"
         << json_escape(exact_string(state.exact_determinant));
  output << ",\"log_determinant\":"
         << static_cast<double>(state.log_determinant);
  output << ",\"minus5_count\":" << level_count(state.levels, -5);
  output << ",\"minus5_edges\":";
  append_level_edges(output, state.levels, -5);
  output << ",\"plus3_count\":" << level_count(state.levels, 3);
  output << ",\"plus3_edges\":";
  append_level_edges(output, state.levels, 3);
  output << ",\"plus7_count\":" << level_count(state.levels, 7);
  output << ",\"plus7_edges\":";
  append_level_edges(output, state.levels, 7);
}

void atomic_write(
    const std::filesystem::path& path, const std::string& contents) {
  if (!path.parent_path().empty()) {
    std::filesystem::create_directories(path.parent_path());
  }
  const std::filesystem::path temporary = path.string() + ".tmp";
  {
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("cannot create checkpoint");
    output.write(
        contents.data(), static_cast<std::streamsize>(contents.size()));
    output.flush();
    if (!output) throw std::runtime_error("cannot flush checkpoint");
  }
  std::filesystem::rename(temporary, path);
}

void write_snapshot(
    const Arguments& arguments,
    std::string_view start_sha256,
    const State& current,
    const State& best,
    const Statistics& statistics,
    const std::vector<SquareHit>& hits,
    double elapsed,
    std::string_view termination) {
  std::ostringstream output;
  output << std::setprecision(18);
  output << "{\"best\":{";
  append_level_state(output, best);
  output << "},\"challenge_id\":\"maxdet-23-v1\"";
  output << ",\"claim_boundary\":"
         << json_escape(
                "Exact Gram screening only; every hit still requires "
                "an exact G=AA^T decomposition with "
                "A in {-1,+1}^{23x23}.");
  output << ",\"complete\":"
         << (termination == "running" ? "false" : "true");
  output << ",\"current\":{";
  append_level_state(output, current);
  output << "},\"elapsed_seconds\":" << elapsed;
  output << ",\"engine\":\"gram-multilevel-tabu\"";
  output << ",\"frontier_root\":"
         << json_escape(std::to_string(kFrontierRoot));
  output << ",\"frontier_squared\":"
         << json_escape(exact_string(frontier_determinant()));
  output << ",\"hits\":[";
  for (std::size_t index = 0; index < hits.size(); ++index) {
    if (index != 0) output << ',';
    const SquareHit& hit = hits[index];
    output << "{\"determinant\":"
           << json_escape(exact_string(hit.determinant));
    output << ",\"divisible_by_2_22\":"
           << (hit.divisible ? "true" : "false");
    output << ",\"minus5_count\":" << level_count(hit.levels, -5);
    output << ",\"minus5_edges\":";
    append_level_edges(output, hit.levels, -5);
    output << ",\"plus3_count\":" << level_count(hit.levels, 3);
    output << ",\"plus3_edges\":";
    append_level_edges(output, hit.levels, 3);
    output << ",\"plus7_count\":" << level_count(hit.levels, 7);
    output << ",\"plus7_edges\":";
    append_level_edges(output, hit.levels, 7);
    output << ",\"positive_definite\":"
           << (hit.positive_definite ? "true" : "false");
    output << ",\"qualified\":"
           << (hit.divisible && hit.positive_definite ? "true" : "false");
    output << ",\"square_root\":"
           << json_escape(exact_string(hit.root)) << '}';
  }
  output << ']';
  output << ",\"level_values\":[-5,-1,3";
  if (arguments.max_plus7 > 0) output << ",7";
  output << ']';
  output << ",\"mode\":" << json_escape(arguments.mode);
  output << ",\"normalization\":"
         << json_escape(
                arguments.max_plus7 > 0
                    ? kFourLevelNormalization
                    : kThreeLevelNormalization);
  output << ",\"parameters\":{";
  output << "\"checkpoint_seconds\":" << arguments.checkpoint_seconds;
  output << ",\"cooling_iterations\":" << arguments.cooling_iterations;
  output << ",\"heartbeat_seconds\":" << arguments.heartbeat_seconds;
  output << ",\"kick_size\":" << arguments.kick_size;
  output << ",\"max_iterations\":" << arguments.max_iterations;
  output << ",\"max_minus5\":" << arguments.max_minus5;
  output << ",\"max_stored_hits\":" << arguments.max_stored_hits;
  output << ",\"min_minus5\":" << arguments.min_minus5;
  output << ",\"max_plus7\":" << arguments.max_plus7;
  output << ",\"min_plus7\":" << arguments.min_plus7;
  output << ",\"restart_iterations\":" << arguments.restart_iterations;
  output << ",\"relocate_minus5\":"
         << (arguments.relocate_minus5 ? "true" : "false");
  output << ",\"relocate_plus7\":"
         << (arguments.relocate_plus7 ? "true" : "false");
  output << ",\"seconds\":" << arguments.seconds;
  output << ",\"tabu_tenure\":" << arguments.tabu_tenure;
  output << ",\"temperature_end\":"
         << static_cast<double>(arguments.temperature_end);
  output << ",\"temperature_start\":"
         << static_cast<double>(arguments.temperature_start) << '}';
  output << ",\"schema_version\":1";
  output << ",\"seed\":" << arguments.seed;
  output << ",\"start\":{";
  output << "\"path\":" << json_escape(arguments.start.string());
  output << ",\"raw_sha256\":" << json_escape(start_sha256) << '}';
  output << ",\"statistics\":{";
  output << "\"above_frontier\":" << statistics.above_frontier;
  output << ",\"accepted_moves\":" << statistics.accepted_moves;
  output << ",\"candidate_moves\":" << statistics.candidate_moves;
  output << ",\"exact_neighbor_screens\":"
         << statistics.exact_neighbor_screens;
  output << ",\"exact_pd_hit_checks\":"
         << statistics.exact_pd_hit_checks;
  output << ",\"exact_pd_state_checks\":"
         << statistics.exact_pd_state_checks;
  output << ",\"exact_squares\":" << statistics.exact_squares;
  output << ",\"iterations\":" << statistics.iterations;
  output << ",\"non_pd_proposals\":" << statistics.non_pd_proposals;
  output << ",\"qualified_survivors\":"
         << statistics.qualified_survivors;
  output << ",\"restarts\":" << statistics.restarts;
  output << ",\"unrecorded_square_observations\":"
         << statistics.unrecorded_square_observations << '}';
  output << ",\"termination\":" << json_escape(termination);
  output << "}\n";
  atomic_write(arguments.output, output.str());
}

struct Candidate {
  int edge_index = -1;
  int new_value = -1;
  int second_edge_index = -1;
  int second_new_value = -1;
  long double projected_log_determinant =
      -std::numeric_limits<long double>::infinity();
  long double priority =
      -std::numeric_limits<long double>::infinity();
  Exact exact_determinant = 0;
};

long double annealing_temperature(
    const Arguments& arguments, std::uint64_t iteration) {
  const long double phase = static_cast<long double>(
      iteration % arguments.cooling_iterations) /
      static_cast<long double>(arguments.cooling_iterations);
  return arguments.temperature_start *
         std::pow(
             arguments.temperature_end / arguments.temperature_start,
             phase);
}

std::array<int, 2> adjacent_values(int value) {
  if (value == -5) return {-1, 99};
  if (value == -1) return {-5, 3};
  if (value == 3) return {-1, 7};
  if (value == 7) return {3, 99};
  throw std::runtime_error("state contains an invalid level");
}

bool counts_allowed(
    int current_minus5,
    int current_plus7,
    int old_value,
    int new_value,
    const Arguments& arguments) {
  int next_minus5 = current_minus5;
  int next_plus7 = current_plus7;
  if (old_value == -5) --next_minus5;
  if (new_value == -5) ++next_minus5;
  if (old_value == 7) --next_plus7;
  if (new_value == 7) ++next_plus7;
  return next_minus5 >= arguments.min_minus5 &&
         next_minus5 <= arguments.max_minus5 &&
         next_plus7 >= arguments.min_plus7 &&
         next_plus7 <= arguments.max_plus7;
}

bool kick_from_best(
    const State& best,
    State& destination,
    const Arguments& arguments,
    std::mt19937_64& randomizer,
    Statistics& statistics,
    std::vector<SquareHit>& hits,
    std::set<LevelKey>& recorded_hits) {
  if (arguments.kick_size == 0) return false;
  std::array<int, kEdgeCount> indices{};
  std::iota(indices.begin(), indices.end(), 0);
  for (int attempt = 0; attempt < 128; ++attempt) {
    destination = best;
    std::shuffle(indices.begin(), indices.end(), randomizer);
    int minus5 = level_count(destination.levels, -5);
    int plus7 = level_count(destination.levels, 7);
    int changed = 0;
    std::array<bool, kEdgeCount> touched{};
    if (arguments.relocate_plus7 && arguments.kick_size >= 2) {
      std::vector<int> old_edges;
      std::vector<int> destinations;
      destinations.reserve(kEdgeCount);
      for (int edge_index = 0; edge_index < kEdgeCount; ++edge_index) {
        if (destination.levels[edge_index] == 7) {
          old_edges.push_back(edge_index);
        } else if (destination.levels[edge_index] == 3) {
          destinations.push_back(edge_index);
        }
      }
      if (old_edges.empty() || destinations.empty()) continue;
      std::uniform_int_distribution<std::size_t> choose_old(
          0, old_edges.size() - 1);
      std::uniform_int_distribution<std::size_t> choose_new(
          0, destinations.size() - 1);
      const int old_index =
          old_edges.size() == 1
              ? old_edges.front()
              : old_edges[choose_old(randomizer)];
      const int new_index = destinations[choose_new(randomizer)];
      destination.levels[old_index] = 3;
      destination.levels[new_index] = 7;
      touched[old_index] = true;
      touched[new_index] = true;
      changed = 2;
    }
    if (arguments.relocate_minus5 &&
        arguments.kick_size - changed >= 2) {
      std::vector<int> old_edges;
      std::vector<int> destinations;
      destinations.reserve(kEdgeCount);
      for (int edge_index = 0; edge_index < kEdgeCount; ++edge_index) {
        if (destination.levels[edge_index] == -5) {
          old_edges.push_back(edge_index);
        } else if (destination.levels[edge_index] == -1) {
          destinations.push_back(edge_index);
        }
      }
      if (old_edges.empty() || destinations.empty()) {
        continue;
      }
      std::uniform_int_distribution<std::size_t> choose_old(
          0, old_edges.size() - 1);
      std::uniform_int_distribution<std::size_t> choose_new(
          0, destinations.size() - 1);
      const int old_index =
          old_edges.size() == 1
              ? old_edges.front()
              : old_edges[choose_old(randomizer)];
      const int new_index = destinations[choose_new(randomizer)];
      destination.levels[old_index] = -1;
      destination.levels[new_index] = -5;
      touched[old_index] = true;
      touched[new_index] = true;
      changed = 2;
    }
    for (int cursor = 0;
         cursor < kEdgeCount && changed < arguments.kick_size;
         ++cursor) {
      const int edge_index = indices[cursor];
      if (touched[edge_index]) continue;
      const int old_value = destination.levels[edge_index];
      std::array<int, 2> choices = adjacent_values(old_value);
      std::shuffle(choices.begin(), choices.end(), randomizer);
      for (const int new_value : choices) {
        if (new_value == 99 ||
            !counts_allowed(
                minus5,
                plus7,
                old_value,
                new_value,
                arguments)) {
          continue;
        }
        destination.levels[edge_index] =
            static_cast<std::int8_t>(new_value);
        touched[edge_index] = true;
        if (old_value == -5) --minus5;
        if (new_value == -5) ++minus5;
        if (old_value == 7) --plus7;
        if (new_value == 7) ++plus7;
        ++changed;
        break;
      }
    }
    if (changed != arguments.kick_size) continue;
    if (!rebuild_numeric(destination)) continue;
    destination.exact_determinant =
        exact_determinant(destination.gram);
    ++statistics.exact_pd_state_checks;
    if (!exact_positive_definite(destination.gram)) continue;
    screen_exact_candidate(
        destination.levels,
        destination.gram,
        destination.exact_determinant,
        statistics,
        hits,
        recorded_hits,
        arguments.max_stored_hits);
    return true;
  }
  return false;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    for (int index = 1; index < argc; ++index) {
      if (std::string_view(argv[index]) == "--help") {
        print_usage(std::cout);
        return 0;
      }
    }
    const Arguments arguments = parse_arguments(argc, argv);
    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);

    {
      std::error_code error;
      const auto status =
          std::filesystem::symlink_status(arguments.output, error);
      const bool missing =
          error == std::errc::no_such_file_or_directory;
      if (missing) {
        error.clear();
      }
      if (error) {
        throw std::runtime_error(
            "cannot inspect output path before starting");
      }
      if (!missing &&
          status.type() != std::filesystem::file_type::not_found) {
        throw std::runtime_error(
            "refusing to overwrite an existing snapshot");
      }
    }

    long double log_crt_modulus = 0;
    for (const std::uint64_t prime : kCrtPrimes) {
      if (!is_prime(prime)) {
        throw std::runtime_error("configured CRT modulus is not prime");
      }
      log_crt_modulus += std::log(static_cast<long double>(prime));
    }
    const long double log_twice_hadamard_bound =
        std::log(2.0L) +
        static_cast<long double>(kOrder) / 2 *
            std::log(
                arguments.max_plus7 > 0 ? 1607.0L : 1079.0L);
    if (!(log_crt_modulus > log_twice_hadamard_bound)) {
      throw std::runtime_error(
          "CRT modulus does not uniquely cover the determinant bound");
    }

    const std::string start_bytes = read_file_bytes(arguments.start);
    const std::string start_sha256 = sha256(start_bytes);
    const SignMatrix start_matrix = read_sign_matrix(arguments.start);
    State current;
    current.levels = levels_from_gram(gram_of(start_matrix));
    const int start_minus5 = level_count(current.levels, -5);
    if (start_minus5 < arguments.min_minus5 ||
        start_minus5 > arguments.max_minus5) {
      throw std::runtime_error(
          "start Gram -5 count is outside requested bounds");
    }
    const int start_plus7 = level_count(current.levels, 7);
    if (start_plus7 < arguments.min_plus7 ||
        start_plus7 > arguments.max_plus7) {
      throw std::runtime_error(
          "start Gram +7 count is outside requested bounds");
    }
    if (!rebuild_numeric(current) ||
        !exact_positive_definite(current.gram)) {
      throw std::runtime_error("start Gram is not exactly positive definite");
    }
    current.exact_determinant = exact_determinant(current.gram);
    State best = current;

    std::mt19937_64 randomizer(arguments.seed);
    std::uniform_real_distribution<long double> unit(0.0L, 1.0L);
    Statistics statistics;
    std::vector<SquareHit> hits;
    std::set<LevelKey> recorded_hits;
    std::array<std::uint64_t, kEdgeCount> tabu_until{};

    const auto started = Clock::now();
    const auto deadline =
        started + std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(arguments.seconds));
    auto next_heartbeat =
        started + std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(
                          arguments.heartbeat_seconds));
    auto next_checkpoint =
        started + std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(
                          arguments.checkpoint_seconds));
    auto elapsed = [&]() {
      return std::chrono::duration<double>(Clock::now() - started).count();
    };

    std::uint64_t last_best_iteration = 0;
    bool local_optimum = false;
    write_snapshot(
        arguments,
        start_sha256,
        current,
        best,
        statistics,
        hits,
        0.0,
        "running");

    while (!stop_requested &&
           Clock::now() < deadline &&
           (arguments.max_iterations == 0 ||
            statistics.iterations < arguments.max_iterations)) {
      const std::size_t hit_count_before = hits.size();
      std::vector<Candidate> candidates;
      candidates.reserve(3 * kEdgeCount);
      const long double temperature =
          annealing_temperature(arguments, statistics.iterations);
      const int current_minus5 = level_count(current.levels, -5);
      const int current_plus7 = level_count(current.levels, 7);

      for (int edge_index = 0; edge_index < kEdgeCount; ++edge_index) {
        const int old_value = current.levels[edge_index];
        for (const int new_value : adjacent_values(old_value)) {
          if (new_value == 99 ||
              !counts_allowed(
                  current_minus5,
                  current_plus7,
                  old_value,
                  new_value,
                  arguments)) {
            continue;
          }
          ++statistics.candidate_moves;
          Levels candidate_levels = current.levels;
          candidate_levels[edge_index] =
              static_cast<std::int8_t>(new_value);
          Gram candidate_gram = current.gram;
          const Edge edge = kEdges[edge_index];
          candidate_gram[edge.first][edge.second] = new_value;
          candidate_gram[edge.second][edge.first] = new_value;
          const Exact determinant = exact_determinant(candidate_gram);
          screen_exact_candidate(
              candidate_levels,
              candidate_gram,
              determinant,
              statistics,
              hits,
              recorded_hits,
              arguments.max_stored_hits);

          const long double projected =
              projected_log_determinant(
                  current, edge_index, new_value);
          if (!std::isfinite(projected)) continue;
          const bool tabu =
              statistics.iterations < tabu_until[edge_index];
          const bool aspiration = determinant > best.exact_determinant;
          if (arguments.mode == "tabu" && tabu && !aspiration) continue;
          if (arguments.mode == "hill" &&
              determinant <= current.exact_determinant) {
            continue;
          }

          long double priority = projected;
          if (arguments.mode == "anneal") {
            const long double sample = std::max(
                unit(randomizer),
                std::numeric_limits<long double>::min());
            const long double gumbel = -std::log(-std::log(sample));
            priority =
                (projected - current.log_determinant) / temperature +
                gumbel;
          }
          candidates.push_back(
              {edge_index,
               new_value,
               -1,
               -1,
               projected,
               priority,
               determinant});
        }
      }

      if (arguments.relocate_minus5) {
        for (int old_index = 0;
             old_index < kEdgeCount;
             ++old_index) {
          if (current.levels[old_index] != -5) continue;
          for (int new_index = 0;
               new_index < kEdgeCount;
               ++new_index) {
            if (current.levels[new_index] != -1) continue;
            ++statistics.candidate_moves;
            Levels candidate_levels = current.levels;
            candidate_levels[old_index] = -1;
            candidate_levels[new_index] = -5;
            Gram candidate_gram = current.gram;
            const Edge old_edge = kEdges[old_index];
            const Edge new_edge = kEdges[new_index];
            candidate_gram[old_edge.first][old_edge.second] = -1;
            candidate_gram[old_edge.second][old_edge.first] = -1;
            candidate_gram[new_edge.first][new_edge.second] = -5;
            candidate_gram[new_edge.second][new_edge.first] = -5;
            const Exact determinant = exact_determinant(candidate_gram);
            screen_exact_candidate(
                candidate_levels,
                candidate_gram,
                determinant,
                statistics,
                hits,
                recorded_hits,
                arguments.max_stored_hits);
            if (determinant <= 0) continue;
            const long double projected =
                std::log(static_cast<long double>(determinant));
            const bool tabu =
                statistics.iterations < tabu_until[old_index] ||
                statistics.iterations < tabu_until[new_index];
            const bool aspiration =
                determinant > best.exact_determinant;
            if (arguments.mode == "tabu" && tabu && !aspiration) {
              continue;
            }
            if (arguments.mode == "hill" &&
                determinant <= current.exact_determinant) {
              continue;
            }
            long double priority = projected;
            if (arguments.mode == "anneal") {
              const long double sample = std::max(
                  unit(randomizer),
                  std::numeric_limits<long double>::min());
              const long double gumbel = -std::log(-std::log(sample));
              priority =
                  (projected - current.log_determinant) / temperature +
                  gumbel;
            }
            candidates.push_back(
                {old_index,
                 -1,
                 new_index,
                 -5,
                 projected,
                 priority,
                 determinant});
          }
        }
      }

      if (arguments.relocate_plus7) {
        for (int old_index = 0;
             old_index < kEdgeCount;
             ++old_index) {
          if (current.levels[old_index] != 7) continue;
          for (int new_index = 0;
               new_index < kEdgeCount;
               ++new_index) {
            if (current.levels[new_index] != 3) continue;
            ++statistics.candidate_moves;
            Levels candidate_levels = current.levels;
            candidate_levels[old_index] = 3;
            candidate_levels[new_index] = 7;
            Gram candidate_gram = current.gram;
            const Edge old_edge = kEdges[old_index];
            const Edge new_edge = kEdges[new_index];
            candidate_gram[old_edge.first][old_edge.second] = 3;
            candidate_gram[old_edge.second][old_edge.first] = 3;
            candidate_gram[new_edge.first][new_edge.second] = 7;
            candidate_gram[new_edge.second][new_edge.first] = 7;
            const Exact determinant = exact_determinant(candidate_gram);
            screen_exact_candidate(
                candidate_levels,
                candidate_gram,
                determinant,
                statistics,
                hits,
                recorded_hits,
                arguments.max_stored_hits);
            if (determinant <= 0) continue;
            const long double projected =
                std::log(static_cast<long double>(determinant));
            const bool tabu =
                statistics.iterations < tabu_until[old_index] ||
                statistics.iterations < tabu_until[new_index];
            const bool aspiration =
                determinant > best.exact_determinant;
            if (arguments.mode == "tabu" && tabu && !aspiration) {
              continue;
            }
            if (arguments.mode == "hill" &&
                determinant <= current.exact_determinant) {
              continue;
            }
            long double priority = projected;
            if (arguments.mode == "anneal") {
              const long double sample = std::max(
                  unit(randomizer),
                  std::numeric_limits<long double>::min());
              const long double gumbel = -std::log(-std::log(sample));
              priority =
                  (projected - current.log_determinant) / temperature +
                  gumbel;
            }
            candidates.push_back(
                {old_index,
                 3,
                 new_index,
                 7,
                 projected,
                 priority,
                 determinant});
          }
        }
      }

      std::sort(
          candidates.begin(),
          candidates.end(),
          [](const Candidate& left, const Candidate& right) {
            if (left.priority != right.priority) {
              return left.priority > right.priority;
            }
            if (left.edge_index != right.edge_index) {
              return left.edge_index < right.edge_index;
            }
            if (left.new_value != right.new_value) {
              return left.new_value < right.new_value;
            }
            return left.second_edge_index < right.second_edge_index;
          });

      bool moved = false;
      Candidate selected;
      State proposed;
      for (const Candidate& candidate : candidates) {
        proposed = current;
        proposed.levels[candidate.edge_index] =
            static_cast<std::int8_t>(candidate.new_value);
        if (candidate.second_edge_index >= 0) {
          proposed.levels[candidate.second_edge_index] =
              static_cast<std::int8_t>(candidate.second_new_value);
        }
        if (!rebuild_numeric(proposed)) {
          ++statistics.non_pd_proposals;
          continue;
        }
        ++statistics.exact_pd_state_checks;
        if (!exact_positive_definite(proposed.gram)) {
          ++statistics.non_pd_proposals;
          continue;
        }
        proposed.exact_determinant = candidate.exact_determinant;
        const long double discrepancy = std::abs(
            proposed.log_determinant -
            candidate.projected_log_determinant);
        if (discrepancy > 1e-8L) {
          throw std::runtime_error(
              "determinant-lemma projection disagrees with Cholesky");
        }
        selected = candidate;
        moved = true;
        break;
      }

      ++statistics.iterations;
      if (moved) {
        current = proposed;
        ++statistics.accepted_moves;
        if (arguments.mode == "tabu") {
          const int jitter_bound = std::max(1, arguments.tabu_tenure / 2);
          std::uniform_int_distribution<int> jitter(0, jitter_bound);
          tabu_until[selected.edge_index] =
              statistics.iterations +
              static_cast<std::uint64_t>(
                  arguments.tabu_tenure + jitter(randomizer));
          if (selected.second_edge_index >= 0) {
            tabu_until[selected.second_edge_index] =
                statistics.iterations +
                static_cast<std::uint64_t>(
                    arguments.tabu_tenure + jitter(randomizer));
          }
        }
        if (current.exact_determinant > best.exact_determinant) {
          best = current;
          last_best_iteration = statistics.iterations;
          write_snapshot(
              arguments,
              start_sha256,
              current,
              best,
              statistics,
              hits,
              elapsed(),
              "running");
          std::cout << "new best determinant="
                    << exact_string(best.exact_determinant)
                    << " minus5=" << level_count(best.levels, -5)
                    << " plus7=" << level_count(best.levels, 7)
                    << " plus3=" << level_count(best.levels, 3)
                    << " iteration=" << statistics.iterations << '\n'
                    << std::flush;
        }
      } else {
        if (!kick_from_best(
                best,
                current,
                arguments,
                randomizer,
                statistics,
                hits,
                recorded_hits)) {
          local_optimum = true;
          break;
        }
        tabu_until.fill(0);
        ++statistics.restarts;
        last_best_iteration = statistics.iterations;
      }

      if (arguments.restart_iterations != 0 &&
          statistics.iterations - last_best_iteration >=
              arguments.restart_iterations) {
        State restarted;
        if (kick_from_best(
                best,
                restarted,
                arguments,
                randomizer,
                statistics,
                hits,
                recorded_hits)) {
          current = restarted;
          tabu_until.fill(0);
          ++statistics.restarts;
          last_best_iteration = statistics.iterations;
        }
      }

      const auto now = Clock::now();
      if (hits.size() != hit_count_before || now >= next_checkpoint) {
        write_snapshot(
            arguments,
            start_sha256,
            current,
            best,
            statistics,
            hits,
            elapsed(),
            "running");
        next_checkpoint =
            now + std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(
                          arguments.checkpoint_seconds));
      }
      if (arguments.heartbeat_seconds > 0 && now >= next_heartbeat) {
        std::cout << "heartbeat iterations=" << statistics.iterations
                  << " screened=" << statistics.exact_neighbor_screens
                  << " above=" << statistics.above_frontier
                  << " squares=" << statistics.exact_squares
                  << " survivors=" << statistics.qualified_survivors
                  << " minus5=" << level_count(current.levels, -5)
                  << " plus7=" << level_count(current.levels, 7)
                  << " best=" << exact_string(best.exact_determinant)
                  << '\n'
                  << std::flush;
        next_heartbeat =
            now + std::chrono::duration_cast<Clock::duration>(
                          std::chrono::duration<double>(
                              arguments.heartbeat_seconds));
      }
    }

    std::string termination;
    if (stop_requested) {
      termination = "signal";
    } else if (local_optimum) {
      termination = "local-optimum";
    } else if (
        arguments.max_iterations != 0 &&
        statistics.iterations >= arguments.max_iterations) {
      termination = "iteration-limit";
    } else {
      termination = "time-limit";
    }
    write_snapshot(
        arguments,
        start_sha256,
        current,
        best,
        statistics,
        hits,
        elapsed(),
        termination);
    std::cout << "finished termination=" << termination
              << " iterations=" << statistics.iterations
              << " screened=" << statistics.exact_neighbor_screens
              << " squares=" << statistics.exact_squares
              << " survivors=" << statistics.qualified_survivors
              << " best=" << exact_string(best.exact_determinant)
              << " elapsed=" << std::fixed << std::setprecision(3)
              << elapsed() << "s\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "gram_multilevel: " << error.what() << '\n';
    return 2;
  }
}
