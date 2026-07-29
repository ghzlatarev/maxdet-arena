#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using Wide = __int128_t;
using UnsignedWide = __uint128_t;

constexpr int kOrder = 23;
constexpr std::uint64_t kAutomorphismGroupOrder = 442'368;
constexpr std::uint64_t kFrontierRoot = 2'779'447'296'000'000ULL;
constexpr std::uint64_t kRequiredDivisor = std::uint64_t{1} << 22;
constexpr std::uint64_t kRadius1Labeled = 9'360;
constexpr std::uint64_t kRadius2Labeled = 21'312'720;
constexpr std::uint64_t kRadius3Labeled = 20'976'452'640ULL;
constexpr std::uint64_t kMaximumRecommendedRepresentatives = 25'000'000;
constexpr double kMaximumRecommendedCpuSeconds = 8.0 * 60.0 * 60.0;
constexpr std::array<std::uint64_t, 4> kCrtPrimes{
    2'147'483'647ULL,
    2'147'483'629ULL,
    2'147'483'587ULL,
    2'147'483'579ULL,
};

using SignMatrix =
    std::array<std::array<std::int8_t, kOrder>, kOrder>;
using Gram = std::array<std::array<int, kOrder>, kOrder>;
using Permutation = std::array<std::uint8_t, kOrder>;

struct Edge {
  std::uint8_t first = 0;
  std::uint8_t second = 0;

  auto operator<=>(const Edge&) const = default;
};

struct Subset {
  std::array<std::uint8_t, 3> values{};
  std::uint8_t size = 0;

  auto operator<=>(const Subset&) const = default;
};

struct OrbitRecord {
  std::uint64_t index = 0;
  std::uint64_t orbit_size = 0;
  Subset removed_global;
  Subset added_global;
};

struct OrbitCatalog {
  int radius = 0;
  std::uint64_t labeled_count = 0;
  std::uint64_t removed_orbit_count = 0;
  std::uint64_t orbit_size_sum = 0;
  std::map<std::uint64_t, std::uint64_t> orbit_size_histogram;
  std::map<std::uint64_t, std::uint64_t> removed_stabilizer_histogram;
  std::vector<OrbitRecord> records;
  double elapsed_seconds = 0.0;
};

struct SquareHit {
  std::uint64_t orbit_index = 0;
  std::uint64_t orbit_size = 0;
  Wide determinant = 0;
  Wide root = 0;
  bool divisible = false;
  bool positive_definite = false;
  Subset removed_global;
  Subset added_global;
  std::vector<Edge> edges;
};

struct ScreenSummary {
  std::uint64_t representative_count = 0;
  std::uint64_t labeled_count = 0;
  std::uint64_t positive_determinant_representatives = 0;
  std::uint64_t positive_determinant_labeled = 0;
  std::uint64_t above_frontier_representatives = 0;
  std::uint64_t above_frontier_labeled = 0;
  std::uint64_t square_representatives = 0;
  std::uint64_t square_labeled = 0;
  std::uint64_t above_frontier_square_representatives = 0;
  std::uint64_t above_frontier_square_labeled = 0;
  std::uint64_t divisible_above_frontier_square_representatives = 0;
  std::uint64_t divisible_above_frontier_square_labeled = 0;
  std::uint64_t positive_definite_route_representatives = 0;
  std::uint64_t positive_definite_route_labeled = 0;
  std::map<std::string, std::uint64_t> square_root_labeled_histogram;
  std::map<std::string, std::uint64_t> square_root_orbit_histogram;
  std::vector<SquareHit> route_hits;
  bool has_maximum = false;
  Wide maximum_determinant = 0;
  std::uint64_t maximum_orbit_index = 0;
  std::uint64_t maximum_orbit_size = 0;
  Subset maximum_removed_global;
  Subset maximum_added_global;
  double elapsed_seconds = 0.0;
};

enum class Mode {
  kNone,
  kSelfTest,
  kCount,
  kCatalog,
  kScreen,
};

struct Arguments {
  Mode mode = Mode::kNone;
  std::filesystem::path reference =
      "references/orrick-et-al-2003/matrix.txt";
  std::filesystem::path catalog;
  std::filesystem::path manifest;
  std::filesystem::path output;
  std::filesystem::path route_snapshot;
  std::string catalog_sha256;
  std::uint64_t shard_count = 1;
  std::uint64_t shard_index = 0;
  double heartbeat_seconds = 10.0;
};

std::string decimal(Wide value) {
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

std::string json_escape(std::string_view value) {
  std::string result{"\""};
  for (const unsigned char character : value) {
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
          constexpr char digits[] = "0123456789abcdef";
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

std::string read_file_bytes(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("cannot open " + path.string());
  }
  std::ostringstream bytes;
  bytes << input.rdbuf();
  if (!input.eof() && input.fail()) {
    throw std::runtime_error("cannot read " + path.string());
  }
  return bytes.str();
}

void atomic_write(
    const std::filesystem::path& path, const std::string& contents) {
  if (path.empty()) {
    throw std::runtime_error("refusing to write an empty path");
  }
  if (!path.parent_path().empty()) {
    std::filesystem::create_directories(path.parent_path());
  }
  const std::filesystem::path temporary =
      path.parent_path() / ("." + path.filename().string() + ".tmp");
  {
    std::ofstream output(
        temporary, std::ios::binary | std::ios::trunc);
    if (!output) {
      throw std::runtime_error(
          "cannot create temporary output " + temporary.string());
    }
    output.write(
        contents.data(), static_cast<std::streamsize>(contents.size()));
    output.flush();
    if (!output) {
      throw std::runtime_error(
          "cannot flush temporary output " + temporary.string());
    }
  }
  std::filesystem::rename(temporary, path);
}

double parse_nonnegative_double(
    const std::string& text, std::string_view option) {
  std::size_t consumed = 0;
  const double value = std::stod(text, &consumed);
  if (consumed != text.size() || value < 0.0 || !std::isfinite(value)) {
    throw std::runtime_error(
        std::string(option) + " must be finite and non-negative");
  }
  return value;
}

std::uint64_t parse_unsigned(
    const std::string& text, std::string_view option) {
  if (text.empty() || text.front() == '-') {
    throw std::runtime_error(
        std::string(option) + " must be an unsigned integer");
  }
  std::size_t consumed = 0;
  const unsigned long long value = std::stoull(text, &consumed);
  if (consumed != text.size()) {
    throw std::runtime_error(
        std::string(option) + " must be an unsigned integer");
  }
  return static_cast<std::uint64_t>(value);
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    const auto value = [&]() -> std::string {
      if (index + 1 >= argc) {
        throw std::runtime_error("missing value after " + option);
      }
      return argv[++index];
    };
    if (option == "--self-test") {
      if (arguments.mode != Mode::kNone) {
        throw std::runtime_error("select exactly one mode");
      }
      arguments.mode = Mode::kSelfTest;
    } else if (option == "--build-catalog") {
      if (arguments.mode != Mode::kNone) {
        throw std::runtime_error("select exactly one mode");
      }
      arguments.mode = Mode::kCatalog;
    } else if (option == "--count-orbits") {
      if (arguments.mode != Mode::kNone) {
        throw std::runtime_error("select exactly one mode");
      }
      arguments.mode = Mode::kCount;
    } else if (option == "--screen") {
      if (arguments.mode != Mode::kNone) {
        throw std::runtime_error("select exactly one mode");
      }
      arguments.mode = Mode::kScreen;
    } else if (option == "--reference") {
      arguments.reference = value();
    } else if (option == "--catalog") {
      arguments.catalog = value();
    } else if (option == "--manifest") {
      arguments.manifest = value();
    } else if (option == "--output") {
      arguments.output = value();
    } else if (option == "--route-snapshot") {
      arguments.route_snapshot = value();
    } else if (option == "--catalog-sha256") {
      arguments.catalog_sha256 = value();
    } else if (option == "--shard-count") {
      arguments.shard_count = parse_unsigned(value(), option);
    } else if (option == "--shard-index") {
      arguments.shard_index = parse_unsigned(value(), option);
    } else if (option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds =
          parse_nonnegative_double(value(), option);
    } else if (option == "--help") {
      std::cout
          << "Usage:\n"
          << "  gram_radius3_orbits --self-test [--reference FILE]\n"
          << "  gram_radius3_orbits --count-orbits [--reference FILE]\n"
          << "  gram_radius3_orbits --build-catalog --catalog FILE "
             "--manifest FILE [--reference FILE]\n"
          << "  gram_radius3_orbits --screen --catalog FILE --output FILE "
             "--catalog-sha256 HEX [--route-snapshot FILE] "
             "[--shard-count N --shard-index I] "
             "[--reference FILE]\n\n"
          << "Build an exact orbit catalog for three present-edge deletions "
             "and three absent-edge additions around the published order-23 "
             "Gram graph, then screen deterministic catalog shards.\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (arguments.mode == Mode::kNone) {
    throw std::runtime_error(
        "select --self-test, --count-orbits, --build-catalog, or --screen");
  }
  if (arguments.shard_count == 0 ||
      arguments.shard_index >= arguments.shard_count) {
    throw std::runtime_error(
        "--shard-index must be in [0, --shard-count)");
  }
  if (arguments.mode == Mode::kCatalog &&
      (arguments.catalog.empty() || arguments.manifest.empty())) {
    throw std::runtime_error(
        "--build-catalog requires --catalog and --manifest");
  }
  if (arguments.mode == Mode::kScreen &&
      (arguments.catalog.empty() || arguments.output.empty() ||
       arguments.catalog_sha256.empty())) {
    throw std::runtime_error(
        "--screen requires --catalog, --catalog-sha256, and --output");
  }
  if (!arguments.catalog_sha256.empty() &&
      (arguments.catalog_sha256.size() != 64 ||
       arguments.catalog_sha256.find_first_not_of("0123456789abcdef") !=
           std::string::npos)) {
    throw std::runtime_error(
        "--catalog-sha256 must be 64 lowercase hexadecimal digits");
  }
  const std::vector<std::filesystem::path> outputs = {
      arguments.catalog, arguments.manifest, arguments.output,
      arguments.route_snapshot};
  for (std::size_t first = 0; first < outputs.size(); ++first) {
    if (outputs[first].empty()) continue;
    for (std::size_t second = first + 1; second < outputs.size(); ++second) {
      if (outputs[second].empty()) continue;
      if (std::filesystem::absolute(outputs[first]).lexically_normal() ==
          std::filesystem::absolute(outputs[second]).lexically_normal()) {
        throw std::runtime_error("all output paths must be distinct");
      }
    }
  }
  return arguments;
}

SignMatrix parse_sign_matrix(std::string_view bytes) {
  SignMatrix matrix{};
  std::istringstream input{std::string(bytes)};
  std::string line;
  int row = 0;
  while (std::getline(input, line)) {
    if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
    if (row >= kOrder) {
      throw std::runtime_error("reference has more than 23 rows");
    }
    std::istringstream tokens(line);
    int value = 0;
    int column = 0;
    while (tokens >> value) {
      if (column >= kOrder || (value != -1 && value != 1)) {
        throw std::runtime_error(
            "reference contains an invalid matrix token");
      }
      matrix[row][column++] = static_cast<std::int8_t>(value);
    }
    if (column != kOrder) {
      throw std::runtime_error("reference row does not have 23 entries");
    }
    ++row;
  }
  if (row != kOrder) {
    throw std::runtime_error("reference does not have 23 rows");
  }
  return matrix;
}

Gram gram_of(const SignMatrix& matrix) {
  Gram gram{};
  for (int row = 0; row < kOrder; ++row) {
    for (int other = 0; other < kOrder; ++other) {
      int inner = 0;
      for (int column = 0; column < kOrder; ++column) {
        inner += matrix[row][column] * matrix[other][column];
      }
      gram[row][other] = inner;
    }
  }
  return gram;
}

std::vector<Edge> all_edges() {
  std::vector<Edge> result;
  result.reserve(kOrder * (kOrder - 1) / 2);
  for (int first = 0; first < kOrder; ++first) {
    for (int second = first + 1; second < kOrder; ++second) {
      result.push_back(
          Edge{static_cast<std::uint8_t>(first),
               static_cast<std::uint8_t>(second)});
    }
  }
  return result;
}

const std::vector<Edge>& global_edges() {
  static const std::vector<Edge> edges = all_edges();
  return edges;
}

std::array<std::array<std::uint8_t, kOrder>, kOrder>
global_edge_numbers() {
  std::array<std::array<std::uint8_t, kOrder>, kOrder> result{};
  const auto& edges = global_edges();
  for (std::size_t index = 0; index < edges.size(); ++index) {
    const Edge edge = edges[index];
    result[edge.first][edge.second] =
        static_cast<std::uint8_t>(index);
    result[edge.second][edge.first] =
        static_cast<std::uint8_t>(index);
  }
  return result;
}

std::vector<Edge> expected_reference_edges() {
  std::vector<Edge> result;
  const std::array<std::vector<int>, 6> blocks{{
      {0, 1, 2},
      {3, 4, 5, 6},
      {7, 8, 9, 10},
      {11, 12, 13, 14},
      {15, 16, 17, 18},
      {19, 20, 21, 22},
  }};
  for (const auto& block : blocks) {
    for (std::size_t first = 0; first < block.size(); ++first) {
      for (std::size_t second = first + 1; second < block.size(); ++second) {
        result.push_back(
            Edge{static_cast<std::uint8_t>(block[first]),
                 static_cast<std::uint8_t>(block[second])});
      }
    }
  }
  const std::array<std::pair<std::array<int, 2>, std::array<int, 2>>, 3>
      connectors{{
          {{{0, 2}}, {{5, 6}}},
          {{{1, 2}}, {{9, 10}}},
          {{{0, 1}}, {{13, 14}}},
      }};
  for (const auto& [left, right] : connectors) {
    for (const int first : left) {
      for (const int second : right) {
        result.push_back(
            Edge{static_cast<std::uint8_t>(std::min(first, second)),
                 static_cast<std::uint8_t>(std::max(first, second))});
      }
    }
  }
  std::sort(result.begin(), result.end());
  if (result.size() != 45 ||
      std::adjacent_find(result.begin(), result.end()) != result.end()) {
    throw std::runtime_error("internal reference edge construction failed");
  }
  return result;
}

std::vector<Edge> validate_reference_gram(const Gram& gram) {
  std::vector<Edge> present;
  std::array<int, kOrder> degrees{};
  for (int row = 0; row < kOrder; ++row) {
    if (gram[row][row] != kOrder) {
      throw std::runtime_error("reference Gram diagonal is not 23");
    }
    for (int column = row + 1; column < kOrder; ++column) {
      if (gram[row][column] != gram[column][row] ||
          (gram[row][column] != -1 && gram[row][column] != 3)) {
        throw std::runtime_error(
            "reference Gram is outside G=24I-J+4A");
      }
      if (gram[row][column] == 3) {
        present.push_back(
            Edge{static_cast<std::uint8_t>(row),
                 static_cast<std::uint8_t>(column)});
        ++degrees[row];
        ++degrees[column];
      }
    }
  }
  if (present != expected_reference_edges()) {
    throw std::runtime_error(
        "reference defect graph is not the published 45-edge graph");
  }
  std::map<int, int> degree_histogram;
  for (const int degree : degrees) ++degree_histogram[degree];
  if (degree_histogram != std::map<int, int>{{3, 14}, {5, 6}, {6, 3}}) {
    throw std::runtime_error(
        "reference defect degree histogram is not 3^14,5^6,6^3");
  }
  return present;
}

std::vector<std::array<int, 4>> permutations_of_four() {
  std::vector<std::array<int, 4>> result;
  std::array<int, 4> permutation{0, 1, 2, 3};
  do {
    result.push_back(permutation);
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  return result;
}

std::vector<std::array<int, 3>> permutations_of_three() {
  std::vector<std::array<int, 3>> result;
  std::array<int, 3> permutation{0, 1, 2};
  do {
    result.push_back(permutation);
  } while (std::next_permutation(permutation.begin(), permutation.end()));
  return result;
}

std::vector<Permutation> generate_automorphism_group() {
  const auto central_permutations = permutations_of_three();
  const auto four_permutations = permutations_of_four();
  // Arms are indexed by the central vertex missing from their connector.
  // Within an arm the first pair has graph degree 3 and the second degree 5.
  constexpr std::array<std::array<int, 4>, 3> arms{{
      {{7, 8, 9, 10}},
      {{3, 4, 5, 6}},
      {{11, 12, 13, 14}},
  }};
  constexpr std::array<std::array<int, 4>, 2> isolated{{
      {{15, 16, 17, 18}},
      {{19, 20, 21, 22}},
  }};

  std::vector<Permutation> connected;
  connected.reserve(384);
  for (const auto& central : central_permutations) {
    for (int first_arm = 0; first_arm < 4; ++first_arm) {
      for (int second_arm = 0; second_arm < 4; ++second_arm) {
        for (int third_arm = 0; third_arm < 4; ++third_arm) {
          const std::array<int, 3> choices{
              first_arm, second_arm, third_arm};
          Permutation permutation{};
          for (int vertex = 0; vertex < 3; ++vertex) {
            permutation[vertex] =
                static_cast<std::uint8_t>(central[vertex]);
          }
          for (int missing = 0; missing < 3; ++missing) {
            const int target_missing = central[missing];
            const int choice = choices[missing];
            const std::array<int, 4> local{
                (choice & 1) ? 1 : 0,
                (choice & 1) ? 0 : 1,
                (choice & 2) ? 3 : 2,
                (choice & 2) ? 2 : 3,
            };
            for (int position = 0; position < 4; ++position) {
              permutation[arms[missing][position]] =
                  static_cast<std::uint8_t>(
                      arms[target_missing][local[position]]);
            }
          }
          connected.push_back(permutation);
        }
      }
    }
  }

  std::vector<Permutation> group;
  group.reserve(kAutomorphismGroupOrder);
  for (const Permutation& connected_part : connected) {
    for (int swap_blocks = 0; swap_blocks < 2; ++swap_blocks) {
      for (const auto& first_local : four_permutations) {
        for (const auto& second_local : four_permutations) {
          const std::array<std::array<int, 4>, 2> locals{
              first_local, second_local};
          Permutation permutation = connected_part;
          for (int source_block = 0; source_block < 2; ++source_block) {
            const int target_block = source_block ^ swap_blocks;
            for (int position = 0; position < 4; ++position) {
              permutation[isolated[source_block][position]] =
                  static_cast<std::uint8_t>(
                      isolated[target_block][locals[source_block][position]]);
            }
          }
          group.push_back(permutation);
        }
      }
    }
  }
  if (group.size() != kAutomorphismGroupOrder) {
    throw std::runtime_error("automorphism group enumeration count mismatch");
  }
  std::sort(group.begin(), group.end());
  if (std::adjacent_find(group.begin(), group.end()) != group.end()) {
    throw std::runtime_error(
        "automorphism group enumeration contains duplicates");
  }
  return group;
}

Edge permute_edge(const Edge& edge, const Permutation& permutation) {
  int first = permutation[edge.first];
  int second = permutation[edge.second];
  if (first > second) std::swap(first, second);
  return Edge{static_cast<std::uint8_t>(first),
              static_cast<std::uint8_t>(second)};
}

void validate_automorphism_group(
    const std::vector<Permutation>& group,
    const std::vector<Edge>& present) {
  const std::set<Edge> expected(present.begin(), present.end());
  for (const Permutation& permutation : group) {
    std::array<bool, kOrder> seen{};
    for (const std::uint8_t image : permutation) {
      if (image >= kOrder || seen[image]) {
        throw std::runtime_error(
            "automorphism enumeration contains a non-permutation");
      }
      seen[image] = true;
    }
    std::set<Edge> image;
    for (const Edge& edge : present) {
      image.insert(permute_edge(edge, permutation));
    }
    if (image != expected) {
      throw std::runtime_error(
          "enumerated permutation does not preserve the defect graph");
    }
  }
  // Structural completeness proof:
  // * the unique 15-vertex component has a degree-6 K3;
  // * permuting that K3 permutes three arms indexed by the missing K3
  //   vertex, and each arm has independent swaps of its degree-3 pair and
  //   degree-5 pair, giving 3! * 4^3 = 384;
  // * the remaining components are two K4s, giving 4!^2 * 2 = 1152.
  // Their direct product is 442368, exactly the distinct automorphisms above.
  if (384ULL * 1'152ULL != group.size()) {
    throw std::runtime_error(
        "structural automorphism-order proof disagrees with enumeration");
  }
}

std::uint64_t choose(std::uint64_t n, int k) {
  if (k < 0 || static_cast<std::uint64_t>(k) > n) return 0;
  if (k > static_cast<int>(n) - k) k = static_cast<int>(n) - k;
  std::uint64_t result = 1;
  for (int index = 1; index <= k; ++index) {
    result =
        result * (n - static_cast<std::uint64_t>(k) + index) /
        static_cast<std::uint64_t>(index);
  }
  return result;
}

std::uint64_t fixed_subsets_from_cycles(
    const std::array<std::uint16_t, 4>& cycle_counts, int radius) {
  const std::uint64_t fixed_points = cycle_counts[1];
  if (radius == 1) return fixed_points;
  if (radius == 2) {
    return choose(fixed_points, 2) + cycle_counts[2];
  }
  if (radius == 3) {
    return choose(fixed_points, 3) +
           fixed_points * cycle_counts[2] + cycle_counts[3];
  }
  throw std::runtime_error("unsupported Burnside radius");
}

std::uint64_t burnside_orbit_count(
    int radius, const std::vector<Permutation>& group,
    const std::vector<Edge>& present) {
  const auto edge_numbers = global_edge_numbers();
  const std::set<Edge> present_set(present.begin(), present.end());
  std::vector<Edge> absent;
  for (const Edge& edge : global_edges()) {
    if (!present_set.contains(edge)) absent.push_back(edge);
  }
  const std::array<const std::vector<Edge>*, 2> colors{
      &present, &absent};
  UnsignedWide fixed_pair_sum = 0;
  for (const Permutation& permutation : group) {
    std::array<std::uint64_t, 2> fixed{};
    for (int color = 0; color < 2; ++color) {
      const std::vector<Edge>& edges = *colors[color];
      std::array<std::int16_t, 253> local{};
      local.fill(-1);
      for (std::size_t index = 0; index < edges.size(); ++index) {
        const Edge edge = edges[index];
        local[edge_numbers[edge.first][edge.second]] =
            static_cast<std::int16_t>(index);
      }
      std::vector<bool> visited(edges.size(), false);
      std::array<std::uint16_t, 4> cycle_counts{};
      for (std::size_t start = 0; start < edges.size(); ++start) {
        if (visited[start]) continue;
        std::size_t current = start;
        int length = 0;
        do {
          if (current >= edges.size() || visited[current]) {
            throw std::runtime_error(
                "induced edge permutation cycle is malformed");
          }
          visited[current] = true;
          ++length;
          const Edge image =
              permute_edge(edges[current], permutation);
          const std::int16_t next =
              local[edge_numbers[image.first][image.second]];
          if (next < 0) {
            throw std::runtime_error(
                "induced edge permutation changed edge color");
          }
          current = static_cast<std::size_t>(next);
        } while (current != start);
        if (length <= radius) ++cycle_counts[length];
      }
      fixed[color] =
          fixed_subsets_from_cycles(cycle_counts, radius);
    }
    fixed_pair_sum +=
        static_cast<UnsignedWide>(fixed[0]) * fixed[1];
  }
  if (fixed_pair_sum % group.size() != 0) {
    throw std::runtime_error("Burnside sum is not divisible by group order");
  }
  const UnsignedWide orbit_count = fixed_pair_sum / group.size();
  if (orbit_count > std::numeric_limits<std::uint64_t>::max()) {
    throw std::runtime_error("Burnside orbit count exceeds uint64");
  }
  return static_cast<std::uint64_t>(orbit_count);
}

std::vector<Subset> combinations(int n, int radius) {
  std::vector<Subset> result;
  result.reserve(static_cast<std::size_t>(choose(n, radius)));
  if (radius == 1) {
    for (int first = 0; first < n; ++first) {
      result.push_back(
          Subset{{static_cast<std::uint8_t>(first), 0, 0}, 1});
    }
  } else if (radius == 2) {
    for (int first = 0; first < n; ++first) {
      for (int second = first + 1; second < n; ++second) {
        result.push_back(
            Subset{{static_cast<std::uint8_t>(first),
                    static_cast<std::uint8_t>(second), 0},
                   2});
      }
    }
  } else if (radius == 3) {
    for (int first = 0; first < n; ++first) {
      for (int second = first + 1; second < n; ++second) {
        for (int third = second + 1; third < n; ++third) {
          result.push_back(
              Subset{{static_cast<std::uint8_t>(first),
                      static_cast<std::uint8_t>(second),
                      static_cast<std::uint8_t>(third)},
                     3});
        }
      }
    }
  } else {
    throw std::runtime_error("only radii 1, 2, and 3 are supported");
  }
  return result;
}

std::uint64_t combination_rank(const Subset& subset, int n) {
  if (subset.size == 1) return subset.values[0];
  const std::uint64_t first = subset.values[0];
  const std::uint64_t second = subset.values[1];
  if (subset.size == 2) {
    return first * (2ULL * n - first - 1ULL) / 2ULL +
           second - first - 1ULL;
  }
  if (subset.size == 3) {
    const std::uint64_t third = subset.values[2];
    const std::uint64_t first_prefix =
        choose(n, 3) - choose(n - first, 3);
    const std::uint64_t gap = second - first - 1ULL;
    const std::uint64_t second_prefix =
        gap * (static_cast<std::uint64_t>(n) - 1ULL) -
        gap * (first + 1ULL + second - 1ULL) / 2ULL;
    return first_prefix + second_prefix + third - second - 1ULL;
  }
  throw std::runtime_error("invalid subset radius");
}

Subset permute_local_subset(
    const Subset& local_subset,
    const std::vector<Edge>& local_edges,
    const std::array<std::int16_t, 253>& global_to_local,
    const std::array<std::array<std::uint8_t, kOrder>, kOrder>&
        edge_numbers,
    const Permutation& permutation) {
  Subset image;
  image.size = local_subset.size;
  for (int index = 0; index < local_subset.size; ++index) {
    const Edge edge = permute_edge(
        local_edges[local_subset.values[index]], permutation);
    const std::uint8_t global = edge_numbers[edge.first][edge.second];
    const std::int16_t local = global_to_local[global];
    if (local < 0) {
      throw std::runtime_error(
          "automorphism moved an edge across present/absent colors");
    }
    image.values[index] = static_cast<std::uint8_t>(local);
  }
  if (image.size >= 2 && image.values[0] > image.values[1]) {
    std::swap(image.values[0], image.values[1]);
  }
  if (image.size == 3) {
    if (image.values[1] > image.values[2]) {
      std::swap(image.values[1], image.values[2]);
    }
    if (image.values[0] > image.values[1]) {
      std::swap(image.values[0], image.values[1]);
    }
  }
  return image;
}

Subset local_to_global(
    const Subset& local, const std::vector<Edge>& local_edges,
    const std::array<std::array<std::uint8_t, kOrder>, kOrder>&
        edge_numbers) {
  Subset global;
  global.size = local.size;
  for (int index = 0; index < local.size; ++index) {
    const Edge edge = local_edges[local.values[index]];
    global.values[index] = edge_numbers[edge.first][edge.second];
  }
  if (global.size >= 2 && global.values[0] > global.values[1]) {
    std::swap(global.values[0], global.values[1]);
  }
  if (global.size == 3) {
    if (global.values[1] > global.values[2]) {
      std::swap(global.values[1], global.values[2]);
    }
    if (global.values[0] > global.values[1]) {
      std::swap(global.values[0], global.values[1]);
    }
  }
  return global;
}

OrbitCatalog build_orbit_catalog(
    int radius, const std::vector<Permutation>& group,
    const std::vector<Edge>& present, double heartbeat_seconds = 0.0) {
  const Clock::time_point started = Clock::now();
  Clock::time_point next_heartbeat =
      started + std::chrono::duration_cast<Clock::duration>(
                    std::chrono::duration<double>(heartbeat_seconds));
  const auto edge_numbers = global_edge_numbers();
  const std::set<Edge> present_set(present.begin(), present.end());
  std::vector<Edge> absent;
  for (const Edge& edge : global_edges()) {
    if (!present_set.contains(edge)) absent.push_back(edge);
  }
  if (present.size() != 45 || absent.size() != 208) {
    throw std::runtime_error("present/absent edge partition mismatch");
  }

  std::array<std::int16_t, 253> present_local{};
  std::array<std::int16_t, 253> absent_local{};
  present_local.fill(-1);
  absent_local.fill(-1);
  for (std::size_t index = 0; index < present.size(); ++index) {
    const Edge edge = present[index];
    present_local[edge_numbers[edge.first][edge.second]] =
        static_cast<std::int16_t>(index);
  }
  for (std::size_t index = 0; index < absent.size(); ++index) {
    const Edge edge = absent[index];
    absent_local[edge_numbers[edge.first][edge.second]] =
        static_cast<std::int16_t>(index);
  }

  const std::vector<Subset> removals =
      combinations(static_cast<int>(present.size()), radius);
  const std::vector<Subset> additions =
      combinations(static_cast<int>(absent.size()), radius);
  for (std::size_t index = 0; index < removals.size(); ++index) {
    if (combination_rank(removals[index], present.size()) != index) {
      throw std::runtime_error("removal combination rank mismatch");
    }
  }
  for (std::size_t index = 0; index < additions.size(); ++index) {
    if (combination_rank(additions[index], absent.size()) != index) {
      throw std::runtime_error("addition combination rank mismatch");
    }
  }

  OrbitCatalog catalog;
  catalog.radius = radius;
  catalog.labeled_count =
      static_cast<std::uint64_t>(removals.size()) * additions.size();
  std::vector<std::int32_t> removal_owner(removals.size(), -1);

  for (std::size_t removal_rank = 0;
       removal_rank < removals.size(); ++removal_rank) {
    if (removal_owner[removal_rank] >= 0) continue;
    const Subset removed = removals[removal_rank];
    const std::int32_t removal_orbit_id =
        static_cast<std::int32_t>(catalog.removed_orbit_count++);
    std::vector<std::uint32_t> stabilizer;
    std::uint64_t removal_orbit_size = 0;
    for (std::uint32_t group_index = 0;
         group_index < group.size(); ++group_index) {
      const Subset image = permute_local_subset(
          removed, present, present_local, edge_numbers,
          group[group_index]);
      const std::uint64_t image_rank =
          combination_rank(image, present.size());
      if (image_rank >= removal_owner.size()) {
        throw std::runtime_error("removal orbit rank overflow");
      }
      if (removal_owner[image_rank] < 0) {
        removal_owner[image_rank] = removal_orbit_id;
        ++removal_orbit_size;
      } else if (removal_owner[image_rank] != removal_orbit_id) {
        throw std::runtime_error("removal orbits overlap");
      }
      if (image == removed) stabilizer.push_back(group_index);
    }
    if (stabilizer.empty() ||
        group.size() % stabilizer.size() != 0 ||
        group.size() / stabilizer.size() != removal_orbit_size) {
      throw std::runtime_error(
          "removal orbit-stabilizer identity failed");
    }
    ++catalog.removed_stabilizer_histogram[stabilizer.size()];

    std::vector<std::int32_t> addition_owner(additions.size(), -1);
    std::int32_t addition_orbit_id = 0;
    for (std::size_t addition_rank = 0;
         addition_rank < additions.size(); ++addition_rank) {
      if (addition_owner[addition_rank] >= 0) continue;
      const Subset added = additions[addition_rank];
      std::uint64_t addition_orbit_size = 0;
      for (const std::uint32_t group_index : stabilizer) {
        const Subset image = permute_local_subset(
            added, absent, absent_local, edge_numbers,
            group[group_index]);
        const std::uint64_t image_rank =
            combination_rank(image, absent.size());
        if (image_rank >= addition_owner.size()) {
          throw std::runtime_error("addition orbit rank overflow");
        }
        if (addition_owner[image_rank] < 0) {
          addition_owner[image_rank] = addition_orbit_id;
          ++addition_orbit_size;
        } else if (addition_owner[image_rank] != addition_orbit_id) {
          throw std::runtime_error("addition orbits overlap");
        }
      }
      if (addition_orbit_size == 0 ||
          stabilizer.size() % addition_orbit_size != 0) {
        throw std::runtime_error(
            "addition orbit-stabilizer identity failed");
      }
      const std::uint64_t pair_orbit_size =
          removal_orbit_size * addition_orbit_size;
      if (pair_orbit_size == 0 ||
          group.size() % pair_orbit_size != 0) {
        throw std::runtime_error(
            "pair orbit size does not divide automorphism group");
      }
      OrbitRecord record;
      record.index = catalog.records.size();
      record.orbit_size = pair_orbit_size;
      record.removed_global =
          local_to_global(removed, present, edge_numbers);
      record.added_global =
          local_to_global(added, absent, edge_numbers);
      catalog.records.push_back(record);
      catalog.orbit_size_sum += pair_orbit_size;
      ++catalog.orbit_size_histogram[pair_orbit_size];
      ++addition_orbit_id;
    }
    const Clock::time_point now = Clock::now();
    if (heartbeat_seconds > 0.0 && now >= next_heartbeat) {
      const double elapsed =
          std::chrono::duration<double>(now - started).count();
      std::cerr
          << "{\"elapsed_seconds\":" << std::fixed
          << std::setprecision(3) << elapsed
          << ",\"event\":\"radius-orbit-catalog-heartbeat\""
          << ",\"radius\":" << radius
          << ",\"removal_rank\":" << removal_rank
          << ",\"removed_orbits\":" << catalog.removed_orbit_count
          << ",\"representatives\":" << catalog.records.size()
          << "}\n";
      next_heartbeat =
          now + std::chrono::duration_cast<Clock::duration>(
                    std::chrono::duration<double>(heartbeat_seconds));
    }
  }
  if (std::find(removal_owner.begin(), removal_owner.end(), -1) !=
      removal_owner.end()) {
    throw std::runtime_error("removal orbit catalog is incomplete");
  }
  if (catalog.orbit_size_sum != catalog.labeled_count) {
    throw std::runtime_error(
        "pair orbit multiplicities do not recover labeled family");
  }
  catalog.elapsed_seconds =
      std::chrono::duration<double>(Clock::now() - started).count();
  return catalog;
}

std::uint64_t modular_power(
    std::uint64_t base, std::uint64_t exponent,
    std::uint64_t modulus) {
  std::uint64_t result = 1;
  while (exponent != 0) {
    if ((exponent & 1U) != 0U) result = (result * base) % modulus;
    base = (base * base) % modulus;
    exponent >>= 1U;
  }
  return result;
}

bool is_prime_by_trial_division(std::uint64_t value) {
  if (value < 2) return false;
  if (value % 2 == 0) return value == 2;
  for (std::uint64_t divisor = 3;
       divisor <= value / divisor; divisor += 2) {
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
      work[row][column] =
          value >= 0
              ? static_cast<std::uint64_t>(value) % prime
              : prime - static_cast<std::uint64_t>(-value) % prime;
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

Wide reconstruct_crt(
    const std::array<std::uint64_t, kCrtPrimes.size()>& residues) {
  // Every supported row has squared norm at most 23^2 + 22*3^2 = 727.
  // Hadamard bounds every supported principal determinant by
  // 727^(23/2) < 2^110. The four-prime CRT modulus exceeds 2^123,
  // so centered reconstruction is unique with ample margin.
  UnsignedWide reconstructed = 0;
  UnsignedWide modulus = 1;
  for (std::size_t index = 0; index < kCrtPrimes.size(); ++index) {
    const std::uint64_t prime = kCrtPrimes[index];
    const std::uint64_t residue = residues[index];
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

Wide exact_determinant(const Gram& matrix, int order = kOrder) {
  std::array<std::uint64_t, kCrtPrimes.size()> residues{};
  for (std::size_t index = 0; index < kCrtPrimes.size(); ++index) {
    residues[index] =
        determinant_modulo(matrix, order, kCrtPrimes[index]);
  }
  return reconstruct_crt(residues);
}

struct ModularBase {
  std::uint64_t prime = 0;
  std::uint64_t determinant = 0;
  std::array<std::array<std::uint64_t, kOrder>, kOrder> inverse{};
};

using ModularBases = std::array<ModularBase, kCrtPrimes.size()>;

ModularBase modular_base(const Gram& base, std::uint64_t prime) {
  ModularBase result;
  result.prime = prime;
  result.determinant = determinant_modulo(base, kOrder, prime);
  if (result.determinant == 0) {
    throw std::runtime_error(
        "base Gram is singular modulo a CRT prime");
  }
  std::array<std::array<std::uint64_t, 2 * kOrder>, kOrder> work{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      const int value = base[row][column];
      work[row][column] =
          value >= 0
              ? static_cast<std::uint64_t>(value) % prime
              : prime - static_cast<std::uint64_t>(-value) % prime;
    }
    work[row][kOrder + row] = 1;
  }
  for (int column = 0; column < kOrder; ++column) {
    int pivot_row = column;
    while (pivot_row < kOrder && work[pivot_row][column] == 0) {
      ++pivot_row;
    }
    if (pivot_row == kOrder) {
      throw std::runtime_error(
          "base Gram modular inverse unexpectedly failed");
    }
    if (pivot_row != column) {
      std::swap(work[pivot_row], work[column]);
    }
    const std::uint64_t inverse_pivot =
        modular_power(work[column][column], prime - 2U, prime);
    for (int inner = 0; inner < 2 * kOrder; ++inner) {
      work[column][inner] =
          (work[column][inner] * inverse_pivot) % prime;
    }
    for (int row = 0; row < kOrder; ++row) {
      if (row == column || work[row][column] == 0) continue;
      const std::uint64_t factor = work[row][column];
      for (int inner = 0; inner < 2 * kOrder; ++inner) {
        const std::uint64_t product =
            (factor * work[column][inner]) % prime;
        work[row][inner] =
            work[row][inner] >= product
                ? work[row][inner] - product
                : work[row][inner] + prime - product;
      }
    }
  }
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      result.inverse[row][column] = work[row][kOrder + column];
    }
  }
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      std::uint64_t product = 0;
      for (int inner = 0; inner < kOrder; ++inner) {
        const int value = base[row][inner];
        const std::uint64_t residue =
            value >= 0
                ? static_cast<std::uint64_t>(value) % prime
                : prime - static_cast<std::uint64_t>(-value) % prime;
        product =
            (product + residue * result.inverse[inner][column]) %
            prime;
      }
      if (product != static_cast<std::uint64_t>(row == column)) {
        throw std::runtime_error(
            "base Gram modular inverse identity check failed");
      }
    }
  }
  return result;
}

ModularBases modular_bases(const Gram& base) {
  ModularBases result;
  for (std::size_t index = 0; index < kCrtPrimes.size(); ++index) {
    if (!is_prime_by_trial_division(kCrtPrimes[index])) {
      throw std::runtime_error("a configured CRT modulus is not prime");
    }
    result[index] = modular_base(base, kCrtPrimes[index]);
  }
  return result;
}

std::uint64_t determinant_modulo_residue_matrix(
    std::array<std::array<std::uint64_t, kOrder>, kOrder> work,
    int order, std::uint64_t prime) {
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
    }
  }
  return determinant;
}

std::array<std::uint64_t, kCrtPrimes.size()>
fast_determinant_residues(
    const OrbitRecord& record, const ModularBases& bases) {
  std::array<bool, kOrder> used{};
  const auto& edges = global_edges();
  for (int index = 0; index < record.removed_global.size; ++index) {
    const Edge edge = edges[record.removed_global.values[index]];
    used[edge.first] = true;
    used[edge.second] = true;
  }
  for (int index = 0; index < record.added_global.size; ++index) {
    const Edge edge = edges[record.added_global.values[index]];
    used[edge.first] = true;
    used[edge.second] = true;
  }
  std::array<int, kOrder> vertices{};
  std::array<int, kOrder> position{};
  position.fill(-1);
  int dimension = 0;
  for (int vertex = 0; vertex < kOrder; ++vertex) {
    if (!used[vertex]) continue;
    vertices[dimension] = vertex;
    position[vertex] = dimension++;
  }
  if (dimension <= 0 || dimension > 12) {
    throw std::runtime_error(
        "rank-update endpoint dimension is outside [1,12]");
  }
  std::array<std::array<int, kOrder>, kOrder> delta{};
  for (int index = 0; index < record.removed_global.size; ++index) {
    const Edge edge = edges[record.removed_global.values[index]];
    const int first = position[edge.first];
    const int second = position[edge.second];
    delta[first][second] -= 4;
    delta[second][first] -= 4;
  }
  for (int index = 0; index < record.added_global.size; ++index) {
    const Edge edge = edges[record.added_global.values[index]];
    const int first = position[edge.first];
    const int second = position[edge.second];
    delta[first][second] += 4;
    delta[second][first] += 4;
  }

  std::array<std::uint64_t, kCrtPrimes.size()> residues{};
  for (std::size_t prime_index = 0;
       prime_index < bases.size(); ++prime_index) {
    const ModularBase& base = bases[prime_index];
    std::array<std::array<std::uint64_t, kOrder>, kOrder> lemma{};
    // det(B + P*Delta*P^T)
    //   = det(B) det(I + Delta * P^T B^{-1} P).
    // This formulation does not require Delta to be invertible.
    for (int row = 0; row < dimension; ++row) {
      for (int column = 0; column < dimension; ++column) {
        std::uint64_t value =
            static_cast<std::uint64_t>(row == column);
        for (int inner = 0; inner < dimension; ++inner) {
          const int coefficient = delta[row][inner];
          if (coefficient == 0) continue;
          const std::uint64_t coefficient_mod =
              coefficient > 0
                  ? static_cast<std::uint64_t>(coefficient)
                  : base.prime -
                        static_cast<std::uint64_t>(-coefficient);
          const std::uint64_t term =
              coefficient_mod *
              base.inverse[vertices[inner]][vertices[column]] %
              base.prime;
          value += term;
          if (value >= base.prime) value -= base.prime;
        }
        lemma[row][column] = value;
      }
    }
    const std::uint64_t correction =
        determinant_modulo_residue_matrix(
            lemma, dimension, base.prime);
    residues[prime_index] =
        base.determinant * correction % base.prime;
  }
  return residues;
}

bool exact_positive_definite(const Gram& matrix) {
  for (int order = 1; order <= kOrder; ++order) {
    if (exact_determinant(matrix, order) <= 0) return false;
  }
  return true;
}

Wide integer_square_root(Wide value) {
  if (value < 0) {
    throw std::runtime_error("square root of negative integer");
  }
  if (value == 0) return 0;
  const UnsignedWide input = static_cast<UnsignedWide>(value);
  unsigned bits = 0;
  for (UnsignedWide copy = input; copy != 0; copy >>= 1U) ++bits;
  UnsignedWide estimate =
      UnsignedWide{1} << ((bits + 1U) / 2U);
  for (;;) {
    const UnsignedWide next = (estimate + input / estimate) >> 1U;
    if (next >= estimate) {
      while ((estimate + 1) <= input / (estimate + 1)) ++estimate;
      while (estimate > input / estimate) --estimate;
      return static_cast<Wide>(estimate);
    }
    estimate = next;
  }
}

Gram candidate_gram(const Gram& base, const OrbitRecord& record) {
  Gram result = base;
  const auto& edges = global_edges();
  for (int index = 0; index < record.removed_global.size; ++index) {
    const Edge edge = edges[record.removed_global.values[index]];
    if (result[edge.first][edge.second] != 3) {
      throw std::runtime_error(
          "catalog removal is not a base defect edge");
    }
    result[edge.first][edge.second] = -1;
    result[edge.second][edge.first] = -1;
  }
  for (int index = 0; index < record.added_global.size; ++index) {
    const Edge edge = edges[record.added_global.values[index]];
    if (result[edge.first][edge.second] != -1) {
      throw std::runtime_error(
          "catalog addition is not a base nonedge");
    }
    result[edge.first][edge.second] = 3;
    result[edge.second][edge.first] = 3;
  }
  return result;
}

std::vector<Edge> defect_edges(const Gram& gram) {
  std::vector<Edge> result;
  for (const Edge& edge : global_edges()) {
    if (gram[edge.first][edge.second] == 3) result.push_back(edge);
  }
  return result;
}

ScreenSummary screen_records(
    const Gram& base, const ModularBases& bases,
    const std::vector<OrbitRecord>& records,
    std::uint64_t shard_count, std::uint64_t shard_index,
    double heartbeat_seconds) {
  const Clock::time_point started = Clock::now();
  Clock::time_point next_heartbeat =
      started + std::chrono::duration_cast<Clock::duration>(
                    std::chrono::duration<double>(heartbeat_seconds));
  const Wide frontier_squared =
      static_cast<Wide>(kFrontierRoot) * kFrontierRoot;
  ScreenSummary summary;

  for (const OrbitRecord& record : records) {
    if (record.index % shard_count != shard_index) continue;
    const auto fast_residues =
        fast_determinant_residues(record, bases);
    const Wide determinant = reconstruct_crt(fast_residues);
    if (!summary.has_maximum ||
        determinant > summary.maximum_determinant) {
      summary.has_maximum = true;
      summary.maximum_determinant = determinant;
      summary.maximum_orbit_index = record.index;
      summary.maximum_orbit_size = record.orbit_size;
      summary.maximum_removed_global = record.removed_global;
      summary.maximum_added_global = record.added_global;
    }
    if ((record.index & 4095U) == 0U) {
      const Gram candidate = candidate_gram(base, record);
      for (std::size_t index = 0; index < kCrtPrimes.size(); ++index) {
        if (fast_residues[index] != determinant_modulo(
                candidate, kOrder, kCrtPrimes[index])) {
          throw std::runtime_error(
              "fast determinant residue spot-check failed");
        }
      }
      if (determinant != exact_determinant(candidate)) {
        throw std::runtime_error(
            "fast determinant CRT spot-check failed");
      }
    }
    ++summary.representative_count;
    summary.labeled_count += record.orbit_size;
    if (determinant > 0) {
      ++summary.positive_determinant_representatives;
      summary.positive_determinant_labeled += record.orbit_size;
    }
    if (determinant > frontier_squared) {
      ++summary.above_frontier_representatives;
      summary.above_frontier_labeled += record.orbit_size;
    }
    if (determinant > 0) {
      const Wide root = integer_square_root(determinant);
      if (root * root == determinant) {
        const Gram candidate = candidate_gram(base, record);
        if (determinant != exact_determinant(candidate)) {
          throw std::runtime_error(
              "square-hit fast determinant cross-check failed");
        }
        ++summary.square_representatives;
        summary.square_labeled += record.orbit_size;
        const std::string root_text = decimal(root);
        summary.square_root_labeled_histogram[root_text] +=
            record.orbit_size;
        ++summary.square_root_orbit_histogram[root_text];
        if (root > static_cast<Wide>(kFrontierRoot)) {
          ++summary.above_frontier_square_representatives;
          summary.above_frontier_square_labeled += record.orbit_size;
          if (root % kRequiredDivisor == 0) {
            ++summary
                 .divisible_above_frontier_square_representatives;
            summary.divisible_above_frontier_square_labeled +=
                record.orbit_size;
            const bool positive_definite =
                exact_positive_definite(candidate);
            if (positive_definite) {
              ++summary.positive_definite_route_representatives;
              summary.positive_definite_route_labeled +=
                  record.orbit_size;
              SquareHit hit;
              hit.orbit_index = record.index;
              hit.orbit_size = record.orbit_size;
              hit.determinant = determinant;
              hit.root = root;
              hit.divisible = true;
              hit.positive_definite = true;
              hit.removed_global = record.removed_global;
              hit.added_global = record.added_global;
              hit.edges = defect_edges(candidate);
              summary.route_hits.push_back(std::move(hit));
            }
          }
        }
      }
    }
    const Clock::time_point now = Clock::now();
    if (heartbeat_seconds > 0.0 && now >= next_heartbeat) {
      const double elapsed =
          std::chrono::duration<double>(now - started).count();
      std::cerr
          << "{\"elapsed_seconds\":" << std::fixed
          << std::setprecision(3) << elapsed
          << ",\"event\":\"radius3-screen-heartbeat\""
          << ",\"representatives\":" << summary.representative_count
          << ",\"route_hits\":" << summary.route_hits.size()
          << ",\"shard_count\":" << shard_count
          << ",\"shard_index\":" << shard_index << "}\n";
      next_heartbeat =
          now + std::chrono::duration_cast<Clock::duration>(
                    std::chrono::duration<double>(heartbeat_seconds));
    }
  }
  summary.elapsed_seconds =
      std::chrono::duration<double>(Clock::now() - started).count();
  return summary;
}

std::string subset_csv(const Subset& subset) {
  std::ostringstream output;
  for (int index = 0; index < subset.size; ++index) {
    if (index != 0) output << ',';
    output << static_cast<unsigned>(subset.values[index]);
  }
  return output.str();
}

Subset parse_subset_csv(const std::string& text, int radius) {
  Subset subset;
  subset.size = static_cast<std::uint8_t>(radius);
  std::istringstream input(text);
  std::string field;
  int count = 0;
  while (std::getline(input, field, ',')) {
    if (count >= radius) {
      throw std::runtime_error("catalog subset has too many entries");
    }
    const std::uint64_t value =
        parse_unsigned(field, "catalog subset entry");
    if (value >= global_edges().size()) {
      throw std::runtime_error("catalog edge number is out of range");
    }
    subset.values[count++] = static_cast<std::uint8_t>(value);
  }
  if (count != radius ||
      !std::is_sorted(
          subset.values.begin(), subset.values.begin() + radius) ||
      std::adjacent_find(
          subset.values.begin(), subset.values.begin() + radius) !=
          subset.values.begin() + radius) {
    throw std::runtime_error(
        "catalog subset must be a sorted distinct radius-tuple");
  }
  return subset;
}

std::string catalog_bytes(
    const OrbitCatalog& catalog, std::string_view reference_sha256,
    std::string_view reference_gram_sha256) {
  std::ostringstream output;
  output << "GRAM_RADIUS3_ORBIT_CATALOG_V1\n";
  output << "reference_raw_sha256 " << reference_sha256 << '\n';
  output << "reference_gram_sha256 " << reference_gram_sha256 << '\n';
  output << "group_order " << kAutomorphismGroupOrder << '\n';
  output << "radius " << catalog.radius << '\n';
  output << "labeled_count " << catalog.labeled_count << '\n';
  output << "representative_count " << catalog.records.size() << '\n';
  output << "edge_numbering lexicographic-zero-based-pairs-on-23-vertices\n";
  output << "records index orbit_size removed_edge_ids added_edge_ids\n";
  for (const OrbitRecord& record : catalog.records) {
    output << record.index << ' ' << record.orbit_size << ' '
           << subset_csv(record.removed_global) << ' '
           << subset_csv(record.added_global) << '\n';
  }
  return output.str();
}

OrbitCatalog parse_catalog(
    std::string_view bytes, std::string_view expected_reference_sha256,
    std::string_view expected_reference_gram_sha256,
    double heartbeat_seconds = 0.0) {
  const Clock::time_point started = Clock::now();
  Clock::time_point next_heartbeat =
      started + std::chrono::duration_cast<Clock::duration>(
                    std::chrono::duration<double>(heartbeat_seconds));
  std::istringstream input{std::string(bytes)};
  std::string line;
  const auto next_line = [&]() -> std::string {
    std::string result;
    if (!std::getline(input, result)) {
      throw std::runtime_error("catalog is truncated");
    }
    return result;
  };
  if (next_line() != "GRAM_RADIUS3_ORBIT_CATALOG_V1") {
    throw std::runtime_error("unsupported catalog schema");
  }
  const auto parse_header = [&](std::string_view key) {
    const std::string header = next_line();
    const std::string prefix = std::string(key) + " ";
    if (!header.starts_with(prefix)) {
      throw std::runtime_error(
          "catalog expected header " + std::string(key));
    }
    return header.substr(prefix.size());
  };
  if (parse_header("reference_raw_sha256") !=
      expected_reference_sha256) {
    throw std::runtime_error(
        "catalog reference raw SHA-256 mismatch");
  }
  if (parse_header("reference_gram_sha256") !=
      expected_reference_gram_sha256) {
    throw std::runtime_error(
        "catalog reference Gram SHA-256 mismatch");
  }
  if (parse_unsigned(
          parse_header("group_order"), "catalog group_order") !=
      kAutomorphismGroupOrder) {
    throw std::runtime_error("catalog group order mismatch");
  }
  OrbitCatalog catalog;
  catalog.radius = static_cast<int>(
      parse_unsigned(parse_header("radius"), "catalog radius"));
  catalog.labeled_count =
      parse_unsigned(parse_header("labeled_count"), "catalog labeled_count");
  const std::uint64_t representative_count = parse_unsigned(
      parse_header("representative_count"),
      "catalog representative_count");
  if (parse_header("edge_numbering") !=
      "lexicographic-zero-based-pairs-on-23-vertices") {
    throw std::runtime_error("catalog edge-numbering convention mismatch");
  }
  if (next_line() !=
      "records index orbit_size removed_edge_ids added_edge_ids") {
    throw std::runtime_error("catalog records header mismatch");
  }
  if (catalog.radius != 3 || catalog.labeled_count != kRadius3Labeled) {
    throw std::runtime_error("catalog is not the expected radius-3 family");
  }
  catalog.records.reserve(representative_count);
  std::pair<Subset, Subset> previous_record;
  bool have_previous_record = false;
  const std::vector<Edge> expected_edges = expected_reference_edges();
  const std::set<Edge> present(
      expected_edges.begin(), expected_edges.end());
  while (std::getline(input, line)) {
    if (line.empty()) {
      throw std::runtime_error("catalog contains a blank record");
    }
    std::istringstream fields(line);
    std::string index_text;
    std::string size_text;
    std::string removed_text;
    std::string added_text;
    std::string trailing;
    if (!(fields >> index_text >> size_text >> removed_text >> added_text) ||
        (fields >> trailing)) {
      throw std::runtime_error("malformed catalog record");
    }
    OrbitRecord record;
    record.index = parse_unsigned(index_text, "catalog record index");
    record.orbit_size =
        parse_unsigned(size_text, "catalog record orbit size");
    record.removed_global =
        parse_subset_csv(removed_text, catalog.radius);
    record.added_global =
        parse_subset_csv(added_text, catalog.radius);
    if (record.index != catalog.records.size() ||
        record.orbit_size == 0 ||
        kAutomorphismGroupOrder % record.orbit_size != 0) {
      throw std::runtime_error(
          "catalog record index or orbit size is invalid");
    }
    const std::pair<Subset, Subset> record_key{
        record.removed_global, record.added_global};
    if (have_previous_record && !(previous_record < record_key)) {
      throw std::runtime_error(
          "catalog representatives are not in strict canonical order");
    }
    previous_record = record_key;
    have_previous_record = true;
    for (int index = 0; index < catalog.radius; ++index) {
      if (!present.contains(
              global_edges()[record.removed_global.values[index]]) ||
          present.contains(
              global_edges()[record.added_global.values[index]])) {
        throw std::runtime_error(
            "catalog record violates base edge colors");
      }
    }
    catalog.orbit_size_sum += record.orbit_size;
    ++catalog.orbit_size_histogram[record.orbit_size];
    catalog.records.push_back(record);
    if ((catalog.records.size() & 262'143U) == 0U) {
      const Clock::time_point now = Clock::now();
      if (heartbeat_seconds > 0.0 && now >= next_heartbeat) {
        const double elapsed =
            std::chrono::duration<double>(now - started).count();
        std::cerr
            << "{\"elapsed_seconds\":" << std::fixed
            << std::setprecision(3) << elapsed
            << ",\"event\":\"radius3-catalog-parse-heartbeat\""
            << ",\"expected_representatives\":"
            << representative_count
            << ",\"parsed_representatives\":"
            << catalog.records.size() << "}\n";
        next_heartbeat =
            now + std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(heartbeat_seconds));
      }
    }
  }
  if (catalog.records.size() != representative_count ||
      catalog.orbit_size_sum != catalog.labeled_count) {
    throw std::runtime_error(
        "catalog record count or orbit-size sum mismatch");
  }
  return catalog;
}

std::string gram_bytes(const Gram& gram) {
  std::ostringstream output;
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      if (column != 0) output << ' ';
      output << gram[row][column];
    }
    output << '\n';
  }
  return output.str();
}

void append_uint_histogram(
    std::ostream& output,
    const std::map<std::uint64_t, std::uint64_t>& histogram) {
  output << '{';
  bool first = true;
  for (const auto& [key, value] : histogram) {
    if (!first) output << ',';
    first = false;
    output << json_escape(std::to_string(key)) << ':' << value;
  }
  output << '}';
}

void append_string_histogram(
    std::ostream& output,
    const std::map<std::string, std::uint64_t>& histogram) {
  output << '{';
  bool first = true;
  for (const auto& [key, value] : histogram) {
    if (!first) output << ',';
    first = false;
    output << json_escape(key) << ':' << value;
  }
  output << '}';
}

void append_edges(std::ostream& output, const std::vector<Edge>& edges) {
  output << '[';
  for (std::size_t index = 0; index < edges.size(); ++index) {
    if (index != 0) output << ',';
    output << '[' << static_cast<unsigned>(edges[index].first) + 1
           << ',' << static_cast<unsigned>(edges[index].second) + 1
           << ']';
  }
  output << ']';
}

void append_subset_edges(std::ostream& output, const Subset& subset) {
  std::vector<Edge> edges;
  for (int index = 0; index < subset.size; ++index) {
    edges.push_back(global_edges()[subset.values[index]]);
  }
  append_edges(output, edges);
}

std::string manifest_json(
    const OrbitCatalog& catalog, std::string_view catalog_sha256,
    std::string_view catalog_path, std::string_view reference_path,
    std::string_view reference_sha256,
    std::string_view reference_gram_sha256,
    std::size_t benchmark_representatives,
    double benchmark_seconds) {
  const double rate =
      benchmark_seconds > 0.0
          ? static_cast<double>(benchmark_representatives) /
                benchmark_seconds
          : 0.0;
  const double projected_seconds =
      rate > 0.0 ? static_cast<double>(catalog.records.size()) / rate
                 : std::numeric_limits<double>::infinity();
  const bool go =
      catalog.records.size() <= kMaximumRecommendedRepresentatives &&
      projected_seconds <= kMaximumRecommendedCpuSeconds;
  std::ostringstream output;
  output << std::setprecision(17);
  output << "{\"catalog\":{";
  output << "\"format\":\"GRAM_RADIUS3_ORBIT_CATALOG_V1\"";
  output << ",\"path\":" << json_escape(std::string(catalog_path));
  output << ",\"sha256\":" << json_escape(std::string(catalog_sha256));
  output << "},\"challenge_id\":\"maxdet-23-v1\"";
  output << ",\"claim_boundary\":"
         << json_escape(
                "Exact quotient under Aut(B0) of the radius-3 "
                "delete/add family; Gram screening remains a necessary "
                "condition and does not construct a sign factor.");
  output << ",\"complete\":true";
  output << ",\"elapsed_seconds\":" << catalog.elapsed_seconds;
  output << ",\"engine\":\"gram-radius3-orbits\"";
  output << ",\"family\":{";
  output << "\"absent_edge_count\":208";
  output << ",\"labeled_count\":" << catalog.labeled_count;
  output << ",\"present_edge_count\":45";
  output << ",\"radius\":3";
  output << ",\"representative_count\":" << catalog.records.size();
  output << "}";
  output << ",\"go_no_go\":{";
  output << "\"benchmark_representatives\":" << benchmark_representatives;
  output << ",\"benchmark_seconds\":" << benchmark_seconds;
  output << ",\"go\":" << (go ? "true" : "false");
  output << ",\"maximum_cpu_seconds\":" << kMaximumRecommendedCpuSeconds;
  output << ",\"maximum_representatives\":"
         << kMaximumRecommendedRepresentatives;
  output << ",\"projected_cpu_seconds\":" << projected_seconds;
  output << ",\"representatives_per_second\":" << rate;
  output << "}";
  output << ",\"normalization\":\"G=24I-J+4A\"";
  output << ",\"orbit_audit\":{";
  output << "\"automorphism_group\":"
         << json_escape(
                "((S2xS2)^3 semidirect S3) x (S4 wr S2)");
  output << ",\"automorphism_group_order\":"
         << kAutomorphismGroupOrder;
  output << ",\"orbit_size_histogram\":";
  append_uint_histogram(output, catalog.orbit_size_histogram);
  output << ",\"orbit_size_sum\":" << catalog.orbit_size_sum;
  output << ",\"removed_orbit_count\":" << catalog.removed_orbit_count;
  output << ",\"removed_stabilizer_histogram\":";
  append_uint_histogram(output, catalog.removed_stabilizer_histogram);
  output << ",\"method\":"
         << json_escape(
                "First quotient deletion subsets by the full group; "
                "for each deletion representative quotient addition "
                "subsets by its exact full stabilizer. Pair orbit size "
                "is deletion-orbit size times stabilizer-orbit size.");
  output << "}";
  output << ",\"reference\":{";
  output << "\"gram_sha256\":"
         << json_escape(std::string(reference_gram_sha256));
  output << ",\"path\":" << json_escape(std::string(reference_path));
  output << ",\"raw_sha256\":"
         << json_escape(std::string(reference_sha256));
  output << "}";
  output << ",\"schema_version\":1";
  output << ",\"termination\":\"completed\"}\n";
  return output.str();
}

std::string screen_report_json(
    const ScreenSummary& summary, std::string_view catalog_path,
    std::string_view catalog_sha256, std::string_view reference_path,
    std::string_view reference_sha256,
    std::string_view reference_gram_sha256,
    std::uint64_t shard_count, std::uint64_t shard_index,
    std::size_t total_representatives) {
  std::ostringstream output;
  output << std::setprecision(17);
  output << "{\"catalog\":{";
  output << "\"path\":" << json_escape(std::string(catalog_path));
  output << ",\"sha256\":" << json_escape(std::string(catalog_sha256));
  output << ",\"total_representatives\":" << total_representatives;
  output << "},\"challenge_id\":\"maxdet-23-v1\"";
  output << ",\"claim_boundary\":"
         << json_escape(
                "Exact Gram-orbit screening only; every routed Gram "
                "still requires an exact {-1,+1} factor.");
  output << ",\"complete\":true";
  output << ",\"completion_scope\":\"assigned-shard\"";
  output << ",\"elapsed_seconds\":" << summary.elapsed_seconds;
  output << ",\"engine\":\"gram-radius3-orbits\"";
  output << ",\"frontier_root\":"
         << json_escape(std::to_string(kFrontierRoot));
  output << ",\"frontier_squared\":"
         << json_escape(decimal(
                static_cast<Wide>(kFrontierRoot) * kFrontierRoot));
  if (summary.has_maximum) {
    const Wide maximum_root =
        summary.maximum_determinant > 0
            ? integer_square_root(summary.maximum_determinant)
            : 0;
    output << ",\"maximum\":{";
    output << "\"added_edges\":";
    append_subset_edges(output, summary.maximum_added_global);
    output << ",\"determinant\":"
           << json_escape(decimal(summary.maximum_determinant));
    output << ",\"is_square\":"
           << (maximum_root * maximum_root ==
                       summary.maximum_determinant
                   ? "true"
                   : "false");
    output << ",\"orbit_index\":" << summary.maximum_orbit_index;
    output << ",\"orbit_size\":" << summary.maximum_orbit_size;
    output << ",\"removed_edges\":";
    append_subset_edges(output, summary.maximum_removed_global);
    output << "}";
  }
  output << ",\"normalization\":\"G=24I-J+4A\"";
  output << ",\"reference\":{";
  output << "\"gram_sha256\":"
         << json_escape(std::string(reference_gram_sha256));
  output << ",\"path\":" << json_escape(std::string(reference_path));
  output << ",\"raw_sha256\":"
         << json_escape(std::string(reference_sha256));
  output << "}";
  output << ",\"schema_version\":1";
  output << ",\"whole_family_complete\":"
         << (shard_count == 1 ? "true" : "false");
  output << ",\"shard_count\":" << shard_count;
  output << ",\"shard_index\":" << shard_index;
  output << ",\"square_root_labeled_histogram\":";
  append_string_histogram(
      output, summary.square_root_labeled_histogram);
  output << ",\"square_root_orbit_histogram\":";
  append_string_histogram(
      output, summary.square_root_orbit_histogram);
  output << ",\"statistics\":{";
  output << "\"above_frontier_labeled\":"
         << summary.above_frontier_labeled;
  output << ",\"above_frontier_representatives\":"
         << summary.above_frontier_representatives;
  output << ",\"above_frontier_square_labeled\":"
         << summary.above_frontier_square_labeled;
  output << ",\"above_frontier_square_representatives\":"
         << summary.above_frontier_square_representatives;
  output << ",\"divisible_above_frontier_square_labeled\":"
         << summary.divisible_above_frontier_square_labeled;
  output << ",\"divisible_above_frontier_square_representatives\":"
         << summary.divisible_above_frontier_square_representatives;
  output << ",\"labeled_count\":" << summary.labeled_count;
  output << ",\"positive_definite_route_labeled\":"
         << summary.positive_definite_route_labeled;
  output << ",\"positive_definite_route_representatives\":"
         << summary.positive_definite_route_representatives;
  output << ",\"positive_determinant_labeled\":"
         << summary.positive_determinant_labeled;
  output << ",\"positive_determinant_representatives\":"
         << summary.positive_determinant_representatives;
  output << ",\"representative_count\":"
         << summary.representative_count;
  output << ",\"square_labeled\":" << summary.square_labeled;
  output << ",\"square_representatives\":"
         << summary.square_representatives;
  output << "}";
  output << ",\"termination\":\"completed\"}\n";
  return output.str();
}

void append_route_hit(std::ostream& output, const SquareHit& hit) {
  output << "{\"added_edges\":";
  append_subset_edges(output, hit.added_global);
  output << ",\"determinant\":" << json_escape(decimal(hit.determinant));
  output << ",\"divisible_by_2_22\":"
         << (hit.divisible ? "true" : "false");
  output << ",\"edge_count\":" << hit.edges.size();
  output << ",\"edges\":";
  append_edges(output, hit.edges);
  output << ",\"orbit_index\":" << hit.orbit_index;
  output << ",\"orbit_size\":" << hit.orbit_size;
  output << ",\"positive_definite\":"
         << (hit.positive_definite ? "true" : "false");
  output << ",\"qualified\":"
         << (hit.divisible && hit.positive_definite ? "true" : "false");
  output << ",\"removed_edges\":";
  append_subset_edges(output, hit.removed_global);
  output << ",\"square_root\":" << json_escape(decimal(hit.root));
  output << '}';
}

std::string route_snapshot_json(
    const ScreenSummary& summary, std::string_view catalog_path,
    std::string_view catalog_sha256, std::uint64_t shard_count,
    std::uint64_t shard_index) {
  std::ostringstream output;
  output << std::setprecision(17);
  output << "{\"catalog_path\":"
         << json_escape(std::string(catalog_path));
  output << ",\"catalog_sha256\":"
         << json_escape(std::string(catalog_sha256));
  output << ",\"challenge_id\":\"maxdet-23-v1\"";
  output << ",\"claim_boundary\":"
         << json_escape(
                "Truthful compatibility snapshot from exact radius-3 "
                "orbit screening; routed hits remain Gram-only.");
  output << ",\"complete\":true";
  output << ",\"completion_scope\":\"assigned-shard\"";
  output << ",\"elapsed_seconds\":" << summary.elapsed_seconds;
  output << ",\"engine\":\"gram-radius3-orbits\"";
  output << ",\"frontier_root\":"
         << json_escape(std::to_string(kFrontierRoot));
  output << ",\"frontier_squared\":"
         << json_escape(decimal(
                static_cast<Wide>(kFrontierRoot) * kFrontierRoot));
  output << ",\"hits\":[";
  for (std::size_t index = 0; index < summary.route_hits.size(); ++index) {
    if (index != 0) output << ',';
    append_route_hit(output, summary.route_hits[index]);
  }
  output << ']';
  output << ",\"mode\":\"reference-radius3-exact-orbits\"";
  output << ",\"normalization\":\"G=24I-J+4A\"";
  output << ",\"parameters\":{";
  output << "\"max_stored_hits\":" << summary.route_hits.size();
  output << ",\"shard_count\":" << shard_count;
  output << ",\"shard_index\":" << shard_index;
  output << "}";
  output << ",\"schema_version\":1";
  output << ",\"whole_family_complete\":"
         << (shard_count == 1 ? "true" : "false");
  output << ",\"statistics\":{";
  output << "\"exact_squares\":" << summary.square_representatives;
  output << ",\"qualified_survivors\":"
         << summary.positive_definite_route_representatives;
  output << ",\"unrecorded_square_observations\":0";
  output << "}";
  output << ",\"termination\":\"completed\"}\n";
  return output.str();
}

void validate_regression(
    int radius, const OrbitCatalog& catalog,
    const ScreenSummary& summary) {
  const std::uint64_t expected_labeled =
      radius == 1 ? kRadius1Labeled : kRadius2Labeled;
  const std::string expected_root =
      radius == 1 ? "2743271424000000" : "2779447296000000";
  if (catalog.labeled_count != expected_labeled ||
      catalog.orbit_size_sum != expected_labeled ||
      summary.labeled_count != expected_labeled ||
      summary.square_labeled != 12 ||
      summary.square_root_labeled_histogram !=
          std::map<std::string, std::uint64_t>{{expected_root, 12}}) {
    throw std::runtime_error(
        "radius-" + std::to_string(radius) +
        " labeled square-root regression failed");
  }
}

void differential_checks(
    const Gram& base, const ModularBases& bases,
    const OrbitCatalog& radius1,
    const OrbitCatalog& radius2) {
  constexpr std::uint64_t independent_prime = 2'147'483'563ULL;
  const std::array<const OrbitCatalog*, 2> catalogs{
      &radius1, &radius2};
  for (const OrbitCatalog* catalog : catalogs) {
    if (catalog->records.empty()) {
      throw std::runtime_error("differential catalog is empty");
    }
    for (std::size_t index = 0;
         index < catalog->records.size(); ++index) {
      const Gram candidate =
          candidate_gram(base, catalog->records[index]);
      const auto fast_residues =
          fast_determinant_residues(catalog->records[index], bases);
      for (std::size_t prime_index = 0;
           prime_index < kCrtPrimes.size(); ++prime_index) {
        if (fast_residues[prime_index] != determinant_modulo(
                candidate, kOrder, kCrtPrimes[prime_index])) {
          throw std::runtime_error(
              "fast and full determinant residues disagree");
        }
      }
      const Wide crt = reconstruct_crt(fast_residues);
      if (crt != exact_determinant(candidate)) {
        throw std::runtime_error(
            "fast and full CRT determinants disagree");
      }
      Wide reduced = crt % static_cast<Wide>(independent_prime);
      if (reduced < 0) reduced += independent_prime;
      if (static_cast<std::uint64_t>(reduced) !=
          determinant_modulo(candidate, kOrder, independent_prime)) {
        throw std::runtime_error(
            "CRT determinant disagrees with independent fifth-prime residue");
      }
    }
  }
}

void differential_radius3_checks(
    const Gram& base, const ModularBases& bases,
    const OrbitCatalog& catalog) {
  if (catalog.records.empty()) {
    throw std::runtime_error("radius-3 differential catalog is empty");
  }
  std::uint64_t state = 0x9e3779b97f4a7c15ULL;
  for (int sample = 0; sample < 256; ++sample) {
    state = state * 6'364'136'223'846'793'005ULL +
            1'442'695'040'888'963'407ULL;
    const std::size_t index =
        static_cast<std::size_t>(state % catalog.records.size());
    const OrbitRecord& record = catalog.records[index];
    const Gram candidate = candidate_gram(base, record);
    const auto fast_residues =
        fast_determinant_residues(record, bases);
    for (std::size_t prime_index = 0;
         prime_index < kCrtPrimes.size(); ++prime_index) {
      if (fast_residues[prime_index] != determinant_modulo(
              candidate, kOrder, kCrtPrimes[prime_index])) {
        throw std::runtime_error(
            "random radius-3 determinant residue check failed");
      }
    }
    if (reconstruct_crt(fast_residues) !=
        exact_determinant(candidate)) {
      throw std::runtime_error(
          "random radius-3 CRT determinant check failed");
    }
  }
}

void run_self_test(
    const Gram& base, const ModularBases& bases,
    const std::vector<Permutation>& group,
    const std::vector<Edge>& present) {
  if (sha256("abc") !=
      "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb"
      "410ff61f20015ad") {
    throw std::runtime_error("SHA-256 self-test failed");
  }
  for (const std::uint64_t prime : kCrtPrimes) {
    if (!is_prime_by_trial_division(prime)) {
      throw std::runtime_error("a configured CRT modulus is not prime");
    }
  }
  const Wide frontier_squared =
      static_cast<Wide>(kFrontierRoot) * kFrontierRoot;
  if (exact_determinant(base) != frontier_squared ||
      !exact_positive_definite(base)) {
    throw std::runtime_error(
        "published Gram exact determinant/PD regression failed");
  }

  const OrbitCatalog radius1 =
      build_orbit_catalog(1, group, present);
  const OrbitCatalog radius2 =
      build_orbit_catalog(2, group, present);
  const ScreenSummary screen1 =
      screen_records(base, bases, radius1.records, 1, 0, 0.0);
  const ScreenSummary screen2 =
      screen_records(base, bases, radius2.records, 1, 0, 0.0);
  validate_regression(1, radius1, screen1);
  validate_regression(2, radius2, screen2);
  if (burnside_orbit_count(1, group, present) !=
          radius1.records.size() ||
      burnside_orbit_count(2, group, present) !=
          radius2.records.size()) {
    throw std::runtime_error(
        "Burnside and stabilizer-fiber orbit counts disagree");
  }
  differential_checks(base, bases, radius1, radius2);

  std::cout
      << "SELF-TEST PASS"
      << " group_order=" << group.size()
      << " radius1_orbits=" << radius1.records.size()
      << " radius1_labeled=" << radius1.labeled_count
      << " radius1_square_labeled=" << screen1.square_labeled
      << " radius2_orbits=" << radius2.records.size()
      << " radius2_labeled=" << radius2.labeled_count
      << " radius2_square_labeled=" << screen2.square_labeled
      << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    const std::string reference_bytes =
        read_file_bytes(arguments.reference);
    const std::string reference_sha256 = sha256(reference_bytes);
    const SignMatrix reference = parse_sign_matrix(reference_bytes);
    const Gram base = gram_of(reference);
    const std::string reference_gram_sha256 = sha256(gram_bytes(base));
    const std::vector<Edge> present = validate_reference_gram(base);
    const std::vector<Permutation> group = generate_automorphism_group();
    validate_automorphism_group(group, present);
    const ModularBases bases = modular_bases(base);

    if (arguments.mode == Mode::kSelfTest) {
      run_self_test(base, bases, group, present);
      return 0;
    }

    if (arguments.mode == Mode::kCount) {
      const Clock::time_point started = Clock::now();
      const std::uint64_t count =
          burnside_orbit_count(3, group, present);
      const double elapsed =
          std::chrono::duration<double>(Clock::now() - started).count();
      std::cout
          << "BURNSIDE COUNT"
          << " radius=3"
          << " representatives=" << count
          << " labeled=" << kRadius3Labeled
          << " go_count="
          << (count <= kMaximumRecommendedRepresentatives ? "true"
                                                          : "false")
          << " elapsed_seconds=" << std::fixed
          << std::setprecision(6) << elapsed << '\n';
      return 0;
    }

    if (arguments.mode == Mode::kCatalog) {
      OrbitCatalog catalog = build_orbit_catalog(
          3, group, present, arguments.heartbeat_seconds);
      if (catalog.labeled_count != kRadius3Labeled) {
        throw std::runtime_error(
            "radius-3 labeled count regression failed");
      }
      const std::uint64_t burnside_count =
          burnside_orbit_count(3, group, present);
      if (catalog.records.size() != burnside_count) {
        throw std::runtime_error(
            "radius-3 stabilizer catalog disagrees with Burnside count");
      }
      differential_radius3_checks(base, bases, catalog);
      const std::string bytes = catalog_bytes(
          catalog, reference_sha256, reference_gram_sha256);
      const std::string catalog_hash = sha256(bytes);

      const std::size_t benchmark_count =
          std::min<std::size_t>(catalog.records.size(), 10'000);
      std::vector<OrbitRecord> benchmark_records(
          catalog.records.begin(),
          catalog.records.begin() + benchmark_count);
      const ScreenSummary benchmark =
          screen_records(base, bases, benchmark_records, 1, 0, 0.0);

      atomic_write(arguments.catalog, bytes);
      const std::string manifest = manifest_json(
          catalog, catalog_hash, arguments.catalog.string(),
          arguments.reference.string(), reference_sha256,
          reference_gram_sha256,
          benchmark.representative_count, benchmark.elapsed_seconds);
      atomic_write(arguments.manifest, manifest);
      std::cout
          << "CATALOG COMPLETE"
          << " representatives=" << catalog.records.size()
          << " labeled=" << catalog.labeled_count
          << " removed_orbits=" << catalog.removed_orbit_count
          << " catalog_sha256=" << catalog_hash
          << " elapsed_seconds=" << std::fixed
          << std::setprecision(6) << catalog.elapsed_seconds
          << '\n';
      return 0;
    }

    const std::string bytes = read_file_bytes(arguments.catalog);
    const std::string catalog_hash = sha256(bytes);
    if (catalog_hash != arguments.catalog_sha256) {
      throw std::runtime_error(
          "catalog SHA-256 does not match --catalog-sha256");
    }
    const OrbitCatalog catalog = parse_catalog(
        bytes, reference_sha256, reference_gram_sha256,
        arguments.heartbeat_seconds);
    differential_radius3_checks(base, bases, catalog);
    const ScreenSummary summary = screen_records(
        base, bases, catalog.records, arguments.shard_count,
        arguments.shard_index, arguments.heartbeat_seconds);
    atomic_write(
        arguments.output,
        screen_report_json(
            summary, arguments.catalog.string(), catalog_hash,
            arguments.reference.string(), reference_sha256,
            reference_gram_sha256, arguments.shard_count,
            arguments.shard_index, catalog.records.size()));
    if (!arguments.route_snapshot.empty()) {
      atomic_write(
          arguments.route_snapshot,
          route_snapshot_json(
              summary, arguments.catalog.string(), catalog_hash,
              arguments.shard_count, arguments.shard_index));
    }
    std::cout
        << "SCREEN COMPLETE"
        << " representatives=" << summary.representative_count
        << " labeled=" << summary.labeled_count
        << " squares=" << summary.square_representatives
        << " above_frontier_squares="
        << summary.above_frontier_square_representatives
        << " route_hits="
        << summary.positive_definite_route_representatives
        << " elapsed_seconds=" << std::fixed
        << std::setprecision(6) << summary.elapsed_seconds
        << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "gram_radius3_orbits: " << error.what() << '\n';
    return 1;
  }
}
