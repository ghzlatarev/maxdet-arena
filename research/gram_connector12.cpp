#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using Wide = __int128_t;
using UnsignedWide = __uint128_t;

constexpr int kOrder = 23;
constexpr int kBlockCount = 6;
constexpr std::array<int, kBlockCount> kBlockSizes{3, 4, 4, 4, 4, 4};
constexpr std::array<int, kBlockCount> kBlockOffsets{0, 3, 7, 11, 15, 19};
constexpr std::uint64_t kExpectedLabeledConfigurations = 12'072'240;
constexpr std::uint64_t kFrontierRoot = 2'779'447'296'000'000ULL;
constexpr std::uint64_t kBaseAutomorphismGroupOrder = 5'733'089'280ULL;

using Matrix = std::array<std::array<int, kOrder>, kOrder>;

struct Endpoint {
  std::uint8_t block = 0;
  std::uint8_t mask = 0;
};

struct Connector {
  Endpoint first;
  Endpoint second;
};

using Configuration = std::array<Connector, 3>;

struct BlockPair {
  std::uint8_t first = 0;
  std::uint8_t second = 0;
};

struct OrbitData {
  std::uint64_t labeled_count = 0;
};

struct SquareHit {
  std::uint64_t key = 0;
  Wide determinant = 0;
  Wide square_root = 0;
  bool divisible_by_2_22 = false;
  bool positive_definite = false;
  std::vector<std::pair<int, int>> edges;
  std::vector<std::pair<int, int>> added_edges;
};

struct Arguments {
  std::filesystem::path output;
  std::filesystem::path route_snapshot;
  double heartbeat_seconds = 10.0;
};

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

double parse_positive_double(
    const std::string& text, std::string_view option) {
  std::size_t consumed = 0;
  const double value = std::stod(text, &consumed);
  if (consumed != text.size() || !(value > 0.0) ||
      !std::isfinite(value)) {
    throw std::runtime_error(
        std::string(option) + " must be a positive finite number");
  }
  return value;
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
    if (option == "--output") {
      arguments.output = value();
    } else if (option == "--route-snapshot") {
      arguments.route_snapshot = value();
    } else if (option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds =
          parse_positive_double(value(), option);
    } else if (option == "--help") {
      std::cout
          << "Usage: gram_connector12 --output FILE "
             "[--route-snapshot FILE] [--heartbeat-seconds S]\n\n"
          << "Exactly enumerate, modulo the ideal block automorphism group, "
             "three K2,2 connectors on three distinct block pairs of "
             "K3 disjoint-union 5K4.\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (arguments.output.empty()) {
    throw std::runtime_error("--output is required");
  }
  if (!arguments.route_snapshot.empty() &&
      std::filesystem::absolute(arguments.output).lexically_normal() ==
          std::filesystem::absolute(arguments.route_snapshot)
              .lexically_normal()) {
    throw std::runtime_error(
        "--route-snapshot must differ from --output");
  }
  return arguments;
}

void atomic_write(
    const std::filesystem::path& path, const std::string& contents) {
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

std::vector<int> valid_masks(int block_size) {
  std::vector<int> result;
  for (int mask = 0; mask < (1 << block_size); ++mask) {
    if (std::popcount(static_cast<unsigned>(mask)) == 2) {
      result.push_back(mask);
    }
  }
  return result;
}

std::uint8_t permute_mask(
    int mask, const std::array<int, 4>& permutation, int size) {
  int result = 0;
  for (int vertex = 0; vertex < size; ++vertex) {
    if ((mask >> vertex) & 1) {
      result |= 1 << permutation[vertex];
    }
  }
  return static_cast<std::uint8_t>(result);
}

struct SequenceCanonicalizer {
  // The low four bits hold the first mask, then the second and third.
  // Tables cover every possible packed input; only sequences of 2-subsets
  // are queried.
  std::array<std::array<std::uint16_t, 4096>, 4> size3{};
  std::array<std::array<std::uint16_t, 4096>, 4> size4{};

  SequenceCanonicalizer() {
    fill_tables(3, size3);
    fill_tables(4, size4);
  }

  void fill_tables(
      int size,
      std::array<std::array<std::uint16_t, 4096>, 4>& tables) {
    const std::vector<int> masks = valid_masks(size);
    std::array<int, 4> permutation{0, 1, 2, 3};
    for (int length = 1; length <= 3; ++length) {
      int count = 1;
      for (int index = 0; index < length; ++index) {
        count *= static_cast<int>(masks.size());
      }
      for (int assignment = 0; assignment < count; ++assignment) {
        int copy = assignment;
        std::array<int, 3> sequence{};
        std::uint16_t packed = 0;
        for (int index = 0; index < length; ++index) {
          sequence[index] =
              masks[copy % static_cast<int>(masks.size())];
          copy /= static_cast<int>(masks.size());
          packed |= static_cast<std::uint16_t>(
              sequence[index] << (4 * index));
        }

        std::uint16_t best = std::numeric_limits<std::uint16_t>::max();
        std::iota(permutation.begin(), permutation.end(), 0);
        do {
          std::uint16_t candidate = 0;
          for (int index = 0; index < length; ++index) {
            candidate |= static_cast<std::uint16_t>(
                permute_mask(sequence[index], permutation, size)
                << (4 * index));
          }
          // Compare sequences lexicographically, not as little-endian
          // integers.
          bool smaller = false;
          bool larger = false;
          for (int index = 0; index < length; ++index) {
            const int left = (candidate >> (4 * index)) & 15;
            const int right = (best >> (4 * index)) & 15;
            if (left < right) {
              smaller = true;
              break;
            }
            if (left > right) {
              larger = true;
              break;
            }
          }
          if (best == std::numeric_limits<std::uint16_t>::max() ||
              (smaller && !larger)) {
            best = candidate;
          }
        } while (std::next_permutation(
            permutation.begin(), permutation.begin() + size));
        tables[length][packed] = best;
      }
    }
  }

  std::uint16_t operator()(
      int block_size, int length, std::uint16_t packed) const {
    return block_size == 3 ? size3[length][packed]
                           : size4[length][packed];
  }
};

constexpr std::array<std::array<int, 3>, 6> kConnectorOrders{{
    {{0, 1, 2}},
    {{0, 2, 1}},
    {{1, 0, 2}},
    {{1, 2, 0}},
    {{2, 0, 1}},
    {{2, 1, 0}},
}};

std::uint64_t canonical_key(
    const Configuration& configuration,
    const SequenceCanonicalizer& sequence_canonicalizer) {
  std::uint64_t best = std::numeric_limits<std::uint64_t>::max();

  for (const auto& order : kConnectorOrders) {
    for (int orientation_mask = 0; orientation_mask < 8;
         ++orientation_mask) {
      std::array<Endpoint, 6> endpoints{};
      for (int position = 0; position < 3; ++position) {
        const Connector& connector = configuration[order[position]];
        const bool reverse = (orientation_mask >> position) & 1;
        endpoints[2 * position] =
            reverse ? connector.second : connector.first;
        endpoints[2 * position + 1] =
            reverse ? connector.first : connector.second;
      }

      std::array<int, kBlockCount> normalized_block{};
      normalized_block.fill(-1);
      normalized_block[0] = 0;
      int next_block = 1;
      for (const Endpoint& endpoint : endpoints) {
        if (normalized_block[endpoint.block] < 0) {
          normalized_block[endpoint.block] = next_block++;
        }
      }

      std::array<std::uint8_t, 6> canonical_masks{};
      for (int block = 0; block < kBlockCount; ++block) {
        std::array<int, 3> positions{};
        int length = 0;
        std::uint16_t packed = 0;
        for (int position = 0; position < 6; ++position) {
          if (endpoints[position].block != block) continue;
          positions[length] = position;
          packed |= static_cast<std::uint16_t>(
              endpoints[position].mask << (4 * length));
          ++length;
        }
        if (length == 0) continue;
        const std::uint16_t canonical = sequence_canonicalizer(
            kBlockSizes[block], length, packed);
        for (int index = 0; index < length; ++index) {
          canonical_masks[positions[index]] =
              static_cast<std::uint8_t>(
                  (canonical >> (4 * index)) & 15);
        }
      }

      std::uint64_t key = 0;
      for (int position = 0; position < 6; ++position) {
        const std::uint64_t field =
            (static_cast<std::uint64_t>(
                 normalized_block[endpoints[position].block])
             << 4) |
            canonical_masks[position];
        key = (key << 7) | field;
      }
      best = std::min(best, key);
    }
  }
  return best;
}

Configuration decode_key(std::uint64_t key) {
  std::array<std::uint8_t, 6> fields{};
  for (int position = 5; position >= 0; --position) {
    fields[position] = static_cast<std::uint8_t>(key & 127U);
    key >>= 7U;
  }
  Configuration result{};
  for (int index = 0; index < 3; ++index) {
    result[index].first = Endpoint{
        static_cast<std::uint8_t>(fields[2 * index] >> 4),
        static_cast<std::uint8_t>(fields[2 * index] & 15)};
    result[index].second = Endpoint{
        static_cast<std::uint8_t>(fields[2 * index + 1] >> 4),
        static_cast<std::uint8_t>(fields[2 * index + 1] & 15)};
  }
  return result;
}

std::vector<std::pair<int, int>> added_edges(
    const Configuration& configuration) {
  std::vector<std::pair<int, int>> result;
  for (const Connector& connector : configuration) {
    for (int left = 0; left < kBlockSizes[connector.first.block]; ++left) {
      if (!((connector.first.mask >> left) & 1)) continue;
      for (int right = 0; right < kBlockSizes[connector.second.block];
           ++right) {
        if (!((connector.second.mask >> right) & 1)) continue;
        int first = kBlockOffsets[connector.first.block] + left;
        int second = kBlockOffsets[connector.second.block] + right;
        if (first > second) std::swap(first, second);
        result.emplace_back(first, second);
      }
    }
  }
  std::sort(result.begin(), result.end());
  if (result.size() != 12 ||
      std::adjacent_find(result.begin(), result.end()) != result.end()) {
    throw std::runtime_error(
        "decoded connector configuration does not have 12 distinct edges");
  }
  return result;
}

std::vector<std::pair<int, int>> all_edges(
    const Configuration& configuration) {
  std::vector<std::pair<int, int>> result;
  for (int block = 0; block < kBlockCount; ++block) {
    for (int left = 0; left < kBlockSizes[block]; ++left) {
      for (int right = left + 1; right < kBlockSizes[block]; ++right) {
        result.emplace_back(
            kBlockOffsets[block] + left,
            kBlockOffsets[block] + right);
      }
    }
  }
  const auto extra = added_edges(configuration);
  result.insert(result.end(), extra.begin(), extra.end());
  std::sort(result.begin(), result.end());
  if (result.size() != 45 ||
      std::adjacent_find(result.begin(), result.end()) != result.end()) {
    throw std::runtime_error(
        "defect graph does not have 45 distinct edges");
  }
  return result;
}

Matrix gram(const Configuration& configuration) {
  Matrix result{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      result[row][column] = row == column ? 23 : -1;
    }
  }
  for (const auto& [first, second] : all_edges(configuration)) {
    result[first][second] = 3;
    result[second][first] = 3;
  }
  return result;
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
    while (pivot_row < order && work[pivot_row][column] == 0) {
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

Wide exact_determinant(
    const Matrix& matrix, int order = kOrder) {
  // Every row of every supported principal submatrix has squared norm at
  // most 23^2 + 22*3^2 = 727. Hadamard therefore bounds every determinant
  // by 727^(23/2) < 2^110. The product of these four primes exceeds 2^124,
  // so symmetric CRT reconstruction is exact and unique.
  constexpr std::array<std::uint64_t, 4> primes{
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
  if (reconstructed > modulus / 2U) {
    return static_cast<Wide>(reconstructed) -
           static_cast<Wide>(modulus);
  }
  return static_cast<Wide>(reconstructed);
}

bool exact_positive_definite(const Matrix& matrix) {
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

Configuration published_configuration() {
  return Configuration{{
      Connector{Endpoint{0, 0b101}, Endpoint{1, 0b1100}},
      Connector{Endpoint{0, 0b110}, Endpoint{2, 0b1100}},
      Connector{Endpoint{0, 0b011}, Endpoint{3, 0b1100}},
  }};
}

std::vector<BlockPair> block_pairs() {
  std::vector<BlockPair> result;
  for (int first = 0; first < kBlockCount; ++first) {
    for (int second = first + 1; second < kBlockCount; ++second) {
      result.push_back(BlockPair{
          static_cast<std::uint8_t>(first),
          static_cast<std::uint8_t>(second)});
    }
  }
  return result;
}

std::vector<Connector> connector_options(const BlockPair& pair) {
  std::vector<Connector> result;
  for (const int first_mask : valid_masks(kBlockSizes[pair.first])) {
    for (const int second_mask : valid_masks(kBlockSizes[pair.second])) {
      result.push_back(Connector{
          Endpoint{pair.first, static_cast<std::uint8_t>(first_mask)},
          Endpoint{pair.second, static_cast<std::uint8_t>(second_mask)}});
    }
  }
  return result;
}

void append_edge_array(
    std::ostream& output,
    const std::vector<std::pair<int, int>>& edges) {
  output << '[';
  for (std::size_t index = 0; index < edges.size(); ++index) {
    if (index != 0) output << ',';
    output << '[' << edges[index].first + 1
           << ',' << edges[index].second + 1 << ']';
  }
  output << ']';
}

void append_hit(std::ostream& output, const SquareHit& hit) {
  output << "{\"added_edges\":";
  append_edge_array(output, hit.added_edges);
  output << ",\"canonical_key\":" << json_escape(std::to_string(hit.key));
  output << ",\"determinant\":" << json_escape(decimal(hit.determinant));
  output << ",\"divisible_by_2_22\":"
         << (hit.divisible_by_2_22 ? "true" : "false");
  output << ",\"edge_count\":" << hit.edges.size();
  output << ",\"edges\":";
  append_edge_array(output, hit.edges);
  output << ",\"positive_definite\":"
         << (hit.positive_definite ? "true" : "false");
  output << ",\"qualified\":"
         << (hit.divisible_by_2_22 && hit.positive_definite ? "true"
                                                            : "false");
  output << ",\"square_root\":"
         << json_escape(decimal(hit.square_root)) << '}';
}

std::string route_snapshot_json(
    const std::vector<SquareHit>& route_hits,
    std::uint64_t exact_squares,
    std::uint64_t qualified_survivors,
    std::uint64_t orbit_count,
    double elapsed_seconds) {
  const Wide frontier_squared =
      static_cast<Wide>(kFrontierRoot) *
      static_cast<Wide>(kFrontierRoot);
  std::ostringstream output;
  output << std::setprecision(17);
  output << "{\"challenge_id\":\"maxdet-23-v1\"";
  output << ",\"claim_boundary\":"
         << json_escape(
                "Compatibility snapshot produced by gram-connector12; "
                "survivors still require exact sign-matrix factorization.");
  output << ",\"complete\":true";
  output << ",\"elapsed_seconds\":" << elapsed_seconds;
  // Existing exact Hasse and shell tools validate this schema label. The
  // producer and mode fields preserve the actual provenance.
  output << ",\"engine\":\"gram-tabu\"";
  output << ",\"frontier_root\":"
         << json_escape(std::to_string(kFrontierRoot));
  output << ",\"frontier_squared\":"
         << json_escape(decimal(frontier_squared));
  output << ",\"hits\":[";
  for (std::size_t index = 0; index < route_hits.size(); ++index) {
    if (index != 0) output << ',';
    append_hit(output, route_hits[index]);
  }
  output << ']';
  output << ",\"mode\":\"connector12-exact-orbit-enumeration\"";
  output << ",\"normalization\":\"G=24I-J+4A\"";
  output << ",\"parameters\":{";
  output << "\"max_stored_hits\":" << route_hits.size();
  output << ",\"orbit_count\":" << orbit_count;
  output << ",\"producer_engine\":\"gram-connector12\"}";
  output << ",\"schema_version\":1";
  output << ",\"statistics\":{";
  output << "\"exact_squares\":" << exact_squares;
  output << ",\"qualified_survivors\":" << qualified_survivors;
  output << ",\"unrecorded_square_observations\":0}";
  output << ",\"termination\":\"completed\"}\n";
  return output.str();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    const Clock::time_point started = Clock::now();
    Clock::time_point next_heartbeat =
        started + std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(
                          arguments.heartbeat_seconds));

    const SequenceCanonicalizer sequence_canonicalizer;
    const Configuration published = published_configuration();
    const std::uint64_t published_key =
        canonical_key(published, sequence_canonicalizer);
    const Wide published_determinant =
        exact_determinant(gram(published));
    const Wide frontier_squared =
        static_cast<Wide>(kFrontierRoot) *
        static_cast<Wide>(kFrontierRoot);
    if (published_determinant != frontier_squared) {
      throw std::runtime_error(
          "published three-K2,2 reproduction determinant mismatch");
    }

    const std::vector<BlockPair> pairs = block_pairs();
    std::array<std::vector<Connector>, 15> options{};
    for (std::size_t index = 0; index < pairs.size(); ++index) {
      options[index] = connector_options(pairs[index]);
    }

    std::unordered_map<std::uint64_t, OrbitData> orbits;
    std::uint64_t labeled = 0;
    for (int first_pair = 0; first_pair < 13; ++first_pair) {
      for (int second_pair = first_pair + 1; second_pair < 14;
           ++second_pair) {
        for (int third_pair = second_pair + 1; third_pair < 15;
             ++third_pair) {
          for (const Connector& first : options[first_pair]) {
            for (const Connector& second : options[second_pair]) {
              for (const Connector& third : options[third_pair]) {
                const Configuration configuration{{first, second, third}};
                const std::uint64_t key =
                    canonical_key(configuration, sequence_canonicalizer);
                ++orbits[key].labeled_count;
                ++labeled;
              }
            }
          }
          const Clock::time_point now = Clock::now();
          if (now >= next_heartbeat) {
            const double elapsed =
                std::chrono::duration<double>(now - started).count();
            std::cerr
                << "{\"elapsed_seconds\":" << std::fixed
                << std::setprecision(3) << elapsed
                << ",\"event\":\"orbit-heartbeat\""
                << ",\"labeled_configurations\":" << labeled
                << ",\"orbit_count\":" << orbits.size() << "}\n";
            next_heartbeat =
                now + std::chrono::duration_cast<Clock::duration>(
                          std::chrono::duration<double>(
                              arguments.heartbeat_seconds));
          }
        }
      }
    }
    if (labeled != kExpectedLabeledConfigurations) {
      throw std::runtime_error(
          "labeled family count mismatch: got " +
          std::to_string(labeled));
    }
    const auto published_orbit = orbits.find(published_key);
    if (published_orbit == orbits.end()) {
      throw std::runtime_error(
          "published connector pattern missing from orbit family");
    }
    const std::uint64_t orbit_multiplicity_sum =
        std::accumulate(
            orbits.begin(), orbits.end(), std::uint64_t{0},
            [](std::uint64_t sum, const auto& entry) {
              return sum + entry.second.labeled_count;
            });
    if (orbit_multiplicity_sum != labeled) {
      throw std::runtime_error("orbit multiplicities do not sum to family");
    }
    std::map<std::uint64_t, std::uint64_t> orbit_size_histogram;
    for (const auto& [key, orbit] : orbits) {
      (void)key;
      if (kBaseAutomorphismGroupOrder % orbit.labeled_count != 0) {
        throw std::runtime_error(
            "orbit size does not divide base automorphism group order");
      }
      ++orbit_size_histogram[orbit.labeled_count];
    }

    std::vector<std::uint64_t> keys;
    keys.reserve(orbits.size());
    for (const auto& [key, ignored] : orbits) {
      (void)ignored;
      keys.push_back(key);
    }
    std::sort(keys.begin(), keys.end());

    Wide maximum_determinant = 0;
    std::uint64_t maximum_key = 0;
    std::uint64_t positive_determinants = 0;
    std::uint64_t exact_squares = 0;
    std::uint64_t frontier_ties = 0;
    std::uint64_t above_frontier_determinants = 0;
    std::uint64_t above_frontier_squares = 0;
    std::uint64_t divisible_above_frontier_squares = 0;
    std::uint64_t qualified_survivors = 0;
    std::vector<SquareHit> square_hits;
    std::vector<SquareHit> route_hits;

    for (std::size_t index = 0; index < keys.size(); ++index) {
      const std::uint64_t key = keys[index];
      const Configuration configuration = decode_key(key);
      if (canonical_key(configuration, sequence_canonicalizer) != key) {
        throw std::runtime_error(
            "canonical key did not survive decode/re-encode");
      }
      const Matrix matrix = gram(configuration);
      const Wide determinant = exact_determinant(matrix);
      if (determinant > 0) ++positive_determinants;
      if (determinant > maximum_determinant) {
        maximum_determinant = determinant;
        maximum_key = key;
      }
      if (determinant > frontier_squared) {
        ++above_frontier_determinants;
      }
      if (determinant <= 0) continue;
      const Wide root = integer_square_root(determinant);
      if (root * root != determinant) continue;
      ++exact_squares;

      SquareHit hit;
      hit.key = key;
      hit.determinant = determinant;
      hit.square_root = root;
      hit.divisible_by_2_22 =
          root % (static_cast<Wide>(1) << 22) == 0;
      hit.edges = all_edges(configuration);
      hit.added_edges = added_edges(configuration);
      if (root >= kFrontierRoot) {
        hit.positive_definite = exact_positive_definite(matrix);
        square_hits.push_back(hit);
      }
      if (root == kFrontierRoot) ++frontier_ties;
      if (root > kFrontierRoot) {
        ++above_frontier_squares;
        if (hit.divisible_by_2_22) {
          ++divisible_above_frontier_squares;
          if (!hit.positive_definite) {
            hit.positive_definite = exact_positive_definite(matrix);
          }
          if (hit.positive_definite) {
            ++qualified_survivors;
            route_hits.push_back(hit);
          }
        }
      }

      const Clock::time_point now = Clock::now();
      if (now >= next_heartbeat) {
        const double elapsed =
            std::chrono::duration<double>(now - started).count();
        std::cerr
            << "{\"elapsed_seconds\":" << std::fixed
            << std::setprecision(3) << elapsed
            << ",\"event\":\"determinant-heartbeat\""
            << ",\"exact_orbits\":" << index + 1
            << ",\"orbit_count\":" << keys.size()
            << ",\"square_count\":" << exact_squares << "}\n";
        next_heartbeat =
            now + std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(
                          arguments.heartbeat_seconds));
      }
    }

    const double elapsed_seconds =
        std::chrono::duration<double>(Clock::now() - started).count();
    const Configuration maximum_configuration = decode_key(maximum_key);
    const std::vector<std::pair<int, int>> maximum_added =
        added_edges(maximum_configuration);

    std::ostringstream report;
    report << std::setprecision(17);
    report << "{\"challenge_id\":\"maxdet-23-v1\"";
    report << ",\"claim_boundary\":"
           << json_escape(
                  "Exact Gram-orbit screening only; a qualified square "
                  "still requires an exact factor in {-1,+1}^{23x23}.");
    report << ",\"complete\":true";
    report << ",\"elapsed_seconds\":" << elapsed_seconds;
    report << ",\"engine\":\"gram-connector12\"";
    report << ",\"family\":{";
    report << "\"base\":\"K3 disjoint-union 5K4\"";
    report << ",\"connector\":\"K2,2\"";
    report << ",\"connector_count\":3";
    report << ",\"distinct_block_pairs\":true";
    report << ",\"exact_added_edge_count\":12";
    report << ",\"labeled_configuration_count\":" << labeled;
    report << ",\"orbit_size_histogram\":{";
    bool first_histogram_entry = true;
    for (const auto& [orbit_size, count] : orbit_size_histogram) {
      if (!first_histogram_entry) report << ',';
      first_histogram_entry = false;
      report << json_escape(std::to_string(orbit_size))
             << ':' << count;
    }
    report << '}';
    report << ",\"symmetry_group\":"
           << json_escape("S3 x (S4 wr S5)");
    report << ",\"symmetry_group_order\":"
           << kBaseAutomorphismGroupOrder;
    report << ",\"symmetry_reduced_orbit_count\":" << orbits.size();
    report << '}';
    report << ",\"frontier_root\":"
           << json_escape(std::to_string(kFrontierRoot));
    report << ",\"frontier_squared\":"
           << json_escape(decimal(frontier_squared));
    report << ",\"maximum\":{";
    report << "\"added_edges\":";
    append_edge_array(report, maximum_added);
    report << ",\"canonical_key\":"
           << json_escape(std::to_string(maximum_key));
    report << ",\"determinant\":"
           << json_escape(decimal(maximum_determinant));
    const Wide maximum_root =
        integer_square_root(maximum_determinant);
    report << ",\"is_square\":"
           << (maximum_root * maximum_root == maximum_determinant
                   ? "true"
                   : "false")
           << '}';
    report << ",\"normalization\":\"G=24I-J+4A\"";
    report << ",\"published_reproduction\":{";
    report << "\"canonical_key\":"
           << json_escape(std::to_string(published_key));
    report << ",\"determinant\":"
           << json_escape(decimal(published_determinant));
    report << ",\"labeled_orbit_size\":"
           << published_orbit->second.labeled_count;
    report << ",\"square_root\":"
           << json_escape(std::to_string(kFrontierRoot));
    report << ",\"verified\":true}";
    report << ",\"schema_version\":1";
    report << ",\"square_hits_at_or_above_frontier\":[";
    for (std::size_t index = 0; index < square_hits.size(); ++index) {
      if (index != 0) report << ',';
      append_hit(report, square_hits[index]);
    }
    report << ']';
    report << ",\"statistics\":{";
    report << "\"above_frontier_determinants\":"
           << above_frontier_determinants;
    report << ",\"above_frontier_squares\":"
           << above_frontier_squares;
    report << ",\"divisible_above_frontier_squares\":"
           << divisible_above_frontier_squares;
    report << ",\"exact_determinants\":" << keys.size();
    report << ",\"exact_squares\":" << exact_squares;
    report << ",\"frontier_ties\":" << frontier_ties;
    report << ",\"orbit_multiplicity_sum\":"
           << orbit_multiplicity_sum;
    report << ",\"positive_determinants\":"
           << positive_determinants;
    report << ",\"qualified_survivors\":"
           << qualified_survivors << '}';
    report << ",\"termination\":\"completed\"}\n";
    atomic_write(arguments.output, report.str());

    if (!arguments.route_snapshot.empty()) {
      atomic_write(
          arguments.route_snapshot,
          route_snapshot_json(
              route_hits, exact_squares, qualified_survivors,
              keys.size(), elapsed_seconds));
    }

    std::cout
        << "complete labeled=" << labeled
        << " orbits=" << keys.size()
        << " exact_squares=" << exact_squares
        << " frontier_ties=" << frontier_ties
        << " above_frontier_squares=" << above_frontier_squares
        << " qualified_survivors=" << qualified_survivors
        << " elapsed=" << std::fixed << std::setprecision(3)
        << elapsed_seconds << "s"
        << " output=" << arguments.output.string() << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "gram_connector12: " << error.what() << '\n';
    return 2;
  }
}
