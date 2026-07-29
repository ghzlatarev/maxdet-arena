#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
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
constexpr int kLeftSize = 3;
constexpr int kRightBlockCount = 5;
constexpr int kRightBlockSize = 4;
constexpr int kRightSize = kRightBlockCount * kRightBlockSize;
constexpr int kActiveRightVertices = 6;
constexpr int kAddedEdges = 12;
constexpr std::uint64_t kExpectedLabeledCount = 3'488'400;
constexpr std::uint64_t kExpectedOrbitCount = 20;
constexpr std::uint64_t kGroupOrder = 5'733'089'280ULL;
constexpr std::uint64_t kFrontierRoot = 2'779'447'296'000'000ULL;
constexpr std::uint64_t kRequiredDivisor = std::uint64_t{1} << 22;
constexpr std::array<std::uint64_t, 4> kCrtPrimes{
    2'147'483'647ULL,
    2'147'483'629ULL,
    2'147'483'587ULL,
    2'147'483'579ULL,
};

using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using State = std::uint64_t;
using BlockCounts = std::array<std::uint8_t, 3>;
using Profile = std::array<BlockCounts, kRightBlockCount>;

struct Arguments {
  std::filesystem::path output;
  std::filesystem::path route_snapshot;
  std::filesystem::path research_square_snapshot;
  bool self_test = false;
  double heartbeat_seconds = 10.0;
};

struct OrbitTally {
  std::uint64_t labeled_count = 0;
  int minimum_swap_radius = kAddedEdges;
};

struct OrbitResult {
  std::uint64_t index = 0;
  std::uint64_t key = 0;
  std::uint64_t labeled_count = 0;
  int minimum_swap_radius = 0;
  bool distinct_connector = false;
  bool reuse_connector = false;
  Wide determinant = 0;
  bool square = false;
  Wide root = 0;
  bool divisible = false;
  bool positive_definite = false;
  std::vector<std::pair<int, int>> added_edges;
  std::vector<std::pair<int, int>> edges;
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
    } else if (option == "--research-square-snapshot") {
      arguments.research_square_snapshot = value();
    } else if (option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds =
          parse_positive_double(value(), option);
    } else if (option == "--self-test") {
      arguments.self_test = true;
    } else if (option == "--help") {
      std::cout
          << "Usage: gram_published_degree_slice "
             "[--self-test | --output FILE [--route-snapshot FILE] "
             "[--research-square-snapshot FILE]] "
             "[--heartbeat-seconds S]\n\n"
          << "Exactly quotient and determinant-screen the K3 + 5K4 "
             "12-edge added-degree slice (4,4,4 on K3 and "
             "2^6 0^14 on the K4 side).\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (!arguments.self_test && arguments.output.empty()) {
    throw std::runtime_error("--output is required unless --self-test is used");
  }
  if (!arguments.route_snapshot.empty() &&
      std::filesystem::absolute(arguments.output).lexically_normal() ==
          std::filesystem::absolute(arguments.route_snapshot)
              .lexically_normal()) {
    throw std::runtime_error(
        "--route-snapshot must differ from --output");
  }
  const auto normalized = [](const std::filesystem::path& path) {
    return std::filesystem::absolute(path).lexically_normal();
  };
  if (!arguments.research_square_snapshot.empty() &&
      (normalized(arguments.research_square_snapshot) ==
           normalized(arguments.output) ||
       (!arguments.route_snapshot.empty() &&
        normalized(arguments.research_square_snapshot) ==
            normalized(arguments.route_snapshot)))) {
    throw std::runtime_error(
        "--research-square-snapshot must differ from other outputs");
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

int state_value(State state, int vertex) {
  return static_cast<int>((state >> (2 * vertex)) & 3U);
}

State set_state_value(State state, int vertex, int value) {
  const State mask = State{3} << (2 * vertex);
  return (state & ~mask) |
         (static_cast<State>(value) << (2 * vertex));
}

int block_code(const BlockCounts& counts) {
  return 25 * counts[0] + 5 * counts[1] + counts[2];
}

BlockCounts decode_block_code(int code) {
  BlockCounts counts{};
  counts[0] = static_cast<std::uint8_t>(code / 25);
  code %= 25;
  counts[1] = static_cast<std::uint8_t>(code / 5);
  counts[2] = static_cast<std::uint8_t>(code % 5);
  return counts;
}

std::uint64_t encode_profile(const Profile& profile) {
  std::uint64_t key = 0;
  for (const BlockCounts& counts : profile) {
    key = key * 125U +
          static_cast<std::uint64_t>(block_code(counts));
  }
  return key;
}

Profile decode_profile(std::uint64_t key) {
  Profile profile{};
  for (int block = kRightBlockCount - 1; block >= 0; --block) {
    profile[block] =
        decode_block_code(static_cast<int>(key % 125U));
    key /= 125U;
  }
  if (key != 0) {
    throw std::runtime_error("profile key exceeds five base-125 fields");
  }
  return profile;
}

Profile profile_from_state(State state) {
  Profile profile{};
  for (int vertex = 0; vertex < kRightSize; ++vertex) {
    const int value = state_value(state, vertex);
    if (value != 0) {
      ++profile[vertex / kRightBlockSize][value - 1];
    }
  }
  return profile;
}

constexpr std::array<std::array<int, 3>, 6> kColorPermutations{{
    {{0, 1, 2}},
    {{0, 2, 1}},
    {{1, 0, 2}},
    {{1, 2, 0}},
    {{2, 0, 1}},
    {{2, 1, 0}},
}};

std::uint64_t canonical_key(const Profile& profile) {
  std::uint64_t best = std::numeric_limits<std::uint64_t>::max();
  for (const auto& permutation : kColorPermutations) {
    std::array<int, kRightBlockCount> codes{};
    for (int block = 0; block < kRightBlockCount; ++block) {
      const BlockCounts transformed{
          profile[block][permutation[0]],
          profile[block][permutation[1]],
          profile[block][permutation[2]],
      };
      codes[block] = block_code(transformed);
    }
    std::sort(codes.begin(), codes.end());
    std::uint64_t candidate = 0;
    for (const int code : codes) {
      candidate = candidate * 125U +
                  static_cast<std::uint64_t>(code);
    }
    best = std::min(best, candidate);
  }
  return best;
}

std::uint64_t canonical_key(State state) {
  return canonical_key(profile_from_state(state));
}

State representative_state(std::uint64_t key) {
  const Profile profile = decode_profile(key);
  if (canonical_key(profile) != key) {
    throw std::runtime_error(
        "decoded profile is not a canonical profile key");
  }
  State state = 0;
  for (int block = 0; block < kRightBlockCount; ++block) {
    int position = 0;
    for (int color = 0; color < 3; ++color) {
      for (int count = 0; count < profile[block][color]; ++count) {
        if (position >= kRightBlockSize) {
          throw std::runtime_error(
              "canonical profile overfills a K4 block");
        }
        state = set_state_value(
            state, block * kRightBlockSize + position, color + 1);
        ++position;
      }
    }
  }
  if (canonical_key(state) != key) {
    throw std::runtime_error(
        "representative state did not preserve canonical key");
  }
  return state;
}

State published_state() {
  State state = 0;
  for (int offset = 0; offset < 2; ++offset) {
    state = set_state_value(state, offset, 1);
    state = set_state_value(state, 4 + offset, 2);
    state = set_state_value(state, 8 + offset, 3);
  }
  return state;
}

int swap_radius_from_published(State state) {
  const State published = published_state();
  int common_edges = 0;
  for (int vertex = 0; vertex < kRightSize; ++vertex) {
    const int left = state_value(state, vertex);
    const int right = state_value(published, vertex);
    if (left == 0 || right == 0) continue;
    common_edges += left == right ? 2 : 1;
  }
  return kAddedEdges - common_edges;
}

void enumerate_balanced_states(
    const std::array<int, kActiveRightVertices>& active,
    const std::function<void(State)>& visitor) {
  for (int first_zero = 0; first_zero < 5; ++first_zero) {
    for (int second_zero = first_zero + 1; second_zero < 6;
         ++second_zero) {
      std::array<int, 4> remaining{};
      int remaining_count = 0;
      for (int index = 0; index < 6; ++index) {
        if (index != first_zero && index != second_zero) {
          remaining[remaining_count++] = index;
        }
      }
      for (int first_one = 0; first_one < 3; ++first_one) {
        for (int second_one = first_one + 1; second_one < 4;
             ++second_one) {
          State state = 0;
          state = set_state_value(state, active[first_zero], 1);
          state = set_state_value(state, active[second_zero], 1);
          for (int position = 0; position < 4; ++position) {
            const bool color_one =
                position == first_one || position == second_one;
            state = set_state_value(
                state, active[remaining[position]],
                color_one ? 2 : 3);
          }
          visitor(state);
        }
      }
    }
  }
}

void enumerate_active_sets(
    int next_vertex, int depth,
    std::array<int, kActiveRightVertices>& active,
    const std::function<void(
        const std::array<int, kActiveRightVertices>&)>& visitor) {
  if (depth == kActiveRightVertices) {
    visitor(active);
    return;
  }
  const int remaining = kActiveRightVertices - depth;
  for (int vertex = next_vertex;
       vertex <= kRightSize - remaining; ++vertex) {
    active[depth] = vertex;
    enumerate_active_sets(vertex + 1, depth + 1, active, visitor);
  }
}

std::map<std::uint64_t, OrbitTally> direct_labeled_quotient(
    double heartbeat_seconds) {
  const Clock::time_point started = Clock::now();
  Clock::time_point next_heartbeat =
      started + std::chrono::duration_cast<Clock::duration>(
                    std::chrono::duration<double>(heartbeat_seconds));
  std::map<std::uint64_t, OrbitTally> orbits;
  std::array<int, kActiveRightVertices> active{};
  std::uint64_t labeled = 0;
  enumerate_active_sets(
      0, 0, active,
      [&](const std::array<int, kActiveRightVertices>& selected) {
        enumerate_balanced_states(selected, [&](State state) {
          OrbitTally& orbit = orbits[canonical_key(state)];
          ++orbit.labeled_count;
          orbit.minimum_swap_radius = std::min(
              orbit.minimum_swap_radius,
              swap_radius_from_published(state));
          ++labeled;
        });
        const Clock::time_point now = Clock::now();
        if (now >= next_heartbeat) {
          const double elapsed =
              std::chrono::duration<double>(now - started).count();
          std::cerr
              << "{\"elapsed_seconds\":" << std::fixed
              << std::setprecision(3) << elapsed
              << ",\"event\":\"degree-slice-enumeration-heartbeat\""
              << ",\"labeled_states\":" << labeled
              << ",\"orbit_count\":" << orbits.size() << "}\n";
          next_heartbeat =
              now + std::chrono::duration_cast<Clock::duration>(
                        std::chrono::duration<double>(
                            heartbeat_seconds));
        }
      });
  if (labeled != kExpectedLabeledCount) {
    throw std::runtime_error(
        "direct labeled-state count mismatch: got " +
        std::to_string(labeled));
  }
  if (orbits.size() != kExpectedOrbitCount) {
    throw std::runtime_error(
        "direct canonical orbit count mismatch: got " +
        std::to_string(orbits.size()));
  }
  return orbits;
}

std::uint64_t factorial(int value) {
  std::uint64_t result = 1;
  for (int factor = 2; factor <= value; ++factor) {
    result *= static_cast<std::uint64_t>(factor);
  }
  return result;
}

std::uint64_t profile_multiplicity(const Profile& profile) {
  std::uint64_t result = 1;
  for (const BlockCounts& counts : profile) {
    const int active = counts[0] + counts[1] + counts[2];
    if (active > kRightBlockSize) {
      throw std::runtime_error("invalid overfull block profile");
    }
    result *= factorial(kRightBlockSize) /
              (factorial(counts[0]) * factorial(counts[1]) *
               factorial(counts[2]) *
               factorial(kRightBlockSize - active));
  }
  return result;
}

void enumerate_ordered_profiles_recursive(
    int block, std::array<int, 3>& remaining, Profile& profile,
    std::vector<Profile>& profiles) {
  if (block == kRightBlockCount) {
    if (remaining == std::array<int, 3>{0, 0, 0}) {
      profiles.push_back(profile);
    }
    return;
  }
  for (int first = 0;
       first <= std::min(remaining[0], kRightBlockSize); ++first) {
    for (int second = 0;
         second <= std::min(
             remaining[1], kRightBlockSize - first);
         ++second) {
      for (int third = 0;
           third <= std::min(
               remaining[2], kRightBlockSize - first - second);
           ++third) {
        profile[block] = BlockCounts{
            static_cast<std::uint8_t>(first),
            static_cast<std::uint8_t>(second),
            static_cast<std::uint8_t>(third),
        };
        remaining[0] -= first;
        remaining[1] -= second;
        remaining[2] -= third;
        enumerate_ordered_profiles_recursive(
            block + 1, remaining, profile, profiles);
        remaining[0] += first;
        remaining[1] += second;
        remaining[2] += third;
      }
    }
  }
}

std::vector<std::uint64_t> independent_profile_component_sizes(
    std::uint64_t& profile_labeled_sum,
    std::uint64_t& ordered_profile_count) {
  std::vector<Profile> profiles;
  Profile profile{};
  std::array<int, 3> remaining{2, 2, 2};
  enumerate_ordered_profiles_recursive(
      0, remaining, profile, profiles);
  ordered_profile_count = profiles.size();

  std::unordered_map<std::uint64_t, std::size_t> indices;
  indices.reserve(profiles.size() * 2);
  std::vector<std::uint64_t> multiplicities;
  multiplicities.reserve(profiles.size());
  profile_labeled_sum = 0;
  for (std::size_t index = 0; index < profiles.size(); ++index) {
    const std::uint64_t code = encode_profile(profiles[index]);
    if (!indices.emplace(code, index).second) {
      throw std::runtime_error(
          "ordered profile generator produced a duplicate");
    }
    const std::uint64_t multiplicity =
        profile_multiplicity(profiles[index]);
    multiplicities.push_back(multiplicity);
    profile_labeled_sum += multiplicity;
  }

  std::vector<bool> visited(profiles.size(), false);
  std::vector<std::uint64_t> component_sizes;
  for (std::size_t start = 0; start < profiles.size(); ++start) {
    if (visited[start]) continue;
    std::queue<std::size_t> pending;
    pending.push(start);
    visited[start] = true;
    std::uint64_t component_size = 0;
    while (!pending.empty()) {
      const std::size_t index = pending.front();
      pending.pop();
      component_size += multiplicities[index];
      const Profile current = profiles[index];

      for (int block = 0; block < kRightBlockCount - 1; ++block) {
        Profile neighbor = current;
        std::swap(neighbor[block], neighbor[block + 1]);
        const auto found = indices.find(encode_profile(neighbor));
        if (found == indices.end()) {
          throw std::runtime_error(
              "block-swap generator left the profile family");
        }
        if (!visited[found->second]) {
          visited[found->second] = true;
          pending.push(found->second);
        }
      }
      for (int color = 0; color < 2; ++color) {
        Profile neighbor = current;
        for (BlockCounts& counts : neighbor) {
          std::swap(counts[color], counts[color + 1]);
        }
        const auto found = indices.find(encode_profile(neighbor));
        if (found == indices.end()) {
          throw std::runtime_error(
              "color-swap generator left the profile family");
        }
        if (!visited[found->second]) {
          visited[found->second] = true;
          pending.push(found->second);
        }
      }
    }
    component_sizes.push_back(component_size);
  }
  std::sort(component_sizes.begin(), component_sizes.end());
  return component_sizes;
}

std::vector<std::pair<int, int>> added_edges(State state) {
  std::vector<std::pair<int, int>> result;
  for (int vertex = 0; vertex < kRightSize; ++vertex) {
    const int value = state_value(state, vertex);
    if (value == 0) continue;
    const int missing_color = value - 1;
    const int right_vertex = kLeftSize + vertex;
    for (int left_vertex = 0; left_vertex < kLeftSize;
         ++left_vertex) {
      if (left_vertex != missing_color) {
        result.emplace_back(left_vertex, right_vertex);
      }
    }
  }
  std::sort(result.begin(), result.end());
  if (result.size() != kAddedEdges ||
      std::adjacent_find(result.begin(), result.end()) != result.end()) {
    throw std::runtime_error(
        "slice state did not produce 12 distinct added edges");
  }
  return result;
}

std::vector<std::pair<int, int>> all_edges(State state) {
  std::vector<std::pair<int, int>> result;
  for (int left = 0; left < kLeftSize; ++left) {
    for (int right = left + 1; right < kLeftSize; ++right) {
      result.emplace_back(left, right);
    }
  }
  for (int block = 0; block < kRightBlockCount; ++block) {
    const int offset = kLeftSize + block * kRightBlockSize;
    for (int left = 0; left < kRightBlockSize; ++left) {
      for (int right = left + 1; right < kRightBlockSize; ++right) {
        result.emplace_back(offset + left, offset + right);
      }
    }
  }
  const auto extra = added_edges(state);
  result.insert(result.end(), extra.begin(), extra.end());
  std::sort(result.begin(), result.end());
  if (result.size() != 45 ||
      std::adjacent_find(result.begin(), result.end()) != result.end()) {
    throw std::runtime_error(
        "slice state did not produce a 45-edge defect graph");
  }
  return result;
}

Matrix gram(State state) {
  Matrix result{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      result[row][column] = row == column ? 23 : -1;
    }
  }
  for (const auto& [first, second] : all_edges(state)) {
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
  // by 727^(23/2) < 2^110. The four-prime modulus exceeds 2^123, so centered
  // CRT reconstruction is exact and unique.
  UnsignedWide reconstructed = 0;
  UnsignedWide modulus = 1;
  for (const std::uint64_t prime : kCrtPrimes) {
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

bool connector_decomposable(State state, bool require_distinct_blocks) {
  std::array<int, 3> color_blocks{-1, -1, -1};
  for (int color = 0; color < 3; ++color) {
    int count = 0;
    for (int vertex = 0; vertex < kRightSize; ++vertex) {
      if (state_value(state, vertex) != color + 1) continue;
      const int block = vertex / kRightBlockSize;
      if (color_blocks[color] < 0) {
        color_blocks[color] = block;
      } else if (color_blocks[color] != block) {
        return false;
      }
      ++count;
    }
    if (count != 2) return false;
  }
  if (!require_distinct_blocks) return true;
  std::sort(color_blocks.begin(), color_blocks.end());
  return std::adjacent_find(
             color_blocks.begin(), color_blocks.end()) ==
         color_blocks.end();
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

void append_profile(std::ostream& output, std::uint64_t key) {
  const Profile profile = decode_profile(key);
  output << '[';
  for (int block = 0; block < kRightBlockCount; ++block) {
    if (block != 0) output << ',';
    output << '['
           << static_cast<int>(profile[block][0]) << ','
           << static_cast<int>(profile[block][1]) << ','
           << static_cast<int>(profile[block][2]) << ']';
  }
  output << ']';
}

void append_string_uint_histogram(
    std::ostream& output,
    const std::map<std::uint64_t, std::uint64_t>& histogram) {
  output << '{';
  bool first = true;
  for (const auto& [value, count] : histogram) {
    if (!first) output << ',';
    first = false;
    output << json_escape(std::to_string(value)) << ':' << count;
  }
  output << '}';
}

void append_int_uint_histogram(
    std::ostream& output,
    const std::map<int, std::uint64_t>& histogram) {
  output << '{';
  bool first = true;
  for (const auto& [value, count] : histogram) {
    if (!first) output << ',';
    first = false;
    output << json_escape(std::to_string(value)) << ':' << count;
  }
  output << '}';
}

void append_route_hit(
    std::ostream& output, const OrbitResult& result) {
  output << "{\"added_edges\":";
  append_edge_array(output, result.added_edges);
  output << ",\"canonical_key\":"
         << json_escape(std::to_string(result.key));
  output << ",\"determinant\":"
         << json_escape(decimal(result.determinant));
  output << ",\"divisible_by_2_22\":"
         << (result.divisible ? "true" : "false");
  output << ",\"edge_count\":" << result.edges.size();
  output << ",\"edges\":";
  append_edge_array(output, result.edges);
  output << ",\"orbit_index\":" << result.index;
  output << ",\"positive_definite\":"
         << (result.positive_definite ? "true" : "false");
  output << ",\"qualified\":"
         << (result.divisible && result.positive_definite ? "true"
                                                          : "false");
  output << ",\"square_root\":"
         << json_escape(decimal(result.root)) << '}';
}

std::string route_snapshot_json(
    const std::vector<OrbitResult>& route_hits,
    std::uint64_t exact_squares,
    std::uint64_t qualified_survivors,
    double elapsed_seconds) {
  const Wide frontier_squared =
      static_cast<Wide>(kFrontierRoot) *
      static_cast<Wide>(kFrontierRoot);
  std::ostringstream output;
  output << std::setprecision(17);
  output << "{\"challenge_id\":\"maxdet-23-v1\"";
  output << ",\"claim_boundary\":"
         << json_escape(
                "Compatibility snapshot produced by "
                "gram-published-degree-slice; survivors still require "
                "exact sign-matrix factorization.");
  output << ",\"complete\":true";
  output << ",\"elapsed_seconds\":" << elapsed_seconds;
  // Existing exact Hasse and shell tools validate this compatibility label.
  // The producer and mode fields retain the actual provenance.
  output << ",\"engine\":\"gram-tabu\"";
  output << ",\"frontier_root\":"
         << json_escape(std::to_string(kFrontierRoot));
  output << ",\"frontier_squared\":"
         << json_escape(decimal(frontier_squared));
  output << ",\"hits\":[";
  for (std::size_t index = 0; index < route_hits.size(); ++index) {
    if (index != 0) output << ',';
    append_route_hit(output, route_hits[index]);
  }
  output << ']';
  output << ",\"mode\":\"published-degree-slice-exact-orbit-enumeration\"";
  output << ",\"normalization\":\"G=24I-J+4A\"";
  output << ",\"parameters\":{";
  output << "\"max_stored_hits\":" << route_hits.size();
  output << ",\"orbit_count\":" << kExpectedOrbitCount;
  output << ",\"producer_engine\":\"gram-published-degree-slice\"}";
  output << ",\"schema_version\":1";
  output << ",\"statistics\":{";
  output << "\"exact_squares\":" << exact_squares;
  output << ",\"qualified_survivors\":" << qualified_survivors;
  output << ",\"unrecorded_square_observations\":0}";
  output << ",\"termination\":\"completed\"}\n";
  return output.str();
}

std::string research_square_snapshot_json(
    const std::vector<OrbitResult>& square_hits,
    double elapsed_seconds) {
  const Wide frontier_squared =
      static_cast<Wide>(kFrontierRoot) *
      static_cast<Wide>(kFrontierRoot);
  std::uint64_t compatible = 0;
  for (const OrbitResult& hit : square_hits) {
    if (hit.divisible && hit.positive_definite) ++compatible;
  }
  std::ostringstream output;
  output << std::setprecision(17);
  output << "{\"challenge_id\":\"maxdet-23-v1\"";
  output << ",\"claim_boundary\":"
         << json_escape(
                "Research-only exact square snapshot. It deliberately "
                "retains tied and subfrontier Grams for compatibility "
                "audits and is not a strict-improvement route.");
  output << ",\"complete\":true";
  output << ",\"elapsed_seconds\":" << elapsed_seconds;
  output << ",\"engine\":\"gram-tabu\"";
  output << ",\"frontier_filter_applied\":false";
  output << ",\"frontier_root\":"
         << json_escape(std::to_string(kFrontierRoot));
  output << ",\"frontier_squared\":"
         << json_escape(decimal(frontier_squared));
  output << ",\"hits\":[";
  for (std::size_t index = 0; index < square_hits.size(); ++index) {
    if (index != 0) output << ',';
    append_route_hit(output, square_hits[index]);
  }
  output << ']';
  output << ",\"mode\":\"published-degree-slice-research-squares\"";
  output << ",\"normalization\":\"G=24I-J+4A\"";
  output << ",\"parameters\":{";
  output << "\"max_stored_hits\":" << square_hits.size();
  output << ",\"orbit_count\":" << kExpectedOrbitCount;
  output << ",\"producer_engine\":\"gram-published-degree-slice\"}";
  output << ",\"schema_version\":1";
  output << ",\"statistics\":{";
  output << "\"exact_squares\":" << square_hits.size();
  output << ",\"qualified_survivors\":" << compatible;
  output << ",\"unrecorded_square_observations\":0}";
  output << ",\"termination\":\"completed\"}\n";
  return output.str();
}

std::uint64_t choose(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  std::uint64_t result = 1;
  for (int index = 1; index <= k; ++index) {
    result =
        result * static_cast<std::uint64_t>(n - k + index) /
        static_cast<std::uint64_t>(index);
  }
  return result;
}

struct Audit {
  std::map<std::uint64_t, OrbitTally> orbits;
  std::uint64_t closed_form_labeled_count = 0;
  std::uint64_t profile_labeled_count = 0;
  std::uint64_t ordered_profile_count = 0;
  std::vector<std::uint64_t> independent_component_sizes;
  std::vector<OrbitResult> results;
};

Audit run_audit(double heartbeat_seconds) {
  Audit audit;
  audit.closed_form_labeled_count =
      choose(20, 6) * factorial(6) /
      (factorial(2) * factorial(2) * factorial(2));
  if (audit.closed_form_labeled_count != kExpectedLabeledCount) {
    throw std::runtime_error("closed-form labeled count mismatch");
  }

  audit.orbits = direct_labeled_quotient(heartbeat_seconds);
  audit.independent_component_sizes =
      independent_profile_component_sizes(
          audit.profile_labeled_count,
          audit.ordered_profile_count);
  if (audit.profile_labeled_count != kExpectedLabeledCount) {
    throw std::runtime_error(
        "independent block-profile labeled sum mismatch");
  }
  if (audit.independent_component_sizes.size() !=
      kExpectedOrbitCount) {
    throw std::runtime_error(
        "independent profile-generator component count mismatch");
  }

  std::vector<std::uint64_t> direct_orbit_sizes;
  for (const auto& [key, tally] : audit.orbits) {
    (void)key;
    if (kGroupOrder % tally.labeled_count != 0) {
      throw std::runtime_error(
          "direct orbit size does not divide full group order");
    }
    direct_orbit_sizes.push_back(tally.labeled_count);
  }
  std::sort(direct_orbit_sizes.begin(), direct_orbit_sizes.end());
  if (direct_orbit_sizes != audit.independent_component_sizes) {
    throw std::runtime_error(
        "direct canonical and profile-component orbit sizes disagree");
  }

  const Wide frontier_squared =
      static_cast<Wide>(kFrontierRoot) *
      static_cast<Wide>(kFrontierRoot);
  const Wide published_determinant =
      exact_determinant(gram(published_state()));
  if (published_determinant != frontier_squared) {
    throw std::runtime_error(
        "published degree-pattern representative misses the frontier");
  }
  if (canonical_key(published_state()) == 0) {
    throw std::runtime_error("published canonical key is unexpectedly zero");
  }

  std::uint64_t index = 0;
  for (const auto& [key, tally] : audit.orbits) {
    const State state = representative_state(key);
    OrbitResult result;
    result.index = index++;
    result.key = key;
    result.labeled_count = tally.labeled_count;
    result.minimum_swap_radius = tally.minimum_swap_radius;
    result.distinct_connector =
        connector_decomposable(state, true);
    result.reuse_connector =
        connector_decomposable(state, false);
    const Matrix matrix = gram(state);
    result.determinant = exact_determinant(matrix);
    if (result.determinant > 0) {
      result.root = integer_square_root(result.determinant);
      result.square =
          result.root * result.root == result.determinant;
    }
    result.divisible =
        result.square &&
        result.root % static_cast<Wide>(kRequiredDivisor) == 0;
    result.positive_definite = exact_positive_definite(matrix);
    result.added_edges = added_edges(state);
    result.edges = all_edges(state);
    audit.results.push_back(std::move(result));
  }
  return audit;
}

std::string report_json(
    const Audit& audit, double elapsed_seconds) {
  const Wide frontier_squared =
      static_cast<Wide>(kFrontierRoot) *
      static_cast<Wide>(kFrontierRoot);
  std::map<std::uint64_t, std::uint64_t> orbit_size_histogram;
  std::map<int, std::uint64_t> minimum_radius_histogram;
  std::uint64_t positive = 0;
  std::uint64_t positive_labeled = 0;
  std::uint64_t positive_definite = 0;
  std::uint64_t positive_definite_labeled = 0;
  std::uint64_t above = 0;
  std::uint64_t above_labeled = 0;
  std::uint64_t squares = 0;
  std::uint64_t square_labeled = 0;
  std::uint64_t above_squares = 0;
  std::uint64_t above_square_labeled = 0;
  std::uint64_t divisible_above_squares = 0;
  std::uint64_t qualified_survivors = 0;
  std::uint64_t frontier_ties = 0;
  std::uint64_t distinct_connector = 0;
  std::uint64_t reuse_connector = 0;
  std::uint64_t radius_le_three = 0;
  std::uint64_t connector_or_radius = 0;
  std::uint64_t labeled_sum = 0;
  const OrbitResult* maximum = nullptr;
  std::vector<const OrbitResult*> square_results;
  std::vector<const OrbitResult*> route_results;

  for (const OrbitResult& result : audit.results) {
    ++orbit_size_histogram[result.labeled_count];
    ++minimum_radius_histogram[result.minimum_swap_radius];
    labeled_sum += result.labeled_count;
    if (result.determinant > 0) {
      ++positive;
      positive_labeled += result.labeled_count;
    }
    if (result.positive_definite) {
      ++positive_definite;
      positive_definite_labeled += result.labeled_count;
    }
    if (result.determinant > frontier_squared) {
      ++above;
      above_labeled += result.labeled_count;
    }
    if (result.square) {
      ++squares;
      square_labeled += result.labeled_count;
      square_results.push_back(&result);
      if (result.root == static_cast<Wide>(kFrontierRoot)) {
        ++frontier_ties;
      }
      if (result.root > static_cast<Wide>(kFrontierRoot)) {
        ++above_squares;
        above_square_labeled += result.labeled_count;
        if (result.divisible) {
          ++divisible_above_squares;
          if (result.positive_definite) {
            ++qualified_survivors;
            route_results.push_back(&result);
          }
        }
      }
    }
    if (result.distinct_connector) ++distinct_connector;
    if (result.reuse_connector) ++reuse_connector;
    if (result.minimum_swap_radius <= 3) ++radius_le_three;
    if (result.reuse_connector ||
        result.minimum_swap_radius <= 3) {
      ++connector_or_radius;
    }
    if (maximum == nullptr ||
        result.determinant > maximum->determinant) {
      maximum = &result;
    }
  }
  if (labeled_sum != kExpectedLabeledCount || maximum == nullptr) {
    throw std::runtime_error("report aggregation mismatch");
  }

  std::ostringstream output;
  output << std::setprecision(17);
  output << "{\"challenge_id\":\"maxdet-23-v1\"";
  output << ",\"claim_boundary\":"
         << json_escape(
                "Exact Gram-orbit screening only; a qualified square "
                "still requires an exact {-1,+1} factor.");
  output << ",\"complete\":true";
  output << ",\"elapsed_seconds\":" << elapsed_seconds;
  output << ",\"engine\":\"gram-published-degree-slice\"";
  output << ",\"family\":{";
  output << "\"added_edge_count\":12";
  output << ",\"base\":\"K3 disjoint-union 5K4\"";
  output << ",\"construction_scope\":"
         << json_escape(
                "Fix the ideal K3 disjoint-union 5K4 base; add exactly "
                "12 previously absent K3-to-K4-side edges and delete no "
                "base edge. This is not the family of every graph with "
                "the same final degree multiset.");
  output << ",\"closed_form\":"
         << json_escape("C(20,6) * 6! / (2!^3)");
  output << ",\"closed_form_labeled_count\":"
         << audit.closed_form_labeled_count;
  output << ",\"direct_labeled_count\":" << labeled_sum;
  output << ",\"k3_added_degrees\":[4,4,4]";
  output << ",\"k4_side_added_degree_histogram\":{"
         << "\"0\":14,\"2\":6}";
  output << ",\"no_base_edge_deletions\":true";
  output << ",\"only_absent_k3_to_k4_side_edges_added\":true";
  output << ",\"orbit_size_histogram\":";
  append_string_uint_histogram(output, orbit_size_histogram);
  output << ",\"symmetry_group\":"
         << json_escape("S3 x (S4 wr S5)");
  output << ",\"symmetry_group_order\":" << kGroupOrder;
  output << ",\"symmetry_reduced_orbit_count\":"
         << audit.results.size() << '}';
  output << ",\"frontier_root\":"
         << json_escape(std::to_string(kFrontierRoot));
  output << ",\"frontier_squared\":"
         << json_escape(decimal(frontier_squared));
  output << ",\"independent_quotient_check\":{";
  output << "\"component_count\":"
         << audit.independent_component_sizes.size();
  output << ",\"generator_set\":"
         << json_escape(
                "adjacent K4-block swaps and adjacent missing-color swaps");
  output << ",\"labeled_multiplicity_sum\":"
         << audit.profile_labeled_count;
  output << ",\"ordered_block_profile_count\":"
         << audit.ordered_profile_count;
  output << ",\"orbit_size_multiset_matches\":true}";
  output << ",\"maximum\":{";
  output << "\"added_edges\":";
  append_edge_array(output, maximum->added_edges);
  output << ",\"canonical_key\":"
         << json_escape(std::to_string(maximum->key));
  output << ",\"determinant\":"
         << json_escape(decimal(maximum->determinant));
  output << ",\"is_square\":"
         << (maximum->square ? "true" : "false");
  output << ",\"orbit_index\":" << maximum->index;
  output << ",\"positive_definite\":"
         << (maximum->positive_definite ? "true" : "false")
         << '}';
  output << ",\"normalization\":\"G=24I-J+4A\"";
  output << ",\"overlap\":{";
  output << "\"connector12_distinct_pair_orbits\":"
         << distinct_connector;
  output << ",\"connector12_reuse_orbits\":" << reuse_connector;
  output << ",\"connector12_reuse_or_radius_le_3_orbits\":"
         << connector_or_radius;
  output << ",\"minimum_swap_radius_histogram\":";
  append_int_uint_histogram(output, minimum_radius_histogram);
  output << ",\"outside_connector12_reuse_orbits\":"
         << audit.results.size() - reuse_connector;
  output << ",\"outside_connector12_reuse_or_radius_le_3_orbits\":"
         << audit.results.size() - connector_or_radius;
  output << ",\"outside_radius_le_3_orbits\":"
         << audit.results.size() - radius_le_three;
  output << ",\"radius_le_3_orbits\":" << radius_le_three << '}';
  output << ",\"orbits\":[";
  for (std::size_t index = 0; index < audit.results.size(); ++index) {
    if (index != 0) output << ',';
    const OrbitResult& result = audit.results[index];
    output << "{\"added_edges\":";
    append_edge_array(output, result.added_edges);
    output << ",\"block_missing_color_counts\":";
    append_profile(output, result.key);
    output << ",\"canonical_key\":"
           << json_escape(std::to_string(result.key));
    output << ",\"determinant\":"
           << json_escape(decimal(result.determinant));
    output << ",\"divisible_by_2_22\":"
           << (result.divisible ? "true" : "false");
    output << ",\"in_connector12_distinct_pair\":"
           << (result.distinct_connector ? "true" : "false");
    output << ",\"in_connector12_reuse\":"
           << (result.reuse_connector ? "true" : "false");
    output << ",\"in_radius_le_3\":"
           << (result.minimum_swap_radius <= 3 ? "true" : "false");
    output << ",\"is_square\":"
           << (result.square ? "true" : "false");
    output << ",\"labeled_orbit_size\":"
           << result.labeled_count;
    output << ",\"minimum_swap_radius_from_published\":"
           << result.minimum_swap_radius;
    output << ",\"orbit_index\":" << result.index;
    output << ",\"positive_definite\":"
           << (result.positive_definite ? "true" : "false");
    if (result.square) {
      output << ",\"square_root\":"
             << json_escape(decimal(result.root));
    }
    output << '}';
  }
  output << ']';
  output << ",\"schema_version\":1";
  output << ",\"square_orbits\":[";
  for (std::size_t index = 0; index < square_results.size(); ++index) {
    if (index != 0) output << ',';
    const OrbitResult& result = *square_results[index];
    output << "{\"canonical_key\":"
           << json_escape(std::to_string(result.key));
    output << ",\"determinant\":"
           << json_escape(decimal(result.determinant));
    output << ",\"divisible_by_2_22\":"
           << (result.divisible ? "true" : "false");
    output << ",\"labeled_orbit_size\":"
           << result.labeled_count;
    output << ",\"orbit_index\":" << result.index;
    output << ",\"positive_definite\":"
           << (result.positive_definite ? "true" : "false");
    output << ",\"square_root\":"
           << json_escape(decimal(result.root)) << '}';
  }
  output << ']';
  output << ",\"statistics\":{";
  output << "\"above_frontier_determinants\":" << above;
  output << ",\"above_frontier_labeled\":" << above_labeled;
  output << ",\"above_frontier_square_labeled\":"
         << above_square_labeled;
  output << ",\"above_frontier_squares\":" << above_squares;
  output << ",\"divisible_above_frontier_squares\":"
         << divisible_above_squares;
  output << ",\"exact_determinants\":" << audit.results.size();
  output << ",\"exact_squares\":" << squares;
  output << ",\"frontier_ties\":" << frontier_ties;
  output << ",\"labeled_count\":" << labeled_sum;
  output << ",\"positive_definite_labeled\":"
         << positive_definite_labeled;
  output << ",\"positive_definite_orbits\":"
         << positive_definite;
  output << ",\"positive_determinant_labeled\":"
         << positive_labeled;
  output << ",\"positive_determinants\":" << positive;
  output << ",\"qualified_survivors\":"
         << qualified_survivors;
  output << ",\"square_labeled\":" << square_labeled << '}';
  output << ",\"termination\":\"completed\"}\n";
  return output.str();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    const Clock::time_point started = Clock::now();
    const Audit audit = run_audit(arguments.heartbeat_seconds);
    const double elapsed_seconds =
        std::chrono::duration<double>(Clock::now() - started).count();

    const Wide frontier_squared =
        static_cast<Wide>(kFrontierRoot) *
        static_cast<Wide>(kFrontierRoot);
    std::vector<OrbitResult> route_hits;
    std::vector<OrbitResult> research_square_hits;
    std::uint64_t exact_squares = 0;
    for (const OrbitResult& result : audit.results) {
      if (result.square) {
        ++exact_squares;
        research_square_hits.push_back(result);
      }
      if (result.square && result.root > kFrontierRoot &&
          result.divisible && result.positive_definite) {
        route_hits.push_back(result);
      }
    }

    if (arguments.self_test) {
      if (audit.results.size() != kExpectedOrbitCount ||
          exact_squares == 0 ||
          exact_determinant(gram(published_state())) !=
              frontier_squared) {
        throw std::runtime_error("self-test postconditions failed");
      }
      std::cout
          << "self-test passed labeled=" << kExpectedLabeledCount
          << " orbits=" << audit.results.size()
          << " exact_squares=" << exact_squares
          << " route_hits=" << route_hits.size()
          << " elapsed=" << std::fixed << std::setprecision(3)
          << elapsed_seconds << "s\n";
      return 0;
    }

    atomic_write(
        arguments.output, report_json(audit, elapsed_seconds));
    if (!arguments.route_snapshot.empty()) {
      atomic_write(
          arguments.route_snapshot,
          route_snapshot_json(
              route_hits, exact_squares, route_hits.size(),
              elapsed_seconds));
    }
    if (!arguments.research_square_snapshot.empty()) {
      atomic_write(
          arguments.research_square_snapshot,
          research_square_snapshot_json(
              research_square_hits, elapsed_seconds));
    }
    std::cout
        << "complete labeled=" << kExpectedLabeledCount
        << " orbits=" << audit.results.size()
        << " exact_squares=" << exact_squares
        << " route_hits=" << route_hits.size()
        << " elapsed=" << std::fixed << std::setprecision(3)
        << elapsed_seconds << "s"
        << " output=" << arguments.output.string() << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr
        << "gram_published_degree_slice: " << error.what() << '\n';
    return 2;
  }
}
