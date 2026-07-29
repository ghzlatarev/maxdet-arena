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
using Clock = std::chrono::steady_clock;
using Exact = __int128_t;
using UnsignedExact = __uint128_t;
using SignMatrix = std::array<std::array<int, kOrder>, kOrder>;
using Gram = std::array<std::array<int, kOrder>, kOrder>;
using FloatingMatrix =
    std::array<std::array<long double, kOrder>, kOrder>;
using Graph = std::array<unsigned char, kEdgeCount>;
using GraphKey = std::array<std::uint64_t, 4>;

constexpr std::uint64_t kFrontierRoot = 2'779'447'296'000'000ULL;
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

std::string exact_string(const Exact& value) {
  if (value == 0) return "0";
  Exact remaining = value;
  const bool negative = remaining < 0;
  if (negative) remaining = -remaining;
  std::string result;
  while (remaining != 0) {
    result.push_back(
        static_cast<char>('0' + static_cast<int>(remaining % 10)));
    remaining /= 10;
  }
  if (negative) result.push_back('-');
  std::reverse(result.begin(), result.end());
  return result;
}

constexpr std::array<std::uint64_t, 4> kCrtPrimes = {
    1'000'000'007ULL,
    1'000'000'009ULL,
    1'000'000'033ULL,
    1'000'000'087ULL,
};

std::uint64_t modular_power(
    std::uint64_t base,
    std::uint64_t exponent,
    std::uint64_t modulus) {
  std::uint64_t result = 1;
  while (exponent != 0) {
    if (exponent & 1U) {
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
    while (pivot_row < order && work[pivot_row][column] == 0) ++pivot_row;
    if (pivot_row == order) return 0;
    if (pivot_row != column) {
      std::swap(work[pivot_row], work[column]);
      determinant =
          determinant == 0 ? 0 : prime - determinant;
    }

    const std::uint64_t pivot = work[column][column];
    determinant = (determinant * pivot) % prime;
    const std::uint64_t inverse =
        modular_power(pivot, prime - 2, prime);
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
        prime - 2,
        prime);
    const std::uint64_t multiplier =
        (difference * inverse) % prime;
    reconstructed += modulus * multiplier;
    modulus *= prime;
  }

  // Every normalized Gram row has squared norm at most
  // 23^2 + 22*3^2 = 727.  Thus Hadamard bounds |det(G)| below
  // 727^(23/2) < 10^33, while this CRT modulus exceeds 10^36.
  // Symmetric reconstruction is consequently unique.
  if (reconstructed > modulus / 2) {
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

Exact integer_square_root(const Exact& value) {
  if (value < 0) {
    throw std::runtime_error("square root of a negative integer");
  }
  if (value == 0) return 0;
  UnsignedExact remaining = static_cast<UnsignedExact>(value);
  unsigned bits = 0;
  for (UnsignedExact copy = remaining; copy != 0; copy >>= 1U) {
    ++bits;
  }
  UnsignedExact estimate =
      UnsignedExact{1} << ((bits + 1) / 2);
  for (;;) {
    const UnsignedExact next =
        (estimate + remaining / estimate) >> 1U;
    if (next >= estimate) return static_cast<Exact>(estimate);
    estimate = next;
  }
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
        gram[row][other] += matrix[row][column] * matrix[other][column];
      }
    }
  }
  return gram;
}

Graph graph_from_gram(const Gram& gram) {
  Graph graph{};
  for (int index = 0; index < kEdgeCount; ++index) {
    const Edge edge = kEdges[index];
    const int value = gram[edge.first][edge.second];
    if (value != -1 && value != 3) {
      throw std::runtime_error(
          "start Gram is not normalized as 24I-J+4A");
    }
    graph[index] = static_cast<unsigned char>(value == 3);
  }
  for (int diagonal = 0; diagonal < kOrder; ++diagonal) {
    if (gram[diagonal][diagonal] != kOrder) {
      throw std::runtime_error("start Gram diagonal is not 23");
    }
  }
  return graph;
}

Gram gram_from_graph(const Graph& graph) {
  Gram gram{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      gram[row][column] =
          row == column ? kOrder : -1;
    }
  }
  for (int index = 0; index < kEdgeCount; ++index) {
    if (!graph[index]) continue;
    const Edge edge = kEdges[index];
    gram[edge.first][edge.second] = 3;
    gram[edge.second][edge.first] = 3;
  }
  return gram;
}

int edge_count(const Graph& graph) {
  return std::accumulate(graph.begin(), graph.end(), 0);
}

GraphKey graph_key(const Graph& graph) {
  GraphKey key{};
  for (int index = 0; index < kEdgeCount; ++index) {
    if (graph[index]) {
      key[index / 64] |= std::uint64_t{1} << (index % 64);
    }
  }
  return key;
}

struct State {
  Graph graph{};
  Gram gram{};
  FloatingMatrix inverse{};
  long double log_determinant =
      -std::numeric_limits<long double>::infinity();
  Exact exact_determinant = 0;
};

bool rebuild_numeric(State& state) {
  state.gram = gram_from_graph(state.graph);
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
  for (int right_hand_side = 0;
       right_hand_side < kOrder;
       ++right_hand_side) {
    std::array<long double, kOrder> forward{};
    for (int row = 0; row < kOrder; ++row) {
      long double value = static_cast<long double>(
          row == right_hand_side);
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
      inverse[row][right_hand_side] = solution[row];
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
    const State& state, int edge_index) {
  const Edge edge = kEdges[edge_index];
  const long double delta = state.graph[edge_index] ? -4.0L : 4.0L;
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
    } else if (option == "--cooling-iterations") {
      arguments.cooling_iterations = parse_unsigned(value(), option);
    } else if (option == "--max-stored-hits") {
      arguments.max_stored_hits = parse_unsigned(value(), option);
    } else if (option == "--temperature-start") {
      arguments.temperature_start =
          parse_double(value(), option);
    } else if (option == "--temperature-end") {
      arguments.temperature_end =
          parse_double(value(), option);
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
  if (arguments.cooling_iterations == 0) {
    throw std::runtime_error("--cooling-iterations must be positive");
  }
  if (!(arguments.temperature_start > 0) ||
      !(arguments.temperature_end > 0) ||
      arguments.temperature_end > arguments.temperature_start) {
    throw std::runtime_error(
        "temperatures must satisfy start >= end > 0");
  }

  const auto start =
      std::filesystem::absolute(arguments.start).lexically_normal();
  const auto output =
      std::filesystem::absolute(arguments.output).lexically_normal();
  if (start == output) {
    throw std::runtime_error("--start and --output must be distinct");
  }
  return arguments;
}

void print_usage(std::ostream& output) {
  output
      << "usage: gram_tabu --start MATRIX --output SNAPSHOT [options]\n\n"
      << "Search normalized order-23 Gram graphs and exact-screen "
         "square determinants.\n\n"
      << "options:\n"
      << "  --mode MODE                 hill, tabu, or anneal "
         "(default tabu)\n"
      << "  --seed N                    PRNG seed (default 23)\n"
      << "  --seconds S                 positive wall limit "
         "(default 60)\n"
      << "  --max-iterations N          0 means wall limit only "
         "(default 0)\n"
      << "  --restart-iterations N      iterations without a best "
         "before a kick\n"
      << "  --kick-size N               toggled graph edges per kick\n"
      << "  --tabu-tenure N             positive tabu tenure\n"
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
  std::uint64_t accepted_moves = 0;
  std::uint64_t exact_neighbor_screens = 0;
  std::uint64_t above_frontier = 0;
  std::uint64_t exact_squares = 0;
  std::uint64_t exact_pd_checks = 0;
  std::uint64_t qualified_survivors = 0;
  std::uint64_t unrecorded_square_observations = 0;
  std::uint64_t non_pd_proposals = 0;
  std::uint64_t restarts = 0;
};

struct SquareHit {
  Graph graph{};
  Exact determinant = 0;
  Exact root = 0;
  bool divisible = false;
  bool positive_definite = false;
};

bool screen_exact_candidate(
    const Graph& graph,
    const Gram& gram,
    const Exact& determinant,
    Statistics& statistics,
    std::vector<SquareHit>& hits,
    std::set<GraphKey>& recorded_hits,
    std::uint64_t max_stored_hits) {
  ++statistics.exact_neighbor_screens;
  if (determinant <= frontier_determinant()) return false;
  ++statistics.above_frontier;

  const Exact root = integer_square_root(determinant);
  if (root * root != determinant) return false;
  ++statistics.exact_squares;
  const Exact divisor = Exact(1) << kRequiredPower;
  const bool divisible = root % divisor == 0;
  ++statistics.exact_pd_checks;
  const bool positive_definite = exact_positive_definite(gram);
  const bool qualified = divisible && positive_definite;
  if (qualified) ++statistics.qualified_survivors;

  const GraphKey key = graph_key(graph);
  if (recorded_hits.find(key) != recorded_hits.end()) return false;
  if (hits.size() >= max_stored_hits) {
    ++statistics.unrecorded_square_observations;
    return false;
  }
  recorded_hits.insert(key);
  hits.push_back(
      {graph, determinant, root, divisible, positive_definite});
  return true;
}

std::string json_escape(std::string_view value) {
  std::string result;
  result.push_back('"');
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

void append_edges(std::ostream& output, const Graph& graph) {
  output << '[';
  bool first_edge = true;
  for (int index = 0; index < kEdgeCount; ++index) {
    if (!graph[index]) continue;
    if (!first_edge) output << ',';
    first_edge = false;
    output << '[' << kEdges[index].first + 1
           << ',' << kEdges[index].second + 1 << ']';
  }
  output << ']';
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
    const State& current,
    const State& best,
    const Statistics& statistics,
    const std::vector<SquareHit>& hits,
    double elapsed,
    std::string_view termination) {
  std::ostringstream output;
  output << std::setprecision(18);
  output << "{\"best\":{";
  output << "\"edge_count\":" << edge_count(best.graph);
  output << ",\"edges\":";
  append_edges(output, best.graph);
  output << ",\"exact_determinant\":"
         << json_escape(exact_string(best.exact_determinant));
  output << ",\"log_determinant\":"
         << static_cast<double>(best.log_determinant);
  output << "},\"challenge_id\":\"maxdet-23-v1\"";
  output << ",\"claim_boundary\":"
         << json_escape(
                "Normalized Gram screening only; a survivor still requires "
                "an exact G=AA^T decomposition with A in {-1,+1}.")
         ;
  output << ",\"complete\":"
         << (termination == "running" ? "false" : "true");
  output << ",\"current\":{";
  output << "\"edge_count\":" << edge_count(current.graph);
  output << ",\"edges\":";
  append_edges(output, current.graph);
  output << ",\"exact_determinant\":"
         << json_escape(exact_string(current.exact_determinant));
  output << ",\"log_determinant\":"
         << static_cast<double>(current.log_determinant);
  output << "},\"elapsed_seconds\":" << elapsed;
  output << ",\"engine\":\"gram-tabu\"";
  output << ",\"frontier_root\":" << json_escape(
      std::to_string(kFrontierRoot));
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
    output << ",\"edge_count\":" << edge_count(hit.graph);
    output << ",\"edges\":";
    append_edges(output, hit.graph);
    output << ",\"positive_definite\":"
           << (hit.positive_definite ? "true" : "false");
    output << ",\"qualified\":"
           << (hit.divisible && hit.positive_definite ? "true" : "false");
    output << ",\"square_root\":"
           << json_escape(exact_string(hit.root)) << '}';
  }
  output << ']';
  output << ",\"mode\":" << json_escape(arguments.mode);
  output << ",\"normalization\":\"G=24I-J+4A\"";
  output << ",\"parameters\":{";
  output << "\"checkpoint_seconds\":" << arguments.checkpoint_seconds;
  output << ",\"cooling_iterations\":" << arguments.cooling_iterations;
  output << ",\"heartbeat_seconds\":" << arguments.heartbeat_seconds;
  output << ",\"kick_size\":" << arguments.kick_size;
  output << ",\"max_iterations\":" << arguments.max_iterations;
  output << ",\"max_stored_hits\":" << arguments.max_stored_hits;
  output << ",\"restart_iterations\":" << arguments.restart_iterations;
  output << ",\"seconds\":" << arguments.seconds;
  output << ",\"tabu_tenure\":" << arguments.tabu_tenure;
  output << ",\"temperature_end\":"
         << static_cast<double>(arguments.temperature_end);
  output << ",\"temperature_start\":"
         << static_cast<double>(arguments.temperature_start) << '}';
  output << ",\"schema_version\":1";
  output << ",\"seed\":" << arguments.seed;
  output << ",\"start\":" << json_escape(arguments.start.string());
  output << ",\"statistics\":{";
  output << "\"above_frontier\":" << statistics.above_frontier;
  output << ",\"accepted_moves\":" << statistics.accepted_moves;
  output << ",\"exact_neighbor_screens\":"
         << statistics.exact_neighbor_screens;
  output << ",\"exact_pd_checks\":" << statistics.exact_pd_checks;
  output << ",\"exact_squares\":" << statistics.exact_squares;
  output << ",\"iterations\":" << statistics.iterations;
  output << ",\"non_pd_proposals\":" << statistics.non_pd_proposals;
  output << ",\"qualified_survivors\":"
         << statistics.qualified_survivors;
  output << ",\"unrecorded_square_observations\":"
         << statistics.unrecorded_square_observations;
  output << ",\"restarts\":" << statistics.restarts << '}';
  output << ",\"termination\":" << json_escape(termination);
  output << "}\n";
  atomic_write(arguments.output, output.str());
}

struct Candidate {
  int edge_index = -1;
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

bool kick_from_best(
    const State& best,
    State& destination,
    int kick_size,
    std::mt19937_64& randomizer,
    Statistics& statistics,
    std::vector<SquareHit>& hits,
    std::set<GraphKey>& recorded_hits,
    std::uint64_t max_stored_hits) {
  if (kick_size == 0) return false;
  std::array<int, kEdgeCount> indices{};
  std::iota(indices.begin(), indices.end(), 0);
  for (int attempt = 0; attempt < 128; ++attempt) {
    destination = best;
    std::shuffle(indices.begin(), indices.end(), randomizer);
    for (int selected = 0; selected < kick_size; ++selected) {
      destination.graph[indices[selected]] ^= 1U;
    }
    if (!rebuild_numeric(destination)) continue;
    destination.exact_determinant =
        exact_determinant(destination.gram);
    screen_exact_candidate(
        destination.graph,
        destination.gram,
        destination.exact_determinant,
        statistics,
        hits,
        recorded_hits,
        max_stored_hits);
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
    long double log_crt_modulus = 0;
    for (const std::uint64_t prime : kCrtPrimes) {
      if (!is_prime(prime)) {
        throw std::runtime_error("configured CRT modulus is not prime");
      }
      log_crt_modulus += std::log(
          static_cast<long double>(prime));
    }
    const long double log_twice_hadamard_bound =
        std::log(2.0L) +
        static_cast<long double>(kOrder) / 2 *
            std::log(727.0L);
    if (!(log_crt_modulus > log_twice_hadamard_bound)) {
      throw std::runtime_error(
          "CRT modulus does not uniquely cover the determinant bound");
    }
    std::mt19937_64 randomizer(arguments.seed);
    std::uniform_real_distribution<long double> unit(0.0L, 1.0L);

    const SignMatrix start_matrix = read_sign_matrix(arguments.start);
    State current;
    current.graph = graph_from_gram(gram_of(start_matrix));
    if (!rebuild_numeric(current)) {
      throw std::runtime_error("start Gram is not positive definite");
    }
    current.exact_determinant = exact_determinant(current.gram);
    State best = current;

    Statistics statistics;
    std::vector<SquareHit> hits;
    std::set<GraphKey> recorded_hits;
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
      candidates.reserve(kEdgeCount);
      const long double temperature =
          annealing_temperature(arguments, statistics.iterations);

      for (int edge_index = 0; edge_index < kEdgeCount; ++edge_index) {
        Graph candidate_graph = current.graph;
        candidate_graph[edge_index] ^= 1U;
        Gram candidate_gram = current.gram;
        const Edge edge = kEdges[edge_index];
        const int updated_value = candidate_graph[edge_index] ? 3 : -1;
        candidate_gram[edge.first][edge.second] = updated_value;
        candidate_gram[edge.second][edge.first] = updated_value;
        const Exact determinant = exact_determinant(candidate_gram);
        screen_exact_candidate(
            candidate_graph,
            candidate_gram,
            determinant,
            statistics,
            hits,
            recorded_hits,
            arguments.max_stored_hits);

        const long double projected =
            projected_log_determinant(current, edge_index);
        if (!std::isfinite(projected)) continue;
        const bool tabu =
            statistics.iterations < tabu_until[edge_index];
        const bool aspiration =
            determinant > best.exact_determinant;
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
            {edge_index, projected, priority, determinant});
      }

      std::sort(
          candidates.begin(),
          candidates.end(),
          [](const Candidate& left, const Candidate& right) {
            if (left.priority != right.priority) {
              return left.priority > right.priority;
            }
            return left.edge_index < right.edge_index;
          });

      bool moved = false;
      Candidate selected;
      State proposed;
      for (const Candidate& candidate : candidates) {
        proposed = current;
        proposed.graph[candidate.edge_index] ^= 1U;
        if (!rebuild_numeric(proposed)) {
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
        }
        if (current.exact_determinant > best.exact_determinant) {
          best = current;
          last_best_iteration = statistics.iterations;
          write_snapshot(
              arguments,
              current,
              best,
              statistics,
              hits,
              elapsed(),
              "running");
          std::cout << "new best determinant="
                    << exact_string(best.exact_determinant)
                    << " edges=" << edge_count(best.graph)
                    << " iteration=" << statistics.iterations << '\n'
                    << std::flush;
        }
      } else {
        if (!kick_from_best(
                best,
                current,
                arguments.kick_size,
                randomizer,
                statistics,
                hits,
                recorded_hits,
                arguments.max_stored_hits)) {
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
                arguments.kick_size,
                randomizer,
                statistics,
                hits,
                recorded_hits,
                arguments.max_stored_hits)) {
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
    std::cerr << "gram_tabu: " << error.what() << '\n';
    return 2;
  }
}
