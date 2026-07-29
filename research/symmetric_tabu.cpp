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
#include <limits>
#include <random>
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
constexpr int kUndirectedEntries = kOrder * (kOrder + 1) / 2;
constexpr std::size_t kVisitTableSize = 1U << 18U;
constexpr int kBaselineTenure = 11;
constexpr int kMaximumTenure = 128;
constexpr int kTargetedEdges = 30;
constexpr int kRandomCompoundMoves = 128;
constexpr std::uint64_t kExactCheckInterval = 256;
constexpr long double kSingularThreshold = 1.0e-18L;
constexpr long double kExactProbeSlack = 1.0e-8L;

using Clock = std::chrono::steady_clock;
using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using Inverse =
    std::array<std::array<long double, kOrder>, kOrder>;
using Wide = __int128_t;

volatile std::sig_atomic_t stop_requested = 0;

void request_stop(int) { stop_requested = 1; }

struct Arguments {
  std::filesystem::path start;
  std::filesystem::path output;
  std::filesystem::path log;
  std::uint64_t seed = 23;
  double seconds = 3600.0;
  double heartbeat_seconds = 60.0;
  std::uint64_t restart_iterations = 32768;
};

struct Edge {
  int first = -1;
  int second = -1;
  int id = -1;
};

struct State {
  Matrix matrix{};
  Inverse inverse{};
  long double log_abs_determinant =
      -std::numeric_limits<long double>::infinity();
};

struct Candidate {
  int first_edge = -1;
  int second_edge = -1;
  long double projected_log_abs_determinant =
      -std::numeric_limits<long double>::infinity();
  bool aspiration = false;
};

struct DirectedUpdate {
  int row = -1;
  int column = -1;
  long double delta = 0.0L;
};

struct Visit {
  std::uint64_t hash = 0;
  std::uint64_t iteration = 0;
  bool occupied = false;
};

struct Statistics {
  std::uint64_t iterations = 0;
  std::uint64_t candidates_evaluated = 0;
  std::uint64_t compound_moves = 0;
  std::uint64_t worsening_moves = 0;
  std::uint64_t exact_checks = 0;
  std::uint64_t inverse_rebuilds = 0;
  std::uint64_t singular_moves_rejected = 0;
  std::uint64_t cycles = 0;
  std::uint64_t tabu_resets = 0;
  std::uint64_t restarts = 0;
};

Wide absolute(Wide value) { return value < 0 ? -value : value; }

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

Matrix read_matrix(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("cannot open start matrix: " + path.string());
  }
  Matrix matrix{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      if (!(input >> matrix[row][column]) ||
          (matrix[row][column] != -1 && matrix[row][column] != 1)) {
        throw std::runtime_error(
            "start matrix must contain exactly 23x23 entries in {-1,+1}");
      }
    }
  }
  std::string extra;
  if (input >> extra) {
    throw std::runtime_error("start matrix contains extra data");
  }
  for (int row = 0; row < kOrder; ++row) {
    for (int column = row + 1; column < kOrder; ++column) {
      if (matrix[row][column] != matrix[column][row]) {
        throw std::runtime_error("start matrix must be symmetric");
      }
    }
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

void atomic_write_matrix(const std::filesystem::path& path,
                         const Matrix& matrix, std::uint64_t nonce) {
  const std::filesystem::path directory =
      path.parent_path().empty() ? std::filesystem::path(".")
                                 : path.parent_path();
  std::filesystem::create_directories(directory);
  const std::string bytes = matrix_bytes(matrix);

  std::filesystem::path temporary;
  int descriptor = -1;
  for (std::uint64_t attempt = 0; attempt < 64; ++attempt) {
    temporary =
        directory /
        ("." + path.filename().string() + ".symmetric-tabu-" +
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
          "cannot create checkpoint temporary file: " +
          std::string(std::strerror(errno)));
    }
  }
  if (descriptor < 0) {
    throw std::runtime_error(
        "cannot allocate a unique checkpoint temporary file");
  }

  bool renamed = false;
  try {
    write_all(descriptor, bytes);
    if (::fsync(descriptor) != 0) {
      throw std::runtime_error(
          "cannot sync checkpoint: " + std::string(std::strerror(errno)));
    }
    if (::close(descriptor) != 0) {
      descriptor = -1;
      throw std::runtime_error(
          "cannot close checkpoint: " + std::string(std::strerror(errno)));
    }
    descriptor = -1;
    if (::rename(temporary.c_str(), path.c_str()) != 0) {
      throw std::runtime_error(
          "cannot atomically install checkpoint: " +
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

  long double log_abs_determinant = 0.0L;
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
      state.log_abs_determinant =
          -std::numeric_limits<long double>::infinity();
      return false;
    }
    if (pivot_row != column) {
      std::swap(augmented[pivot_row], augmented[column]);
    }
    const long double pivot = augmented[column][column];
    log_abs_determinant += std::log(std::fabs(pivot));
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
  state.log_abs_determinant = log_abs_determinant;
  return true;
}

std::vector<Edge> make_edges() {
  std::vector<Edge> edges;
  edges.reserve(kUndirectedEntries);
  for (int first = 0; first < kOrder; ++first) {
    for (int second = first; second < kOrder; ++second) {
      edges.push_back(
          Edge{first, second, static_cast<int>(edges.size())});
    }
  }
  return edges;
}

bool disjoint(const Edge& first, const Edge& second) {
  return first.first != second.first &&
         first.first != second.second &&
         first.second != second.first &&
         first.second != second.second;
}

void append_updates(const State& state, const Edge& edge,
                    std::array<DirectedUpdate, 4>& updates, int& count) {
  const long double delta =
      -2.0L *
      static_cast<long double>(state.matrix[edge.first][edge.second]);
  updates[count++] =
      DirectedUpdate{edge.first, edge.second, delta};
  if (edge.first != edge.second) {
    updates[count++] =
        DirectedUpdate{edge.second, edge.first, delta};
  }
}

long double small_determinant(
    std::array<std::array<long double, 4>, 4> matrix, int size) {
  long double determinant = 1.0L;
  int sign = 1;
  for (int column = 0; column < size; ++column) {
    int pivot_row = column;
    for (int row = column + 1; row < size; ++row) {
      if (std::fabs(matrix[row][column]) >
          std::fabs(matrix[pivot_row][column])) {
        pivot_row = row;
      }
    }
    if (std::fabs(matrix[pivot_row][column]) <
        kSingularThreshold) {
      return 0.0L;
    }
    if (pivot_row != column) {
      std::swap(matrix[pivot_row], matrix[column]);
      sign = -sign;
    }
    const long double pivot = matrix[column][column];
    determinant *= pivot;
    for (int row = column + 1; row < size; ++row) {
      const long double factor = matrix[row][column] / pivot;
      for (int inner = column + 1; inner < size; ++inner) {
        matrix[row][inner] -= factor * matrix[column][inner];
      }
    }
  }
  return static_cast<long double>(sign) * determinant;
}

long double projected_log_score(const State& state, const Edge& first,
                                const Edge* second) {
  std::array<DirectedUpdate, 4> updates{};
  int count = 0;
  append_updates(state, first, updates, count);
  if (second != nullptr) append_updates(state, *second, updates, count);

  std::array<std::array<long double, 4>, 4> lemma{};
  for (int row = 0; row < count; ++row) {
    for (int column = 0; column < count; ++column) {
      lemma[row][column] =
          (row == column ? 1.0L : 0.0L) +
          updates[column].delta *
              state.inverse[updates[row].column][updates[column].row];
    }
  }
  const long double ratio = small_determinant(lemma, count);
  if (std::fabs(ratio) < kSingularThreshold) {
    return -std::numeric_limits<long double>::infinity();
  }
  return state.log_abs_determinant + std::log(std::fabs(ratio));
}

void flip_edge(Matrix& matrix, const Edge& edge) {
  matrix[edge.first][edge.second] *= -1;
  if (edge.first != edge.second) {
    matrix[edge.second][edge.first] *= -1;
  }
}

std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

std::array<std::uint64_t, kUndirectedEntries> make_zobrist(
    std::uint64_t seed) {
  std::array<std::uint64_t, kUndirectedEntries> values{};
  std::uint64_t state = seed ^ 0x53796d6d65747279ULL;
  for (std::uint64_t& value : values) {
    state = splitmix64(state);
    value = state;
  }
  return values;
}

std::uint64_t matrix_hash(
    const Matrix& matrix, const std::vector<Edge>& edges,
    const std::array<std::uint64_t, kUndirectedEntries>& zobrist) {
  std::uint64_t hash = 0;
  for (const Edge& edge : edges) {
    if (matrix[edge.first][edge.second] == 1) {
      hash ^= zobrist[edge.id];
    }
  }
  return hash;
}

long double temperature(std::uint64_t iteration,
                        std::uint64_t last_improvement) {
  constexpr std::uint64_t period = 4096;
  const long double phase =
      static_cast<long double>(iteration % period) /
      static_cast<long double>(period);
  const long double triangle =
      phase <= 0.5L ? 2.0L * phase : 2.0L * (1.0L - phase);
  constexpr long double minimum = 0.018L;
  constexpr long double maximum = 0.30L;
  long double result =
      std::exp(std::log(minimum) +
               triangle * std::log(maximum / minimum));
  const std::uint64_t stagnant = iteration - last_improvement;
  const long double reheat =
      1.0L +
      std::min(1.0L, static_cast<long double>(stagnant) / 16384.0L);
  result *= reheat;
  return std::min(result, 0.45L);
}

std::vector<Candidate> build_candidates(
    const State& state, const std::vector<Edge>& edges,
    const std::array<std::uint64_t, kUndirectedEntries>& tabu_until,
    std::uint64_t iteration, long double best_exact_log,
    std::mt19937_64& randomizer, Statistics& statistics,
    bool ignore_tabu) {
  std::vector<Candidate> candidates;
  candidates.reserve(kUndirectedEntries + 512);
  std::vector<std::pair<long double, int>> ranked_off_diagonal;
  ranked_off_diagonal.reserve(kUndirectedEntries - kOrder);

  auto add_candidate = [&](int first_id, int second_id,
                           long double projected) {
    ++statistics.candidates_evaluated;
    if (!std::isfinite(projected)) return;
    const bool tabu =
        !ignore_tabu &&
        (iteration < tabu_until[first_id] ||
         (second_id >= 0 && iteration < tabu_until[second_id]));
    const bool aspiration = projected > best_exact_log + 1.0e-12L;
    if (tabu && !aspiration) return;
    candidates.push_back(
        Candidate{first_id, second_id, projected, aspiration});
  };

  for (const Edge& edge : edges) {
    const long double projected =
        projected_log_score(state, edge, nullptr);
    add_candidate(edge.id, -1, projected);
    if (edge.first != edge.second && std::isfinite(projected)) {
      ranked_off_diagonal.emplace_back(projected, edge.id);
    }
  }

  const std::size_t target_count =
      std::min<std::size_t>(kTargetedEdges,
                            ranked_off_diagonal.size());
  std::partial_sort(
      ranked_off_diagonal.begin(),
      ranked_off_diagonal.begin() +
          static_cast<std::ptrdiff_t>(target_count),
      ranked_off_diagonal.end(),
      [](const auto& first, const auto& second) {
        return first.first > second.first;
      });
  ranked_off_diagonal.resize(target_count);
  for (std::size_t first = 0; first < target_count; ++first) {
    for (std::size_t second = first + 1; second < target_count;
         ++second) {
      const Edge& first_edge = edges[ranked_off_diagonal[first].second];
      const Edge& second_edge =
          edges[ranked_off_diagonal[second].second];
      if (!disjoint(first_edge, second_edge)) continue;
      add_candidate(
          first_edge.id, second_edge.id,
          projected_log_score(state, first_edge, &second_edge));
    }
  }

  std::array<int, kOrder> vertices{};
  for (int index = 0; index < kOrder; ++index) vertices[index] = index;
  for (int sample = 0; sample < kRandomCompoundMoves; ++sample) {
    for (int index = 0; index < 4; ++index) {
      const int remaining = kOrder - index;
      const int selected =
          index + static_cast<int>(randomizer() %
                                   static_cast<std::uint64_t>(remaining));
      std::swap(vertices[index], vertices[selected]);
    }
    int first_vertex = vertices[0];
    int second_vertex = vertices[1];
    int third_vertex = vertices[2];
    int fourth_vertex = vertices[3];
    if (first_vertex > second_vertex) {
      std::swap(first_vertex, second_vertex);
    }
    if (third_vertex > fourth_vertex) {
      std::swap(third_vertex, fourth_vertex);
    }
    const auto find_edge = [&](int first, int second) -> const Edge& {
      const auto found = std::find_if(
          edges.begin(), edges.end(), [&](const Edge& edge) {
            return edge.first == first && edge.second == second;
          });
      if (found == edges.end()) {
        throw std::runtime_error("internal edge lookup failed");
      }
      return *found;
    };
    const Edge& first_edge = find_edge(first_vertex, second_vertex);
    const Edge& second_edge = find_edge(third_vertex, fourth_vertex);
    add_candidate(
        first_edge.id, second_edge.id,
        projected_log_score(state, first_edge, &second_edge));
  }
  return candidates;
}

Candidate choose_candidate(const std::vector<Candidate>& candidates,
                           long double current_log,
                           long double current_temperature,
                           std::mt19937_64& randomizer) {
  if (candidates.empty()) {
    throw std::runtime_error("no admissible symmetric move");
  }
  const long double maximum =
      std::max_element(
          candidates.begin(), candidates.end(),
          [](const Candidate& first, const Candidate& second) {
            return first.projected_log_abs_determinant <
                   second.projected_log_abs_determinant;
          })
          ->projected_log_abs_determinant;

  std::vector<long double> weights;
  weights.reserve(candidates.size());
  long double total = 0.0L;
  for (const Candidate& candidate : candidates) {
    const long double exponent =
        (candidate.projected_log_abs_determinant - maximum) /
        current_temperature;
    const long double weight =
        exponent < -80.0L ? 0.0L : std::exp(exponent);
    weights.push_back(weight);
    total += weight;
  }
  if (!(total > 0.0L) || !std::isfinite(total)) {
    throw std::runtime_error("invalid simulated-tempering weights");
  }

  const long double unit =
      std::generate_canonical<long double, 64>(randomizer);
  long double threshold = unit * total;
  for (std::size_t index = 0; index < candidates.size(); ++index) {
    threshold -= weights[index];
    if (threshold <= 0.0L) return candidates[index];
  }
  (void)current_log;
  return candidates.back();
}

std::uint64_t strict_unsigned(std::string_view text,
                              std::string_view option) {
  if (text.empty() ||
      !std::all_of(text.begin(), text.end(),
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
    } else if (option == "--log") {
      arguments.log = value();
    } else if (option == "--seed") {
      arguments.seed = strict_unsigned(value(), option);
    } else if (option == "--seconds") {
      arguments.seconds = strict_double(value(), option, false);
    } else if (option == "--heartbeat-seconds" ||
               option == "--heartbeat") {
      arguments.heartbeat_seconds =
          strict_double(value(), option, true);
    } else if (option == "--restart-iterations") {
      arguments.restart_iterations = strict_unsigned(value(), option);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (arguments.start.empty()) {
    throw std::runtime_error("--start is required");
  }
  if (arguments.output.empty()) {
    throw std::runtime_error("--output is required");
  }
  if (arguments.log.empty()) {
    throw std::runtime_error("--log is required");
  }
  const auto normalized = [](const std::filesystem::path& path) {
    return std::filesystem::absolute(path).lexically_normal();
  };
  const auto start = normalized(arguments.start);
  const auto output = normalized(arguments.output);
  const auto log = normalized(arguments.log);
  if (start == output || start == log || output == log) {
    throw std::runtime_error(
        "--start, --output, and --log must be distinct");
  }
  return arguments;
}

void log_record(std::ofstream& log, const Arguments& arguments,
                const Statistics& statistics, const char* event,
                double elapsed_seconds, Wide best_score,
                Wide last_exact_score, int tenure,
                long double current_temperature) {
  log << "{\"absolute_determinant\":\""
      << wide_to_string(best_score)
      << "\",\"candidates_evaluated\":"
      << statistics.candidates_evaluated
      << ",\"compound_moves\":" << statistics.compound_moves
      << ",\"cycles\":" << statistics.cycles
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(3) << elapsed_seconds
      << ",\"event\":\"" << event
      << "\",\"exact_checks\":" << statistics.exact_checks
      << ",\"inverse_rebuilds\":" << statistics.inverse_rebuilds
      << ",\"iterations\":" << statistics.iterations
      << ",\"last_exact_determinant\":\""
      << wide_to_string(last_exact_score)
      << "\",\"restarts\":" << statistics.restarts
      << ",\"seed\":" << arguments.seed
      << ",\"singular_moves_rejected\":"
      << statistics.singular_moves_rejected
      << ",\"tabu_resets\":" << statistics.tabu_resets
      << ",\"temperature\":" << std::setprecision(6)
      << static_cast<double>(current_temperature)
      << ",\"tenure\":" << tenure
      << ",\"worsening_moves\":" << statistics.worsening_moves
      << "}\n";
  log.flush();
  if (!log) throw std::runtime_error("cannot append research log");
}

bool exact_check_and_promote(
    const State& state, Matrix& best_matrix, Wide& best_score,
    Wide& last_exact_score,
    const Arguments& arguments, Statistics& statistics,
    std::ofstream& log, const Clock::time_point& started, int tenure,
    long double current_temperature, std::uint64_t checkpoint_nonce) {
  const Wide exact_score = absolute(exact_determinant(state.matrix));
  ++statistics.exact_checks;
  last_exact_score = exact_score;
  if (exact_score <= best_score) return false;

  best_score = exact_score;
  best_matrix = state.matrix;
  atomic_write_matrix(arguments.output, state.matrix, checkpoint_nonce);
  const double elapsed =
      std::chrono::duration<double>(Clock::now() - started).count();
  log_record(log, arguments, statistics, "new_best", elapsed, best_score,
             last_exact_score, tenure, current_temperature);
  std::cout << "new best |det|=" << wide_to_string(best_score)
            << " iteration=" << statistics.iterations << '\n'
            << std::flush;
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    if (!arguments.log.parent_path().empty()) {
      std::filesystem::create_directories(arguments.log.parent_path());
    }
    std::ofstream log(arguments.log, std::ios::app);
    if (!log) throw std::runtime_error("cannot open research log");

    const std::vector<Edge> edges = make_edges();
    State state;
    state.matrix = read_matrix(arguments.start);
    if (!rebuild_inverse(state)) {
      throw std::runtime_error("start matrix must be nonsingular");
    }

    Statistics statistics;
    statistics.inverse_rebuilds = 1;
    Matrix best_matrix = state.matrix;
    Wide best_score = absolute(exact_determinant(state.matrix));
    ++statistics.exact_checks;
    if (best_score == 0) {
      throw std::runtime_error("start matrix must have nonzero determinant");
    }
    Wide last_exact_score = best_score;
    std::uint64_t checkpoint_nonce = 0;
    atomic_write_matrix(
        arguments.output, state.matrix, checkpoint_nonce++);

    std::mt19937_64 randomizer(arguments.seed);
    const auto zobrist = make_zobrist(arguments.seed);
    std::uint64_t current_hash =
        matrix_hash(state.matrix, edges, zobrist);
    std::vector<Visit> visits(kVisitTableSize);
    visits[current_hash & (kVisitTableSize - 1U)] =
        Visit{current_hash, 0, true};
    std::array<std::uint64_t, kUndirectedEntries> tabu_until{};
    int tenure =
        kBaselineTenure + static_cast<int>(arguments.seed % 5U);
    const int baseline_tenure = tenure;
    std::uint64_t last_cycle_iteration = 0;
    std::uint64_t last_improvement_iteration = 0;

    const auto started = Clock::now();
    const auto deadline =
        started + std::chrono::duration<double>(arguments.seconds);
    auto next_heartbeat =
        started +
        std::chrono::duration<double>(arguments.heartbeat_seconds);

    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);
    log_record(log, arguments, statistics, "start", 0.0, best_score,
               last_exact_score, tenure,
               temperature(0, last_improvement_iteration));

    while (!stop_requested && Clock::now() < deadline) {
      ++statistics.iterations;
      if (arguments.restart_iterations != 0 &&
          statistics.iterations % arguments.restart_iterations == 0) {
        state.matrix = best_matrix;
        std::vector<int> kick_ids(kUndirectedEntries);
        for (int id = 0; id < kUndirectedEntries; ++id) {
          kick_ids[static_cast<std::size_t>(id)] = id;
        }
        std::shuffle(kick_ids.begin(), kick_ids.end(), randomizer);
        const int kick_size =
            2 + static_cast<int>(randomizer() % 17U);
        for (int kick = 0; kick < kick_size; ++kick) {
          flip_edge(
              state.matrix,
              edges[kick_ids[static_cast<std::size_t>(kick)]]);
        }
        if (!rebuild_inverse(state)) {
          state.matrix = best_matrix;
          if (!rebuild_inverse(state)) {
            throw std::runtime_error(
                "cannot rebuild the nonsingular incumbent at restart");
          }
        }
        ++statistics.inverse_rebuilds;
        ++statistics.restarts;
        tabu_until.fill(0);
        std::fill(visits.begin(), visits.end(), Visit{});
        current_hash = matrix_hash(state.matrix, edges, zobrist);
        visits[current_hash & (kVisitTableSize - 1U)] =
            Visit{current_hash, statistics.iterations, true};
        tenure = baseline_tenure;
        last_cycle_iteration = statistics.iterations;
        last_improvement_iteration = statistics.iterations;
        const double elapsed =
            std::chrono::duration<double>(
                Clock::now() - started).count();
        log_record(
            log, arguments, statistics, "restart", elapsed, best_score,
            last_exact_score, tenure,
            temperature(
                statistics.iterations, last_improvement_iteration));
      }
      const long double current_temperature =
          temperature(statistics.iterations,
                      last_improvement_iteration);
      const long double best_exact_log =
          std::log(static_cast<long double>(best_score));
      std::vector<Candidate> candidates = build_candidates(
          state, edges, tabu_until, statistics.iterations,
          best_exact_log, randomizer, statistics, false);
      if (candidates.empty()) {
        tabu_until.fill(0);
        ++statistics.tabu_resets;
        candidates = build_candidates(
            state, edges, tabu_until, statistics.iterations,
            best_exact_log, randomizer, statistics, true);
      }

      const Candidate selected =
          choose_candidate(candidates, state.log_abs_determinant,
                           current_temperature, randomizer);
      const long double previous_log = state.log_abs_determinant;
      const Matrix before = state.matrix;
      flip_edge(state.matrix, edges[selected.first_edge]);
      if (selected.second_edge >= 0) {
        flip_edge(state.matrix, edges[selected.second_edge]);
      }
      if (!rebuild_inverse(state)) {
        state.matrix = before;
        if (!rebuild_inverse(state)) {
          throw std::runtime_error(
              "cannot recover inverse after a singular move");
        }
        ++statistics.inverse_rebuilds;
        ++statistics.singular_moves_rejected;
        tabu_until[selected.first_edge] =
            statistics.iterations + kMaximumTenure;
        if (selected.second_edge >= 0) {
          tabu_until[selected.second_edge] =
              statistics.iterations + kMaximumTenure;
        }
        continue;
      }
      ++statistics.inverse_rebuilds;
      if (selected.second_edge >= 0) ++statistics.compound_moves;
      if (state.log_abs_determinant < previous_log) {
        ++statistics.worsening_moves;
      }

      const int jitter = static_cast<int>(randomizer() % 7U);
      const std::uint64_t expires =
          statistics.iterations +
          static_cast<std::uint64_t>(tenure + jitter + 1);
      tabu_until[selected.first_edge] = expires;
      current_hash ^= zobrist[selected.first_edge];
      if (selected.second_edge >= 0) {
        tabu_until[selected.second_edge] = expires;
        current_hash ^= zobrist[selected.second_edge];
      }

      Visit& visit =
          visits[current_hash & (kVisitTableSize - 1U)];
      if (visit.occupied && visit.hash == current_hash &&
          statistics.iterations > visit.iteration) {
        const std::uint64_t cycle_length =
            statistics.iterations - visit.iteration;
        if (cycle_length <= 512U) {
          ++statistics.cycles;
          last_cycle_iteration = statistics.iterations;
          tenure = std::min(
              kMaximumTenure,
              tenure + 2 +
                  static_cast<int>(
                      std::min<std::uint64_t>(cycle_length / 8U, 10U)));
        }
      }
      visit = Visit{current_hash, statistics.iterations, true};
      if (statistics.iterations - last_cycle_iteration >= 1024U &&
          (statistics.iterations & 127U) == 0U &&
          tenure > baseline_tenure) {
        --tenure;
      }

      const bool near_best =
          selected.aspiration ||
          state.log_abs_determinant + kExactProbeSlack >= best_exact_log;
      const bool periodic =
          statistics.iterations % kExactCheckInterval == 0U;
      if (near_best || periodic) {
        if (exact_check_and_promote(
                state, best_matrix, best_score, last_exact_score, arguments,
                statistics, log, started, tenure, current_temperature,
                checkpoint_nonce++)) {
          last_improvement_iteration = statistics.iterations;
        }
      }

      const auto now = Clock::now();
      if (arguments.heartbeat_seconds > 0.0 &&
          now >= next_heartbeat) {
        const double elapsed =
            std::chrono::duration<double>(now - started).count();
        log_record(log, arguments, statistics, "heartbeat", elapsed,
                   best_score, last_exact_score, tenure,
                   current_temperature);
        next_heartbeat =
            now +
            std::chrono::duration<double>(
                arguments.heartbeat_seconds);
      }
    }

    const long double final_temperature =
        temperature(statistics.iterations,
                    last_improvement_iteration);
    exact_check_and_promote(
        state, best_matrix, best_score, last_exact_score, arguments, statistics,
        log, started, tenure, final_temperature, checkpoint_nonce++);
    const double elapsed =
        std::chrono::duration<double>(Clock::now() - started).count();
    log_record(log, arguments, statistics,
               stop_requested ? "stopped" : "finished", elapsed,
               best_score, last_exact_score, tenure, final_temperature);
    const double moves_per_second =
        elapsed > 0.0
            ? static_cast<double>(statistics.iterations) / elapsed
            : 0.0;
    std::cout << "finished |det|=" << wide_to_string(best_score)
              << " iterations=" << statistics.iterations
              << " compound=" << statistics.compound_moves
              << " moves_per_second=" << std::fixed
              << std::setprecision(1) << moves_per_second << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "symmetric_tabu: " << error.what() << '\n';
    return 2;
  }
}
