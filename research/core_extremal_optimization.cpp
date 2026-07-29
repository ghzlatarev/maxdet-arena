// Exact extremal-optimization search for the order-23 maximal determinant
// problem.
//
// The underlying state and determinant arithmetic come from the independently
// audited dephased-core engine in core_adjugate_tabu.cpp.  A 23 by 23 sign
// matrix is represented by a 22 by 22 binary core B, with
//
//     |det(H)| = 2^22 |det(B)|.
//
// All candidate determinants and all accepted rank-one updates are exact
// integers.  Unlike hill climbing or tabu search, this engine always takes a
// nonsingular move.  It ranks the 484 bit flips by their exact resulting
// determinant (minus an adaptive integer recency penalty), samples a rank from
// a truncated power law, and applies that move even when it is downhill.
//
// Periodic secondary pools expose all 45 dephasing-complement directions and
// balanced same-row/same-column two-bit swaps.  Exact multi-bit kicks rebuild
// determinant and adjugate from scratch.  The global incumbent is kept
// separately and checkpointed atomically.
//
// This is a stochastic research solver, not an exhaustive search.

#include <optional>

#define main core_adjugate_tabu_embedded_main
#include "core_adjugate_tabu.cpp"
#undef main

namespace {

constexpr std::uint64_t kEoFrontier =
    UINT64_C(2779447296000000);
constexpr std::uint64_t kNeverMoved =
    std::numeric_limits<std::uint64_t>::max();
constexpr std::uint64_t kEoIdentityCheckInterval = UINT64_C(8192);
constexpr std::uint64_t kEoDeterminantCheckInterval = UINT64_C(65536);
constexpr int kComplementMoveCount = 2 * kCoreOrder + 1;
constexpr int kMaximumPairCandidates =
    2 * kCoreOrder * (kCoreOrder * (kCoreOrder - 1) / 2);

struct EoArguments {
  fs::path start;
  fs::path output;
  fs::path log;
  fs::path summary;
  std::uint64_t seed = UINT64_C(34001);
  std::uint64_t frontier = kEoFrontier;
  double seconds = 300.0;
  double heartbeat_seconds = 15.0;
  double tau = 1.45;
  std::uint64_t max_iterations = 0;
  std::uint64_t recency_window = 32;
  std::uint64_t recency_penalty = UINT64_C(250000);
  std::uint64_t recency_penalty_max = UINT64_C(4000000);
  std::uint64_t complement_interval = 64;
  std::uint64_t pair_interval = 257;
  std::uint64_t kick_interval = UINT64_C(400000);
  int kick_min_flips = 10;
  int kick_max_flips = 18;
  int kick_attempts = 24;
  int self_test_rounds = 0;
};

enum class EoMoveKind : std::uint8_t {
  kBit,
  kRowComplement,
  kColumnComplement,
  kWholeComplement,
  kRowSwap,
  kColumnSwap,
};

struct EoMove {
  EoMoveKind kind = EoMoveKind::kBit;
  int first = -1;
  int second = -1;
  int third = -1;
  int recency_id = -1;
  std::int64_t determinant = 0;
};

struct RankedMove {
  EoMove move{};
  std::int64_t adjusted_fitness = 0;
  std::uint64_t determinant_magnitude = 0;
  std::uint64_t tie_breaker = 0;
};

struct EoStatistics {
  std::uint64_t iterations = 0;
  std::uint64_t exact_bit_directions = 0;
  std::uint64_t exact_complement_directions = 0;
  std::uint64_t exact_pair_directions = 0;
  std::uint64_t singular_directions = 0;
  std::uint64_t bit_moves = 0;
  std::uint64_t row_complements = 0;
  std::uint64_t column_complements = 0;
  std::uint64_t whole_complements = 0;
  std::uint64_t row_swaps = 0;
  std::uint64_t column_swaps = 0;
  std::uint64_t uphill_moves = 0;
  std::uint64_t downhill_moves = 0;
  std::uint64_t level_moves = 0;
  std::uint64_t sampled_rank_sum = 0;
  std::uint64_t sampled_rank_max = 0;
  std::uint64_t potential_cycles = 0;
  std::uint64_t immediate_refips = 0;
  std::uint64_t penalty_increases = 0;
  std::uint64_t penalty_decays = 0;
  std::uint64_t exact_kicks = 0;
  std::uint64_t kick_attempts = 0;
  std::uint64_t kick_singular = 0;
  std::uint64_t promotions = 0;
  std::uint64_t equal_best_hits = 0;
  std::uint64_t distinct_tie_checkpoints = 0;
  std::uint64_t identity_checks = 0;
  std::uint64_t determinant_checks = 0;
  std::uint64_t sampled_hamming_max = 0;
  std::uint64_t best_strictly_below = 0;
};

std::uint64_t eo_strict_positive(
    std::string_view text, std::string_view option,
    std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max()) {
  const std::uint64_t value = strict_unsigned(text, option);
  if (value == 0 || value > maximum) {
    throw std::runtime_error(
        std::string(option) + " is outside its valid positive range");
  }
  return value;
}

int eo_strict_positive_int(
    std::string_view text, std::string_view option, int maximum) {
  return static_cast<int>(
      eo_strict_positive(
          text, option, static_cast<std::uint64_t>(maximum)));
}

[[maybe_unused]] EoArguments parse_eo_arguments(
    int argc, char** argv) {
  EoArguments arguments;
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
    } else if (option == "--summary") {
      arguments.summary = value();
    } else if (option == "--seed") {
      arguments.seed = strict_unsigned(value(), option);
    } else if (option == "--frontier") {
      arguments.frontier = strict_unsigned(value(), option);
    } else if (option == "--seconds") {
      arguments.seconds = strict_double(value(), option, false);
    } else if (option == "--heartbeat" ||
               option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds =
          strict_double(value(), option, true);
    } else if (option == "--tau") {
      arguments.tau = strict_double(value(), option, false);
      if (arguments.tau < 0.1 || arguments.tau > 10.0) {
        throw std::runtime_error("--tau must be in [0.1,10]");
      }
    } else if (option == "--max-iterations") {
      arguments.max_iterations = strict_unsigned(value(), option);
    } else if (option == "--recency-window") {
      arguments.recency_window =
          eo_strict_positive(value(), option, UINT64_C(1000000));
    } else if (option == "--recency-penalty") {
      arguments.recency_penalty =
          eo_strict_positive(value(), option, UINT64_C(1000000000));
    } else if (option == "--recency-penalty-max") {
      arguments.recency_penalty_max =
          eo_strict_positive(value(), option, UINT64_C(1000000000));
    } else if (option == "--complement-interval") {
      arguments.complement_interval = strict_unsigned(value(), option);
    } else if (option == "--pair-interval") {
      arguments.pair_interval = strict_unsigned(value(), option);
    } else if (option == "--kick-interval") {
      arguments.kick_interval = strict_unsigned(value(), option);
    } else if (option == "--kick-min-flips") {
      arguments.kick_min_flips =
          eo_strict_positive_int(value(), option, kCoreEntries);
    } else if (option == "--kick-max-flips") {
      arguments.kick_max_flips =
          eo_strict_positive_int(value(), option, kCoreEntries);
    } else if (option == "--kick-attempts") {
      arguments.kick_attempts =
          eo_strict_positive_int(
              value(), option, std::numeric_limits<int>::max());
    } else if (option == "--self-test") {
      arguments.self_test_rounds =
          eo_strict_positive_int(value(), option, 100000);
    } else if (option == "--help") {
      std::cout
          << "usage: core_extremal_optimization --start MATRIX "
             "--output MATRIX --log JSONL --summary JSON [options]\n"
          << "  --seed N --frontier N --seconds S "
             "--heartbeat-seconds S --max-iterations N\n"
          << "  --tau X --recency-window N --recency-penalty N "
             "--recency-penalty-max N\n"
          << "  --complement-interval N --pair-interval N\n"
          << "  --kick-interval N --kick-min-flips N "
             "--kick-max-flips N --kick-attempts N\n"
          << "audit mode: --start MATRIX --self-test ROUNDS\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }

  if (arguments.start.empty()) {
    throw std::runtime_error("--start is required");
  }
  if (arguments.kick_min_flips > arguments.kick_max_flips) {
    throw std::runtime_error(
        "--kick-min-flips must not exceed --kick-max-flips");
  }
  if (arguments.recency_penalty > arguments.recency_penalty_max) {
    throw std::runtime_error(
        "--recency-penalty must not exceed --recency-penalty-max");
  }
  if (arguments.self_test_rounds == 0) {
    if (arguments.output.empty() || arguments.log.empty() ||
        arguments.summary.empty()) {
      throw std::runtime_error(
          "--output, --log, and --summary are required in search mode");
    }
    const std::array<fs::path, 4> paths{
        arguments.start, arguments.output, arguments.log,
        arguments.summary};
    for (std::size_t first = 0; first < paths.size(); ++first) {
      for (std::size_t second = first + 1;
           second < paths.size(); ++second) {
        if (paths[first] == paths[second]) {
          throw std::runtime_error(
              "input and output paths must be distinct");
        }
      }
    }
    for (const fs::path& path :
         {arguments.output, arguments.log, arguments.summary}) {
      if (fs::exists(path)) {
        throw std::runtime_error(
            "refusing to overwrite existing output: " + path.string());
      }
    }
  }
  return arguments;
}

State eo_initial_state(const fs::path& path) {
  const SignMatrix input = read_sign_matrix(path);
  State state;
  state.core = dephase_to_core(input);
  state.determinant = exact_core_determinant(state.core);
  if (state.determinant == 0) {
    throw std::runtime_error("start matrix must be nonsingular");
  }
  state.adjugate = exact_adjugate(state.core);
  check_adjugate_identity(state);

  std::vector<std::vector<Wide>> sign_wide(
      kSignOrder, std::vector<Wide>(kSignOrder));
  for (int row = 0; row < kSignOrder; ++row) {
    for (int column = 0; column < kSignOrder; ++column) {
      sign_wide[row][column] = input[row][column];
    }
  }
  const Wide sign_determinant = bareiss(sign_wide);
  const Wide scaled =
      static_cast<Wide>(state.determinant) * kScale;
  if (sign_determinant != scaled &&
      sign_determinant != -scaled) {
    throw std::runtime_error(
        "dephased core determinant disagrees with sign matrix");
  }
  return state;
}

class RankSampler {
 public:
  RankSampler(int maximum_rank, double tau) {
    cumulative_.reserve(static_cast<std::size_t>(maximum_rank));
    long double total = 0.0L;
    for (int rank = 1; rank <= maximum_rank; ++rank) {
      total += std::pow(
          static_cast<long double>(rank),
          -static_cast<long double>(tau));
      cumulative_.push_back(total);
    }
  }

  std::size_t sample(
      std::size_t available, std::mt19937_64& randomizer) const {
    if (available == 0 || available > cumulative_.size()) {
      throw std::runtime_error("invalid EO rank pool size");
    }
    const long double unit =
        std::generate_canonical<long double, 64>(randomizer);
    const long double target =
        unit * cumulative_[available - 1];
    const auto found = std::lower_bound(
        cumulative_.begin(),
        cumulative_.begin() +
            static_cast<std::ptrdiff_t>(available),
        target);
    return static_cast<std::size_t>(
        std::distance(cumulative_.begin(), found));
  }

 private:
  std::vector<long double> cumulative_;
};

std::uint64_t recency_urgency(
    std::uint64_t last_moved, std::uint64_t iteration,
    std::uint64_t window) {
  if (last_moved == kNeverMoved || iteration <= last_moved) {
    return 0;
  }
  const std::uint64_t age = iteration - last_moved;
  return age <= window ? window + 1 - age : 0;
}

std::int64_t adjusted_fitness(
    const EoMove& move,
    const std::array<std::uint64_t, kCoreEntries>& last_moved,
    std::uint64_t iteration, std::uint64_t recency_window,
    std::uint64_t penalty) {
  std::uint64_t urgency = 0;
  if (move.kind == EoMoveKind::kBit) {
    urgency = recency_urgency(
        last_moved[static_cast<std::size_t>(move.recency_id)],
        iteration, recency_window);
  } else if (move.kind == EoMoveKind::kRowSwap) {
    const int row = move.first;
    const int first_id = row * kCoreOrder + move.second;
    const int second_id = row * kCoreOrder + move.third;
    urgency =
        recency_urgency(
            last_moved[static_cast<std::size_t>(first_id)],
            iteration, recency_window) +
        recency_urgency(
            last_moved[static_cast<std::size_t>(second_id)],
            iteration, recency_window);
  } else if (move.kind == EoMoveKind::kColumnSwap) {
    const int column = move.first;
    const int first_id = move.second * kCoreOrder + column;
    const int second_id = move.third * kCoreOrder + column;
    urgency =
        recency_urgency(
            last_moved[static_cast<std::size_t>(first_id)],
            iteration, recency_window) +
        recency_urgency(
            last_moved[static_cast<std::size_t>(second_id)],
            iteration, recency_window);
  }
  const Wide adjusted =
      static_cast<Wide>(magnitude(move.determinant)) -
      static_cast<Wide>(penalty) * urgency;
  return checked_narrow(
      adjusted, INT64_MAX, "EO adjusted fitness");
}

std::uint64_t eo_move_key(const EoMove& move) {
  std::uint64_t key =
      static_cast<std::uint64_t>(move.kind);
  for (const int coordinate :
       {move.first, move.second, move.third}) {
    key = key * UINT64_C(1021) +
          static_cast<std::uint64_t>(coordinate + 1);
  }
  return key;
}

EoMove select_ranked_move(
    std::vector<RankedMove>& candidates,
    const RankSampler& sampler, std::mt19937_64& randomizer,
    EoStatistics& statistics) {
  if (candidates.empty()) {
    throw std::runtime_error("no nonsingular EO move");
  }
  const std::size_t selected_rank =
      sampler.sample(candidates.size(), randomizer);
  auto better = [](const RankedMove& left,
                   const RankedMove& right) {
    if (left.adjusted_fitness != right.adjusted_fitness) {
      return left.adjusted_fitness > right.adjusted_fitness;
    }
    if (left.determinant_magnitude !=
        right.determinant_magnitude) {
      return left.determinant_magnitude >
             right.determinant_magnitude;
    }
    return left.tie_breaker > right.tie_breaker;
  };
  std::nth_element(
      candidates.begin(),
      candidates.begin() +
          static_cast<std::ptrdiff_t>(selected_rank),
      candidates.end(), better);
  statistics.sampled_rank_sum += selected_rank + 1;
  statistics.sampled_rank_max = std::max(
      statistics.sampled_rank_max,
      static_cast<std::uint64_t>(selected_rank + 1));
  return candidates[selected_rank].move;
}

void append_ranked(
    std::vector<RankedMove>& candidates, const EoMove& move,
    const std::array<std::uint64_t, kCoreEntries>& last_moved,
    std::uint64_t iteration, std::uint64_t recency_window,
    std::uint64_t recency_penalty,
    std::uint64_t tie_salt, EoStatistics& statistics) {
  if (move.determinant == 0) {
    ++statistics.singular_directions;
    return;
  }
  candidates.push_back(RankedMove{
      move,
      adjusted_fitness(
          move, last_moved, iteration, recency_window,
          recency_penalty),
      magnitude(move.determinant),
      splitmix64(eo_move_key(move) ^ tie_salt)});
}

EoMove choose_bit_move(
    const State& state,
    std::vector<RankedMove>& candidates,
    const std::array<std::uint64_t, kCoreEntries>& last_moved,
    std::uint64_t iteration, std::uint64_t recency_window,
    std::uint64_t recency_penalty,
    const RankSampler& sampler, std::mt19937_64& randomizer,
    EoStatistics& statistics) {
  candidates.clear();
  const std::uint64_t tie_salt = randomizer();
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      ++statistics.exact_bit_directions;
      const int id = row * kCoreOrder + column;
      append_ranked(
          candidates,
          EoMove{
              EoMoveKind::kBit, row, column, -1, id,
              bit_candidate_determinant(state, row, column)},
          last_moved, iteration, recency_window,
          recency_penalty, tie_salt, statistics);
    }
  }
  return select_ranked_move(
      candidates, sampler, randomizer, statistics);
}

EoMove choose_complement_move(
    const State& state,
    std::vector<RankedMove>& candidates,
    const std::array<std::uint64_t, kCoreEntries>& last_moved,
    std::uint64_t iteration, std::uint64_t recency_window,
    std::uint64_t recency_penalty,
    const RankSampler& sampler, std::mt19937_64& randomizer,
    EoStatistics& statistics) {
  candidates.clear();
  const std::uint64_t tie_salt = randomizer();
  for (int row = 0; row < kCoreOrder; ++row) {
    ++statistics.exact_complement_directions;
    append_ranked(
        candidates,
        EoMove{
            EoMoveKind::kRowComplement, row, -1, -1, -1,
            row_complement_determinant(state, row)},
        last_moved, iteration, recency_window,
        recency_penalty, tie_salt, statistics);
  }
  for (int column = 0; column < kCoreOrder; ++column) {
    ++statistics.exact_complement_directions;
    append_ranked(
        candidates,
        EoMove{
            EoMoveKind::kColumnComplement, column, -1, -1,
            -1, column_complement_determinant(state, column)},
        last_moved, iteration, recency_window,
        recency_penalty, tie_salt, statistics);
  }
  ++statistics.exact_complement_directions;
  append_ranked(
      candidates,
      EoMove{
          EoMoveKind::kWholeComplement, -1, -1, -1, -1,
          whole_complement_determinant(state)},
      last_moved, iteration, recency_window,
      recency_penalty, tie_salt, statistics);
  return select_ranked_move(
      candidates, sampler, randomizer, statistics);
}

std::int64_t row_swap_determinant(
    const State& state, int row, int first_column,
    int second_column) {
  Wide result = state.determinant;
  const std::int64_t first_delta =
      state.core[row][first_column] == 0U ? 1 : -1;
  const std::int64_t second_delta =
      state.core[row][second_column] == 0U ? 1 : -1;
  result +=
      static_cast<Wide>(first_delta) *
          state.adjugate[first_column][row] +
      static_cast<Wide>(second_delta) *
          state.adjugate[second_column][row];
  return checked_narrow(
      result, kDeterminantBound, "row-swap determinant");
}

std::int64_t column_swap_determinant(
    const State& state, int column, int first_row,
    int second_row) {
  Wide result = state.determinant;
  const std::int64_t first_delta =
      state.core[first_row][column] == 0U ? 1 : -1;
  const std::int64_t second_delta =
      state.core[second_row][column] == 0U ? 1 : -1;
  result +=
      static_cast<Wide>(first_delta) *
          state.adjugate[column][first_row] +
      static_cast<Wide>(second_delta) *
          state.adjugate[column][second_row];
  return checked_narrow(
      result, kDeterminantBound, "column-swap determinant");
}

EoMove choose_balanced_pair_move(
    const State& state,
    std::vector<RankedMove>& candidates,
    const std::array<std::uint64_t, kCoreEntries>& last_moved,
    std::uint64_t iteration, std::uint64_t recency_window,
    std::uint64_t recency_penalty,
    const RankSampler& sampler, std::mt19937_64& randomizer,
    EoStatistics& statistics) {
  candidates.clear();
  const std::uint64_t tie_salt = randomizer();
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int first = 0; first < kCoreOrder; ++first) {
      for (int second = first + 1; second < kCoreOrder; ++second) {
        if (state.core[row][first] == state.core[row][second]) {
          continue;
        }
        ++statistics.exact_pair_directions;
        append_ranked(
            candidates,
            EoMove{
                EoMoveKind::kRowSwap, row, first, second, -1,
                row_swap_determinant(state, row, first, second)},
            last_moved, iteration, recency_window,
            recency_penalty, tie_salt, statistics);
      }
    }
  }
  for (int column = 0; column < kCoreOrder; ++column) {
    for (int first = 0; first < kCoreOrder; ++first) {
      for (int second = first + 1; second < kCoreOrder; ++second) {
        if (state.core[first][column] ==
            state.core[second][column]) {
          continue;
        }
        ++statistics.exact_pair_directions;
        append_ranked(
            candidates,
            EoMove{
                EoMoveKind::kColumnSwap, column, first, second,
                -1,
                column_swap_determinant(
                    state, column, first, second)},
            last_moved, iteration, recency_window,
            recency_penalty, tie_salt, statistics);
      }
    }
  }
  return select_ranked_move(
      candidates, sampler, randomizer, statistics);
}

void mark_recent_bits(
    const EoMove& move, std::uint64_t iteration,
    std::array<std::uint64_t, kCoreEntries>& last_moved,
    EoStatistics& statistics) {
  auto mark = [&](int id) {
    std::uint64_t& last =
        last_moved[static_cast<std::size_t>(id)];
    if (last != kNeverMoved && iteration - last <= 2) {
      ++statistics.immediate_refips;
    }
    last = iteration;
  };
  if (move.kind == EoMoveKind::kBit) {
    mark(move.recency_id);
  } else if (move.kind == EoMoveKind::kRowSwap) {
    mark(move.first * kCoreOrder + move.second);
    mark(move.first * kCoreOrder + move.third);
  } else if (move.kind == EoMoveKind::kColumnSwap) {
    mark(move.second * kCoreOrder + move.first);
    mark(move.third * kCoreOrder + move.first);
  }
}

void apply_eo_move(
    State& state, const EoMove& move, std::uint64_t& hash,
    const std::array<std::uint64_t, kCoreEntries>& zobrist,
    EoStatistics& statistics) {
  Vector u{};
  Vector v{};
  if (move.kind == EoMoveKind::kBit) {
    const int row = move.first;
    const int column = move.second;
    u[row] = state.core[row][column] == 0U ? 1 : -1;
    v[column] = 1;
    apply_rank_one(state, u, v, move.determinant);
    state.core[row][column] ^= 1U;
    hash ^= zobrist[row * kCoreOrder + column];
    ++statistics.bit_moves;
    return;
  }
  if (move.kind == EoMoveKind::kRowSwap) {
    const int row = move.first;
    u[row] = 1;
    for (const int column : {move.second, move.third}) {
      v[column] =
          state.core[row][column] == 0U ? 1 : -1;
    }
    apply_rank_one(state, u, v, move.determinant);
    for (const int column : {move.second, move.third}) {
      state.core[row][column] ^= 1U;
      hash ^= zobrist[row * kCoreOrder + column];
    }
    ++statistics.row_swaps;
    return;
  }
  if (move.kind == EoMoveKind::kColumnSwap) {
    const int column = move.first;
    v[column] = 1;
    for (const int row : {move.second, move.third}) {
      u[row] =
          state.core[row][column] == 0U ? 1 : -1;
    }
    apply_rank_one(state, u, v, move.determinant);
    for (const int row : {move.second, move.third}) {
      state.core[row][column] ^= 1U;
      hash ^= zobrist[row * kCoreOrder + column];
    }
    ++statistics.column_swaps;
    return;
  }
  if (move.kind == EoMoveKind::kRowComplement) {
    const int row = move.first;
    u[row] = 1;
    for (int column = 0; column < kCoreOrder; ++column) {
      v[column] =
          state.core[row][column] == 0U ? 1 : -1;
    }
    apply_rank_one(state, u, v, move.determinant);
    for (int column = 0; column < kCoreOrder; ++column) {
      state.core[row][column] ^= 1U;
      hash ^= zobrist[row * kCoreOrder + column];
    }
    ++statistics.row_complements;
    return;
  }
  if (move.kind == EoMoveKind::kColumnComplement) {
    const int column = move.first;
    v[column] = 1;
    for (int row = 0; row < kCoreOrder; ++row) {
      u[row] =
          state.core[row][column] == 0U ? 1 : -1;
    }
    apply_rank_one(state, u, v, move.determinant);
    for (int row = 0; row < kCoreOrder; ++row) {
      state.core[row][column] ^= 1U;
      hash ^= zobrist[row * kCoreOrder + column];
    }
    ++statistics.column_complements;
    return;
  }

  u.fill(-1);
  v.fill(1);
  apply_rank_one(state, u, v, move.determinant);
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      state.core[row][column] ^= 1U;
      hash ^= zobrist[row * kCoreOrder + column];
      state.adjugate[row][column] =
          -state.adjugate[row][column];
    }
  }
  ++statistics.whole_complements;
}

bool apply_eo_kick(
    State& state, const CoreMatrix& best_core,
    std::uint64_t kick_index, std::uint64_t& hash,
    const std::array<std::uint64_t, kCoreEntries>& zobrist,
    const EoArguments& arguments, std::mt19937_64& randomizer,
    EoStatistics& statistics) {
  std::array<int, kCoreEntries> entries{};
  std::iota(entries.begin(), entries.end(), 0);
  std::uniform_int_distribution<int> flip_count(
      arguments.kick_min_flips, arguments.kick_max_flips);
  std::optional<CoreMatrix> retained;
  std::uint64_t retained_magnitude = 0;
  const CoreMatrix base =
      (kick_index & UINT64_C(1)) == 0 ? best_core : state.core;
  for (int attempt = 0; attempt < arguments.kick_attempts;
       ++attempt) {
    ++statistics.kick_attempts;
    std::shuffle(entries.begin(), entries.end(), randomizer);
    CoreMatrix candidate = base;
    const int flips = flip_count(randomizer);
    for (int index = 0; index < flips; ++index) {
      const int flat = entries[static_cast<std::size_t>(index)];
      candidate[flat / kCoreOrder][flat % kCoreOrder] ^= 1U;
    }
    const std::int64_t determinant =
        exact_core_determinant(candidate);
    if (determinant == 0) {
      ++statistics.kick_singular;
      continue;
    }
    const std::uint64_t candidate_magnitude =
        magnitude(determinant);
    if (!retained.has_value() ||
        candidate_magnitude > retained_magnitude) {
      retained = candidate;
      retained_magnitude = candidate_magnitude;
    }
  }
  if (!retained.has_value()) {
    return false;
  }
  state.core = *retained;
  state.determinant = exact_core_determinant(state.core);
  state.adjugate = exact_adjugate(state.core);
  check_adjugate_identity(state);
  hash = core_hash(state.core, zobrist);
  ++statistics.exact_kicks;
  return true;
}

std::uint64_t core_hamming(
    const CoreMatrix& first, const CoreMatrix& second) {
  std::uint64_t distance = 0;
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      distance +=
          first[row][column] == second[row][column] ? 0U : 1U;
    }
  }
  return distance;
}

void eo_log_record(
    std::ofstream& log, const EoArguments& arguments,
    const EoStatistics& statistics, std::string_view event,
    double elapsed, const State& state,
    std::uint64_t best_magnitude,
    std::uint64_t current_penalty) {
  const std::uint64_t exact_directions =
      statistics.exact_bit_directions +
      statistics.exact_complement_directions +
      statistics.exact_pair_directions;
  const long double average_rank =
      statistics.iterations == 0
          ? 0.0L
          : static_cast<long double>(
                statistics.sampled_rank_sum) /
                static_cast<long double>(statistics.iterations);
  log << "{\"best_absolute_determinant\":\""
      << wide_to_string(
             static_cast<Wide>(best_magnitude) * kScale)
      << "\",\"best_core_quotient\":" << best_magnitude
      << ",\"best_strictly_below\":"
      << statistics.best_strictly_below
      << ",\"bit_moves\":" << statistics.bit_moves
      << ",\"column_complements\":"
      << statistics.column_complements
      << ",\"column_swaps\":" << statistics.column_swaps
      << ",\"current_recency_penalty\":" << current_penalty
      << ",\"distinct_tie_checkpoints\":"
      << statistics.distinct_tie_checkpoints
      << ",\"downhill_moves\":" << statistics.downhill_moves
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(6) << elapsed
      << ",\"engine\":\"exact-core-extremal-optimization-v1\""
      << ",\"equal_best_hits\":" << statistics.equal_best_hits
      << ",\"event\":\"" << json_escape(event) << "\""
      << ",\"exact_bit_directions\":"
      << statistics.exact_bit_directions
      << ",\"exact_complement_directions\":"
      << statistics.exact_complement_directions
      << ",\"exact_directions\":" << exact_directions
      << ",\"exact_kicks\":" << statistics.exact_kicks
      << ",\"exact_pair_directions\":"
      << statistics.exact_pair_directions
      << ",\"frontier\":\"" << arguments.frontier << "\""
      << ",\"identity_checks\":" << statistics.identity_checks
      << ",\"immediate_refips\":"
      << statistics.immediate_refips
      << ",\"iterations\":" << statistics.iterations
      << ",\"kick_attempts\":" << statistics.kick_attempts
      << ",\"level_moves\":" << statistics.level_moves
      << ",\"mean_sampled_rank\":" << std::setprecision(6)
      << static_cast<double>(average_rank)
      << ",\"output_path\":\""
      << json_escape(arguments.output.string()) << "\""
      << ",\"penalty_decays\":" << statistics.penalty_decays
      << ",\"penalty_increases\":"
      << statistics.penalty_increases
      << ",\"potential_cycles\":"
      << statistics.potential_cycles
      << ",\"promotions\":" << statistics.promotions
      << ",\"row_complements\":" << statistics.row_complements
      << ",\"row_swaps\":" << statistics.row_swaps
      << ",\"sampled_hamming_max\":"
      << statistics.sampled_hamming_max
      << ",\"sampled_rank_max\":"
      << statistics.sampled_rank_max
      << ",\"seed\":" << arguments.seed
      << ",\"singular_directions\":"
      << statistics.singular_directions
      << ",\"tau\":" << arguments.tau
      << ",\"uphill_moves\":" << statistics.uphill_moves
      << ",\"whole_complements\":"
      << statistics.whole_complements
      << ",\"working_core_determinant\":"
      << state.determinant << "}\n";
  log.flush();
  if (!log) {
    throw std::runtime_error("failed writing EO research log");
  }
}

std::string eo_summary_text(
    const EoArguments& arguments,
    const EoStatistics& statistics, double elapsed,
    const State& state, std::uint64_t best_magnitude,
    std::uint64_t current_penalty,
    std::string_view finish_reason) {
  const std::uint64_t exact_directions =
      statistics.exact_bit_directions +
      statistics.exact_complement_directions +
      statistics.exact_pair_directions;
  const long double average_rank =
      statistics.iterations == 0
          ? 0.0L
          : static_cast<long double>(
                statistics.sampled_rank_sum) /
                static_cast<long double>(statistics.iterations);
  std::ostringstream output;
  output << "{\n"
         << "  \"schema_version\": 1,\n"
         << "  \"engine\": "
            "\"exact-core-extremal-optimization-v1\",\n"
         << "  \"start\": \""
         << json_escape(arguments.start.string()) << "\",\n"
         << "  \"output\": \""
         << json_escape(arguments.output.string()) << "\",\n"
         << "  \"seed\": " << arguments.seed << ",\n"
         << "  \"seconds_limit\": " << arguments.seconds << ",\n"
         << "  \"elapsed_seconds\": " << std::fixed
         << std::setprecision(6) << elapsed << ",\n"
         << "  \"finish_reason\": \""
         << json_escape(finish_reason) << "\",\n"
         << "  \"tau\": " << arguments.tau << ",\n"
         << "  \"recency_window\": "
         << arguments.recency_window << ",\n"
         << "  \"recency_penalty_initial\": "
         << arguments.recency_penalty << ",\n"
         << "  \"recency_penalty_final\": "
         << current_penalty << ",\n"
         << "  \"recency_penalty_max\": "
         << arguments.recency_penalty_max << ",\n"
         << "  \"complement_interval\": "
         << arguments.complement_interval << ",\n"
         << "  \"pair_interval\": "
         << arguments.pair_interval << ",\n"
         << "  \"kick_interval\": "
         << arguments.kick_interval << ",\n"
         << "  \"kick_flip_range\": ["
         << arguments.kick_min_flips << ", "
         << arguments.kick_max_flips << "],\n"
         << "  \"frontier\": \"" << arguments.frontier << "\",\n"
         << "  \"best_absolute_determinant\": \""
         << wide_to_string(
                static_cast<Wide>(best_magnitude) * kScale)
         << "\",\n"
         << "  \"best_core_quotient\": "
         << best_magnitude << ",\n"
         << "  \"working_core_determinant\": "
         << state.determinant << ",\n"
         << "  \"strict_frontier_improvement\": "
         << (static_cast<Wide>(best_magnitude) * kScale >
                     static_cast<Wide>(arguments.frontier)
                 ? "true"
                 : "false")
         << ",\n"
         << "  \"iterations\": " << statistics.iterations
         << ",\n"
         << "  \"exact_directions\": " << exact_directions
         << ",\n"
         << "  \"exact_bit_directions\": "
         << statistics.exact_bit_directions << ",\n"
         << "  \"exact_complement_directions\": "
         << statistics.exact_complement_directions << ",\n"
         << "  \"exact_pair_directions\": "
         << statistics.exact_pair_directions << ",\n"
         << "  \"singular_directions\": "
         << statistics.singular_directions << ",\n"
         << "  \"moves\": {\"bit\": "
         << statistics.bit_moves << ", \"row_complement\": "
         << statistics.row_complements
         << ", \"column_complement\": "
         << statistics.column_complements
         << ", \"whole_complement\": "
         << statistics.whole_complements
         << ", \"row_swap\": " << statistics.row_swaps
         << ", \"column_swap\": " << statistics.column_swaps
         << "},\n"
         << "  \"move_direction\": {\"uphill\": "
         << statistics.uphill_moves << ", \"downhill\": "
         << statistics.downhill_moves << ", \"level\": "
         << statistics.level_moves << "},\n"
         << "  \"mean_sampled_rank\": "
         << std::setprecision(6)
         << static_cast<double>(average_rank) << ",\n"
         << "  \"sampled_rank_max\": "
         << statistics.sampled_rank_max << ",\n"
         << "  \"potential_cycles\": "
         << statistics.potential_cycles << ",\n"
         << "  \"immediate_refips\": "
         << statistics.immediate_refips << ",\n"
         << "  \"penalty_increases\": "
         << statistics.penalty_increases << ",\n"
         << "  \"penalty_decays\": "
         << statistics.penalty_decays << ",\n"
         << "  \"exact_kicks\": " << statistics.exact_kicks
         << ",\n"
         << "  \"kick_attempts\": "
         << statistics.kick_attempts << ",\n"
         << "  \"kick_singular\": "
         << statistics.kick_singular << ",\n"
         << "  \"promotions\": " << statistics.promotions
         << ",\n"
         << "  \"equal_best_hits\": "
         << statistics.equal_best_hits << ",\n"
         << "  \"distinct_tie_checkpoints\": "
         << statistics.distinct_tie_checkpoints << ",\n"
         << "  \"best_strictly_below\": "
         << statistics.best_strictly_below << ",\n"
         << "  \"sampled_hamming_max\": "
         << statistics.sampled_hamming_max << ",\n"
         << "  \"identity_checks\": "
         << statistics.identity_checks << ",\n"
         << "  \"determinant_checks\": "
         << statistics.determinant_checks << "\n"
         << "}\n";
  return output.str();
}

CoreMatrix transformed_core(
    const CoreMatrix& input, const EoMove& move) {
  CoreMatrix result = input;
  if (move.kind == EoMoveKind::kBit) {
    result[move.first][move.second] ^= 1U;
  } else if (move.kind == EoMoveKind::kRowSwap) {
    result[move.first][move.second] ^= 1U;
    result[move.first][move.third] ^= 1U;
  } else if (move.kind == EoMoveKind::kColumnSwap) {
    result[move.second][move.first] ^= 1U;
    result[move.third][move.first] ^= 1U;
  } else if (move.kind == EoMoveKind::kRowComplement) {
    for (int column = 0; column < kCoreOrder; ++column) {
      result[move.first][column] ^= 1U;
    }
  } else if (move.kind == EoMoveKind::kColumnComplement) {
    for (int row = 0; row < kCoreOrder; ++row) {
      result[row][move.first] ^= 1U;
    }
  } else {
    for (auto& row : result) {
      for (std::uint8_t& entry : row) entry ^= 1U;
    }
  }
  return result;
}

void audit_candidate(
    const State& state, const EoMove& move,
    std::string_view label) {
  const CoreMatrix transformed =
      transformed_core(state.core, move);
  const std::int64_t exact =
      exact_core_determinant(transformed);
  if (exact != move.determinant) {
    throw std::runtime_error(
        std::string("self-test candidate mismatch for ") +
        std::string(label) + ": predicted " +
        std::to_string(move.determinant) + ", exact " +
        std::to_string(exact));
  }
}

int run_eo_self_test(
    const EoArguments& arguments, State state) {
  std::mt19937_64 randomizer(arguments.seed);
  const auto zobrist = make_zobrist(arguments.seed);
  std::uint64_t hash = core_hash(state.core, zobrist);
  EoStatistics statistics;

  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      audit_candidate(
          state,
          EoMove{
              EoMoveKind::kBit, row, column, -1,
              row * kCoreOrder + column,
              bit_candidate_determinant(state, row, column)},
          "all initial bit directions");
    }
  }
  for (int row = 0; row < kCoreOrder; ++row) {
    audit_candidate(
        state,
        EoMove{
            EoMoveKind::kRowComplement, row, -1, -1, -1,
            row_complement_determinant(state, row)},
        "all initial row complements");
  }
  for (int column = 0; column < kCoreOrder; ++column) {
    audit_candidate(
        state,
        EoMove{
            EoMoveKind::kColumnComplement, column, -1, -1, -1,
            column_complement_determinant(state, column)},
        "all initial column complements");
  }
  audit_candidate(
      state,
      EoMove{
          EoMoveKind::kWholeComplement, -1, -1, -1, -1,
          whole_complement_determinant(state)},
      "initial whole complement");

  for (int round = 0; round < arguments.self_test_rounds;
       ++round) {
    const int row =
        static_cast<int>(randomizer() % kCoreOrder);
    const int column =
        static_cast<int>(randomizer() % kCoreOrder);
    int other_column =
        static_cast<int>(randomizer() % (kCoreOrder - 1));
    if (other_column >= column) ++other_column;
    int other_row =
        static_cast<int>(randomizer() % (kCoreOrder - 1));
    if (other_row >= row) ++other_row;

    const std::array<EoMove, 6> probes{
        EoMove{
            EoMoveKind::kBit, row, column, -1,
            row * kCoreOrder + column,
            bit_candidate_determinant(state, row, column)},
        EoMove{
            EoMoveKind::kRowComplement, row, -1, -1, -1,
            row_complement_determinant(state, row)},
        EoMove{
            EoMoveKind::kColumnComplement, column, -1, -1,
            -1, column_complement_determinant(state, column)},
        EoMove{
            EoMoveKind::kWholeComplement, -1, -1, -1, -1,
            whole_complement_determinant(state)},
        EoMove{
            EoMoveKind::kRowSwap, row, column, other_column,
            -1,
            row_swap_determinant(
                state, row, column, other_column)},
        EoMove{
            EoMoveKind::kColumnSwap, column, row, other_row,
            -1,
            column_swap_determinant(
                state, column, row, other_row)}};
    for (const EoMove& probe : probes) {
      audit_candidate(state, probe, "random mixed direction");
    }

    std::vector<EoMove> nonsingular;
    for (const EoMove& probe : probes) {
      if (probe.determinant != 0) nonsingular.push_back(probe);
    }
    if (nonsingular.empty()) {
      throw std::runtime_error(
          "self-test found no nonsingular probe");
    }
    const EoMove chosen =
        nonsingular[static_cast<std::size_t>(
            round) % nonsingular.size()];
    apply_eo_move(
        state, chosen, hash, zobrist, statistics);
    const std::int64_t exact =
        exact_core_determinant(state.core);
    if (exact != state.determinant) {
      throw std::runtime_error(
          "self-test incremental determinant mismatch");
    }
    check_adjugate_identity(state);
    if ((round & 3) == 0 &&
        exact_adjugate(state.core) != state.adjugate) {
      throw std::runtime_error(
          "self-test exact adjugate mismatch");
    }
    if (hash != core_hash(state.core, zobrist)) {
      throw std::runtime_error(
          "self-test incremental hash mismatch");
    }
  }

  RankSampler sampler(kCoreEntries, arguments.tau);
  for (int trial = 0; trial < 10000; ++trial) {
    const std::size_t rank =
        sampler.sample(kCoreEntries, randomizer);
    if (rank >= static_cast<std::size_t>(kCoreEntries)) {
      throw std::runtime_error("rank sampler escaped pool");
    }
  }
  std::cout
      << "SELF_TEST_OK rounds=" << arguments.self_test_rounds
      << " initial_directions=" << kMoveCount
      << " random_direction_checks="
      << static_cast<std::uint64_t>(
             arguments.self_test_rounds) *
             UINT64_C(6)
      << " applied_updates=" << arguments.self_test_rounds
      << '\n';
  return 0;
}

[[maybe_unused]] int run_eo_search(
    const EoArguments& arguments, State state) {
  if (!arguments.log.parent_path().empty()) {
    fs::create_directories(arguments.log.parent_path());
  }
  if (!arguments.summary.parent_path().empty()) {
    fs::create_directories(arguments.summary.parent_path());
  }
  std::ofstream log(arguments.log, std::ios::out | std::ios::trunc);
  if (!log) {
    throw std::runtime_error("cannot create EO research log");
  }

  CoreMatrix best_core = state.core;
  std::uint64_t best_magnitude = magnitude(state.determinant);
  std::uint64_t checkpoint_nonce = 0;
  atomic_write(
      arguments.output,
      sign_matrix_text(core_to_sign(best_core)),
      checkpoint_nonce++);

  std::mt19937_64 randomizer(arguments.seed);
  const auto zobrist = make_zobrist(arguments.seed);
  std::uint64_t hash = core_hash(state.core, zobrist);
  std::vector<Visit> visits(kVisitTableSize);
  visits[hash & (kVisitTableSize - 1U)] =
      Visit{hash, 0, true};
  std::array<std::uint64_t, kCoreEntries> last_moved{};
  last_moved.fill(kNeverMoved);
  std::uint64_t current_penalty = arguments.recency_penalty;
  std::uint64_t last_cycle_iteration = 0;
  std::uint64_t last_kick_iteration = 0;
  std::uint64_t kick_index = 0;
  EoStatistics statistics;

  RankSampler bit_sampler(kCoreEntries, arguments.tau);
  RankSampler complement_sampler(kComplementMoveCount, arguments.tau);
  RankSampler pair_sampler(kMaximumPairCandidates, arguments.tau);
  std::vector<RankedMove> ranked_moves;
  ranked_moves.reserve(kMaximumPairCandidates);

  const auto started = Clock::now();
  const auto deadline =
      started + std::chrono::duration<double>(arguments.seconds);
  auto next_heartbeat =
      started +
      std::chrono::duration<double>(
          arguments.heartbeat_seconds);
  std::signal(SIGINT, request_stop);
  std::signal(SIGTERM, request_stop);

  auto elapsed_now = [&]() {
    return std::chrono::duration<double>(
               Clock::now() - started)
        .count();
  };
  auto checkpoint_if_best = [&]() {
    const std::uint64_t current_magnitude =
        magnitude(state.determinant);
    if (current_magnitude < best_magnitude) {
      statistics.best_strictly_below = std::max(
          statistics.best_strictly_below, current_magnitude);
      return;
    }
    if (current_magnitude == best_magnitude) {
      ++statistics.equal_best_hits;
      if (state.core != best_core) {
        best_core = state.core;
        ++statistics.distinct_tie_checkpoints;
        atomic_write(
            arguments.output,
            sign_matrix_text(core_to_sign(best_core)),
            checkpoint_nonce++);
        eo_log_record(
            log, arguments, statistics, "equal_best",
            elapsed_now(), state, best_magnitude,
            current_penalty);
      }
      return;
    }

    best_magnitude = current_magnitude;
    best_core = state.core;
    ++statistics.promotions;
    atomic_write(
        arguments.output,
        sign_matrix_text(core_to_sign(best_core)),
        checkpoint_nonce++);
    eo_log_record(
        log, arguments, statistics, "new_best",
        elapsed_now(), state, best_magnitude,
        current_penalty);
    const Wide full_score =
        static_cast<Wide>(best_magnitude) * kScale;
    std::cout
        << (full_score >
                    static_cast<Wide>(arguments.frontier)
                ? "STRICT_FRONTIER_PROMOTION "
                : "new_best ")
        << "|det|=" << wide_to_string(full_score)
        << " quotient=" << best_magnitude
        << " iteration=" << statistics.iterations << '\n'
        << std::flush;
  };

  eo_log_record(
      log, arguments, statistics, "start", 0.0, state,
      best_magnitude, current_penalty);

  while (!stop_requested && Clock::now() < deadline &&
         (arguments.max_iterations == 0 ||
          statistics.iterations < arguments.max_iterations)) {
    ++statistics.iterations;
    const std::uint64_t old_magnitude =
        magnitude(state.determinant);

    EoMove move;
    if (arguments.pair_interval != 0 &&
        statistics.iterations % arguments.pair_interval == 0) {
      move = choose_balanced_pair_move(
          state, ranked_moves, last_moved, statistics.iterations,
          arguments.recency_window, current_penalty,
          pair_sampler, randomizer, statistics);
    } else if (
        arguments.complement_interval != 0 &&
        statistics.iterations %
                arguments.complement_interval ==
            0) {
      move = choose_complement_move(
          state, ranked_moves, last_moved, statistics.iterations,
          arguments.recency_window, current_penalty,
          complement_sampler, randomizer, statistics);
    } else {
      move = choose_bit_move(
          state, ranked_moves, last_moved, statistics.iterations,
          arguments.recency_window, current_penalty,
          bit_sampler, randomizer, statistics);
    }
    mark_recent_bits(
        move, statistics.iterations, last_moved, statistics);
    apply_eo_move(state, move, hash, zobrist, statistics);

    const std::uint64_t new_magnitude =
        magnitude(state.determinant);
    if (new_magnitude > old_magnitude) {
      ++statistics.uphill_moves;
    } else if (new_magnitude < old_magnitude) {
      ++statistics.downhill_moves;
    } else {
      ++statistics.level_moves;
    }

    Visit& visit =
        visits[hash & (kVisitTableSize - 1U)];
    if (visit.occupied && visit.hash == hash &&
        statistics.iterations > visit.iteration) {
      const std::uint64_t cycle_length =
          statistics.iterations - visit.iteration;
      if (cycle_length <=
          4 * arguments.recency_window) {
        ++statistics.potential_cycles;
        last_cycle_iteration = statistics.iterations;
        const std::uint64_t increase =
            std::max<std::uint64_t>(
                UINT64_C(1), current_penalty / UINT64_C(8));
        current_penalty = std::min(
            arguments.recency_penalty_max,
            current_penalty + increase);
        ++statistics.penalty_increases;
      }
    }
    visit = Visit{hash, statistics.iterations, true};
    if (statistics.iterations - last_cycle_iteration >=
            UINT64_C(512) &&
        (statistics.iterations & UINT64_C(127)) == 0 &&
        current_penalty > arguments.recency_penalty) {
      current_penalty = std::max(
          arguments.recency_penalty,
          current_penalty -
              std::max<std::uint64_t>(
                  UINT64_C(1),
                  current_penalty / UINT64_C(16)));
      ++statistics.penalty_decays;
    }

    checkpoint_if_best();

    if ((statistics.iterations & UINT64_C(127)) == 0) {
      statistics.sampled_hamming_max = std::max(
          statistics.sampled_hamming_max,
          core_hamming(state.core, best_core));
    }

    if (arguments.kick_interval != 0 &&
        statistics.iterations - last_kick_iteration >=
            arguments.kick_interval) {
      ++kick_index;
      const bool kicked = apply_eo_kick(
          state, best_core, kick_index, hash, zobrist,
          arguments, randomizer, statistics);
      last_moved.fill(kNeverMoved);
      std::fill(visits.begin(), visits.end(), Visit{});
      visits[hash & (kVisitTableSize - 1U)] =
          Visit{hash, statistics.iterations, true};
      current_penalty = arguments.recency_penalty;
      last_cycle_iteration = statistics.iterations;
      last_kick_iteration = statistics.iterations;
      eo_log_record(
          log, arguments, statistics,
          kicked ? "exact_kick" : "kick_restart_failed",
          elapsed_now(), state, best_magnitude,
          current_penalty);
      checkpoint_if_best();
    }

    if (statistics.iterations %
            kEoIdentityCheckInterval ==
        0) {
      check_adjugate_identity(state);
      ++statistics.identity_checks;
    }
    if (statistics.iterations %
            kEoDeterminantCheckInterval ==
        0) {
      if (exact_core_determinant(state.core) !=
          state.determinant) {
        throw std::runtime_error(
            "incremental EO determinant invariant failed");
      }
      ++statistics.determinant_checks;
    }

    const auto now = Clock::now();
    if (arguments.heartbeat_seconds > 0.0 &&
        now >= next_heartbeat) {
      eo_log_record(
          log, arguments, statistics, "heartbeat",
          elapsed_now(), state, best_magnitude,
          current_penalty);
      next_heartbeat =
          now +
          std::chrono::duration<double>(
              arguments.heartbeat_seconds);
    }
  }

  check_adjugate_identity(state);
  ++statistics.identity_checks;
  if (exact_core_determinant(state.core) !=
      state.determinant) {
    throw std::runtime_error(
        "final EO determinant invariant failed");
  }
  ++statistics.determinant_checks;
  atomic_write(
      arguments.output,
      sign_matrix_text(core_to_sign(best_core)),
      checkpoint_nonce++);

  const std::string finish_reason =
      stop_requested
          ? "signal"
          : (arguments.max_iterations != 0 &&
                     statistics.iterations >=
                         arguments.max_iterations
                 ? "iteration_limit"
                 : "time_limit");
  const double elapsed = elapsed_now();
  eo_log_record(
      log, arguments, statistics, "finished", elapsed,
      state, best_magnitude, current_penalty);
  atomic_write(
      arguments.summary,
      eo_summary_text(
          arguments, statistics, elapsed, state,
          best_magnitude, current_penalty, finish_reason),
      checkpoint_nonce++);

  const std::uint64_t exact_directions =
      statistics.exact_bit_directions +
      statistics.exact_complement_directions +
      statistics.exact_pair_directions;
  std::cout
      << "finished |det|="
      << wide_to_string(
             static_cast<Wide>(best_magnitude) * kScale)
      << " quotient=" << best_magnitude
      << " iterations=" << statistics.iterations
      << " exact_directions=" << exact_directions
      << " downhill_moves=" << statistics.downhill_moves
      << " kicks=" << statistics.exact_kicks << '\n';
  return 0;
}

}  // namespace

#ifndef CORE_EXTREMAL_OPTIMIZATION_NO_MAIN
int main(int argc, char** argv) {
  try {
    const EoArguments arguments =
        parse_eo_arguments(argc, argv);
    State state = eo_initial_state(arguments.start);
    if (arguments.self_test_rounds != 0) {
      return run_eo_self_test(arguments, std::move(state));
    }
    return run_eo_search(arguments, std::move(state));
  } catch (const std::exception& error) {
    std::cerr
        << "core_extremal_optimization: " << error.what()
        << '\n';
    return 2;
  }
}
#endif
