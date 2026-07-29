// Exact replica-exchange / parallel-tempering pilot for order-23 maximal
// determinant search.
//
// The objective distribution at inverse temperature beta is
//
//                     pi_beta(B) proportional to |det(B)|^beta
//
// over nonsingular dephased 22x22 binary cores.  Local Metropolis decisions
// therefore use the ratio of two exact integer determinants.  Adjacent
// replica swaps use
//
//   log alpha = (beta_i-beta_j) (log |det(B_j)|-log |det(B_i)|),
//
// which is the detailed-balance ratio for exchanging the configurations at
// fixed temperatures.
//
// The symmetric local kernel mixes single-bit flips, balanced two-bit flips
// in one row or column, row/column complements, and whole-core complements.
// All candidate determinants and accepted adjugate updates are exact.  A
// short acceptance-feedback pass contracts beta gaps at swap bottlenecks and
// expands them where exchange is easy.  This is inspired by feedback-
// optimized parallel tempering, but is not a full equilibrium/FOPT study.

#include <optional>

#define main core_adjugate_tabu_embedded_main
#include "core_adjugate_tabu.cpp"
#undef main

namespace {

constexpr std::uint64_t kPTFrontier =
    UINT64_C(2779447296000000);
constexpr std::uint64_t kPTAuditPrime = UINT64_C(1000000007);
constexpr int kPTMoveKinds = 6;

struct PTArguments {
  std::vector<fs::path> seeds;
  fs::path output;
  fs::path log;
  fs::path summary;
  std::uint64_t seed = 34001;
  std::uint64_t frontier = kPTFrontier;
  double seconds = 300.0;
  double heartbeat_seconds = 15.0;
  double feedback_seconds = 20.0;
  double cold_downhill_target = 0.02;
  double hot_downhill_target = 0.80;
  int expected_seeds = 3;
  int replica_count = 16;
  int differential_samples = 72;
  std::uint64_t reseed_sweeps = UINT64_C(200000);
  std::uint64_t identity_interval = UINT64_C(16384);
  std::uint64_t exact_check_interval = UINT64_C(65536);
};

enum class PTMoveKind : std::uint8_t {
  kBit = 0,
  kBalancedRowPair = 1,
  kBalancedColumnPair = 2,
  kRowComplement = 3,
  kColumnComplement = 4,
  kWholeComplement = 5,
};

struct PTMove {
  PTMoveKind kind = PTMoveKind::kBit;
  int first = -1;
  int second = -1;
  int third = -1;
  std::int64_t determinant = 0;
  bool valid = false;
};

struct Walker {
  int id = -1;
  int last_endpoint = -1;
  bool armed_from_cold = false;
};

struct Replica {
  State state{};
  std::uint64_t hash = 0;
  Walker walker{};
};

struct TemperatureCalibration {
  std::size_t downhill_samples = 0;
  long double median_log_drop = 0.0L;
  long double beta_cold = 0.0L;
  long double beta_hot = 0.0L;
};

struct PTStatistics {
  std::uint64_t sweeps = 0;
  std::uint64_t local_proposals = 0;
  std::uint64_t local_accepts = 0;
  std::uint64_t local_rejects = 0;
  std::uint64_t singular_rejects = 0;
  std::uint64_t null_proposals = 0;
  std::array<std::uint64_t, kPTMoveKinds> proposed_by_kind{};
  std::array<std::uint64_t, kPTMoveKinds> accepted_by_kind{};
  std::vector<std::uint64_t> swap_attempts;
  std::vector<std::uint64_t> swap_accepts;
  std::vector<std::uint64_t> calibration_swap_attempts;
  std::vector<std::uint64_t> calibration_swap_accepts;
  std::vector<std::uint64_t> post_swap_attempts;
  std::vector<std::uint64_t> post_swap_accepts;
  std::uint64_t endpoint_crossings = 0;
  std::uint64_t completed_round_trips = 0;
  std::uint64_t post_endpoint_crossings = 0;
  std::uint64_t post_completed_round_trips = 0;
  std::uint64_t promotions = 0;
  std::uint64_t elite_reseed_events = 0;
  std::uint64_t elite_reseeded_replicas = 0;
  std::uint64_t walker_tracking_resets = 0;
  std::uint64_t identity_checks = 0;
  std::uint64_t exact_checks = 0;
  std::uint64_t differential_checks = 0;
  bool feedback_applied = false;
  double feedback_elapsed_seconds = 0.0;
  Statistics core_updates{};
};

struct PTLogger {
  std::ofstream stream;

  explicit PTLogger(const fs::path& path) {
    if (fs::exists(path)) {
      throw std::runtime_error(
          "refusing to overwrite log: " + path.string());
    }
    if (!path.parent_path().empty()) {
      fs::create_directories(path.parent_path());
    }
    stream.open(path, std::ios::out | std::ios::trunc);
    if (!stream) {
      throw std::runtime_error(
          "cannot create parallel-tempering log: " +
          path.string());
    }
  }

  void line(const std::string& record) {
    stream << record << '\n';
    stream.flush();
    if (!stream) {
      throw std::runtime_error(
          "failed writing parallel-tempering JSONL log");
    }
  }
};

int pt_positive_int(std::string_view text,
                    std::string_view option, int maximum) {
  const std::uint64_t parsed = strict_unsigned(text, option);
  if (parsed == 0 ||
      parsed > static_cast<std::uint64_t>(maximum)) {
    throw std::runtime_error(
        std::string(option) + " is outside its valid range");
  }
  return static_cast<int>(parsed);
}

double pt_probability(std::string_view text,
                      std::string_view option) {
  const double parsed = strict_double(text, option, false);
  if (parsed >= 1.0) {
    throw std::runtime_error(
        std::string(option) + " must be less than one");
  }
  return parsed;
}

PTArguments parse_pt_arguments(int argc, char** argv) {
  PTArguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    auto value = [&]() -> std::string_view {
      ++index;
      if (index >= argc) {
        throw std::runtime_error("missing value for " + option);
      }
      return argv[index];
    };
    if (option == "--seed-matrix") {
      arguments.seeds.emplace_back(value());
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
    } else if (option == "--feedback-seconds") {
      arguments.feedback_seconds =
          strict_double(value(), option, true);
    } else if (option == "--cold-downhill-target") {
      arguments.cold_downhill_target =
          pt_probability(value(), option);
    } else if (option == "--hot-downhill-target") {
      arguments.hot_downhill_target =
          pt_probability(value(), option);
    } else if (option == "--expected-seeds") {
      arguments.expected_seeds =
          pt_positive_int(value(), option, 64);
    } else if (option == "--replicas") {
      arguments.replica_count =
          pt_positive_int(value(), option, 64);
    } else if (option == "--differential-samples") {
      const std::uint64_t parsed = strict_unsigned(value(), option);
      if (parsed >
          static_cast<std::uint64_t>(
              std::numeric_limits<int>::max())) {
        throw std::runtime_error(
            "--differential-samples is too large");
      }
      arguments.differential_samples = static_cast<int>(parsed);
    } else if (option == "--reseed-sweeps") {
      arguments.reseed_sweeps = strict_unsigned(value(), option);
    } else if (option == "--identity-interval") {
      arguments.identity_interval = strict_unsigned(value(), option);
    } else if (option == "--exact-check-interval") {
      arguments.exact_check_interval =
          strict_unsigned(value(), option);
    } else if (option == "--help") {
      std::cout
          << "usage: core_parallel_tempering "
             "--seed-matrix H0 --seed-matrix H1 "
             "--seed-matrix H2 --output MATRIX --log JSONL "
             "--summary JSON [options]\n"
          << "  --seed N --frontier N --seconds S --heartbeat S\n"
          << "  --replicas N --feedback-seconds S\n"
          << "  --cold-downhill-target P "
             "--hot-downhill-target P\n"
          << "  --differential-samples N "
             "--reseed-sweeps N "
             "--identity-interval N --exact-check-interval N\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }

  if (arguments.output.empty() || arguments.log.empty() ||
      arguments.summary.empty()) {
    throw std::runtime_error(
        "--output, --log, and --summary are required");
  }
  if (static_cast<int>(arguments.seeds.size()) !=
      arguments.expected_seeds) {
    throw std::runtime_error(
        "expected " + std::to_string(arguments.expected_seeds) +
        " --seed-matrix inputs but received " +
        std::to_string(arguments.seeds.size()));
  }
  if (arguments.replica_count < 4) {
    throw std::runtime_error(
        "--replicas must be at least four");
  }
  if (arguments.cold_downhill_target >=
      arguments.hot_downhill_target) {
    throw std::runtime_error(
        "cold downhill target must be below hot downhill target");
  }
  if (arguments.feedback_seconds >= arguments.seconds &&
      arguments.feedback_seconds != 0.0) {
    throw std::runtime_error(
        "--feedback-seconds must be zero or less than --seconds");
  }
  const std::array<fs::path, 3> outputs{
      arguments.output, arguments.log, arguments.summary};
  for (std::size_t first = 0; first < outputs.size(); ++first) {
    for (std::size_t second = first + 1;
         second < outputs.size(); ++second) {
      if (outputs[first] == outputs[second]) {
        throw std::runtime_error(
            "output, log, and summary paths must differ");
      }
    }
  }
  for (const fs::path& seed : arguments.seeds) {
    for (const fs::path& output : outputs) {
      if (seed == output) {
        throw std::runtime_error(
            "seed and output paths must differ");
      }
    }
  }
  return arguments;
}

std::uint64_t pt_modular_power(std::uint64_t base,
                               std::uint64_t exponent) {
  std::uint64_t result = 1;
  while (exponent != 0) {
    if ((exponent & 1U) != 0) {
      result = static_cast<std::uint64_t>(
          (static_cast<__uint128_t>(result) * base) %
          kPTAuditPrime);
    }
    base = static_cast<std::uint64_t>(
        (static_cast<__uint128_t>(base) * base) %
        kPTAuditPrime);
    exponent >>= 1U;
  }
  return result;
}

std::uint64_t pt_core_determinant_mod(const CoreMatrix& core) {
  std::array<std::array<std::uint64_t, kCoreOrder>, kCoreOrder>
      work{};
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      work[row][column] = core[row][column];
    }
  }
  std::uint64_t determinant = 1;
  bool negative = false;
  for (int column = 0; column < kCoreOrder; ++column) {
    int pivot = column;
    while (pivot < kCoreOrder && work[pivot][column] == 0) {
      ++pivot;
    }
    if (pivot == kCoreOrder) {
      return 0;
    }
    if (pivot != column) {
      std::swap(work[pivot], work[column]);
      negative = !negative;
    }
    const std::uint64_t pivot_value = work[column][column];
    determinant = static_cast<std::uint64_t>(
        (static_cast<__uint128_t>(determinant) *
         pivot_value) %
        kPTAuditPrime);
    const std::uint64_t inverse =
        pt_modular_power(pivot_value, kPTAuditPrime - 2);
    for (int row = column + 1; row < kCoreOrder; ++row) {
      if (work[row][column] == 0) {
        continue;
      }
      const std::uint64_t factor =
          static_cast<std::uint64_t>(
              (static_cast<__uint128_t>(
                   work[row][column]) *
               inverse) %
              kPTAuditPrime);
      for (int inner = column; inner < kCoreOrder; ++inner) {
        const std::uint64_t product =
            static_cast<std::uint64_t>(
                (static_cast<__uint128_t>(factor) *
                 work[column][inner]) %
                kPTAuditPrime);
        work[row][inner] =
            work[row][inner] >= product
                ? work[row][inner] - product
                : work[row][inner] + kPTAuditPrime - product;
      }
    }
  }
  if (negative && determinant != 0) {
    determinant = kPTAuditPrime - determinant;
  }
  return determinant;
}

std::uint64_t pt_signed_mod(std::int64_t value) {
  if (value >= 0) {
    return static_cast<std::uint64_t>(value) % kPTAuditPrime;
  }
  const std::uint64_t absolute =
      static_cast<std::uint64_t>(-value) % kPTAuditPrime;
  return absolute == 0 ? 0 : kPTAuditPrime - absolute;
}

std::vector<std::vector<Wide>> pt_sign_as_wide(
    const SignMatrix& matrix) {
  std::vector<std::vector<Wide>> result(
      kSignOrder, std::vector<Wide>(kSignOrder));
  for (int row = 0; row < kSignOrder; ++row) {
    for (int column = 0; column < kSignOrder; ++column) {
      result[row][column] = matrix[row][column];
    }
  }
  return result;
}

void pt_differential_check(
    const State& state, std::string_view context,
    PTStatistics& statistics) {
  const std::int64_t exact =
      exact_core_determinant(state.core);
  if (exact != state.determinant) {
    throw std::runtime_error(
        std::string(context) +
        ": cached and Bareiss determinants differ");
  }
  if (pt_core_determinant_mod(state.core) !=
      pt_signed_mod(state.determinant)) {
    throw std::runtime_error(
        std::string(context) +
        ": modular and Bareiss determinants differ");
  }
  const Wide sign_determinant =
      bareiss(pt_sign_as_wide(core_to_sign(state.core)));
  if (sign_determinant !=
      static_cast<Wide>(state.determinant) * kScale) {
    throw std::runtime_error(
        std::string(context) +
        ": sign/core determinant quotient differs");
  }
  ++statistics.differential_checks;
}

std::size_t pt_kind_index(PTMoveKind kind) {
  return static_cast<std::size_t>(kind);
}

const char* pt_kind_name(std::size_t index) {
  constexpr std::array<const char*, kPTMoveKinds> names{
      "bit", "balanced_row_pair", "balanced_column_pair",
      "row_complement", "column_complement",
      "whole_complement"};
  return names.at(index);
}

PTMove pt_propose_kind(const State& state, PTMoveKind kind,
                       std::mt19937_64& randomizer) {
  std::uniform_int_distribution<int> coordinate(
      0, kCoreOrder - 1);
  if (kind == PTMoveKind::kBit) {
    const int row = coordinate(randomizer);
    const int column = coordinate(randomizer);
    return PTMove{
        kind, row, column, -1,
        bit_candidate_determinant(state, row, column), true};
  }
  if (kind == PTMoveKind::kBalancedRowPair) {
    const int row = coordinate(randomizer);
    std::array<int, kCoreOrder> zeros{};
    std::array<int, kCoreOrder> ones{};
    int zero_count = 0;
    int one_count = 0;
    for (int column = 0; column < kCoreOrder; ++column) {
      if (state.core[row][column] == 0U) {
        zeros[static_cast<std::size_t>(zero_count++)] = column;
      } else {
        ones[static_cast<std::size_t>(one_count++)] = column;
      }
    }
    if (zero_count == 0 || one_count == 0) {
      return PTMove{kind, row, -1, -1, 0, false};
    }
    std::uniform_int_distribution<int> zero_choice(0, zero_count - 1);
    std::uniform_int_distribution<int> one_choice(0, one_count - 1);
    const int zero_column =
        zeros[static_cast<std::size_t>(zero_choice(randomizer))];
    const int one_column =
        ones[static_cast<std::size_t>(one_choice(randomizer))];
    const Wide candidate =
        static_cast<Wide>(state.determinant) +
        state.adjugate[zero_column][row] -
        state.adjugate[one_column][row];
    return PTMove{
        kind, row, zero_column, one_column,
        checked_narrow(
            candidate, kDeterminantBound,
            "balanced-row-pair determinant"),
        true};
  }
  if (kind == PTMoveKind::kBalancedColumnPair) {
    const int column = coordinate(randomizer);
    std::array<int, kCoreOrder> zeros{};
    std::array<int, kCoreOrder> ones{};
    int zero_count = 0;
    int one_count = 0;
    for (int row = 0; row < kCoreOrder; ++row) {
      if (state.core[row][column] == 0U) {
        zeros[static_cast<std::size_t>(zero_count++)] = row;
      } else {
        ones[static_cast<std::size_t>(one_count++)] = row;
      }
    }
    if (zero_count == 0 || one_count == 0) {
      return PTMove{kind, column, -1, -1, 0, false};
    }
    std::uniform_int_distribution<int> zero_choice(0, zero_count - 1);
    std::uniform_int_distribution<int> one_choice(0, one_count - 1);
    const int zero_row =
        zeros[static_cast<std::size_t>(zero_choice(randomizer))];
    const int one_row =
        ones[static_cast<std::size_t>(one_choice(randomizer))];
    const Wide candidate =
        static_cast<Wide>(state.determinant) +
        state.adjugate[column][zero_row] -
        state.adjugate[column][one_row];
    return PTMove{
        kind, column, zero_row, one_row,
        checked_narrow(
            candidate, kDeterminantBound,
            "balanced-column-pair determinant"),
        true};
  }
  if (kind == PTMoveKind::kRowComplement) {
    const int row = coordinate(randomizer);
    return PTMove{
        kind, row, -1, -1,
        row_complement_determinant(state, row), true};
  }
  if (kind == PTMoveKind::kColumnComplement) {
    const int column = coordinate(randomizer);
    return PTMove{
        kind, column, -1, -1,
        column_complement_determinant(state, column), true};
  }
  return PTMove{
      kind, -1, -1, -1,
      whole_complement_determinant(state), true};
}

PTMove pt_propose(const State& state,
                  std::mt19937_64& randomizer) {
  std::uniform_int_distribution<int> mixture(0, 99);
  const int draw = mixture(randomizer);
  PTMoveKind kind = PTMoveKind::kBit;
  if (draw < 70) {
    kind = PTMoveKind::kBit;
  } else if (draw < 80) {
    kind = PTMoveKind::kBalancedRowPair;
  } else if (draw < 90) {
    kind = PTMoveKind::kBalancedColumnPair;
  } else if (draw < 94) {
    kind = PTMoveKind::kRowComplement;
  } else if (draw < 98) {
    kind = PTMoveKind::kColumnComplement;
  } else {
    kind = PTMoveKind::kWholeComplement;
  }
  return pt_propose_kind(state, kind, randomizer);
}

CoreMatrix pt_materialize(const CoreMatrix& source,
                          const PTMove& move) {
  CoreMatrix result = source;
  if (move.kind == PTMoveKind::kBit) {
    result[move.first][move.second] ^= 1U;
  } else if (move.kind == PTMoveKind::kBalancedRowPair) {
    result[move.first][move.second] ^= 1U;
    result[move.first][move.third] ^= 1U;
  } else if (move.kind == PTMoveKind::kBalancedColumnPair) {
    result[move.second][move.first] ^= 1U;
    result[move.third][move.first] ^= 1U;
  } else if (move.kind == PTMoveKind::kRowComplement) {
    for (int column = 0; column < kCoreOrder; ++column) {
      result[move.first][column] ^= 1U;
    }
  } else if (move.kind == PTMoveKind::kColumnComplement) {
    for (int row = 0; row < kCoreOrder; ++row) {
      result[row][move.first] ^= 1U;
    }
  } else {
    for (auto& row : result) {
      for (std::uint8_t& entry : row) {
        entry ^= 1U;
      }
    }
  }
  return result;
}

void pt_apply_pair(Replica& replica, const PTMove& move) {
  Vector u{};
  Vector v{};
  if (move.kind == PTMoveKind::kBalancedRowPair) {
    const int row = move.first;
    const int zero_column = move.second;
    const int one_column = move.third;
    if (replica.state.core[row][zero_column] != 0U ||
        replica.state.core[row][one_column] != 1U) {
      throw std::runtime_error(
          "balanced row-pair orientation invariant failed");
    }
    u[row] = 1;
    v[zero_column] = 1;
    v[one_column] = -1;
    apply_rank_one(
        replica.state, u, v, move.determinant);
    replica.state.core[row][zero_column] ^= 1U;
    replica.state.core[row][one_column] ^= 1U;
    return;
  }
  const int column = move.first;
  const int zero_row = move.second;
  const int one_row = move.third;
  if (replica.state.core[zero_row][column] != 0U ||
      replica.state.core[one_row][column] != 1U) {
    throw std::runtime_error(
        "balanced column-pair orientation invariant failed");
  }
  u[zero_row] = 1;
  u[one_row] = -1;
  v[column] = 1;
  apply_rank_one(
      replica.state, u, v, move.determinant);
  replica.state.core[zero_row][column] ^= 1U;
  replica.state.core[one_row][column] ^= 1U;
}

void pt_apply(Replica& replica, const PTMove& move,
              const std::array<std::uint64_t, kCoreEntries>&
                  zobrist,
              PTStatistics& statistics) {
  if (move.kind == PTMoveKind::kBalancedRowPair ||
      move.kind == PTMoveKind::kBalancedColumnPair) {
    pt_apply_pair(replica, move);
    if (move.kind == PTMoveKind::kBalancedRowPair) {
      replica.hash ^=
          zobrist[move.first * kCoreOrder + move.second];
      replica.hash ^=
          zobrist[move.first * kCoreOrder + move.third];
    } else {
      replica.hash ^=
          zobrist[move.second * kCoreOrder + move.first];
      replica.hash ^=
          zobrist[move.third * kCoreOrder + move.first];
    }
    return;
  }

  Move core_move;
  core_move.determinant = move.determinant;
  core_move.aspiration = false;
  if (move.kind == PTMoveKind::kBit) {
    core_move.kind = MoveKind::kBit;
    core_move.first = move.first;
    core_move.second = move.second;
    core_move.id = move.first * kCoreOrder + move.second;
  } else if (move.kind == PTMoveKind::kRowComplement) {
    core_move.kind = MoveKind::kRowComplement;
    core_move.first = move.first;
    core_move.second = -1;
    core_move.id = kCoreEntries + move.first;
  } else if (move.kind == PTMoveKind::kColumnComplement) {
    core_move.kind = MoveKind::kColumnComplement;
    core_move.first = move.first;
    core_move.second = -1;
    core_move.id =
        kCoreEntries + kCoreOrder + move.first;
  } else {
    core_move.kind = MoveKind::kWholeComplement;
    core_move.first = -1;
    core_move.second = -1;
    core_move.id = kMoveCount - 1;
  }
  apply_move(
      replica.state, core_move, replica.hash, zobrist,
      statistics.core_updates);
}

bool pt_metropolis_accept(
    std::uint64_t current, std::uint64_t candidate,
    long double beta, std::mt19937_64& randomizer) {
  if (candidate >= current) {
    return true;
  }
  const long double log_acceptance =
      beta *
      (std::log(static_cast<long double>(candidate)) -
       std::log(static_cast<long double>(current)));
  std::uniform_real_distribution<long double> uniform(0.0L, 1.0L);
  const long double draw = uniform(randomizer);
  return draw == 0.0L || std::log(draw) < log_acceptance;
}

long double pt_swap_log_ratio(
    std::uint64_t cold_score, std::uint64_t hot_score,
    long double beta_cold, long double beta_hot) {
  return
      (beta_cold - beta_hot) *
      (std::log(static_cast<long double>(hot_score)) -
       std::log(static_cast<long double>(cold_score)));
}

void pt_audit_swap_sign() {
  constexpr std::uint64_t high = UINT64_C(3);
  constexpr std::uint64_t low = UINT64_C(2);
  constexpr long double beta_cold = 5.0L;
  constexpr long double beta_hot = 1.0L;
  const long double unfavorable =
      pt_swap_log_ratio(
          high, low, beta_cold, beta_hot);
  const long double favorable =
      pt_swap_log_ratio(
          low, high, beta_cold, beta_hot);
  if (!(unfavorable < 0.0L) || !(favorable > 0.0L) ||
      std::abs(unfavorable + favorable) >
          std::numeric_limits<long double>::epsilon() * 16.0L) {
    throw std::runtime_error(
        "adjacent-swap detailed-balance sign audit failed");
  }
}

bool pt_swap_accept(
    const Replica& cold, const Replica& hot,
    long double beta_cold, long double beta_hot,
    std::mt19937_64& randomizer) {
  const long double log_acceptance =
      pt_swap_log_ratio(
          magnitude(cold.state.determinant),
          magnitude(hot.state.determinant),
          beta_cold, beta_hot);
  if (log_acceptance >= 0.0L) {
    return true;
  }
  std::uniform_real_distribution<long double> uniform(0.0L, 1.0L);
  const long double draw = uniform(randomizer);
  return draw == 0.0L || std::log(draw) < log_acceptance;
}

std::vector<State> pt_load_seeds(
    const PTArguments& arguments, PTStatistics& statistics) {
  std::vector<State> seeds;
  std::vector<std::string> unique;
  for (const fs::path& path : arguments.seeds) {
    State state;
    state.core = dephase_to_core(read_sign_matrix(path));
    state.determinant =
        exact_core_determinant(state.core);
    if (state.determinant == 0) {
      throw std::runtime_error(
          "seed is singular: " + path.string());
    }
    const Wide score =
        static_cast<Wide>(magnitude(state.determinant)) * kScale;
    if (score != static_cast<Wide>(arguments.frontier)) {
      throw std::runtime_error(
          "seed is not at the configured frontier: " +
          path.string());
    }
    state.adjugate = exact_adjugate(state.core);
    check_adjugate_identity(state);
    pt_differential_check(
        state, path.string(), statistics);
    std::ostringstream key;
    for (const auto& row : state.core) {
      for (const std::uint8_t entry : row) {
        key << static_cast<unsigned>(entry);
      }
    }
    if (std::find(unique.begin(), unique.end(), key.str()) !=
        unique.end()) {
      throw std::runtime_error(
          "duplicate dephased seed core: " + path.string());
    }
    unique.push_back(key.str());
    seeds.push_back(std::move(state));
  }
  return seeds;
}

TemperatureCalibration pt_calibrate_temperatures(
    const std::vector<State>& seeds,
    const PTArguments& arguments) {
  std::vector<long double> drops;
  for (const State& state : seeds) {
    const std::uint64_t current = magnitude(state.determinant);
    for (int row = 0; row < kCoreOrder; ++row) {
      for (int column = 0; column < kCoreOrder; ++column) {
        const std::uint64_t candidate =
            magnitude(bit_candidate_determinant(
                state, row, column));
        if (candidate != 0 && candidate < current) {
          drops.push_back(
              std::log(static_cast<long double>(current)) -
              std::log(static_cast<long double>(candidate)));
        }
      }
    }
    for (int row = 0; row < kCoreOrder; ++row) {
      const std::uint64_t candidate =
          magnitude(row_complement_determinant(state, row));
      if (candidate != 0 && candidate < current) {
        drops.push_back(
            std::log(static_cast<long double>(current)) -
            std::log(static_cast<long double>(candidate)));
      }
    }
    for (int column = 0; column < kCoreOrder; ++column) {
      const std::uint64_t candidate =
          magnitude(column_complement_determinant(
              state, column));
      if (candidate != 0 && candidate < current) {
        drops.push_back(
            std::log(static_cast<long double>(current)) -
            std::log(static_cast<long double>(candidate)));
      }
    }
    const std::uint64_t whole =
        magnitude(whole_complement_determinant(state));
    if (whole != 0 && whole < current) {
      drops.push_back(
          std::log(static_cast<long double>(current)) -
          std::log(static_cast<long double>(whole)));
    }
  }
  if (drops.empty()) {
    throw std::runtime_error(
        "temperature calibration found no downhill moves");
  }
  std::sort(drops.begin(), drops.end());
  const long double median =
      drops[drops.size() / 2];
  if (!(median > 0.0L) || !std::isfinite(median)) {
    throw std::runtime_error(
        "temperature calibration produced invalid median gap");
  }
  const long double beta_cold =
      -std::log(
          static_cast<long double>(
              arguments.cold_downhill_target)) /
      median;
  const long double beta_hot =
      -std::log(
          static_cast<long double>(
              arguments.hot_downhill_target)) /
      median;
  if (!(beta_cold > beta_hot) || !(beta_hot > 0.0L) ||
      !std::isfinite(beta_cold) ||
      !std::isfinite(beta_hot)) {
    throw std::runtime_error(
        "temperature calibration produced invalid endpoints");
  }
  return TemperatureCalibration{
      drops.size(), median, beta_cold, beta_hot};
}

std::vector<long double> pt_geometric_ladder(
    int count, long double beta_cold,
    long double beta_hot) {
  std::vector<long double> result(
      static_cast<std::size_t>(count));
  const long double ratio =
      std::pow(
          beta_hot / beta_cold,
          1.0L / static_cast<long double>(count - 1));
  result.front() = beta_cold;
  for (int index = 1; index < count - 1; ++index) {
    result[static_cast<std::size_t>(index)] =
        result[static_cast<std::size_t>(index - 1)] * ratio;
  }
  result.back() = beta_hot;
  return result;
}

std::vector<long double> pt_feedback_ladder(
    const std::vector<long double>& beta,
    const std::vector<std::uint64_t>& attempts,
    const std::vector<std::uint64_t>& accepts) {
  if (beta.size() < 2 ||
      attempts.size() + 1 != beta.size() ||
      accepts.size() != attempts.size()) {
    throw std::runtime_error(
        "invalid temperature feedback dimensions");
  }
  std::vector<long double> weighted_gaps(attempts.size());
  long double weighted_total = 0.0L;
  for (std::size_t index = 0; index < attempts.size(); ++index) {
    const long double gap = beta[index] - beta[index + 1];
    long double acceptance = 0.02L;
    if (attempts[index] != 0) {
      acceptance =
          static_cast<long double>(accepts[index]) /
          static_cast<long double>(attempts[index]);
      acceptance =
          std::clamp(acceptance, 0.02L, 1.0L);
    }
    // Contract hard intervals and expand easy ones.  Repeating this
    // conservative update lets a phase-boundary bottleneck move across
    // several neighboring rungs without sacrificing the hot end of the
    // ladder in one aggressive quantile remap.
    weighted_gaps[index] = gap * std::sqrt(acceptance);
    weighted_total += weighted_gaps[index];
  }
  if (!(weighted_total > 0.0L) ||
      !std::isfinite(weighted_total)) {
    throw std::runtime_error(
        "temperature feedback produced invalid total gap");
  }
  const long double total = beta.front() - beta.back();
  const long double scale = total / weighted_total;
  std::vector<long double> result(beta.size());
  result.front() = beta.front();
  for (std::size_t index = 0; index < attempts.size(); ++index) {
    result[index + 1] =
        result[index] - scale * weighted_gaps[index];
  }
  result.back() = beta.back();
  for (std::size_t index = 0; index + 1 < result.size(); ++index) {
    if (!(result[index] > result[index + 1])) {
      throw std::runtime_error(
          "temperature feedback lost strict ordering");
    }
  }
  return result;
}

std::vector<Replica> pt_initialize_replicas(
    const std::vector<State>& seeds, int count,
    const std::array<std::uint64_t, kCoreEntries>& zobrist) {
  std::vector<Replica> replicas;
  replicas.reserve(static_cast<std::size_t>(count));
  for (int index = 0; index < count; ++index) {
    Replica replica;
    replica.state =
        seeds[static_cast<std::size_t>(index) % seeds.size()];
    replica.hash = core_hash(replica.state.core, zobrist);
    replica.walker.id = index;
    if (index == 0) {
      replica.walker.last_endpoint = 0;
      replica.walker.armed_from_cold = true;
    } else if (index == count - 1) {
      replica.walker.last_endpoint = 1;
    }
    replicas.push_back(std::move(replica));
  }
  return replicas;
}

void pt_reset_walker_tracking(
    std::vector<Replica>& replicas) {
  for (std::size_t index = 0; index < replicas.size(); ++index) {
    Walker& walker = replicas[index].walker;
    walker.last_endpoint = -1;
    walker.armed_from_cold = false;
    if (index == 0) {
      walker.last_endpoint = 0;
      walker.armed_from_cold = true;
    } else if (index + 1 == replicas.size()) {
      walker.last_endpoint = 1;
    }
  }
}

void pt_elite_reseed(
    std::vector<Replica>& replicas,
    const std::vector<State>& class_seeds,
    const State& retained_best,
    const std::array<std::uint64_t, kCoreEntries>& zobrist,
    PTStatistics& statistics) {
  const std::size_t count =
      std::min(replicas.size(), class_seeds.size());
  if (count == 0) {
    throw std::runtime_error(
        "elite reseed requires at least one replica and seed");
  }
  for (std::size_t index = 0; index < count; ++index) {
    replicas[index].state =
        index == 0 ? retained_best : class_seeds[index];
    replicas[index].hash =
        core_hash(replicas[index].state.core, zobrist);
  }
  ++statistics.elite_reseed_events;
  statistics.elite_reseeded_replicas += count;
  pt_reset_walker_tracking(replicas);
  ++statistics.walker_tracking_resets;
}

void pt_visit_endpoint(
    Walker& walker, int endpoint, bool post_feedback,
    PTStatistics& statistics) {
  if (walker.last_endpoint == endpoint) {
    return;
  }
  if (walker.last_endpoint == -1) {
    walker.last_endpoint = endpoint;
    if (endpoint == 0) {
      walker.armed_from_cold = true;
    }
    return;
  }
  ++statistics.endpoint_crossings;
  if (post_feedback) {
    ++statistics.post_endpoint_crossings;
  }
  if (endpoint == 0) {
    if (walker.armed_from_cold &&
        walker.last_endpoint == 1) {
      ++statistics.completed_round_trips;
      if (post_feedback) {
        ++statistics.post_completed_round_trips;
      }
    }
    walker.armed_from_cold = true;
  }
  walker.last_endpoint = endpoint;
}

void pt_track_endpoints(
    std::vector<Replica>& replicas,
    PTStatistics& statistics) {
  pt_visit_endpoint(
      replicas.front().walker, 0,
      statistics.feedback_applied, statistics);
  pt_visit_endpoint(
      replicas.back().walker, 1,
      statistics.feedback_applied, statistics);
}

std::string pt_path_array_json(
    const std::vector<fs::path>& paths) {
  std::ostringstream output;
  output << '[';
  for (std::size_t index = 0; index < paths.size(); ++index) {
    if (index != 0) {
      output << ',';
    }
    output << '"' << json_escape(paths[index].string()) << '"';
  }
  output << ']';
  return output.str();
}

std::string pt_u64_array_json(
    const std::vector<std::uint64_t>& values) {
  std::ostringstream output;
  output << '[';
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index != 0) {
      output << ',';
    }
    output << values[index];
  }
  output << ']';
  return output.str();
}

std::string pt_move_array_json(
    const std::array<std::uint64_t, kPTMoveKinds>& values) {
  std::ostringstream output;
  output << '{';
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index != 0) {
      output << ',';
    }
    output << '"' << pt_kind_name(index) << "\":"
           << values[index];
  }
  output << '}';
  return output.str();
}

std::string pt_real_array_json(
    const std::vector<long double>& values,
    bool reciprocal = false) {
  std::ostringstream output;
  output << '[' << std::setprecision(17);
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index != 0) {
      output << ',';
    }
    const long double value =
        reciprocal ? 1.0L / values[index] : values[index];
    output << static_cast<double>(value);
  }
  output << ']';
  return output.str();
}

std::string pt_rate_array_json(
    const std::vector<std::uint64_t>& attempts,
    const std::vector<std::uint64_t>& accepts) {
  std::ostringstream output;
  output << '[' << std::setprecision(17);
  for (std::size_t index = 0; index < attempts.size(); ++index) {
    if (index != 0) {
      output << ',';
    }
    const double rate =
        attempts[index] == 0
            ? 0.0
            : static_cast<double>(accepts[index]) /
                  static_cast<double>(attempts[index]);
    output << rate;
  }
  output << ']';
  return output.str();
}

std::string pt_score_string(std::uint64_t quotient) {
  return wide_to_string(
      static_cast<Wide>(quotient) * kScale);
}

std::pair<std::uint64_t, std::uint64_t>
pt_replica_score_range(const std::vector<Replica>& replicas) {
  std::uint64_t minimum =
      std::numeric_limits<std::uint64_t>::max();
  std::uint64_t maximum = 0;
  for (const Replica& replica : replicas) {
    minimum =
        std::min(minimum, magnitude(replica.state.determinant));
    maximum =
        std::max(maximum, magnitude(replica.state.determinant));
  }
  return {minimum, maximum};
}

void pt_log_state(
    PTLogger& logger, const PTArguments& arguments,
    const PTStatistics& statistics, std::string_view event,
    double elapsed, const std::vector<long double>& beta,
    const std::vector<Replica>& replicas,
    std::uint64_t best_quotient) {
  const auto [minimum, maximum] =
      pt_replica_score_range(replicas);
  std::ostringstream record;
  record
      << "{\"event\":\"" << event << "\""
      << ",\"engine\":\"exact-core-parallel-tempering-v1\""
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(6) << elapsed
      << ",\"seed\":" << arguments.seed
      << ",\"sweeps\":" << statistics.sweeps
      << ",\"local_proposals\":"
      << statistics.local_proposals
      << ",\"local_accepts\":" << statistics.local_accepts
      << ",\"local_rejects\":" << statistics.local_rejects
      << ",\"singular_rejects\":"
      << statistics.singular_rejects
      << ",\"null_proposals\":"
      << statistics.null_proposals
      << ",\"proposed_by_kind\":"
      << pt_move_array_json(statistics.proposed_by_kind)
      << ",\"accepted_by_kind\":"
      << pt_move_array_json(statistics.accepted_by_kind)
      << ",\"swap_attempts\":"
      << pt_u64_array_json(statistics.swap_attempts)
      << ",\"swap_accepts\":"
      << pt_u64_array_json(statistics.swap_accepts)
      << ",\"swap_acceptance\":"
      << pt_rate_array_json(
             statistics.swap_attempts,
             statistics.swap_accepts)
      << ",\"post_swap_acceptance\":"
      << pt_rate_array_json(
             statistics.post_swap_attempts,
             statistics.post_swap_accepts)
      << ",\"endpoint_crossings\":"
      << statistics.endpoint_crossings
      << ",\"completed_round_trips\":"
      << statistics.completed_round_trips
      << ",\"post_endpoint_crossings\":"
      << statistics.post_endpoint_crossings
      << ",\"post_completed_round_trips\":"
      << statistics.post_completed_round_trips
      << ",\"elite_reseed_events\":"
      << statistics.elite_reseed_events
      << ",\"elite_reseeded_replicas\":"
      << statistics.elite_reseeded_replicas
      << ",\"walker_tracking_resets\":"
      << statistics.walker_tracking_resets
      << ",\"feedback_applied\":"
      << (statistics.feedback_applied ? "true" : "false")
      << ",\"beta\":" << pt_real_array_json(beta)
      << ",\"temperature\":"
      << pt_real_array_json(beta, true)
      << ",\"replica_min_quotient\":" << minimum
      << ",\"replica_max_quotient\":" << maximum
      << ",\"promotions\":" << statistics.promotions
      << ",\"identity_checks\":" << statistics.identity_checks
      << ",\"exact_checks\":" << statistics.exact_checks
      << ",\"differential_checks\":"
      << statistics.differential_checks
      << ",\"best_core_quotient\":" << best_quotient
      << ",\"best_absolute_determinant\":\""
      << pt_score_string(best_quotient) << "\""
      << ",\"above_frontier_unverified\":"
      << (static_cast<Wide>(best_quotient) * kScale >
                  static_cast<Wide>(arguments.frontier)
              ? "true"
              : "false")
      << '}';
  logger.line(record.str());
}

void pt_run_move_differentials(
    const std::vector<State>& seeds,
    const PTArguments& arguments,
    const std::array<std::uint64_t, kCoreEntries>& zobrist,
    std::mt19937_64& randomizer,
    PTStatistics& statistics) {
  for (int sample = 0;
       sample < arguments.differential_samples; ++sample) {
    const State& source =
        seeds[static_cast<std::size_t>(sample) % seeds.size()];
    const PTMoveKind kind =
        static_cast<PTMoveKind>(sample % kPTMoveKinds);
    const PTMove move =
        pt_propose_kind(source, kind, randomizer);
    if (!move.valid) {
      --sample;
      continue;
    }
    const CoreMatrix materialized =
        pt_materialize(source.core, move);
    const std::int64_t exact =
        exact_core_determinant(materialized);
    if (exact != move.determinant) {
      throw std::runtime_error(
          std::string("proposal determinant audit failed for ") +
          pt_kind_name(pt_kind_index(kind)));
    }
    if (exact == 0) {
      continue;
    }
    Replica replica;
    replica.state = source;
    replica.hash = core_hash(replica.state.core, zobrist);
    pt_apply(replica, move, zobrist, statistics);
    if (replica.state.core != materialized ||
        replica.state.determinant != exact) {
      throw std::runtime_error(
          std::string("move application audit failed for ") +
          pt_kind_name(pt_kind_index(kind)));
    }
    check_adjugate_identity(replica.state);
    pt_differential_check(
        replica.state,
        std::string("differential move sample ") +
            std::to_string(sample),
        statistics);
  }
}

std::string pt_final_summary(
    const PTArguments& arguments,
    const TemperatureCalibration& calibration,
    const PTStatistics& statistics,
    const std::vector<long double>& initial_beta,
    const std::vector<long double>& final_beta,
    const std::vector<Replica>& replicas,
    std::uint64_t best_quotient, double elapsed,
    bool stopped_by_signal) {
  const auto [minimum, maximum] =
      pt_replica_score_range(replicas);
  std::ostringstream output;
  output
      << "{\"schema_version\":1"
      << ",\"engine\":\"exact-core-parallel-tempering-v1\""
      << ",\"complete\":"
      << (stopped_by_signal ? "false" : "true")
      << ",\"stopped_by_signal\":"
      << (stopped_by_signal ? "true" : "false")
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(6) << elapsed
      << ",\"seed\":" << arguments.seed
      << ",\"input_seeds\":"
      << pt_path_array_json(arguments.seeds)
      << ",\"replicas\":" << replicas.size()
      << ",\"objective_distribution\":"
      << "\"pi_beta(B) proportional to |det(B)|^beta\""
      << ",\"swap_sign_audit\":true"
      << ",\"downhill_calibration_samples\":"
      << calibration.downhill_samples
      << ",\"median_log_determinant_drop\":"
      << static_cast<double>(calibration.median_log_drop)
      << ",\"cold_downhill_target\":"
      << arguments.cold_downhill_target
      << ",\"hot_downhill_target\":"
      << arguments.hot_downhill_target
      << ",\"initial_beta\":"
      << pt_real_array_json(initial_beta)
      << ",\"initial_temperature\":"
      << pt_real_array_json(initial_beta, true)
      << ",\"final_beta\":"
      << pt_real_array_json(final_beta)
      << ",\"final_temperature\":"
      << pt_real_array_json(final_beta, true)
      << ",\"feedback_applied\":"
      << (statistics.feedback_applied ? "true" : "false")
      << ",\"feedback_elapsed_seconds\":"
      << statistics.feedback_elapsed_seconds
      << ",\"calibration_swap_attempts\":"
      << pt_u64_array_json(
             statistics.calibration_swap_attempts)
      << ",\"calibration_swap_accepts\":"
      << pt_u64_array_json(
             statistics.calibration_swap_accepts)
      << ",\"calibration_swap_acceptance\":"
      << pt_rate_array_json(
             statistics.calibration_swap_attempts,
             statistics.calibration_swap_accepts)
      << ",\"post_swap_attempts\":"
      << pt_u64_array_json(statistics.post_swap_attempts)
      << ",\"post_swap_accepts\":"
      << pt_u64_array_json(statistics.post_swap_accepts)
      << ",\"post_swap_acceptance\":"
      << pt_rate_array_json(
             statistics.post_swap_attempts,
             statistics.post_swap_accepts)
      << ",\"sweeps\":" << statistics.sweeps
      << ",\"local_proposals\":"
      << statistics.local_proposals
      << ",\"local_accepts\":" << statistics.local_accepts
      << ",\"local_rejects\":" << statistics.local_rejects
      << ",\"singular_rejects\":"
      << statistics.singular_rejects
      << ",\"null_proposals\":"
      << statistics.null_proposals
      << ",\"proposed_by_kind\":"
      << pt_move_array_json(statistics.proposed_by_kind)
      << ",\"accepted_by_kind\":"
      << pt_move_array_json(statistics.accepted_by_kind)
      << ",\"swap_attempts\":"
      << pt_u64_array_json(statistics.swap_attempts)
      << ",\"swap_accepts\":"
      << pt_u64_array_json(statistics.swap_accepts)
      << ",\"swap_acceptance\":"
      << pt_rate_array_json(
             statistics.swap_attempts,
             statistics.swap_accepts)
      << ",\"endpoint_crossings\":"
      << statistics.endpoint_crossings
      << ",\"completed_round_trips\":"
      << statistics.completed_round_trips
      << ",\"post_endpoint_crossings\":"
      << statistics.post_endpoint_crossings
      << ",\"post_completed_round_trips\":"
      << statistics.post_completed_round_trips
      << ",\"reseed_sweeps\":"
      << arguments.reseed_sweeps
      << ",\"elite_reseed_events\":"
      << statistics.elite_reseed_events
      << ",\"elite_reseeded_replicas\":"
      << statistics.elite_reseeded_replicas
      << ",\"walker_tracking_resets\":"
      << statistics.walker_tracking_resets
      << ",\"replica_min_quotient\":" << minimum
      << ",\"replica_max_quotient\":" << maximum
      << ",\"promotions\":" << statistics.promotions
      << ",\"identity_checks\":" << statistics.identity_checks
      << ",\"exact_checks\":" << statistics.exact_checks
      << ",\"differential_checks\":"
      << statistics.differential_checks
      << ",\"best_core_quotient\":" << best_quotient
      << ",\"best_absolute_determinant\":\""
      << pt_score_string(best_quotient) << "\""
      << ",\"frontier\":\"" << arguments.frontier << "\""
      << ",\"above_frontier_unverified\":"
      << (static_cast<Wide>(best_quotient) * kScale >
                  static_cast<Wide>(arguments.frontier)
              ? "true"
              : "false")
      << ",\"claim_boundary\":["
      << "\"All determinant objectives and retained scores are exact integers.\","
      << "\"Metropolis probabilities use logs of ratios of exact determinants.\","
      << "\"This pilot makes no equilibrium, optimality, or novelty claim.\","
      << "\"Any strict promotion still requires arena verification and a separate H/HT audit.\""
      << "]}";
  return output.str() + "\n";
}

int run_parallel_tempering(const PTArguments& arguments) {
  if (fs::exists(arguments.output) ||
      fs::exists(arguments.log) ||
      fs::exists(arguments.summary)) {
    throw std::runtime_error(
        "output, log, and summary paths must be fresh");
  }
  PTLogger logger(arguments.log);
  std::signal(SIGINT, request_stop);
  std::signal(SIGTERM, request_stop);
  PTStatistics statistics;
  pt_audit_swap_sign();
  const std::size_t interval_count =
      static_cast<std::size_t>(arguments.replica_count - 1);
  statistics.swap_attempts.assign(interval_count, 0);
  statistics.swap_accepts.assign(interval_count, 0);
  statistics.calibration_swap_attempts.assign(interval_count, 0);
  statistics.calibration_swap_accepts.assign(interval_count, 0);
  statistics.post_swap_attempts.assign(interval_count, 0);
  statistics.post_swap_accepts.assign(interval_count, 0);

  std::mt19937_64 randomizer(arguments.seed);
  const auto zobrist = make_zobrist(
      arguments.seed ^ UINT64_C(0x50545245584d4158));
  const std::vector<State> seeds =
      pt_load_seeds(arguments, statistics);
  const TemperatureCalibration calibration =
      pt_calibrate_temperatures(seeds, arguments);
  const std::vector<long double> initial_beta =
      pt_geometric_ladder(
          arguments.replica_count, calibration.beta_cold,
          calibration.beta_hot);
  std::vector<long double> beta = initial_beta;
  std::vector<Replica> replicas =
      pt_initialize_replicas(
          seeds, arguments.replica_count, zobrist);

  pt_run_move_differentials(
      seeds, arguments, zobrist, randomizer, statistics);

  State retained_best = seeds.front();
  for (const State& seed : seeds) {
    if (magnitude(seed.determinant) >
        magnitude(retained_best.determinant)) {
      retained_best = seed;
    }
  }
  std::uint64_t nonce = 0;
  atomic_write(
      arguments.output,
      sign_matrix_text(core_to_sign(retained_best.core)), nonce++);

  const auto started = Clock::now();
  const auto deadline =
      started + std::chrono::duration<double>(arguments.seconds);
  auto next_heartbeat =
      started +
      std::chrono::duration<double>(
          arguments.heartbeat_seconds);
  const auto feedback_time =
      started +
      std::chrono::duration<double>(
          arguments.feedback_seconds);
  std::uint64_t next_reseed_sweep =
      arguments.feedback_seconds == 0.0
          ? arguments.reseed_sweeps
          : 0;

  {
    std::ostringstream record;
    record
        << "{\"event\":\"started\""
        << ",\"engine\":\"exact-core-parallel-tempering-v1\""
        << ",\"seed\":" << arguments.seed
        << ",\"seconds\":" << arguments.seconds
        << ",\"frontier\":\"" << arguments.frontier << "\""
        << ",\"input_seeds\":"
        << pt_path_array_json(arguments.seeds)
        << ",\"replicas\":" << arguments.replica_count
        << ",\"downhill_calibration_samples\":"
        << calibration.downhill_samples
        << ",\"median_log_determinant_drop\":"
        << static_cast<double>(calibration.median_log_drop)
        << ",\"initial_beta\":"
        << pt_real_array_json(initial_beta)
        << ",\"initial_temperature\":"
        << pt_real_array_json(initial_beta, true)
        << ",\"feedback_seconds\":"
        << arguments.feedback_seconds
        << ",\"reseed_sweeps\":"
        << arguments.reseed_sweeps
        << ",\"differential_checks\":"
        << statistics.differential_checks
        << ",\"best_absolute_determinant\":\""
        << pt_score_string(magnitude(retained_best.determinant))
        << "\"}";
    logger.line(record.str());
  }

  while (!stop_requested && Clock::now() < deadline) {
    for (std::size_t index = 0;
         index < replicas.size() && !stop_requested; ++index) {
      Replica& replica = replicas[index];
      const PTMove move =
          pt_propose(replica.state, randomizer);
      ++statistics.local_proposals;
      ++statistics.proposed_by_kind[
          pt_kind_index(move.kind)];
      if (!move.valid) {
        ++statistics.null_proposals;
        continue;
      }
      if (move.determinant == 0) {
        ++statistics.singular_rejects;
        ++statistics.local_rejects;
        continue;
      }
      const std::uint64_t current =
          magnitude(replica.state.determinant);
      const std::uint64_t candidate =
          magnitude(move.determinant);
      if (!pt_metropolis_accept(
              current, candidate, beta[index], randomizer)) {
        ++statistics.local_rejects;
        continue;
      }
      pt_apply(replica, move, zobrist, statistics);
      ++statistics.local_accepts;
      ++statistics.accepted_by_kind[
          pt_kind_index(move.kind)];

      if (candidate > magnitude(retained_best.determinant)) {
        retained_best = replica.state;
        ++statistics.promotions;
        atomic_write(
            arguments.output,
            sign_matrix_text(core_to_sign(retained_best.core)),
            nonce++);
        const double elapsed =
            std::chrono::duration<double>(
                Clock::now() - started)
                .count();
        pt_log_state(
            logger, arguments, statistics, "new_best",
            elapsed, beta, replicas,
            magnitude(retained_best.determinant));
        std::cout
            << "new best |det|="
            << pt_score_string(
                   magnitude(retained_best.determinant))
            << " sweep=" << statistics.sweeps << '\n'
            << std::flush;
      }

      if (arguments.identity_interval != 0 &&
          statistics.local_accepts %
                  arguments.identity_interval ==
              0) {
        check_adjugate_identity(replica.state);
        ++statistics.identity_checks;
      }
    }

    const std::size_t parity =
        static_cast<std::size_t>(statistics.sweeps & UINT64_C(1));
    for (std::size_t index = parity;
         index + 1 < replicas.size(); index += 2) {
      ++statistics.swap_attempts[index];
      if (statistics.feedback_applied) {
        ++statistics.post_swap_attempts[index];
      }
      if (!pt_swap_accept(
              replicas[index], replicas[index + 1],
              beta[index], beta[index + 1], randomizer)) {
        continue;
      }
      ++statistics.swap_accepts[index];
      if (statistics.feedback_applied) {
        ++statistics.post_swap_accepts[index];
      }
      std::swap(replicas[index], replicas[index + 1]);
    }
    pt_track_endpoints(replicas, statistics);
    ++statistics.sweeps;

    if (arguments.reseed_sweeps != 0 &&
        next_reseed_sweep != 0 &&
        statistics.sweeps >= next_reseed_sweep) {
      pt_elite_reseed(
          replicas, seeds, retained_best, zobrist, statistics);
      next_reseed_sweep =
          statistics.sweeps + arguments.reseed_sweeps;
    }

    if (arguments.exact_check_interval != 0 &&
        statistics.local_proposals %
                arguments.exact_check_interval <
            static_cast<std::uint64_t>(
                arguments.replica_count)) {
      const std::size_t audit_index =
          static_cast<std::size_t>(
              statistics.sweeps %
              static_cast<std::uint64_t>(replicas.size()));
      pt_differential_check(
          replicas[audit_index].state,
          "periodic replica audit", statistics);
      ++statistics.exact_checks;
    }

    const auto now = Clock::now();
    if (!statistics.feedback_applied &&
        arguments.feedback_seconds > 0.0 &&
        now >= feedback_time) {
      statistics.calibration_swap_attempts =
          statistics.swap_attempts;
      statistics.calibration_swap_accepts =
          statistics.swap_accepts;
      beta = pt_feedback_ladder(
          beta, statistics.calibration_swap_attempts,
          statistics.calibration_swap_accepts);
      statistics.feedback_applied = true;
      statistics.feedback_elapsed_seconds =
          std::chrono::duration<double>(now - started).count();
      pt_reset_walker_tracking(replicas);
      ++statistics.walker_tracking_resets;
      if (arguments.reseed_sweeps != 0) {
        next_reseed_sweep =
            statistics.sweeps + arguments.reseed_sweeps;
      }
      pt_log_state(
          logger, arguments, statistics, "feedback",
          statistics.feedback_elapsed_seconds, beta, replicas,
          magnitude(retained_best.determinant));
    }

    if (arguments.heartbeat_seconds > 0.0 &&
        now >= next_heartbeat) {
      const double elapsed =
          std::chrono::duration<double>(now - started).count();
      pt_log_state(
          logger, arguments, statistics, "heartbeat",
          elapsed, beta, replicas,
          magnitude(retained_best.determinant));
      next_heartbeat =
          now +
          std::chrono::duration<double>(
              arguments.heartbeat_seconds);
    }
  }

  if (!statistics.feedback_applied) {
    statistics.calibration_swap_attempts =
        statistics.swap_attempts;
    statistics.calibration_swap_accepts =
        statistics.swap_accepts;
  }
  for (const Replica& replica : replicas) {
    check_adjugate_identity(replica.state);
    ++statistics.identity_checks;
  }
  check_adjugate_identity(retained_best);
  ++statistics.identity_checks;
  pt_differential_check(
      retained_best, "final global best", statistics);
  atomic_write(
      arguments.output,
      sign_matrix_text(core_to_sign(retained_best.core)), nonce++);

  const double elapsed =
      std::chrono::duration<double>(Clock::now() - started).count();
  atomic_write(
      arguments.summary,
      pt_final_summary(
          arguments, calibration, statistics, initial_beta, beta,
          replicas, magnitude(retained_best.determinant), elapsed,
          stop_requested != 0),
      nonce++);
  pt_log_state(
      logger, arguments, statistics,
      stop_requested ? "stopped" : "finished",
      elapsed, beta, replicas,
      magnitude(retained_best.determinant));

  std::cout
      << "finished |det|="
      << pt_score_string(magnitude(retained_best.determinant))
      << " sweeps=" << statistics.sweeps
      << " local_proposals=" << statistics.local_proposals
      << " round_trips=" << statistics.completed_round_trips
      << '\n';
  if (static_cast<Wide>(
          magnitude(retained_best.determinant)) *
          kScale >
      static_cast<Wide>(arguments.frontier)) {
    std::cout
        << "UNVERIFIED FRONTIER CANDIDATE: run ./arena verify "
        << arguments.output << '\n';
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    return run_parallel_tempering(
        parse_pt_arguments(argc, argv));
  } catch (const std::exception& error) {
    std::cerr
        << "core_parallel_tempering: "
        << error.what() << '\n';
    return 1;
  }
}
