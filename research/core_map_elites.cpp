// Exact quality-diversity / MAP-Elites pilot for the order-23 maximal
// determinant search.
//
// The engine reuses the independently differential-tested dephased-core
// arithmetic and extremal-optimization move kernels.  It differs from the
// single-walker EO search by retaining the best exact determinant found in
// each structural niche.  Niches are defined by:
//
//   * the nearest supplied frontier seed;
//   * Hamming distance from that seed;
//   * the count of nearly balanced rows and columns; and
//   * off-diagonal row-Gram energy.
//
// Parents are sampled from small tournaments that alternate low-use niches
// and high-quality niches.  Each parent receives an optional exact multi-bit
// destroy step followed by a bounded EO repair burst.  This is an
// optimization heuristic, not an exhaustive search or a novelty classifier.
// H/HT-equivalence claims still require the separate pinned classifier.

#include <map>
#include <set>
#include <tuple>
#include <unordered_set>

#define CORE_EXTREMAL_OPTIMIZATION_NO_MAIN
#include "core_extremal_optimization.cpp"
#undef CORE_EXTREMAL_OPTIMIZATION_NO_MAIN

namespace {

constexpr std::uint64_t kQdFrontier =
    UINT64_C(2779447296000000);

struct QdArguments {
  std::vector<fs::path> seeds;
  fs::path output;
  fs::path archive_directory;
  fs::path log;
  fs::path summary;
  std::uint64_t seed = UINT64_C(35001);
  std::uint64_t frontier = kQdFrontier;
  double seconds = 300.0;
  double heartbeat_seconds = 15.0;
  double tau = 1.45;
  int expected_seeds = 3;
  int burst_moves = 512;
  int archive_stride = 64;
  int tournament_size = 8;
  int maximum_cells = 20000;
  int archive_outputs = 32;
  int kick_min_flips = 4;
  int kick_max_flips = 24;
  int kick_attempts = 12;
  int kick_numerator = 3;
  int kick_denominator = 4;
  int pair_interval = 127;
  int complement_interval = 61;
  int distance_bin_width = 24;
  int gram_bin_width = 256;
  int balance_bin_width = 5;
  std::uint64_t recency_window = 32;
  std::uint64_t recency_penalty = UINT64_C(250000);
  int differential_rounds = 16;
};

struct Descriptor {
  int nearest_seed = 0;
  int distance_bin = 0;
  int balance_bin = 0;
  int gram_bin = 0;

  auto as_tuple() const {
    return std::tie(
        nearest_seed, distance_bin, balance_bin, gram_bin);
  }

  bool operator<(const Descriptor& other) const {
    return as_tuple() < other.as_tuple();
  }

  bool operator==(const Descriptor& other) const {
    return as_tuple() == other.as_tuple();
  }
};

struct DescriptorMetrics {
  Descriptor cell{};
  std::uint64_t nearest_distance = 0;
  std::uint64_t balanced_lines = 0;
  std::uint64_t gram_energy = 0;
};

struct Elite {
  CoreMatrix core{};
  std::int64_t determinant = 0;
  std::uint64_t hash = 0;
  std::uint64_t selections = 0;
  std::uint64_t inserted_at_move = 0;
  DescriptorMetrics metrics{};
};

struct QdStatistics {
  std::uint64_t bursts = 0;
  std::uint64_t moves = 0;
  std::uint64_t parent_tournaments = 0;
  std::uint64_t low_use_parent_choices = 0;
  std::uint64_t quality_parent_choices = 0;
  std::uint64_t exact_kick_evaluations = 0;
  std::uint64_t exact_kicks = 0;
  std::uint64_t kick_failures = 0;
  std::uint64_t exact_state_rebuilds = 0;
  std::uint64_t descriptor_evaluations = 0;
  std::uint64_t archive_insertions = 0;
  std::uint64_t archive_replacements = 0;
  std::uint64_t archive_rejections = 0;
  std::uint64_t archive_capacity_rejections = 0;
  std::uint64_t frontier_state_visits = 0;
  std::uint64_t distinct_frontier_cores = 0;
  std::uint64_t promotions = 0;
  std::uint64_t exact_determinant_checks = 0;
  std::uint64_t identity_checks = 0;
  std::uint64_t hash_checks = 0;
  std::uint64_t maximum_nearest_seed_distance = 0;
  std::uint64_t minimum_observed_gram_energy =
      std::numeric_limits<std::uint64_t>::max();
  std::uint64_t maximum_observed_gram_energy = 0;
  EoStatistics eo{};
};

int qd_positive_int(
    std::string_view text, std::string_view option,
    int maximum = std::numeric_limits<int>::max()) {
  const std::uint64_t value = strict_unsigned(text, option);
  if (value == 0 ||
      value > static_cast<std::uint64_t>(maximum)) {
    throw std::runtime_error(
        std::string(option) + " is outside its valid range");
  }
  return static_cast<int>(value);
}

QdArguments parse_qd_arguments(int argc, char** argv) {
  QdArguments arguments;
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
    } else if (option == "--archive-dir") {
      arguments.archive_directory = value();
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
    } else if (
        option == "--heartbeat" ||
        option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds =
          strict_double(value(), option, true);
    } else if (option == "--tau") {
      arguments.tau = strict_double(value(), option, false);
      if (arguments.tau < 0.1 || arguments.tau > 10.0) {
        throw std::runtime_error("--tau must be in [0.1,10]");
      }
    } else if (option == "--expected-seeds") {
      arguments.expected_seeds =
          qd_positive_int(value(), option, 64);
    } else if (option == "--burst-moves") {
      arguments.burst_moves =
          qd_positive_int(value(), option, 1000000);
    } else if (option == "--archive-stride") {
      arguments.archive_stride =
          qd_positive_int(value(), option, 1000000);
    } else if (option == "--tournament-size") {
      arguments.tournament_size =
          qd_positive_int(value(), option, 1024);
    } else if (option == "--maximum-cells") {
      arguments.maximum_cells =
          qd_positive_int(value(), option, 1000000);
    } else if (option == "--archive-outputs") {
      arguments.archive_outputs =
          qd_positive_int(value(), option, 100000);
    } else if (option == "--kick-min-flips") {
      arguments.kick_min_flips =
          qd_positive_int(value(), option, kCoreEntries);
    } else if (option == "--kick-max-flips") {
      arguments.kick_max_flips =
          qd_positive_int(value(), option, kCoreEntries);
    } else if (option == "--kick-attempts") {
      arguments.kick_attempts =
          qd_positive_int(value(), option, 1000000);
    } else if (option == "--kick-numerator") {
      arguments.kick_numerator =
          qd_positive_int(value(), option, 1000000);
    } else if (option == "--kick-denominator") {
      arguments.kick_denominator =
          qd_positive_int(value(), option, 1000000);
    } else if (option == "--pair-interval") {
      arguments.pair_interval =
          qd_positive_int(value(), option, 1000000);
    } else if (option == "--complement-interval") {
      arguments.complement_interval =
          qd_positive_int(value(), option, 1000000);
    } else if (option == "--distance-bin-width") {
      arguments.distance_bin_width =
          qd_positive_int(value(), option, kCoreEntries);
    } else if (option == "--gram-bin-width") {
      arguments.gram_bin_width =
          qd_positive_int(value(), option, 1000000000);
    } else if (option == "--balance-bin-width") {
      arguments.balance_bin_width =
          qd_positive_int(value(), option, 2 * kCoreOrder);
    } else if (option == "--recency-window") {
      arguments.recency_window =
          eo_strict_positive(value(), option, UINT64_C(1000000));
    } else if (option == "--recency-penalty") {
      arguments.recency_penalty =
          eo_strict_positive(value(), option, UINT64_C(1000000000));
    } else if (option == "--differential-rounds") {
      arguments.differential_rounds =
          qd_positive_int(value(), option, 100000);
    } else if (option == "--help") {
      std::cout
          << "usage: core_map_elites "
             "--seed-matrix H0 --seed-matrix H1 "
             "--seed-matrix H2 --output MATRIX "
             "--archive-dir DIR --log JSONL --summary JSON "
             "[options]\n"
          << "  --seed N --frontier N --seconds S --heartbeat S\n"
          << "  --burst-moves N --archive-stride N "
             "--tournament-size N --maximum-cells N\n"
          << "  --kick-min-flips N --kick-max-flips N "
             "--kick-attempts N\n"
          << "  --kick-numerator N --kick-denominator N\n"
          << "  --pair-interval N --complement-interval N\n"
          << "  --distance-bin-width N --gram-bin-width N "
             "--balance-bin-width N\n"
          << "  --tau X --recency-window N "
             "--recency-penalty N --differential-rounds N\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }

  if (static_cast<int>(arguments.seeds.size()) !=
      arguments.expected_seeds) {
    throw std::runtime_error(
        "expected " + std::to_string(arguments.expected_seeds) +
        " --seed-matrix inputs but received " +
        std::to_string(arguments.seeds.size()));
  }
  if (
      arguments.output.empty() ||
      arguments.archive_directory.empty() ||
      arguments.log.empty() || arguments.summary.empty()) {
    throw std::runtime_error(
        "--output, --archive-dir, --log, and --summary are required");
  }
  if (arguments.kick_min_flips > arguments.kick_max_flips) {
    throw std::runtime_error(
        "--kick-min-flips must not exceed --kick-max-flips");
  }
  if (arguments.kick_numerator > arguments.kick_denominator) {
    throw std::runtime_error(
        "--kick-numerator must not exceed --kick-denominator");
  }
  const std::array<fs::path, 3> files{
      arguments.output, arguments.log, arguments.summary};
  for (std::size_t first = 0; first < files.size(); ++first) {
    if (fs::exists(files[first])) {
      throw std::runtime_error(
          "refusing to overwrite output: " +
          files[first].string());
    }
    for (std::size_t second = first + 1;
         second < files.size(); ++second) {
      if (files[first] == files[second]) {
        throw std::runtime_error(
            "output, log, and summary paths must differ");
      }
    }
  }
  if (fs::exists(arguments.archive_directory)) {
    throw std::runtime_error(
        "refusing to overwrite archive directory: " +
        arguments.archive_directory.string());
  }
  return arguments;
}

std::string qd_core_key(const CoreMatrix& core) {
  std::string result;
  result.reserve(kCoreEntries);
  for (const auto& row : core) {
    for (const std::uint8_t value : row) {
      result.push_back(static_cast<char>('0' + value));
    }
  }
  return result;
}

DescriptorMetrics describe_core(
    const CoreMatrix& core,
    const std::vector<CoreMatrix>& seeds,
    const QdArguments& arguments) {
  DescriptorMetrics result;
  result.nearest_distance =
      std::numeric_limits<std::uint64_t>::max();
  for (std::size_t index = 0; index < seeds.size(); ++index) {
    const std::uint64_t distance =
        core_hamming(core, seeds[index]);
    if (distance < result.nearest_distance) {
      result.nearest_distance = distance;
      result.cell.nearest_seed = static_cast<int>(index);
    }
  }

  const SignMatrix sign = core_to_sign(core);
  for (int row = 1; row < kSignOrder; ++row) {
    int sum = 0;
    for (int column = 0; column < kSignOrder; ++column) {
      sum += sign[row][column];
    }
    if (std::abs(sum) <= 3) ++result.balanced_lines;
  }
  for (int column = 1; column < kSignOrder; ++column) {
    int sum = 0;
    for (int row = 0; row < kSignOrder; ++row) {
      sum += sign[row][column];
    }
    if (std::abs(sum) <= 3) ++result.balanced_lines;
  }

  for (int first = 0; first < kSignOrder; ++first) {
    for (int second = first + 1;
         second < kSignOrder; ++second) {
      int inner = 0;
      for (int column = 0; column < kSignOrder; ++column) {
        inner += sign[first][column] * sign[second][column];
      }
      result.gram_energy +=
          static_cast<std::uint64_t>(inner * inner);
    }
  }

  result.cell.distance_bin = std::min(
      20, static_cast<int>(
              result.nearest_distance /
              static_cast<std::uint64_t>(
                  arguments.distance_bin_width)));
  result.cell.balance_bin = static_cast<int>(
      result.balanced_lines /
      static_cast<std::uint64_t>(
          arguments.balance_bin_width));
  result.cell.gram_bin = std::min(
      63, static_cast<int>(
              result.gram_energy /
              static_cast<std::uint64_t>(
                  arguments.gram_bin_width)));
  return result;
}

State state_from_core(
    const CoreMatrix& core, QdStatistics& statistics) {
  State state;
  state.core = core;
  state.determinant = exact_core_determinant(state.core);
  if (state.determinant == 0) {
    throw std::runtime_error(
        "cannot rebuild a singular archive elite");
  }
  state.adjugate = exact_adjugate(state.core);
  check_adjugate_identity(state);
  ++statistics.exact_state_rebuilds;
  ++statistics.identity_checks;
  return state;
}

bool destroy_state(
    State& state, const QdArguments& arguments,
    std::mt19937_64& randomizer, QdStatistics& statistics) {
  std::array<int, kCoreEntries> coordinates{};
  std::iota(coordinates.begin(), coordinates.end(), 0);
  std::uniform_int_distribution<int> flip_count(
      arguments.kick_min_flips, arguments.kick_max_flips);

  std::optional<CoreMatrix> retained;
  std::int64_t retained_determinant = 0;
  std::uint64_t retained_magnitude = 0;
  for (int attempt = 0; attempt < arguments.kick_attempts;
       ++attempt) {
    std::shuffle(
        coordinates.begin(), coordinates.end(), randomizer);
    CoreMatrix candidate = state.core;
    const int flips = flip_count(randomizer);
    for (int index = 0; index < flips; ++index) {
      const int flat =
          coordinates[static_cast<std::size_t>(index)];
      candidate[flat / kCoreOrder][flat % kCoreOrder] ^= 1U;
    }
    const std::int64_t determinant =
        exact_core_determinant(candidate);
    ++statistics.exact_kick_evaluations;
    if (determinant == 0) continue;
    const std::uint64_t candidate_magnitude =
        magnitude(determinant);
    if (!retained.has_value() ||
        candidate_magnitude > retained_magnitude) {
      retained = candidate;
      retained_determinant = determinant;
      retained_magnitude = candidate_magnitude;
    }
  }
  if (!retained.has_value()) {
    ++statistics.kick_failures;
    return false;
  }

  state.core = *retained;
  state.determinant = retained_determinant;
  state.adjugate = exact_adjugate(state.core);
  check_adjugate_identity(state);
  ++statistics.exact_kicks;
  ++statistics.identity_checks;
  return true;
}

const Descriptor& choose_parent_cell(
    std::map<Descriptor, Elite>& archive,
    const std::vector<Descriptor>& cells,
    const QdArguments& arguments,
    std::mt19937_64& randomizer, QdStatistics& statistics) {
  if (cells.empty()) {
    throw std::runtime_error("quality-diversity archive is empty");
  }
  ++statistics.parent_tournaments;
  const bool choose_quality = (randomizer() & UINT64_C(3)) == 0;
  std::size_t chosen =
      static_cast<std::size_t>(randomizer() % cells.size());
  for (int trial = 1; trial < arguments.tournament_size; ++trial) {
    const std::size_t candidate =
        static_cast<std::size_t>(randomizer() % cells.size());
    const Elite& current = archive.at(cells[chosen]);
    const Elite& proposed = archive.at(cells[candidate]);
    if (choose_quality) {
      if (magnitude(proposed.determinant) >
          magnitude(current.determinant)) {
        chosen = candidate;
      }
    } else if (proposed.selections < current.selections) {
      chosen = candidate;
    }
  }
  if (choose_quality) {
    ++statistics.quality_parent_choices;
  } else {
    ++statistics.low_use_parent_choices;
  }
  return cells[chosen];
}

bool update_archive(
    const State& state, std::uint64_t hash,
    const std::vector<CoreMatrix>& seed_cores,
    const QdArguments& arguments,
    std::map<Descriptor, Elite>& archive,
    std::vector<Descriptor>& cells,
    QdStatistics& statistics) {
  DescriptorMetrics metrics =
      describe_core(state.core, seed_cores, arguments);
  ++statistics.descriptor_evaluations;
  statistics.maximum_nearest_seed_distance = std::max(
      statistics.maximum_nearest_seed_distance,
      metrics.nearest_distance);
  statistics.minimum_observed_gram_energy = std::min(
      statistics.minimum_observed_gram_energy,
      metrics.gram_energy);
  statistics.maximum_observed_gram_energy = std::max(
      statistics.maximum_observed_gram_energy,
      metrics.gram_energy);

  const auto found = archive.find(metrics.cell);
  if (found == archive.end()) {
    if (static_cast<int>(archive.size()) >=
        arguments.maximum_cells) {
      ++statistics.archive_capacity_rejections;
      return false;
    }
    archive.emplace(
        metrics.cell,
        Elite{
            state.core, state.determinant, hash, 0,
            statistics.moves, metrics});
    cells.push_back(metrics.cell);
    ++statistics.archive_insertions;
    return true;
  }
  if (magnitude(state.determinant) >
      magnitude(found->second.determinant)) {
    const std::uint64_t selections = found->second.selections;
    found->second = Elite{
        state.core, state.determinant, hash, selections,
        statistics.moves, metrics};
    ++statistics.archive_replacements;
    return true;
  }
  ++statistics.archive_rejections;
  return false;
}

std::string qd_log_json(
    std::string_view event, double elapsed,
    const QdArguments& arguments, const QdStatistics& statistics,
    const std::map<Descriptor, Elite>& archive,
    std::uint64_t best_magnitude) {
  std::ostringstream output;
  output
      << "{\"event\":\"" << event << "\""
      << ",\"engine\":\"exact-core-map-elites-v1\""
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(6) << elapsed
      << ",\"seed\":" << arguments.seed
      << ",\"bursts\":" << statistics.bursts
      << ",\"moves\":" << statistics.moves
      << ",\"archive_cells\":" << archive.size()
      << ",\"archive_insertions\":"
      << statistics.archive_insertions
      << ",\"archive_replacements\":"
      << statistics.archive_replacements
      << ",\"exact_kicks\":" << statistics.exact_kicks
      << ",\"exact_kick_evaluations\":"
      << statistics.exact_kick_evaluations
      << ",\"descriptor_evaluations\":"
      << statistics.descriptor_evaluations
      << ",\"exact_bit_directions\":"
      << statistics.eo.exact_bit_directions
      << ",\"exact_pair_directions\":"
      << statistics.eo.exact_pair_directions
      << ",\"exact_complement_directions\":"
      << statistics.eo.exact_complement_directions
      << ",\"frontier_state_visits\":"
      << statistics.frontier_state_visits
      << ",\"distinct_frontier_cores\":"
      << statistics.distinct_frontier_cores
      << ",\"promotions\":" << statistics.promotions
      << ",\"best_core_quotient\":" << best_magnitude
      << ",\"best_absolute_determinant\":\""
      << wide_to_string(
             static_cast<Wide>(best_magnitude) * kScale)
      << "\"}\n";
  return output.str();
}

void write_qd_summary(
    const QdArguments& arguments, const QdStatistics& statistics,
    const std::map<Descriptor, Elite>& archive,
    double elapsed, std::uint64_t best_magnitude,
    std::uint64_t checkpoint_nonce, bool stopped_by_signal,
    const std::vector<fs::path>& archive_outputs) {
  const std::uint64_t exact_directions =
      statistics.eo.exact_bit_directions +
      statistics.eo.exact_pair_directions +
      statistics.eo.exact_complement_directions;
  std::ostringstream output;
  output
      << "{\n"
      << "  \"schema_version\": 1,\n"
      << "  \"engine\": \"exact-core-map-elites-v1\",\n"
      << "  \"complete\": "
      << (stopped_by_signal ? "false" : "true") << ",\n"
      << "  \"reason\": \""
      << (stopped_by_signal ? "signal" : "time_limit")
      << "\",\n"
      << "  \"stopped_by_signal\": "
      << (stopped_by_signal ? "true" : "false") << ",\n"
      << "  \"elapsed_seconds\": " << std::fixed
      << std::setprecision(6) << elapsed << ",\n"
      << "  \"seed\": " << arguments.seed << ",\n"
      << "  \"frontier\": \"" << arguments.frontier << "\",\n"
      << "  \"best_core_quotient\": " << best_magnitude
      << ",\n"
      << "  \"best_absolute_determinant\": \""
      << wide_to_string(
             static_cast<Wide>(best_magnitude) * kScale)
      << "\",\n"
      << "  \"above_frontier_unverified\": "
      << (static_cast<Wide>(best_magnitude) * kScale >
                  static_cast<Wide>(arguments.frontier)
              ? "true"
              : "false")
      << ",\n"
      << "  \"bursts\": " << statistics.bursts << ",\n"
      << "  \"moves\": " << statistics.moves << ",\n"
      << "  \"archive_cells\": " << archive.size() << ",\n"
      << "  \"archive_insertions\": "
      << statistics.archive_insertions << ",\n"
      << "  \"archive_replacements\": "
      << statistics.archive_replacements << ",\n"
      << "  \"archive_rejections\": "
      << statistics.archive_rejections << ",\n"
      << "  \"archive_capacity_rejections\": "
      << statistics.archive_capacity_rejections << ",\n"
      << "  \"parent_tournaments\": "
      << statistics.parent_tournaments << ",\n"
      << "  \"low_use_parent_choices\": "
      << statistics.low_use_parent_choices << ",\n"
      << "  \"quality_parent_choices\": "
      << statistics.quality_parent_choices << ",\n"
      << "  \"exact_kicks\": " << statistics.exact_kicks
      << ",\n"
      << "  \"exact_kick_evaluations\": "
      << statistics.exact_kick_evaluations << ",\n"
      << "  \"kick_failures\": "
      << statistics.kick_failures << ",\n"
      << "  \"exact_state_rebuilds\": "
      << statistics.exact_state_rebuilds << ",\n"
      << "  \"descriptor_evaluations\": "
      << statistics.descriptor_evaluations << ",\n"
      << "  \"exact_direction_evaluations\": "
      << exact_directions << ",\n"
      << "  \"exact_bit_directions\": "
      << statistics.eo.exact_bit_directions << ",\n"
      << "  \"exact_pair_directions\": "
      << statistics.eo.exact_pair_directions << ",\n"
      << "  \"exact_complement_directions\": "
      << statistics.eo.exact_complement_directions << ",\n"
      << "  \"frontier_state_visits\": "
      << statistics.frontier_state_visits << ",\n"
      << "  \"distinct_frontier_cores\": "
      << statistics.distinct_frontier_cores << ",\n"
      << "  \"promotions\": " << statistics.promotions
      << ",\n"
      << "  \"identity_checks\": "
      << statistics.identity_checks << ",\n"
      << "  \"exact_determinant_checks\": "
      << statistics.exact_determinant_checks << ",\n"
      << "  \"hash_checks\": " << statistics.hash_checks
      << ",\n"
      << "  \"maximum_nearest_seed_distance\": "
      << statistics.maximum_nearest_seed_distance << ",\n"
      << "  \"minimum_observed_gram_energy\": "
      << (statistics.minimum_observed_gram_energy ==
                  std::numeric_limits<std::uint64_t>::max()
              ? 0
              : statistics.minimum_observed_gram_energy)
      << ",\n"
      << "  \"maximum_observed_gram_energy\": "
      << statistics.maximum_observed_gram_energy << ",\n"
      << "  \"archive_output_count\": "
      << archive_outputs.size() << ",\n"
      << "  \"checkpoint_nonce\": " << checkpoint_nonce
      << ",\n"
      << "  \"descriptor\": {\n"
      << "    \"nearest_seed\": true,\n"
      << "    \"distance_bin_width\": "
      << arguments.distance_bin_width << ",\n"
      << "    \"balance_bin_width\": "
      << arguments.balance_bin_width << ",\n"
      << "    \"gram_bin_width\": "
      << arguments.gram_bin_width << "\n"
      << "  },\n"
      << "  \"search\": {\n"
      << "    \"tau\": " << arguments.tau << ",\n"
      << "    \"burst_moves\": " << arguments.burst_moves
      << ",\n"
      << "    \"archive_stride\": "
      << arguments.archive_stride << ",\n"
      << "    \"tournament_size\": "
      << arguments.tournament_size << ",\n"
      << "    \"kick_flip_range\": ["
      << arguments.kick_min_flips << ", "
      << arguments.kick_max_flips << "],\n"
      << "    \"kick_probability\": \""
      << arguments.kick_numerator << "/"
      << arguments.kick_denominator << "\",\n"
      << "    \"pair_interval\": "
      << arguments.pair_interval << ",\n"
      << "    \"complement_interval\": "
      << arguments.complement_interval << "\n"
      << "  },\n"
      << "  \"claim_boundary\": [\n"
      << "    \"All retained objectives are exact integer determinants.\",\n"
      << "    \"Direction evaluations repeat across states and are not unique matrices.\",\n"
      << "    \"Archive cells are structural search niches, not H/HT-equivalence classes.\",\n"
      << "    \"This stochastic pilot makes no optimality or literature-novelty claim.\"\n"
      << "  ]\n"
      << "}\n";
  atomic_write(
      arguments.summary, output.str(), checkpoint_nonce);
}

int run_qd(const QdArguments& arguments) {
  std::vector<State> seed_states;
  std::vector<CoreMatrix> seed_cores;
  seed_states.reserve(arguments.seeds.size());
  seed_cores.reserve(arguments.seeds.size());
  for (const fs::path& path : arguments.seeds) {
    seed_states.push_back(eo_initial_state(path));
    seed_cores.push_back(seed_states.back().core);
  }

  EoArguments differential;
  differential.start = arguments.seeds.front();
  differential.seed =
      arguments.seed ^ UINT64_C(0x51444d4150454c49);
  differential.tau = arguments.tau;
  differential.self_test_rounds =
      arguments.differential_rounds;
  run_eo_self_test(differential, seed_states.front());

  if (!arguments.log.parent_path().empty()) {
    fs::create_directories(arguments.log.parent_path());
  }
  if (!arguments.summary.parent_path().empty()) {
    fs::create_directories(arguments.summary.parent_path());
  }
  if (!arguments.output.parent_path().empty()) {
    fs::create_directories(arguments.output.parent_path());
  }
  fs::create_directories(arguments.archive_directory);
  std::ofstream log(arguments.log, std::ios::out | std::ios::trunc);
  if (!log) {
    throw std::runtime_error(
        "cannot create quality-diversity log");
  }

  std::mt19937_64 randomizer(arguments.seed);
  const auto zobrist = make_zobrist(arguments.seed);
  std::map<Descriptor, Elite> archive;
  std::vector<Descriptor> cells;
  QdStatistics statistics;
  std::unordered_set<std::string> frontier_cores;
  CoreMatrix best_core = seed_states.front().core;
  std::uint64_t best_magnitude =
      magnitude(seed_states.front().determinant);
  std::uint64_t checkpoint_nonce = 0;

  auto observe_frontier = [&](const State& state) {
    const Wide absolute =
        static_cast<Wide>(magnitude(state.determinant)) * kScale;
    if (absolute == static_cast<Wide>(arguments.frontier)) {
      ++statistics.frontier_state_visits;
      frontier_cores.insert(qd_core_key(state.core));
      statistics.distinct_frontier_cores =
          frontier_cores.size();
    }
  };

  for (std::size_t index = 0;
       index < seed_states.size(); ++index) {
    State& state = seed_states[index];
    const std::uint64_t hash =
        core_hash(state.core, zobrist);
    update_archive(
        state, hash, seed_cores, arguments,
        archive, cells, statistics);
    observe_frontier(state);
    if (magnitude(state.determinant) > best_magnitude) {
      best_magnitude = magnitude(state.determinant);
      best_core = state.core;
    }
  }
  atomic_write(
      arguments.output,
      sign_matrix_text(core_to_sign(best_core)),
      checkpoint_nonce++);

  RankSampler bit_sampler(kCoreEntries, arguments.tau);
  RankSampler pair_sampler(
      kMaximumPairCandidates, arguments.tau);
  RankSampler complement_sampler(
      kComplementMoveCount, arguments.tau);
  std::vector<RankedMove> ranked_moves;
  ranked_moves.reserve(kMaximumPairCandidates);

  const auto started = Clock::now();
  const auto deadline =
      started + std::chrono::duration<double>(
                    arguments.seconds);
  auto next_heartbeat =
      started + std::chrono::duration<double>(
                    arguments.heartbeat_seconds);
  std::signal(SIGINT, request_stop);
  std::signal(SIGTERM, request_stop);

  auto elapsed_now = [&]() {
    return std::chrono::duration<double>(
               Clock::now() - started)
        .count();
  };
  log << qd_log_json(
      "started", 0.0, arguments, statistics,
      archive, best_magnitude);
  log.flush();

  while (
      stop_requested == 0 && Clock::now() < deadline) {
    const Descriptor parent_cell = choose_parent_cell(
        archive, cells, arguments, randomizer, statistics);
    Elite parent = archive.at(parent_cell);
    ++archive.at(parent_cell).selections;
    State state =
        state_from_core(parent.core, statistics);
    std::uint64_t hash =
        core_hash(state.core, zobrist);

    if (
        randomizer() %
                static_cast<std::uint64_t>(
                    arguments.kick_denominator) <
            static_cast<std::uint64_t>(
                arguments.kick_numerator)) {
      const bool destroyed = destroy_state(
          state, arguments, randomizer, statistics);
      hash = core_hash(state.core, zobrist);
      if (destroyed) {
        observe_frontier(state);
        const std::uint64_t kicked_magnitude =
            magnitude(state.determinant);
        if (kicked_magnitude > best_magnitude) {
          best_magnitude = kicked_magnitude;
          best_core = state.core;
          ++statistics.promotions;
          atomic_write(
              arguments.output,
              sign_matrix_text(core_to_sign(best_core)),
              checkpoint_nonce++);
        }
        update_archive(
            state, hash, seed_cores, arguments,
            archive, cells, statistics);
      }
    }

    std::array<std::uint64_t, kCoreEntries> last_moved{};
    last_moved.fill(kNeverMoved);
    for (int within = 0;
         within < arguments.burst_moves; ++within) {
      if (stop_requested != 0 || Clock::now() >= deadline) {
        break;
      }
      ++statistics.moves;
      EoMove move;
      if (
          statistics.moves %
                  static_cast<std::uint64_t>(
                      arguments.pair_interval) ==
              0) {
        move = choose_balanced_pair_move(
            state, ranked_moves, last_moved,
            statistics.moves, arguments.recency_window,
            arguments.recency_penalty, pair_sampler,
            randomizer, statistics.eo);
      } else if (
          statistics.moves %
                  static_cast<std::uint64_t>(
                      arguments.complement_interval) ==
              0) {
        move = choose_complement_move(
            state, ranked_moves, last_moved,
            statistics.moves, arguments.recency_window,
            arguments.recency_penalty,
            complement_sampler, randomizer, statistics.eo);
      } else {
        move = choose_bit_move(
            state, ranked_moves, last_moved,
            statistics.moves, arguments.recency_window,
            arguments.recency_penalty, bit_sampler,
            randomizer, statistics.eo);
      }
      const std::uint64_t before =
          magnitude(state.determinant);
      mark_recent_bits(
          move, statistics.moves, last_moved,
          statistics.eo);
      apply_eo_move(
          state, move, hash, zobrist, statistics.eo);
      const std::uint64_t after =
          magnitude(state.determinant);
      if (after > before) {
        ++statistics.eo.uphill_moves;
      } else if (after < before) {
        ++statistics.eo.downhill_moves;
      } else {
        ++statistics.eo.level_moves;
      }

      observe_frontier(state);
      if (after > best_magnitude) {
        best_magnitude = after;
        best_core = state.core;
        ++statistics.promotions;
        atomic_write(
            arguments.output,
            sign_matrix_text(core_to_sign(best_core)),
            checkpoint_nonce++);
      }

      if (
          statistics.moves %
                  static_cast<std::uint64_t>(
                      arguments.archive_stride) ==
              0) {
        update_archive(
            state, hash, seed_cores, arguments,
            archive, cells, statistics);
      }
      if ((statistics.moves & UINT64_C(8191)) == 0) {
        const std::int64_t exact =
            exact_core_determinant(state.core);
        ++statistics.exact_determinant_checks;
        if (exact != state.determinant) {
          throw std::runtime_error(
              "quality-diversity determinant invariant failed");
        }
        check_adjugate_identity(state);
        ++statistics.identity_checks;
        if (hash != core_hash(state.core, zobrist)) {
          throw std::runtime_error(
              "quality-diversity hash invariant failed");
        }
        ++statistics.hash_checks;
      }
    }
    ++statistics.bursts;
    update_archive(
        state, hash, seed_cores, arguments,
        archive, cells, statistics);

    const auto now = Clock::now();
    if (now >= next_heartbeat) {
      log << qd_log_json(
          "heartbeat", elapsed_now(), arguments,
          statistics, archive, best_magnitude);
      log.flush();
      next_heartbeat =
          now + std::chrono::duration<double>(
                    arguments.heartbeat_seconds);
    }
  }

  atomic_write(
      arguments.output,
      sign_matrix_text(core_to_sign(best_core)),
      checkpoint_nonce++);

  std::vector<const Elite*> ordered;
  ordered.reserve(archive.size());
  for (const auto& [cell, elite] : archive) {
    (void)cell;
    ordered.push_back(&elite);
  }
  std::sort(
      ordered.begin(), ordered.end(),
      [](const Elite* first, const Elite* second) {
        if (magnitude(first->determinant) !=
            magnitude(second->determinant)) {
          return magnitude(first->determinant) >
                 magnitude(second->determinant);
        }
        return first->metrics.cell.as_tuple() <
               second->metrics.cell.as_tuple();
      });
  const std::size_t output_count = std::min(
      ordered.size(),
      static_cast<std::size_t>(arguments.archive_outputs));
  std::vector<fs::path> archive_outputs;
  archive_outputs.reserve(output_count);
  for (std::size_t index = 0; index < output_count; ++index) {
    const Elite& elite = *ordered[index];
    std::ostringstream name;
    name << "elite-" << std::setw(3) << std::setfill('0')
         << index << "-q" << magnitude(elite.determinant)
         << "-n" << elite.metrics.cell.nearest_seed
         << "-d" << elite.metrics.cell.distance_bin
         << "-b" << elite.metrics.cell.balance_bin
         << "-g" << elite.metrics.cell.gram_bin
         << ".matrix.txt";
    const fs::path path =
        arguments.archive_directory / name.str();
    atomic_write(
        path,
        sign_matrix_text(core_to_sign(elite.core)),
        static_cast<std::uint64_t>(index));
    archive_outputs.push_back(path);
  }

  const double elapsed = elapsed_now();
  const bool stopped_by_signal = stop_requested != 0;
  write_qd_summary(
      arguments, statistics, archive, elapsed,
      best_magnitude, checkpoint_nonce, stopped_by_signal,
      archive_outputs);
  log << qd_log_json(
      stopped_by_signal ? "stopped" : "finished",
      elapsed, arguments, statistics,
      archive, best_magnitude);
  log.flush();

  std::cout
      << "finished elapsed=" << std::fixed
      << std::setprecision(6) << elapsed
      << " moves=" << statistics.moves
      << " cells=" << archive.size()
      << " exact_directions="
      << statistics.eo.exact_bit_directions +
             statistics.eo.exact_pair_directions +
             statistics.eo.exact_complement_directions
      << " best="
      << wide_to_string(
             static_cast<Wide>(best_magnitude) * kScale)
      << " promotions=" << statistics.promotions << '\n';
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    return run_qd(parse_qd_arguments(argc, argv));
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
