// Exact linkage-learning / gene-pool optimal-mixing pilot for order-23
// maximal determinant search.
//
// The 23x23 sign matrix is represented by its dephased 22x22 binary core.
// A population is initialized from the 24 known neutral-network frontier
// factors, one explicitly aligned additional H-class frontier factor, and
// exact-polished perturbations.  Each generation learns a normalized-mutual-
// information linkage tree, augments it with the twelve known intra-class
// neutral supports, the exact transported reference-to-H2 12-entry support,
// and row/column blocks, and performs donor mixing over that family.
//
// Every proposed offspring is scored by exact Bareiss determinant arithmetic.
// A mix is accepted only if absolute determinant is nondecreasing.  Exact
// single-core-bit polishing follows each offspring.  A bounded archive keeps
// high-scoring, mutually distant cores available as donors.
//
// This is an experimental standalone solver.  Any retained matrix still
// requires ./arena verify; equivalence or novelty claims require a separate
// pinned H/HT audit.

#include <map>
#include <queue>
#include <set>
#include <tuple>
#include <unordered_set>

#define main core_adjugate_tabu_embedded_main
#include "core_adjugate_tabu.cpp"
#undef main

namespace {

constexpr std::uint64_t kFrontier =
    UINT64_C(2779447296000000);
constexpr std::uint64_t kAuditPrime = UINT64_C(1000000007);

struct GomeaArguments {
  fs::path seed_directory;
  std::vector<fs::path> extra_seeds;
  fs::path bridge_reference;
  fs::path bridge_seed;
  fs::path alignment_metadata;
  fs::path output;
  fs::path archive_directory;
  fs::path log;
  fs::path summary;
  std::uint64_t seed = 33001;
  std::uint64_t frontier = kFrontier;
  double seconds = 300.0;
  double heartbeat_seconds = 15.0;
  int expected_seed_count = 25;
  int population_size = 32;
  int archive_size = 128;
  int initial_kick_min = 8;
  int initial_kick_max = 24;
  int linkage_refresh = 4;
  int maximum_linkage_size = 128;
  int polish_moves = 4;
  int forced_improvement_generations = 8;
  int differential_samples = 64;
};

struct Individual {
  CoreMatrix core{};
  std::int64_t determinant = 0;
  int stagnant_generations = 0;
};

struct ArchiveEntry {
  CoreMatrix core{};
  std::int64_t determinant = 0;
};

struct LinkageSubset {
  std::vector<int> variables;
  std::string source;
};

struct BridgeSupportInfo {
  LinkageSubset linkage;
  std::vector<std::pair<int, int>> raw_coordinates;
  CoreMatrix reference_core{};
  CoreMatrix bridge_core{};
  std::int64_t reference_determinant = 0;
  std::int64_t bridge_determinant = 0;
};

struct SearchStatistics {
  std::uint64_t generations = 0;
  std::uint64_t linkage_rebuilds = 0;
  std::uint64_t exact_evaluations = 0;
  std::uint64_t mix_attempts = 0;
  std::uint64_t mix_changed = 0;
  std::uint64_t accepted_improvements = 0;
  std::uint64_t accepted_equals = 0;
  std::uint64_t rejected_equal_duplicates = 0;
  std::uint64_t rejected_improving_duplicates = 0;
  std::uint64_t rejected_polish_duplicates = 0;
  std::uint64_t rejected_decreases = 0;
  std::uint64_t rejected_singular = 0;
  std::uint64_t polish_attempts = 0;
  std::uint64_t polish_moves = 0;
  std::uint64_t forced_improvements = 0;
  std::uint64_t forced_replacements = 0;
  std::uint64_t diversity_resets = 0;
  std::uint64_t bridge_endpoint_donations = 0;
  std::uint64_t archive_insertions = 0;
  std::uint64_t archive_replacements = 0;
  std::uint64_t promotions = 0;
  std::uint64_t differential_checks = 0;
  std::map<std::string, std::uint64_t> attempted_by_source;
  std::map<std::string, std::uint64_t> accepted_by_source;
};

struct GomeaLogger {
  std::ofstream stream;

  explicit GomeaLogger(const fs::path& path) {
    if (fs::exists(path)) {
      throw std::runtime_error(
          "refusing to overwrite log: " + path.string());
    }
    if (!path.parent_path().empty()) {
      fs::create_directories(path.parent_path());
    }
    stream.open(path, std::ios::out | std::ios::trunc);
    if (!stream) {
      throw std::runtime_error("cannot create log: " + path.string());
    }
  }

  void line(const std::string& record) {
    stream << record << '\n';
    stream.flush();
    if (!stream) {
      throw std::runtime_error("failed writing GOMEA JSONL log");
    }
  }
};

int strict_positive_int(std::string_view text,
                        std::string_view option,
                        int maximum) {
  const std::uint64_t parsed = strict_unsigned(text, option);
  if (parsed == 0 ||
      parsed > static_cast<std::uint64_t>(maximum)) {
    throw std::runtime_error(
        std::string(option) + " is outside its valid range");
  }
  return static_cast<int>(parsed);
}

GomeaArguments parse_gomea_arguments(int argc, char** argv) {
  GomeaArguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    auto value = [&]() -> std::string_view {
      ++index;
      if (index >= argc) {
        throw std::runtime_error("missing value for " + option);
      }
      return argv[index];
    };
    if (option == "--seed-dir") {
      arguments.seed_directory = value();
    } else if (option == "--seed-matrix") {
      arguments.extra_seeds.emplace_back(value());
    } else if (option == "--bridge-reference") {
      arguments.bridge_reference = value();
    } else if (option == "--bridge-seed") {
      arguments.bridge_seed = value();
    } else if (option == "--alignment-metadata") {
      arguments.alignment_metadata = value();
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
    } else if (option == "--heartbeat-seconds" ||
               option == "--heartbeat") {
      arguments.heartbeat_seconds =
          strict_double(value(), option, true);
    } else if (option == "--expected-seeds") {
      arguments.expected_seed_count =
          strict_positive_int(value(), option, 512);
    } else if (option == "--population") {
      arguments.population_size =
          strict_positive_int(value(), option, 512);
    } else if (option == "--archive-size") {
      arguments.archive_size =
          strict_positive_int(value(), option, 4096);
    } else if (option == "--initial-kick-min") {
      arguments.initial_kick_min =
          strict_positive_int(value(), option, kCoreEntries);
    } else if (option == "--initial-kick-max") {
      arguments.initial_kick_max =
          strict_positive_int(value(), option, kCoreEntries);
    } else if (option == "--linkage-refresh") {
      arguments.linkage_refresh =
          strict_positive_int(value(), option, 1000000);
    } else if (option == "--max-linkage-size") {
      arguments.maximum_linkage_size =
          strict_positive_int(value(), option, kCoreEntries);
    } else if (option == "--polish-moves") {
      arguments.polish_moves =
          strict_positive_int(value(), option, 1000);
    } else if (option == "--forced-improvement-generations") {
      arguments.forced_improvement_generations =
          strict_positive_int(value(), option, 1000000);
    } else if (option == "--differential-samples") {
      const std::uint64_t parsed = strict_unsigned(value(), option);
      if (parsed >
          static_cast<std::uint64_t>(
              std::numeric_limits<int>::max())) {
        throw std::runtime_error(
            "--differential-samples is too large");
      }
      arguments.differential_samples = static_cast<int>(parsed);
    } else if (option == "--help") {
      std::cout
          << "usage: core_linkage_gomea --seed-dir DIR "
             "--seed-matrix MATRIX --output MATRIX "
             "--archive-dir DIR --log JSONL --summary JSON [options]\n"
          << "  --bridge-reference MATRIX --bridge-seed MATRIX\n"
          << "  --alignment-metadata JSON\n"
          << "  --seed N --frontier N --seconds S "
             "--heartbeat-seconds S\n"
          << "  --expected-seeds N --population N --archive-size N\n"
          << "  --initial-kick-min N --initial-kick-max N\n"
          << "  --linkage-refresh N --max-linkage-size N\n"
          << "  --polish-moves N "
             "--forced-improvement-generations N\n"
          << "  --differential-samples N\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }

  if (arguments.seed_directory.empty() ||
      arguments.bridge_reference.empty() ||
      arguments.bridge_seed.empty() ||
      arguments.alignment_metadata.empty() ||
      arguments.output.empty() ||
      arguments.archive_directory.empty() ||
      arguments.log.empty() || arguments.summary.empty()) {
    throw std::runtime_error(
        "--seed-dir, --bridge-reference, --bridge-seed, "
        "--alignment-metadata, --output, --archive-dir, --log, "
        "and --summary are required");
  }
  if (!fs::is_regular_file(arguments.alignment_metadata)) {
    throw std::runtime_error(
        "alignment metadata does not exist: " +
        arguments.alignment_metadata.string());
  }
  const std::array<fs::path, 4> outputs{
      arguments.output, arguments.archive_directory,
      arguments.log, arguments.summary};
  for (std::size_t first = 0; first < outputs.size(); ++first) {
    for (std::size_t second = first + 1;
         second < outputs.size(); ++second) {
      if (outputs[first] == outputs[second]) {
        throw std::runtime_error(
            "output, archive, log, and summary paths must differ");
      }
    }
  }
  if (arguments.population_size <
      arguments.expected_seed_count) {
    throw std::runtime_error(
        "population cannot be smaller than expected seed count");
  }
  if (arguments.archive_size <
      arguments.expected_seed_count) {
    throw std::runtime_error(
        "archive cannot be smaller than expected seed count");
  }
  if (arguments.initial_kick_min >
      arguments.initial_kick_max) {
    throw std::runtime_error(
        "initial kick minimum exceeds maximum");
  }
  return arguments;
}

std::vector<fs::path> collect_seed_paths(
    const GomeaArguments& arguments) {
  if (!fs::is_directory(arguments.seed_directory)) {
    throw std::runtime_error(
        "seed directory does not exist: " +
        arguments.seed_directory.string());
  }
  std::vector<fs::path> paths;
  for (const fs::directory_entry& entry :
       fs::directory_iterator(arguments.seed_directory)) {
    if (!entry.is_regular_file()) {
      continue;
    }
    const std::string filename =
        entry.path().filename().string();
    if (filename.starts_with("tie-") &&
        filename.ends_with(".matrix.txt")) {
      paths.push_back(entry.path());
    }
  }
  std::sort(paths.begin(), paths.end());
  paths.insert(paths.end(), arguments.extra_seeds.begin(),
               arguments.extra_seeds.end());
  if (static_cast<int>(paths.size()) !=
      arguments.expected_seed_count) {
    throw std::runtime_error(
        "expected " +
        std::to_string(arguments.expected_seed_count) +
        " seed matrices but collected " +
        std::to_string(paths.size()));
  }
  return paths;
}

std::string core_key(const CoreMatrix& core) {
  std::string key;
  key.reserve(kCoreEntries);
  for (const auto& row : core) {
    for (const std::uint8_t entry : row) {
      key.push_back(entry == 0U ? '0' : '1');
    }
  }
  return key;
}

std::uint64_t core_fingerprint(const CoreMatrix& core) {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  for (const auto& row : core) {
    for (const std::uint8_t entry : row) {
      hash ^= static_cast<std::uint64_t>(entry + 1U);
      hash *= UINT64_C(1099511628211);
    }
  }
  return hash;
}

std::string hex_u64(std::uint64_t value) {
  std::ostringstream output;
  output << std::hex << std::setw(16) << std::setfill('0') << value;
  return output.str();
}

int hamming_distance(const CoreMatrix& left,
                     const CoreMatrix& right) {
  int distance = 0;
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      distance += left[row][column] != right[row][column] ? 1 : 0;
    }
  }
  return distance;
}

std::vector<int> difference_subset(const CoreMatrix& left,
                                   const CoreMatrix& right) {
  std::vector<int> result;
  for (int flat = 0; flat < kCoreEntries; ++flat) {
    if (left[flat / kCoreOrder][flat % kCoreOrder] !=
        right[flat / kCoreOrder][flat % kCoreOrder]) {
      result.push_back(flat);
    }
  }
  return result;
}

std::string subset_key(const std::vector<int>& variables) {
  std::string key(
      static_cast<std::size_t>((kCoreEntries + 7) / 8), '\0');
  for (const int variable : variables) {
    const std::size_t byte =
        static_cast<std::size_t>(variable / 8);
    const unsigned shift =
        static_cast<unsigned>(variable % 8);
    key[byte] = static_cast<char>(
        static_cast<unsigned char>(key[byte]) |
        static_cast<unsigned char>(1U << shift));
  }
  return key;
}

std::uint64_t modular_power(std::uint64_t base,
                            std::uint64_t exponent) {
  std::uint64_t result = 1;
  while (exponent != 0) {
    if ((exponent & 1U) != 0) {
      result =
          static_cast<std::uint64_t>(
              (static_cast<__uint128_t>(result) * base) %
              kAuditPrime);
    }
    base =
        static_cast<std::uint64_t>(
            (static_cast<__uint128_t>(base) * base) %
            kAuditPrime);
    exponent >>= 1;
  }
  return result;
}

std::uint64_t core_determinant_mod(const CoreMatrix& core) {
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
    determinant =
        static_cast<std::uint64_t>(
            (static_cast<__uint128_t>(determinant) *
             pivot_value) %
            kAuditPrime);
    const std::uint64_t inverse =
        modular_power(pivot_value, kAuditPrime - 2);
    for (int row = column + 1; row < kCoreOrder; ++row) {
      if (work[row][column] == 0) {
        continue;
      }
      const std::uint64_t factor =
          static_cast<std::uint64_t>(
              (static_cast<__uint128_t>(
                   work[row][column]) *
               inverse) %
              kAuditPrime);
      for (int inner = column; inner < kCoreOrder; ++inner) {
        const std::uint64_t product =
            static_cast<std::uint64_t>(
                (static_cast<__uint128_t>(factor) *
                 work[column][inner]) %
                kAuditPrime);
        work[row][inner] =
            work[row][inner] >= product
                ? work[row][inner] - product
                : work[row][inner] + kAuditPrime - product;
      }
    }
  }
  if (negative && determinant != 0) {
    determinant = kAuditPrime - determinant;
  }
  return determinant;
}

std::uint64_t signed_mod_prime(std::int64_t value) {
  if (value >= 0) {
    return static_cast<std::uint64_t>(value) % kAuditPrime;
  }
  const std::uint64_t magnitude_value =
      static_cast<std::uint64_t>(-value) % kAuditPrime;
  return magnitude_value == 0 ? 0 : kAuditPrime - magnitude_value;
}

std::vector<std::vector<Wide>> sign_as_wide(
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

void differential_check_core(const CoreMatrix& core,
                             std::int64_t determinant,
                             std::string_view context,
                             SearchStatistics& statistics) {
  const std::int64_t bareiss_determinant =
      exact_core_determinant(core);
  if (bareiss_determinant != determinant) {
    throw std::runtime_error(
        std::string(context) +
        ": cached and Bareiss core determinants differ");
  }
  if (core_determinant_mod(core) !=
      signed_mod_prime(determinant)) {
    throw std::runtime_error(
        std::string(context) +
        ": modular and Bareiss core determinants differ");
  }
  const SignMatrix sign = core_to_sign(core);
  const Wide sign_determinant = bareiss(sign_as_wide(sign));
  if (sign_determinant !=
      static_cast<Wide>(determinant) * kScale) {
    throw std::runtime_error(
        std::string(context) +
        ": sign/core determinant quotient differs");
  }
  ++statistics.differential_checks;
}

std::vector<LinkageSubset> derive_known_neutral_supports(
    const std::vector<Individual>& seeds) {
  if (seeds.size() < 24) {
    throw std::runtime_error(
        "at least 24 neutral seeds are required");
  }
  std::map<std::string, std::vector<int>> unique;
  std::size_t core_12 = 0;
  std::size_t core_31 = 0;
  for (std::size_t first = 0; first < 24; ++first) {
    for (std::size_t second = first + 1; second < 24; ++second) {
      std::vector<int> subset =
          difference_subset(seeds[first].core, seeds[second].core);
      // All twelve original sign-matrix generators flip 12 entries.
      // Dephasing preserves six as 12 core bits.  The six generators
      // that touch the gauge border induce 31-bit core supports.
      if (subset.size() == 12 || subset.size() == 31) {
        unique.emplace(subset_key(subset), std::move(subset));
      }
    }
  }
  if (unique.size() != 12) {
    throw std::runtime_error(
        "expected exactly 12 distinct induced supports from the "
        "known raw 12-entry neutral generators, found " +
        std::to_string(unique.size()));
  }
  std::vector<LinkageSubset> result;
  for (auto& [key, variables] : unique) {
    static_cast<void>(key);
    if (variables.size() == 12) {
      ++core_12;
    } else if (variables.size() == 31) {
      ++core_31;
    }
    result.push_back(
        LinkageSubset{
            std::move(variables),
            "known-neutral-raw12-induced-core"});
  }
  if (core_12 != 6 || core_31 != 6) {
    throw std::runtime_error(
        "expected the known neutral generators to induce six "
        "12-bit and six 31-bit core supports");
  }
  return result;
}

BridgeSupportInfo derive_bridge_support(
    const GomeaArguments& arguments,
    const std::vector<Individual>& seeds,
    SearchStatistics& statistics) {
  const SignMatrix reference =
      read_sign_matrix(arguments.bridge_reference);
  const SignMatrix bridge_seed =
      read_sign_matrix(arguments.bridge_seed);
  std::vector<std::pair<int, int>> raw_coordinates;
  for (int row = 0; row < kSignOrder; ++row) {
    for (int column = 0; column < kSignOrder; ++column) {
      if (reference[row][column] != bridge_seed[row][column]) {
        raw_coordinates.emplace_back(row, column);
      }
    }
  }
  if (raw_coordinates.size() != 12) {
    throw std::runtime_error(
        "transported reference and aligned bridge seed must differ "
        "in exactly 12 sign entries, found " +
        std::to_string(raw_coordinates.size()));
  }

  const CoreMatrix reference_core = dephase_to_core(reference);
  const CoreMatrix bridge_core = dephase_to_core(bridge_seed);
  const std::int64_t reference_determinant =
      exact_core_determinant(reference_core);
  const std::int64_t bridge_determinant =
      exact_core_determinant(bridge_core);
  differential_check_core(
      reference_core, reference_determinant,
      "transported bridge reference", statistics);
  differential_check_core(
      bridge_core, bridge_determinant,
      "aligned bridge seed", statistics);
  if (magnitude(reference_determinant) *
              static_cast<std::uint64_t>(kScale) !=
          arguments.frontier ||
      magnitude(bridge_determinant) *
              static_cast<std::uint64_t>(kScale) !=
          arguments.frontier) {
    throw std::runtime_error(
        "bridge endpoints are not both at the configured frontier");
  }
  bool seed_is_present = false;
  for (const Individual& seed : seeds) {
    if (seed.core == bridge_core) {
      seed_is_present = true;
      break;
    }
  }
  if (!seed_is_present) {
    throw std::runtime_error(
        "--bridge-seed is not present among the population seeds");
  }

  std::vector<int> core_variables =
      difference_subset(reference_core, bridge_core);
  if (core_variables.size() != 12) {
    throw std::runtime_error(
        "transported bridge must remain a 12-variable support after "
        "dephasing, found " +
        std::to_string(core_variables.size()));
  }
  return BridgeSupportInfo{
      LinkageSubset{
          std::move(core_variables),
          "transported-reference-h2-neutral-12"},
      std::move(raw_coordinates),
      reference_core,
      bridge_core,
      reference_determinant,
      bridge_determinant};
}

std::string coordinate_pairs_json(
    const std::vector<std::pair<int, int>>& coordinates) {
  std::ostringstream output;
  output << '[';
  for (std::size_t index = 0; index < coordinates.size(); ++index) {
    if (index != 0) {
      output << ',';
    }
    output << '[' << coordinates[index].first + 1 << ','
           << coordinates[index].second + 1 << ']';
  }
  output << ']';
  return output.str();
}

double binary_entropy(double probability) {
  if (probability <= 0.0 || probability >= 1.0) {
    return 0.0;
  }
  return -probability * std::log(probability) -
         (1.0 - probability) * std::log(1.0 - probability);
}

double normalized_mutual_information(
    const std::vector<Individual>& population,
    int first, int second) {
  std::array<std::uint64_t, 4> counts{};
  for (const Individual& individual : population) {
    const int left =
        individual.core[first / kCoreOrder]
                       [first % kCoreOrder] != 0U
            ? 1
            : 0;
    const int right =
        individual.core[second / kCoreOrder]
                       [second % kCoreOrder] != 0U
            ? 1
            : 0;
    ++counts[static_cast<std::size_t>(2 * left + right)];
  }
  const double total = static_cast<double>(population.size());
  const double first_one =
      static_cast<double>(counts[2] + counts[3]) / total;
  const double second_one =
      static_cast<double>(counts[1] + counts[3]) / total;
  const double first_entropy = binary_entropy(first_one);
  const double second_entropy = binary_entropy(second_one);
  if (first_entropy == 0.0 || second_entropy == 0.0) {
    return 0.0;
  }
  double information = 0.0;
  for (int left = 0; left <= 1; ++left) {
    for (int right = 0; right <= 1; ++right) {
      const double joint =
          static_cast<double>(
              counts[static_cast<std::size_t>(2 * left + right)]) /
          total;
      if (joint == 0.0) {
        continue;
      }
      const double left_probability =
          left == 1 ? first_one : 1.0 - first_one;
      const double right_probability =
          right == 1 ? second_one : 1.0 - second_one;
      information +=
          joint * std::log(
                      joint /
                      (left_probability * right_probability));
    }
  }
  return information /
         std::sqrt(first_entropy * second_entropy);
}

struct ClusterEdge {
  double similarity = 0.0;
  int first = -1;
  int second = -1;
};

struct ClusterEdgeLess {
  bool operator()(const ClusterEdge& left,
                  const ClusterEdge& right) const {
    if (left.similarity != right.similarity) {
      return left.similarity < right.similarity;
    }
    if (left.first != right.first) {
      return left.first > right.first;
    }
    return left.second > right.second;
  }
};

void add_unique_subset(
    std::vector<LinkageSubset>& subsets,
    std::unordered_set<std::string>& seen,
    std::vector<int> variables, std::string source) {
  std::sort(variables.begin(), variables.end());
  variables.erase(
      std::unique(variables.begin(), variables.end()),
      variables.end());
  if (variables.empty()) {
    return;
  }
  const std::string key = subset_key(variables);
  if (seen.insert(key).second) {
    subsets.push_back(
        LinkageSubset{std::move(variables), std::move(source)});
  }
}

std::vector<LinkageSubset> learn_linkage(
    const std::vector<Individual>& population,
    const std::vector<LinkageSubset>& known,
    int maximum_linkage_size) {
  const int maximum_clusters = 2 * kCoreEntries - 1;
  std::vector<std::vector<int>> clusters(
      static_cast<std::size_t>(maximum_clusters));
  std::vector<int> cluster_sizes(
      static_cast<std::size_t>(maximum_clusters), 0);
  std::vector<bool> active(
      static_cast<std::size_t>(maximum_clusters), false);
  std::vector<std::vector<double>> similarity(
      static_cast<std::size_t>(maximum_clusters),
      std::vector<double>(
          static_cast<std::size_t>(maximum_clusters), 0.0));
  std::priority_queue<ClusterEdge,
                      std::vector<ClusterEdge>,
                      ClusterEdgeLess>
      queue;

  for (int variable = 0; variable < kCoreEntries; ++variable) {
    clusters[static_cast<std::size_t>(variable)] = {variable};
    cluster_sizes[static_cast<std::size_t>(variable)] = 1;
    active[static_cast<std::size_t>(variable)] = true;
  }
  for (int first = 0; first < kCoreEntries; ++first) {
    for (int second = first + 1; second < kCoreEntries; ++second) {
      const double value = normalized_mutual_information(
          population, first, second);
      similarity[static_cast<std::size_t>(first)]
                [static_cast<std::size_t>(second)] = value;
      similarity[static_cast<std::size_t>(second)]
                [static_cast<std::size_t>(first)] = value;
      queue.push(ClusterEdge{value, first, second});
    }
  }

  std::vector<LinkageSubset> result;
  std::unordered_set<std::string> seen;
  result.reserve(static_cast<std::size_t>(2 * kCoreEntries + 64));
  for (int variable = 0; variable < kCoreEntries; ++variable) {
    add_unique_subset(result, seen, {variable}, "singleton");
  }
  for (const LinkageSubset& subset : known) {
    add_unique_subset(result, seen, subset.variables, subset.source);
  }
  for (int row = 0; row < kCoreOrder; ++row) {
    std::vector<int> variables;
    for (int column = 0; column < kCoreOrder; ++column) {
      variables.push_back(row * kCoreOrder + column);
    }
    add_unique_subset(
        result, seen, std::move(variables), "core-row");
  }
  for (int column = 0; column < kCoreOrder; ++column) {
    std::vector<int> variables;
    for (int row = 0; row < kCoreOrder; ++row) {
      variables.push_back(row * kCoreOrder + column);
    }
    add_unique_subset(
        result, seen, std::move(variables), "core-column");
  }

  int next_cluster = kCoreEntries;
  int active_count = kCoreEntries;
  while (active_count > 1) {
    ClusterEdge edge;
    bool found = false;
    while (!queue.empty()) {
      edge = queue.top();
      queue.pop();
      if (active[static_cast<std::size_t>(edge.first)] &&
          active[static_cast<std::size_t>(edge.second)]) {
        found = true;
        break;
      }
    }
    if (!found) {
      throw std::runtime_error("linkage priority queue exhausted");
    }
    const int merged = next_cluster++;
    std::vector<int>& merged_variables =
        clusters[static_cast<std::size_t>(merged)];
    merged_variables =
        clusters[static_cast<std::size_t>(edge.first)];
    merged_variables.insert(
        merged_variables.end(),
        clusters[static_cast<std::size_t>(edge.second)].begin(),
        clusters[static_cast<std::size_t>(edge.second)].end());
    std::sort(merged_variables.begin(), merged_variables.end());
    const int first_size =
        cluster_sizes[static_cast<std::size_t>(edge.first)];
    const int second_size =
        cluster_sizes[static_cast<std::size_t>(edge.second)];
    cluster_sizes[static_cast<std::size_t>(merged)] =
        first_size + second_size;
    active[static_cast<std::size_t>(edge.first)] = false;
    active[static_cast<std::size_t>(edge.second)] = false;
    active[static_cast<std::size_t>(merged)] = true;
    --active_count;

    if (static_cast<int>(merged_variables.size()) <=
            maximum_linkage_size &&
        static_cast<int>(merged_variables.size()) <
            kCoreEntries) {
      add_unique_subset(
          result, seen, merged_variables, "learned-linkage");
    }
    for (int other = 0; other < merged; ++other) {
      if (!active[static_cast<std::size_t>(other)]) {
        continue;
      }
      const double first_similarity =
          similarity[static_cast<std::size_t>(edge.first)]
                    [static_cast<std::size_t>(other)];
      const double second_similarity =
          similarity[static_cast<std::size_t>(edge.second)]
                    [static_cast<std::size_t>(other)];
      const double merged_similarity =
          (static_cast<double>(first_size) * first_similarity +
           static_cast<double>(second_size) * second_similarity) /
          static_cast<double>(first_size + second_size);
      similarity[static_cast<std::size_t>(merged)]
                [static_cast<std::size_t>(other)] =
          merged_similarity;
      similarity[static_cast<std::size_t>(other)]
                [static_cast<std::size_t>(merged)] =
          merged_similarity;
      queue.push(ClusterEdge{
          merged_similarity, std::min(merged, other),
          std::max(merged, other)});
    }
  }
  return result;
}

bool copy_subset(CoreMatrix& destination,
                 const CoreMatrix& donor,
                 const std::vector<int>& variables) {
  bool changed = false;
  for (const int flat : variables) {
    std::uint8_t& entry =
        destination[flat / kCoreOrder][flat % kCoreOrder];
    const std::uint8_t donor_entry =
        donor[flat / kCoreOrder][flat % kCoreOrder];
    if (entry != donor_entry) {
      entry = donor_entry;
      changed = true;
    }
  }
  return changed;
}

bool exact_polish(Individual& individual, int maximum_moves,
                  SearchStatistics& statistics) {
  if (individual.determinant == 0 || maximum_moves <= 0) {
    return false;
  }
  State state;
  state.core = individual.core;
  state.determinant = individual.determinant;
  state.adjugate = exact_adjugate(state.core);
  const auto zobrist = make_zobrist(UINT64_C(0x474f4d4541504f4c));
  std::uint64_t hash = core_hash(state.core, zobrist);
  Statistics core_statistics;
  bool improved = false;
  for (int accepted = 0; accepted < maximum_moves; ++accepted) {
    Move best;
    std::uint64_t best_score = magnitude(state.determinant);
    for (int row = 0; row < kCoreOrder; ++row) {
      for (int column = 0; column < kCoreOrder; ++column) {
        ++statistics.polish_attempts;
        const std::int64_t candidate =
            bit_candidate_determinant(state, row, column);
        const std::uint64_t candidate_score = magnitude(candidate);
        if (candidate != 0 && candidate_score > best_score) {
          best = Move{
              MoveKind::kBit, row, column,
              row * kCoreOrder + column, candidate, false};
          best_score = candidate_score;
        }
      }
    }
    if (best.id < 0) {
      break;
    }
    apply_move(state, best, hash, zobrist, core_statistics);
    ++statistics.polish_moves;
    improved = true;
  }
  if (exact_core_determinant(state.core) != state.determinant) {
    throw std::runtime_error(
        "exact polish determinant invariant failed");
  }
  individual.core = state.core;
  individual.determinant = state.determinant;
  return improved;
}

int minimum_archive_distance(
    const CoreMatrix& candidate,
    const std::vector<ArchiveEntry>& archive,
    std::optional<std::size_t> excluded = std::nullopt) {
  int minimum = kCoreEntries + 1;
  for (std::size_t index = 0; index < archive.size(); ++index) {
    if (excluded && index == *excluded) {
      continue;
    }
    minimum =
        std::min(minimum,
                 hamming_distance(candidate, archive[index].core));
  }
  return minimum == kCoreEntries + 1 ? kCoreEntries : minimum;
}

bool archive_insert(
    const Individual& individual,
    std::vector<ArchiveEntry>& archive,
    std::unordered_set<std::string>& archive_keys,
    int capacity, SearchStatistics& statistics) {
  const std::string key = core_key(individual.core);
  if (archive_keys.contains(key)) {
    return false;
  }
  if (static_cast<int>(archive.size()) < capacity) {
    archive.push_back(
        ArchiveEntry{individual.core, individual.determinant});
    archive_keys.insert(key);
    ++statistics.archive_insertions;
    return true;
  }

  std::uint64_t worst_score =
      std::numeric_limits<std::uint64_t>::max();
  for (const ArchiveEntry& entry : archive) {
    worst_score =
        std::min(worst_score, magnitude(entry.determinant));
  }
  const std::uint64_t candidate_score =
      magnitude(individual.determinant);
  if (candidate_score < worst_score) {
    return false;
  }

  std::size_t replacement = archive.size();
  int replacement_novelty = kCoreEntries + 1;
  for (std::size_t index = 0; index < archive.size(); ++index) {
    if (magnitude(archive[index].determinant) != worst_score) {
      continue;
    }
    const int novelty =
        minimum_archive_distance(
            archive[index].core, archive, index);
    if (novelty < replacement_novelty) {
      replacement_novelty = novelty;
      replacement = index;
    }
  }
  if (replacement == archive.size()) {
    return false;
  }
  const int candidate_novelty =
      minimum_archive_distance(individual.core, archive);
  if (candidate_score == worst_score &&
      candidate_novelty <= replacement_novelty) {
    return false;
  }
  archive_keys.erase(core_key(archive[replacement].core));
  archive[replacement] =
      ArchiveEntry{individual.core, individual.determinant};
  archive_keys.insert(key);
  ++statistics.archive_replacements;
  return true;
}

std::pair<int, double> archive_distance_summary(
    const std::vector<ArchiveEntry>& archive) {
  if (archive.size() < 2) {
    return {0, 0.0};
  }
  int minimum = kCoreEntries;
  std::uint64_t sum = 0;
  std::uint64_t pairs = 0;
  for (std::size_t first = 0; first < archive.size(); ++first) {
    for (std::size_t second = first + 1;
         second < archive.size(); ++second) {
      const int distance =
          hamming_distance(
              archive[first].core, archive[second].core);
      minimum = std::min(minimum, distance);
      sum += static_cast<std::uint64_t>(distance);
      ++pairs;
    }
  }
  return {
      minimum,
      static_cast<double>(sum) / static_cast<double>(pairs)};
}

std::vector<Individual> load_seed_individuals(
    const std::vector<fs::path>& paths,
    const GomeaArguments& arguments,
    SearchStatistics& statistics) {
  std::vector<Individual> seeds;
  std::unordered_set<std::string> unique;
  for (const fs::path& path : paths) {
    const SignMatrix sign = read_sign_matrix(path);
    const CoreMatrix core = dephase_to_core(sign);
    const std::int64_t determinant =
        exact_core_determinant(core);
    differential_check_core(
        core, determinant, path.string(), statistics);
    const std::uint64_t score =
        magnitude(determinant) *
        static_cast<std::uint64_t>(kScale);
    if (score != arguments.frontier) {
      throw std::runtime_error(
          "seed is not at configured exact frontier: " +
          path.string());
    }
    const std::string key = core_key(core);
    if (!unique.insert(key).second) {
      throw std::runtime_error(
          "duplicate dephased seed core: " + path.string());
    }
    seeds.push_back(Individual{core, determinant, 0});
  }
  return seeds;
}

CoreMatrix random_kick(const CoreMatrix& source, int cardinality,
                       std::mt19937_64& randomizer) {
  std::array<int, kCoreEntries> variables{};
  std::iota(variables.begin(), variables.end(), 0);
  std::shuffle(variables.begin(), variables.end(), randomizer);
  CoreMatrix result = source;
  for (int index = 0; index < cardinality; ++index) {
    const int flat = variables[static_cast<std::size_t>(index)];
    result[flat / kCoreOrder][flat % kCoreOrder] ^= 1U;
  }
  return result;
}

std::vector<Individual> initialize_population(
    const std::vector<Individual>& seeds,
    const GomeaArguments& arguments,
    std::mt19937_64& randomizer,
    SearchStatistics& statistics) {
  std::vector<Individual> population = seeds;
  std::unordered_set<std::string> unique;
  for (const Individual& seed : seeds) {
    unique.insert(core_key(seed.core));
  }
  std::uniform_int_distribution<int> kick_size(
      arguments.initial_kick_min, arguments.initial_kick_max);
  std::uniform_int_distribution<std::size_t> seed_choice(
      0, seeds.size() - 1);
  int attempts = 0;
  while (static_cast<int>(population.size()) <
         arguments.population_size) {
    if (++attempts > 100000) {
      throw std::runtime_error(
          "could not initialize unique GOMEA population");
    }
    const Individual& parent = seeds[seed_choice(randomizer)];
    Individual candidate;
    candidate.core =
        random_kick(parent.core, kick_size(randomizer), randomizer);
    candidate.determinant =
        exact_core_determinant(candidate.core);
    ++statistics.exact_evaluations;
    if (candidate.determinant == 0) {
      continue;
    }
    static_cast<void>(
        exact_polish(candidate, arguments.polish_moves, statistics));
    const std::string key = core_key(candidate.core);
    if (unique.insert(key).second) {
      population.push_back(std::move(candidate));
    }
  }
  return population;
}

void run_differential_samples(
    const std::vector<Individual>& seeds,
    const GomeaArguments& arguments,
    std::mt19937_64& randomizer,
    SearchStatistics& statistics) {
  if (arguments.differential_samples == 0) {
    return;
  }
  std::uniform_int_distribution<std::size_t> seed_choice(
      0, seeds.size() - 1);
  std::uniform_int_distribution<int> subset_size(1, 128);
  std::array<int, kCoreEntries> variables{};
  std::iota(variables.begin(), variables.end(), 0);
  for (int sample = 0; sample < arguments.differential_samples;
       ++sample) {
    const Individual& recipient = seeds[seed_choice(randomizer)];
    const Individual& donor = seeds[seed_choice(randomizer)];
    std::shuffle(variables.begin(), variables.end(), randomizer);
    CoreMatrix candidate = recipient.core;
    const int cardinality = subset_size(randomizer);
    for (int index = 0; index < cardinality; ++index) {
      const int flat = variables[static_cast<std::size_t>(index)];
      candidate[flat / kCoreOrder][flat % kCoreOrder] =
          donor.core[flat / kCoreOrder][flat % kCoreOrder];
    }
    const std::int64_t determinant =
        exact_core_determinant(candidate);
    ++statistics.exact_evaluations;
    differential_check_core(
        candidate, determinant,
        "random donor subset " + std::to_string(sample),
        statistics);
  }
}

const CoreMatrix& choose_donor(
    const std::vector<Individual>& population,
    const std::vector<ArchiveEntry>& archive,
    std::size_t recipient,
    std::mt19937_64& randomizer) {
  const bool use_archive =
      !archive.empty() &&
      (randomizer() % UINT64_C(4) == 0U);
  if (use_archive) {
    return archive[
        static_cast<std::size_t>(
            randomizer() %
            static_cast<std::uint64_t>(archive.size()))]
        .core;
  }
  if (population.size() == 1) {
    return population.front().core;
  }
  std::size_t donor = recipient;
  while (donor == recipient) {
    donor =
        static_cast<std::size_t>(
            randomizer() %
            static_cast<std::uint64_t>(population.size()));
  }
  return population[donor].core;
}

bool population_contains_core(
    const CoreMatrix& core,
    const std::vector<Individual>& population,
    std::size_t excluded) {
  for (std::size_t index = 0; index < population.size(); ++index) {
    if (index != excluded && population[index].core == core) {
      return true;
    }
  }
  return false;
}

bool force_improvement(
    Individual& individual, const Individual& global_best,
    const std::vector<LinkageSubset>& linkage,
    const std::vector<Individual>& population,
    std::size_t recipient,
    std::mt19937_64& randomizer,
    SearchStatistics& statistics) {
  std::vector<std::size_t> order(linkage.size());
  std::iota(order.begin(), order.end(), 0);
  std::shuffle(order.begin(), order.end(), randomizer);
  const std::uint64_t original_score =
      magnitude(individual.determinant);
  for (const std::size_t subset_index : order) {
    CoreMatrix candidate = individual.core;
    if (!copy_subset(
            candidate, global_best.core,
            linkage[subset_index].variables)) {
      continue;
    }
    ++statistics.mix_attempts;
    ++statistics.attempted_by_source[
        linkage[subset_index].source];
    const std::int64_t determinant =
        exact_core_determinant(candidate);
    ++statistics.exact_evaluations;
    if (determinant == 0) {
      ++statistics.rejected_singular;
      continue;
    }
    if (magnitude(determinant) > original_score &&
        !population_contains_core(
            candidate, population, recipient)) {
      individual.core = candidate;
      individual.determinant = determinant;
      individual.stagnant_generations = 0;
      ++statistics.forced_improvements;
      ++statistics.accepted_by_source[
          linkage[subset_index].source];
      return true;
    }
  }
  if (magnitude(global_best.determinant) >
          magnitude(individual.determinant) &&
      !population_contains_core(
          global_best.core, population, recipient)) {
    individual = global_best;
    individual.stagnant_generations = 0;
    ++statistics.forced_replacements;
    return true;
  }
  return false;
}

bool force_diverse_archive_reset(
    Individual& individual,
    const std::vector<Individual>& population,
    const std::vector<ArchiveEntry>& archive,
    std::size_t recipient, SearchStatistics& statistics) {
  const std::uint64_t current_score =
      magnitude(individual.determinant);
  const ArchiveEntry* selected = nullptr;
  int selected_novelty = -1;
  std::uint64_t selected_score = 0;
  for (const ArchiveEntry& entry : archive) {
    const std::uint64_t entry_score = magnitude(entry.determinant);
    if (entry_score < current_score ||
        population_contains_core(
            entry.core, population, recipient)) {
      continue;
    }
    int novelty = kCoreEntries;
    for (std::size_t index = 0; index < population.size(); ++index) {
      if (index == recipient) {
        continue;
      }
      novelty = std::min(
          novelty,
          hamming_distance(entry.core, population[index].core));
    }
    if (selected == nullptr || novelty > selected_novelty ||
        (novelty == selected_novelty &&
         entry_score > selected_score)) {
      selected = &entry;
      selected_novelty = novelty;
      selected_score = entry_score;
    }
  }
  if (selected == nullptr) {
    return false;
  }
  individual.core = selected->core;
  individual.determinant = selected->determinant;
  individual.stagnant_generations = 0;
  ++statistics.diversity_resets;
  return true;
}

std::string seed_paths_json(const std::vector<fs::path>& paths) {
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

std::string score_string(std::int64_t determinant) {
  return wide_to_string(
      static_cast<Wide>(magnitude(determinant)) * kScale);
}

std::string counters_json(
    const std::map<std::string, std::uint64_t>& counters) {
  std::ostringstream output;
  output << '{';
  bool first = true;
  for (const auto& [name, count] : counters) {
    if (!first) {
      output << ',';
    }
    first = false;
    output << '"' << json_escape(name) << "\":" << count;
  }
  output << '}';
  return output.str();
}

void log_state(
    GomeaLogger& logger, const GomeaArguments& arguments,
    const SearchStatistics& statistics, std::string_view event,
    double elapsed, const Individual& global_best,
    const std::vector<Individual>& population,
    const std::vector<ArchiveEntry>& archive,
    std::size_t linkage_size, double linkage_seconds,
    double mixing_seconds) {
  std::unordered_set<std::string> population_keys;
  std::uint64_t population_min =
      std::numeric_limits<std::uint64_t>::max();
  std::uint64_t population_max = 0;
  for (const Individual& individual : population) {
    population_keys.insert(core_key(individual.core));
    population_min =
        std::min(population_min, magnitude(individual.determinant));
    population_max =
        std::max(population_max, magnitude(individual.determinant));
  }
  const auto [archive_min_distance, archive_mean_distance] =
      archive_distance_summary(archive);
  std::ostringstream record;
  record
      << "{\"event\":\"" << event << "\""
      << ",\"engine\":\"exact-core-linkage-gomea-v1\""
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(6) << elapsed
      << ",\"seed\":" << arguments.seed
      << ",\"generations\":" << statistics.generations
      << ",\"linkage_rebuilds\":" << statistics.linkage_rebuilds
      << ",\"linkage_subsets\":" << linkage_size
      << ",\"linkage_seconds\":" << linkage_seconds
      << ",\"mixing_seconds\":" << mixing_seconds
      << ",\"exact_evaluations\":"
      << statistics.exact_evaluations
      << ",\"mix_attempts\":" << statistics.mix_attempts
      << ",\"mix_changed\":" << statistics.mix_changed
      << ",\"mix_attempted_by_source\":"
      << counters_json(statistics.attempted_by_source)
      << ",\"accepted_improvements\":"
      << statistics.accepted_improvements
      << ",\"accepted_equals\":"
      << statistics.accepted_equals
      << ",\"rejected_equal_duplicates\":"
      << statistics.rejected_equal_duplicates
      << ",\"rejected_improving_duplicates\":"
      << statistics.rejected_improving_duplicates
      << ",\"rejected_polish_duplicates\":"
      << statistics.rejected_polish_duplicates
      << ",\"accepted_by_source\":"
      << counters_json(statistics.accepted_by_source)
      << ",\"rejected_decreases\":"
      << statistics.rejected_decreases
      << ",\"rejected_singular\":"
      << statistics.rejected_singular
      << ",\"polish_attempts\":" << statistics.polish_attempts
      << ",\"polish_moves\":" << statistics.polish_moves
      << ",\"forced_improvements\":"
      << statistics.forced_improvements
      << ",\"forced_replacements\":"
      << statistics.forced_replacements
      << ",\"diversity_resets\":"
      << statistics.diversity_resets
      << ",\"bridge_endpoint_donations\":"
      << statistics.bridge_endpoint_donations
      << ",\"archive_insertions\":"
      << statistics.archive_insertions
      << ",\"archive_replacements\":"
      << statistics.archive_replacements
      << ",\"archive_size\":" << archive.size()
      << ",\"archive_min_hamming\":" << archive_min_distance
      << ",\"archive_mean_hamming\":"
      << archive_mean_distance
      << ",\"population_size\":" << population.size()
      << ",\"population_unique\":" << population_keys.size()
      << ",\"population_min_quotient\":" << population_min
      << ",\"population_max_quotient\":" << population_max
      << ",\"promotions\":" << statistics.promotions
      << ",\"differential_checks\":"
      << statistics.differential_checks
      << ",\"best_core_quotient\":"
      << magnitude(global_best.determinant)
      << ",\"best_absolute_determinant\":\""
      << score_string(global_best.determinant) << "\""
      << ",\"above_frontier_unverified\":"
      << (magnitude(global_best.determinant) *
                      static_cast<std::uint64_t>(kScale) >
                  arguments.frontier
              ? "true"
              : "false")
      << '}';
  logger.line(record.str());
}

void write_archive(
    const fs::path& directory,
    std::vector<ArchiveEntry> archive,
    std::uint64_t& nonce) {
  std::sort(
      archive.begin(), archive.end(),
      [](const ArchiveEntry& left, const ArchiveEntry& right) {
        const std::uint64_t left_score = magnitude(left.determinant);
        const std::uint64_t right_score = magnitude(right.determinant);
        if (left_score != right_score) {
          return left_score > right_score;
        }
        return core_key(left.core) < core_key(right.core);
      });
  fs::create_directories(directory);
  std::ostringstream manifest;
  manifest << "{\"schema_version\":1,\"members\":[";
  for (std::size_t index = 0; index < archive.size(); ++index) {
    std::ostringstream filename;
    filename << "member-" << std::setw(3) << std::setfill('0')
             << index << "-q" << magnitude(archive[index].determinant)
             << "-" << hex_u64(core_fingerprint(archive[index].core))
             << ".matrix.txt";
    atomic_write(
        directory / filename.str(),
        sign_matrix_text(core_to_sign(archive[index].core)),
        nonce++);
    if (index != 0) {
      manifest << ',';
    }
    manifest
        << "{\"path\":\"" << json_escape(filename.str())
        << "\",\"core_quotient\":"
        << archive[index].determinant
        << ",\"absolute_determinant\":\""
        << score_string(archive[index].determinant)
        << "\",\"core_fingerprint\":\""
        << hex_u64(core_fingerprint(archive[index].core))
        << "\"}";
  }
  manifest << "]}";
  atomic_write(
      directory / "manifest.json", manifest.str() + "\n", nonce++);
}

std::string final_summary_json(
    const GomeaArguments& arguments,
    const SearchStatistics& statistics,
    const Individual& global_best,
    const std::vector<Individual>& population,
    const std::vector<ArchiveEntry>& archive,
    const std::vector<fs::path>& seeds,
    const BridgeSupportInfo& bridge,
    std::size_t linkage_size, double elapsed,
    bool stopped_by_signal) {
  const auto [minimum_distance, mean_distance] =
      archive_distance_summary(archive);
  std::unordered_set<std::string> population_keys;
  for (const Individual& individual : population) {
    population_keys.insert(core_key(individual.core));
  }
  std::ostringstream output;
  output
      << "{\"schema_version\":1"
      << ",\"engine\":\"exact-core-linkage-gomea-v1\""
      << ",\"complete\":" << (stopped_by_signal ? "false" : "true")
      << ",\"stopped_by_signal\":"
      << (stopped_by_signal ? "true" : "false")
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(6) << elapsed
      << ",\"seed\":" << arguments.seed
      << ",\"input_seeds\":" << seed_paths_json(seeds)
      << ",\"alignment_metadata\":\""
      << json_escape(arguments.alignment_metadata.string()) << "\""
      << ",\"bridge_reference\":\""
      << json_escape(arguments.bridge_reference.string()) << "\""
      << ",\"bridge_seed\":\""
      << json_escape(arguments.bridge_seed.string()) << "\""
      << ",\"bridge_raw_support_size\":"
      << bridge.raw_coordinates.size()
      << ",\"bridge_core_support_size\":"
      << bridge.linkage.variables.size()
      << ",\"bridge_raw_coordinates_one_based\":"
      << coordinate_pairs_json(bridge.raw_coordinates)
      << ",\"population_size\":" << population.size()
      << ",\"population_unique\":" << population_keys.size()
      << ",\"archive_size\":" << archive.size()
      << ",\"archive_min_hamming\":" << minimum_distance
      << ",\"archive_mean_hamming\":" << mean_distance
      << ",\"linkage_subsets\":" << linkage_size
      << ",\"known_neutral_raw_12_supports\":12"
      << ",\"known_neutral_induced_core_12_supports\":6"
      << ",\"known_neutral_induced_core_31_supports\":6"
      << ",\"transported_bridge_12_supports\":1"
      << ",\"generations\":" << statistics.generations
      << ",\"linkage_rebuilds\":" << statistics.linkage_rebuilds
      << ",\"exact_evaluations\":"
      << statistics.exact_evaluations
      << ",\"mix_attempts\":" << statistics.mix_attempts
      << ",\"mix_changed\":" << statistics.mix_changed
      << ",\"mix_attempted_by_source\":"
      << counters_json(statistics.attempted_by_source)
      << ",\"accepted_improvements\":"
      << statistics.accepted_improvements
      << ",\"accepted_equals\":"
      << statistics.accepted_equals
      << ",\"rejected_equal_duplicates\":"
      << statistics.rejected_equal_duplicates
      << ",\"rejected_improving_duplicates\":"
      << statistics.rejected_improving_duplicates
      << ",\"rejected_polish_duplicates\":"
      << statistics.rejected_polish_duplicates
      << ",\"accepted_by_source\":"
      << counters_json(statistics.accepted_by_source)
      << ",\"rejected_decreases\":"
      << statistics.rejected_decreases
      << ",\"rejected_singular\":"
      << statistics.rejected_singular
      << ",\"polish_attempts\":" << statistics.polish_attempts
      << ",\"polish_moves\":" << statistics.polish_moves
      << ",\"forced_improvements\":"
      << statistics.forced_improvements
      << ",\"forced_replacements\":"
      << statistics.forced_replacements
      << ",\"diversity_resets\":"
      << statistics.diversity_resets
      << ",\"bridge_endpoint_donations\":"
      << statistics.bridge_endpoint_donations
      << ",\"archive_insertions\":"
      << statistics.archive_insertions
      << ",\"archive_replacements\":"
      << statistics.archive_replacements
      << ",\"promotions\":" << statistics.promotions
      << ",\"differential_checks\":"
      << statistics.differential_checks
      << ",\"best_core_quotient\":"
      << magnitude(global_best.determinant)
      << ",\"best_absolute_determinant\":\""
      << score_string(global_best.determinant) << "\""
      << ",\"frontier\":\"" << arguments.frontier << "\""
      << ",\"above_frontier_unverified\":"
      << (magnitude(global_best.determinant) *
                      static_cast<std::uint64_t>(kScale) >
                  arguments.frontier
              ? "true"
              : "false")
      << ",\"claim_boundary\":["
      << "\"All objective comparisons were exact integer comparisons.\","
      << "\"The retained matrix still requires arena verification.\","
      << "\"Archive diversity is not an H-equivalence or novelty claim.\""
      << "]}";
  return output.str() + "\n";
}

int run_gomea(const GomeaArguments& arguments) {
  if (fs::exists(arguments.output) ||
      fs::exists(arguments.archive_directory) ||
      fs::exists(arguments.log) ||
      fs::exists(arguments.summary)) {
    throw std::runtime_error(
        "output, archive, log, and summary paths must be fresh");
  }
  GomeaLogger logger(arguments.log);
  std::signal(SIGINT, request_stop);
  std::signal(SIGTERM, request_stop);
  SearchStatistics statistics;
  std::mt19937_64 randomizer(arguments.seed);
  const std::vector<fs::path> seed_paths =
      collect_seed_paths(arguments);
  std::vector<Individual> seed_individuals =
      load_seed_individuals(seed_paths, arguments, statistics);
  std::vector<LinkageSubset> known_supports =
      derive_known_neutral_supports(seed_individuals);
  const BridgeSupportInfo bridge =
      derive_bridge_support(
          arguments, seed_individuals, statistics);
  known_supports.push_back(bridge.linkage);
  run_differential_samples(
      seed_individuals, arguments, randomizer, statistics);

  std::vector<Individual> population =
      initialize_population(
          seed_individuals, arguments, randomizer, statistics);
  Individual global_best = population.front();
  for (const Individual& individual : population) {
    if (magnitude(individual.determinant) >
        magnitude(global_best.determinant)) {
      global_best = individual;
    }
  }
  std::vector<ArchiveEntry> archive;
  std::unordered_set<std::string> archive_keys;
  for (const Individual& individual : population) {
    static_cast<void>(archive_insert(
        individual, archive, archive_keys,
        arguments.archive_size, statistics));
  }
  const Individual bridge_reference{
      bridge.reference_core, bridge.reference_determinant, 0};
  static_cast<void>(archive_insert(
      bridge_reference, archive, archive_keys,
      arguments.archive_size, statistics));

  std::uint64_t nonce = 0;
  atomic_write(
      arguments.output,
      sign_matrix_text(core_to_sign(global_best.core)), nonce++);
  const auto started = Clock::now();
  const auto deadline =
      started + std::chrono::duration<double>(arguments.seconds);
  auto next_heartbeat =
      started +
      std::chrono::duration<double>(
          arguments.heartbeat_seconds);

  {
    std::ostringstream record;
    record
        << "{\"event\":\"started\""
        << ",\"engine\":\"exact-core-linkage-gomea-v1\""
        << ",\"seed\":" << arguments.seed
        << ",\"seconds\":" << arguments.seconds
        << ",\"frontier\":\"" << arguments.frontier << "\""
        << ",\"seed_paths\":" << seed_paths_json(seed_paths)
        << ",\"input_seed_count\":" << seed_individuals.size()
        << ",\"population_size\":" << population.size()
        << ",\"archive_capacity\":" << arguments.archive_size
        << ",\"known_neutral_raw_12_supports\":"
        << known_supports.size() - 1
        << ",\"known_neutral_induced_core_12_supports\":6"
        << ",\"known_neutral_induced_core_31_supports\":6"
        << ",\"transported_bridge_12_supports\":1"
        << ",\"bridge_raw_coordinates_one_based\":"
        << coordinate_pairs_json(bridge.raw_coordinates)
        << ",\"alignment_metadata\":\""
        << json_escape(arguments.alignment_metadata.string()) << "\""
        << ",\"linkage_refresh\":"
        << arguments.linkage_refresh
        << ",\"maximum_linkage_size\":"
        << arguments.maximum_linkage_size
        << ",\"polish_moves\":" << arguments.polish_moves
        << ",\"differential_checks\":"
        << statistics.differential_checks
        << ",\"best_absolute_determinant\":\""
        << score_string(global_best.determinant) << "\"}";
    logger.line(record.str());
  }

  std::vector<LinkageSubset> linkage;
  double latest_linkage_seconds = 0.0;
  double latest_mixing_seconds = 0.0;
  while (!stop_requested && Clock::now() < deadline) {
    if (linkage.empty() ||
        statistics.generations %
                static_cast<std::uint64_t>(
                    arguments.linkage_refresh) ==
            0) {
      const auto linkage_started = Clock::now();
      linkage = learn_linkage(
          population, known_supports,
          arguments.maximum_linkage_size);
      latest_linkage_seconds =
          std::chrono::duration<double>(
              Clock::now() - linkage_started)
              .count();
      ++statistics.linkage_rebuilds;
    } else {
      latest_linkage_seconds = 0.0;
    }

    const auto mixing_started = Clock::now();
    std::vector<std::size_t> subset_order(linkage.size());
    std::iota(
        subset_order.begin(), subset_order.end(),
        static_cast<std::size_t>(0));
    for (std::size_t recipient = 0;
         recipient < population.size() && !stop_requested;
         ++recipient) {
      Individual& individual = population[recipient];
      const std::uint64_t starting_score =
          magnitude(individual.determinant);
      std::shuffle(
          subset_order.begin(), subset_order.end(), randomizer);
      std::uint64_t subset_counter = 0;
      for (const std::size_t subset_index : subset_order) {
        if ((++subset_counter & UINT64_C(63)) == 0U &&
            Clock::now() >= deadline) {
          break;
        }
        const bool is_bridge_subset =
            linkage[subset_index].source ==
            bridge.linkage.source;
        const CoreMatrix& donor =
            is_bridge_subset
                ? (individual.core == bridge.bridge_core
                       ? bridge.reference_core
                       : individual.core == bridge.reference_core
                             ? bridge.bridge_core
                             : (randomizer() % UINT64_C(2) == 0U
                                    ? bridge.reference_core
                                    : bridge.bridge_core))
                : choose_donor(
                      population, archive, recipient, randomizer);
        if (is_bridge_subset) {
          ++statistics.bridge_endpoint_donations;
        }
        CoreMatrix candidate = individual.core;
        if (!copy_subset(
                candidate, donor,
                linkage[subset_index].variables)) {
          continue;
        }
        ++statistics.mix_attempts;
        ++statistics.mix_changed;
        ++statistics.attempted_by_source[
            linkage[subset_index].source];
        const std::int64_t determinant =
            exact_core_determinant(candidate);
        ++statistics.exact_evaluations;
        if (determinant == 0) {
          ++statistics.rejected_singular;
          continue;
        }
        const std::uint64_t candidate_score =
            magnitude(determinant);
        const std::uint64_t current_score =
            magnitude(individual.determinant);
        if (candidate_score < current_score) {
          ++statistics.rejected_decreases;
          continue;
        }
        if (population_contains_core(
                candidate, population, recipient)) {
          if (candidate_score == current_score) {
            ++statistics.rejected_equal_duplicates;
          } else {
            ++statistics.rejected_improving_duplicates;
          }
          continue;
        }
        individual.core = candidate;
        individual.determinant = determinant;
        if (candidate_score > current_score) {
          ++statistics.accepted_improvements;
        } else {
          ++statistics.accepted_equals;
        }
        ++statistics.accepted_by_source[
            linkage[subset_index].source];
        // Preserve transient exact ties/improvements before later linkage
        // subsets overwrite them during the same optimal-mixing pass.
        static_cast<void>(archive_insert(
            individual, archive, archive_keys,
            arguments.archive_size, statistics));
        if (candidate_score >
            magnitude(global_best.determinant)) {
          global_best = individual;
          ++statistics.promotions;
          atomic_write(
              arguments.output,
              sign_matrix_text(core_to_sign(global_best.core)),
              nonce++);
          const double elapsed =
              std::chrono::duration<double>(
                  Clock::now() - started)
                  .count();
          log_state(
              logger, arguments, statistics, "new_best",
              elapsed, global_best, population, archive,
              linkage.size(), latest_linkage_seconds,
              std::chrono::duration<double>(
                  Clock::now() - mixing_started)
                  .count());
          std::cout
              << "new best |det|="
              << score_string(global_best.determinant)
              << " generation=" << statistics.generations
              << '\n'
              << std::flush;
        }
      }

      Individual polished = individual;
      if (exact_polish(
              polished, arguments.polish_moves, statistics)) {
        if (population_contains_core(
                polished.core, population, recipient)) {
          ++statistics.rejected_polish_duplicates;
        } else {
          individual = std::move(polished);
        }
      }
      if (magnitude(individual.determinant) >
          magnitude(global_best.determinant)) {
        global_best = individual;
        ++statistics.promotions;
        atomic_write(
            arguments.output,
            sign_matrix_text(core_to_sign(global_best.core)),
            nonce++);
      }
      if (magnitude(individual.determinant) > starting_score) {
        individual.stagnant_generations = 0;
      } else {
        ++individual.stagnant_generations;
      }
      if (individual.stagnant_generations >=
          arguments.forced_improvement_generations) {
        static_cast<void>(force_improvement(
            individual, global_best, linkage, population,
            recipient, randomizer, statistics));
        if (individual.stagnant_generations >=
            arguments.forced_improvement_generations) {
          static_cast<void>(force_diverse_archive_reset(
              individual, population, archive, recipient,
              statistics));
        }
      }
      static_cast<void>(archive_insert(
          individual, archive, archive_keys,
          arguments.archive_size, statistics));
      if (Clock::now() >= deadline) {
        break;
      }
    }
    latest_mixing_seconds =
        std::chrono::duration<double>(
            Clock::now() - mixing_started)
            .count();
    ++statistics.generations;

    const auto now = Clock::now();
    if (arguments.heartbeat_seconds > 0.0 &&
        now >= next_heartbeat) {
      const double elapsed =
          std::chrono::duration<double>(now - started).count();
      log_state(
          logger, arguments, statistics, "heartbeat", elapsed,
          global_best, population, archive, linkage.size(),
          latest_linkage_seconds, latest_mixing_seconds);
      next_heartbeat =
          now +
          std::chrono::duration<double>(
              arguments.heartbeat_seconds);
    }
  }

  differential_check_core(
      global_best.core, global_best.determinant,
      "final global best", statistics);
  atomic_write(
      arguments.output,
      sign_matrix_text(core_to_sign(global_best.core)), nonce++);
  write_archive(arguments.archive_directory, archive, nonce);
  const double elapsed =
      std::chrono::duration<double>(Clock::now() - started).count();
  atomic_write(
      arguments.summary,
      final_summary_json(
          arguments, statistics, global_best, population, archive,
          seed_paths, bridge, linkage.size(), elapsed,
          stop_requested != 0),
      nonce++);
  log_state(
      logger, arguments, statistics,
      stop_requested ? "stopped" : "finished", elapsed,
      global_best, population, archive, linkage.size(),
      latest_linkage_seconds, latest_mixing_seconds);

  std::cout
      << "finished |det|=" << score_string(global_best.determinant)
      << " generations=" << statistics.generations
      << " exact_evaluations=" << statistics.exact_evaluations
      << " archive=" << archive.size() << '\n';
  if (magnitude(global_best.determinant) *
          static_cast<std::uint64_t>(kScale) >
      arguments.frontier) {
    std::cout
        << "UNVERIFIED FRONTIER CANDIDATE: run ./arena verify "
        << arguments.output << '\n';
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    return run_gomea(parse_gomea_arguments(argc, argv));
  } catch (const std::exception& error) {
    std::cerr << "core_linkage_gomea: " << error.what() << '\n';
    return 1;
  }
}
