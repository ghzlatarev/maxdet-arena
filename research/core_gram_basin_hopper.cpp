// Exact Gram-sketch basin hopping for the order-23 maximal determinant
// search.
//
// Each basin starts with the best of N exact random k-bit kicks from a
// retained parent.  The kick endpoint is observed before repair.  A
// 10,000--40,000 move reactive-tabu epoch then crosses the basin, retaining
// its exact peak and gated transient/polished samples.
//
// The diversity key is deliberately called a Gram sketch, not a canonical
// form.  For odd order, multiplying a signed line by the product of its
// entries gives an orientation invariant under line negation.  Its normalized
// off-diagonal Gram entries are all 3 modulo 4.  For both rows and columns we
// retain:
//
//   * the global normalized signed Gram-label histogram; and
//   * the lexicographically sorted per-vertex incident-label histograms.
//
// The row and column records are stored as an unordered pair.  The sketch is
// invariant under signed row/column permutations and transpose, but distinct
// matrices can collide.  H/HT classification therefore remains the job of the
// separate exact classifier.
//
// Optional --coronal-pareto selection also forms the two exact oriented
// coronals.  For the row orientation, D contains the odd-order row products,
// G = D H H^T D has every off-diagonal entry congruent to 3 modulo 4,
//
//   W = (G - (24 I - J)) / 4,  M = 6 I + W = (G + J) / 4,
//   q^2 = det(M) (4 - kappa),  kappa = 1^T M^-1 1.
//
// We compute det(M) by exact Bareiss elimination and derive kappa as the
// reduced rational 4 - q^2/det(M), without a floating or exact inverse.  Row
// and column points are sorted as an unordered pair so this search feature is
// transpose invariant.  They guide archive survival and parent selection
// only; the promoted objective remains the exact determinant quotient q.
//
// All determinant objectives are exact integers in the dephased 22 by 22
// binary core.  A strict order-23 frontier improvement is exactly a core
// quotient of at least 662,671,876.

#include <map>
#include <optional>
#include <set>
#include <tuple>
#include <unordered_set>

#define CORE_ADJUGATE_TABU_NO_MAIN
#include "core_adjugate_tabu.cpp"
#undef CORE_ADJUGATE_TABU_NO_MAIN

namespace {

constexpr std::uint64_t kKnownCoreFrontier = UINT64_C(662671875);
constexpr std::uint64_t kStrictCoreTarget = UINT64_C(662671876);
constexpr int kGramLabelCount = 12;
constexpr std::size_t kOneSidedSketchLength =
    static_cast<std::size_t>(
        kGramLabelCount * (kSignOrder + 1));
constexpr std::uint64_t kEpochIdentityCheckInterval = UINT64_C(4096);
constexpr std::uint64_t kEpochDeterminantCheckInterval = UINT64_C(16384);

// M is positive definite with diagonal 6 and entries in [-5, 6], so
// det(M) <= 6^23.  For this positive-definite M, compound-matrix
// Cauchy--Schwarz bounds every minor by 6^k.  Bareiss numerator products are
// therefore below 6^46 < 10^36, inside signed 128-bit range.  This gives
// dependency-free exact arithmetic.
using BigInt = Wide;

struct HopperArguments {
  std::vector<fs::path> seeds;
  fs::path output;
  fs::path archive_directory;
  fs::path log;
  fs::path summary;
  std::uint64_t seed = UINT64_C(37001);
  double seconds = 900.0;
  double heartbeat_seconds = 15.0;
  std::uint64_t maximum_epochs = 0;
  std::uint64_t quotient_gate = UINT64_C(600000000);
  bool gate_all = false;
  int epoch_min_moves = 10000;
  int epoch_max_moves = 40000;
  int kick_flips = 32;
  int kick_attempts = 16;
  int transient_stride = 2048;
  int polish_moves = 96;
  int archive_capacity = 256;
  int archive_outputs = 64;
  int quality_tournament = 8;
  int self_test_rounds = 0;
  bool coronal_pareto = false;
};

struct GramSketch {
  std::vector<std::uint16_t> labels;

  bool operator<(const GramSketch& other) const {
    return labels < other.labels;
  }

  bool operator==(const GramSketch& other) const {
    return labels == other.labels;
  }
};

struct SketchHash {
  std::size_t operator()(const GramSketch& sketch) const noexcept {
    std::uint64_t hash = UINT64_C(1469598103934665603);
    for (const std::uint16_t label : sketch.labels) {
      hash ^= static_cast<std::uint64_t>(label & UINT16_C(0xff));
      hash *= UINT64_C(1099511628211);
      hash ^= static_cast<std::uint64_t>(label >> 8U);
      hash *= UINT64_C(1099511628211);
    }
    return static_cast<std::size_t>(hash);
  }
};

struct CoronalPoint {
  BigInt determinant_m = 0;
  BigInt kappa_numerator = 0;
  BigInt kappa_denominator = 1;

  bool operator==(const CoronalPoint& other) const {
    return determinant_m == other.determinant_m &&
        kappa_numerator == other.kappa_numerator &&
        kappa_denominator == other.kappa_denominator;
  }
};

struct CoronalPair {
  std::array<CoronalPoint, 2> orientations{};

  bool operator==(const CoronalPair& other) const {
    return orientations == other.orientations;
  }
};

struct SketchElite {
  CoreMatrix core{};
  std::int64_t determinant = 0;
  std::uint64_t selections = 0;
  std::uint64_t inserted_at_epoch = 0;
  std::uint64_t last_updated_epoch = 0;
  std::uint64_t sketch_hash = 0;
  std::string origin;
  std::optional<CoronalPair> coronal;
};

struct HopperStatistics {
  std::uint64_t epochs_started = 0;
  std::uint64_t epochs_completed = 0;
  std::uint64_t tabu_moves = 0;
  std::uint64_t polish_moves = 0;
  std::uint64_t kick_evaluations = 0;
  std::uint64_t kick_singular = 0;
  std::uint64_t kick_endpoints = 0;
  std::uint64_t parent_quality = 0;
  std::uint64_t parent_coronal_det_m = 0;
  std::uint64_t parent_coronal_kappa = 0;
  std::uint64_t parent_coronal_pareto = 0;
  std::uint64_t parent_least_used = 0;
  std::uint64_t parent_uniform = 0;
  std::uint64_t observations = 0;
  std::uint64_t seed_observations = 0;
  std::uint64_t kick_observations = 0;
  std::uint64_t transient_observations = 0;
  std::uint64_t peak_observations = 0;
  std::uint64_t polished_observations = 0;
  std::uint64_t gate_rejections = 0;
  std::uint64_t sketch_evaluations = 0;
  std::uint64_t sketch_discoveries = 0;
  std::uint64_t initialization_sketch_discoveries = 0;
  std::uint64_t initialization_promotions = 0;
  std::uint64_t initialization_archive_size = 0;
  bool initialization_complete = false;
  std::uint64_t archive_insertions = 0;
  std::uint64_t archive_replacements = 0;
  std::uint64_t archive_rejections = 0;
  std::uint64_t archive_capacity_rejections = 0;
  std::uint64_t archive_evictions = 0;
  std::uint64_t coronal_capacity_rejections = 0;
  std::uint64_t coronal_pareto_evictions = 0;
  std::uint64_t promotions = 0;
  std::uint64_t strict_target_states = 0;
  std::uint64_t identity_checks = 0;
  std::uint64_t determinant_checks = 0;
  std::uint64_t descriptor_invariance_checks = 0;
  std::uint64_t coronal_pair_evaluations = 0;
  std::uint64_t coronal_formula_checks = 0;
  std::uint64_t coronal_invariance_checks = 0;
  Statistics kernel{};
};

enum class ParentMode : std::uint8_t {
  kQuality,
  kCoronalDetM,
  kCoronalKappa,
  kCoronalPareto,
  kLeastUsed,
  kUniform,
};

enum class ObservationKind : std::uint8_t {
  kSeed,
  kKick,
  kTransient,
  kPeak,
  kPolished,
};

const char* observation_name(ObservationKind kind) {
  switch (kind) {
    case ObservationKind::kSeed:
      return "seed";
    case ObservationKind::kKick:
      return "kick";
    case ObservationKind::kTransient:
      return "transient";
    case ObservationKind::kPeak:
      return "epoch_peak";
    case ObservationKind::kPolished:
      return "polished";
  }
  throw std::runtime_error("unknown observation kind");
}

bool archive_gate_applies(
    const HopperArguments& arguments, ObservationKind kind,
    bool force_seed) {
  if (force_seed || kind == ObservationKind::kSeed) return false;
  if (arguments.gate_all) return true;
  return kind == ObservationKind::kTransient ||
      kind == ObservationKind::kPolished;
}

const char* gate_scope_name(const HopperArguments& arguments) {
  return arguments.gate_all
      ? "all_nonseed_archive_observations"
      : "transient_and_polished_only";
}

std::uint64_t search_added(
    std::uint64_t total, std::uint64_t initialization,
    bool initialization_complete) {
  if (!initialization_complete) return 0;
  if (total < initialization) {
    throw std::runtime_error(
        "initialization baseline exceeded its monotone counter");
  }
  return total - initialization;
}

bool path_entry_exists(const fs::path& path) {
  return fs::symlink_status(path).type() !=
      fs::file_type::not_found;
}

fs::path normalized_absolute_path(const fs::path& path) {
  return fs::weakly_canonical(fs::absolute(path));
}

bool path_is_within_or_equal(
    const fs::path& candidate, const fs::path& directory) {
  auto candidate_part = candidate.begin();
  for (auto directory_part = directory.begin();
       directory_part != directory.end();
       ++directory_part, ++candidate_part) {
    if (
        candidate_part == candidate.end() ||
        *candidate_part != *directory_part) {
      return false;
    }
  }
  return true;
}

int hopper_positive_int(
    std::string_view text, std::string_view option, int maximum) {
  const std::uint64_t parsed = strict_unsigned(text, option);
  if (parsed == 0 ||
      parsed > static_cast<std::uint64_t>(maximum)) {
    throw std::runtime_error(
        std::string(option) + " is outside its valid positive range");
  }
  return static_cast<int>(parsed);
}

HopperArguments parse_hopper_arguments(int argc, char** argv) {
  HopperArguments arguments;
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
    } else if (option == "--seconds") {
      arguments.seconds = strict_double(value(), option, false);
    } else if (
        option == "--heartbeat" ||
        option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds =
          strict_double(value(), option, true);
    } else if (option == "--max-epochs") {
      arguments.maximum_epochs = strict_unsigned(value(), option);
    } else if (option == "--quotient-gate") {
      arguments.quotient_gate = strict_unsigned(value(), option);
    } else if (option == "--gate-all") {
      arguments.gate_all = true;
    } else if (option == "--coronal-pareto") {
      arguments.coronal_pareto = true;
    } else if (option == "--epoch-min-moves") {
      arguments.epoch_min_moves =
          hopper_positive_int(value(), option, 40000);
    } else if (option == "--epoch-max-moves") {
      arguments.epoch_max_moves =
          hopper_positive_int(value(), option, 40000);
    } else if (option == "--kick-flips") {
      arguments.kick_flips =
          hopper_positive_int(value(), option, kCoreEntries);
    } else if (option == "--kick-attempts") {
      arguments.kick_attempts = hopper_positive_int(
          value(), option, std::numeric_limits<int>::max());
    } else if (option == "--transient-stride") {
      arguments.transient_stride = hopper_positive_int(
          value(), option, std::numeric_limits<int>::max());
    } else if (option == "--polish-moves") {
      arguments.polish_moves = hopper_positive_int(
          value(), option, std::numeric_limits<int>::max());
    } else if (option == "--archive-capacity") {
      arguments.archive_capacity =
          hopper_positive_int(value(), option, 4096);
    } else if (option == "--archive-outputs") {
      arguments.archive_outputs =
          hopper_positive_int(value(), option, 4096);
    } else if (option == "--quality-tournament") {
      arguments.quality_tournament =
          hopper_positive_int(value(), option, 4096);
    } else if (option == "--self-test") {
      arguments.self_test_rounds =
          hopper_positive_int(value(), option, 1000000);
    } else if (option == "--help") {
      std::cout
          << "usage: core_gram_basin_hopper "
             "--seed-matrix H0 [--seed-matrix H1 ...] "
             "--output MATRIX --archive-dir DIR --log JSONL "
             "--summary JSON [options]\n"
          << "  --seed N --seconds S --heartbeat-seconds S "
             "--max-epochs N\n"
          << "  --quotient-gate Q "
             "[--gate-all] "
             "--epoch-min-moves N --epoch-max-moves N\n"
          << "  --kick-flips K --kick-attempts N "
             "--transient-stride N --polish-moves N\n"
          << "  --archive-capacity N --archive-outputs N "
             "--quality-tournament N [--coronal-pareto]\n"
          << "audit mode: --self-test ROUNDS [--seed N]\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }

  if (arguments.self_test_rounds != 0) {
    return arguments;
  }
  if (arguments.seeds.empty()) {
    throw std::runtime_error(
        "at least one --seed-matrix is required");
  }
  if (
      arguments.output.empty() ||
      arguments.archive_directory.empty() ||
      arguments.log.empty() || arguments.summary.empty()) {
    throw std::runtime_error(
        "--output, --archive-dir, --log, and --summary are required");
  }
  if (
      arguments.epoch_min_moves < 10000 ||
      arguments.epoch_max_moves < 10000) {
    throw std::runtime_error(
        "reactive-tabu epochs must contain at least 10000 moves");
  }
  if (arguments.epoch_min_moves > arguments.epoch_max_moves) {
    throw std::runtime_error(
        "--epoch-min-moves must not exceed --epoch-max-moves");
  }
  if (arguments.archive_outputs > arguments.archive_capacity) {
    throw std::runtime_error(
        "--archive-outputs must not exceed --archive-capacity");
  }
  if (
      arguments.coronal_pareto &&
      arguments.archive_capacity < 3) {
    throw std::runtime_error(
        "--coronal-pareto requires --archive-capacity at least 3 "
        "to preserve the three exact objective extremes");
  }

  const std::array<fs::path, 3> files{
      arguments.output, arguments.log, arguments.summary};
  std::array<fs::path, 3> normalized_files{};
  for (std::size_t first = 0; first < files.size(); ++first) {
    if (path_entry_exists(files[first])) {
      throw std::runtime_error(
          "refusing to overwrite existing output: " +
          files[first].string());
    }
    normalized_files[first] =
        normalized_absolute_path(files[first]);
    for (std::size_t second = first + 1;
         second < files.size(); ++second) {
      const fs::path normalized_second =
          normalized_absolute_path(files[second]);
      if (normalized_files[first] == normalized_second) {
        throw std::runtime_error(
            "output, log, and summary paths must differ");
      }
    }
  }
  if (path_entry_exists(arguments.archive_directory)) {
    throw std::runtime_error(
        "refusing to overwrite archive directory: " +
        arguments.archive_directory.string());
  }
  const fs::path normalized_archive =
      normalized_absolute_path(arguments.archive_directory);
  for (std::size_t index = 0; index < files.size(); ++index) {
    if (
        path_is_within_or_equal(
            normalized_files[index], normalized_archive) ||
        path_is_within_or_equal(
            normalized_archive, normalized_files[index])) {
      throw std::runtime_error(
          "output, log, and summary paths must be outside and "
          "disjoint from the archive directory (including its "
          "reserved manifest.json)");
    }
  }
  for (const fs::path& seed : arguments.seeds) {
    const fs::path normalized_seed =
        normalized_absolute_path(seed);
    for (const fs::path& output : normalized_files) {
      if (normalized_seed == output) {
        throw std::runtime_error(
            "seed and output paths must differ");
      }
    }
  }
  return arguments;
}

State exact_state_from_core(const CoreMatrix& core) {
  State state;
  state.core = core;
  state.determinant = exact_core_determinant(core);
  if (state.determinant == 0) {
    throw std::runtime_error(
        "cannot construct state from singular core");
  }
  state.adjugate = exact_adjugate(core);
  check_adjugate_identity(state);
  return state;
}

State read_exact_state(const fs::path& path) {
  const SignMatrix input = read_sign_matrix(path);
  State state = exact_state_from_core(dephase_to_core(input));
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

using GramHistogram =
    std::array<std::uint16_t, kGramLabelCount>;

std::vector<std::uint16_t> one_sided_gram_sketch(
    const SignMatrix& matrix, bool use_columns) {
  std::array<int, kSignOrder> line_products{};
  for (int line = 0; line < kSignOrder; ++line) {
    int product = 1;
    for (int coordinate = 0; coordinate < kSignOrder;
         ++coordinate) {
      product *= use_columns
          ? matrix[coordinate][line]
          : matrix[line][coordinate];
    }
    line_products[line] = product;
  }

  GramHistogram global{};
  std::array<GramHistogram, kSignOrder> vertices{};
  for (int first = 0; first < kSignOrder; ++first) {
    for (int second = first + 1;
         second < kSignOrder; ++second) {
      int inner_product = 0;
      for (int coordinate = 0; coordinate < kSignOrder;
           ++coordinate) {
        const int first_value = use_columns
            ? matrix[coordinate][first]
            : matrix[first][coordinate];
        const int second_value = use_columns
            ? matrix[coordinate][second]
            : matrix[second][coordinate];
        inner_product += first_value * second_value;
      }
      const int normalized =
          line_products[first] * line_products[second] *
          inner_product;
      const int numerator = kSignOrder - normalized;
      if (
          normalized < -21 || normalized > 23 ||
          numerator < 0 || numerator % 4 != 0) {
        throw std::runtime_error(
            "normalized signed Gram label violated mod-4 invariant");
      }
      const int label = numerator / 4;
      if (label < 0 || label >= kGramLabelCount) {
        throw std::runtime_error(
            "normalized signed Gram label is out of range");
      }
      ++global[static_cast<std::size_t>(label)];
      ++vertices[first][static_cast<std::size_t>(label)];
      ++vertices[second][static_cast<std::size_t>(label)];
    }
  }
  std::sort(vertices.begin(), vertices.end());

  std::vector<std::uint16_t> result;
  result.reserve(kOneSidedSketchLength);
  result.insert(result.end(), global.begin(), global.end());
  for (const GramHistogram& vertex : vertices) {
    result.insert(result.end(), vertex.begin(), vertex.end());
  }
  if (result.size() != kOneSidedSketchLength) {
    throw std::runtime_error("internal Gram-sketch length mismatch");
  }
  return result;
}

GramSketch gram_sketch(const SignMatrix& matrix) {
  std::vector<std::uint16_t> rows =
      one_sided_gram_sketch(matrix, false);
  std::vector<std::uint16_t> columns =
      one_sided_gram_sketch(matrix, true);
  if (columns < rows) std::swap(rows, columns);
  GramSketch result;
  result.labels.reserve(2 * kOneSidedSketchLength);
  result.labels.insert(
      result.labels.end(), rows.begin(), rows.end());
  result.labels.insert(
      result.labels.end(), columns.begin(), columns.end());
  return result;
}

GramSketch gram_sketch(const CoreMatrix& core) {
  return gram_sketch(core_to_sign(core));
}

BigInt big_absolute(BigInt value) {
  return value < 0 ? -value : value;
}

BigInt big_gcd(BigInt left, BigInt right) {
  left = big_absolute(std::move(left));
  right = big_absolute(std::move(right));
  while (right != 0) {
    BigInt remainder = left % right;
    left = std::move(right);
    right = std::move(remainder);
  }
  return left;
}

std::string big_to_string(const BigInt& value) {
  return wide_to_string(value);
}

BigInt decimal_big(std::string_view text) {
  if (text.empty()) {
    throw std::runtime_error("empty exact decimal integer");
  }
  BigInt result = 0;
  for (const unsigned char character : text) {
    if (character < static_cast<unsigned char>('0') ||
        character > static_cast<unsigned char>('9')) {
      throw std::runtime_error("invalid exact decimal integer");
    }
    result =
        10 * result +
        static_cast<int>(
            character - static_cast<unsigned char>('0'));
  }
  return result;
}

BigInt big_bareiss(
    const std::vector<std::vector<BigInt>>& source) {
  if (source.empty()) return 1;
  const std::size_t order = source.size();
  for (const auto& row : source) {
    if (row.size() != order) {
      throw std::runtime_error(
          "big-integer Bareiss matrix must be square");
    }
  }
  std::vector<std::vector<BigInt>> work = source;
  BigInt previous = 1;
  int sign = 1;
  for (std::size_t column = 0; column + 1 < order;
       ++column) {
    std::size_t pivot = column;
    while (pivot < order && work[pivot][column] == 0) ++pivot;
    if (pivot == order) return 0;
    if (pivot != column) {
      std::swap(work[pivot], work[column]);
      sign = -sign;
    }
    const BigInt pivot_value = work[column][column];
    for (std::size_t row = column + 1; row < order; ++row) {
      for (std::size_t inner = column + 1; inner < order;
           ++inner) {
        const BigInt numerator =
            work[row][inner] * pivot_value -
            work[row][column] * work[column][inner];
        if (numerator % previous != 0) {
          throw std::runtime_error(
              "non-exact big-integer Bareiss division");
        }
        work[row][inner] = numerator / previous;
      }
      work[row][column] = 0;
    }
    previous = pivot_value;
  }
  return sign < 0 ? -work.back().back() : work.back().back();
}

int compare_kappa(
    const CoronalPoint& first, const CoronalPoint& second) {
  const BigInt left =
      first.kappa_numerator * second.kappa_denominator;
  const BigInt right =
      second.kappa_numerator * first.kappa_denominator;
  if (left < right) return -1;
  if (left > right) return 1;
  return 0;
}

bool coronal_point_less(
    const CoronalPoint& first, const CoronalPoint& second) {
  if (first.determinant_m != second.determinant_m) {
    return first.determinant_m < second.determinant_m;
  }
  return compare_kappa(first, second) < 0;
}

CoronalPoint oriented_coronal_point(
    const SignMatrix& matrix, bool use_columns,
    std::uint64_t core_quotient) {
  std::array<int, kSignOrder> line_products{};
  for (int line = 0; line < kSignOrder; ++line) {
    int product = 1;
    for (int coordinate = 0; coordinate < kSignOrder;
         ++coordinate) {
      product *= use_columns
          ? matrix[coordinate][line]
          : matrix[line][coordinate];
    }
    line_products[line] = product;
  }

  std::vector<std::vector<BigInt>> matrix_m(
      kSignOrder, std::vector<BigInt>(kSignOrder));
  for (int first = 0; first < kSignOrder; ++first) {
    for (int second = 0; second < kSignOrder; ++second) {
      int inner_product = 0;
      for (int coordinate = 0; coordinate < kSignOrder;
           ++coordinate) {
        const int first_value = use_columns
            ? matrix[coordinate][first]
            : matrix[first][coordinate];
        const int second_value = use_columns
            ? matrix[coordinate][second]
            : matrix[second][coordinate];
        inner_product += first_value * second_value;
      }
      const int normalized =
          line_products[first] * line_products[second] *
          inner_product;
      if (first == second) {
        if (normalized != kSignOrder) {
          throw std::runtime_error(
              "normalized Gram diagonal is not 23");
        }
      } else if (
          normalized < -21 || normalized > 23 ||
          (normalized - 3) % 4 != 0) {
        throw std::runtime_error(
            "coronal Gram entry violated the 3 mod 4 invariant");
      }
      const int shifted = normalized + 1;
      if (shifted % 4 != 0) {
        throw std::runtime_error(
            "coronal M entry is not integral");
      }
      matrix_m[first][second] = shifted / 4;
    }
  }

  CoronalPoint point;
  point.determinant_m = big_bareiss(matrix_m);
  if (point.determinant_m <= 0) {
    throw std::runtime_error(
        "nonsingular sign matrix produced non-positive det(M)");
  }
  const BigInt quotient = core_quotient;
  const BigInt quotient_squared = quotient * quotient;
  BigInt numerator =
      4 * point.determinant_m - quotient_squared;
  if (numerator <= 0 || numerator >= 4 * point.determinant_m) {
    throw std::runtime_error(
        "derived coronal kappa is outside the exact interval (0,4)");
  }
  const BigInt divisor = big_gcd(numerator, point.determinant_m);
  point.kappa_numerator = numerator / divisor;
  point.kappa_denominator = point.determinant_m / divisor;
  if (
      point.determinant_m *
              (4 * point.kappa_denominator -
               point.kappa_numerator) !=
          quotient_squared * point.kappa_denominator) {
    throw std::runtime_error(
        "exact q^2 = det(M) * (4-kappa) check failed");
  }
  return point;
}

CoronalPair coronal_pair(
    const SignMatrix& matrix, std::uint64_t core_quotient) {
  CoronalPair pair{
      {oriented_coronal_point(
           matrix, false, core_quotient),
       oriented_coronal_point(
           matrix, true, core_quotient)}};
  if (coronal_point_less(
          pair.orientations[1], pair.orientations[0])) {
    std::swap(pair.orientations[0], pair.orientations[1]);
  }
  return pair;
}

CoronalPair coronal_pair(
    const CoreMatrix& core, std::int64_t determinant) {
  return coronal_pair(core_to_sign(core), magnitude(determinant));
}

const CoronalPoint& maximum_coronal_det_m(
    const CoronalPair& pair) {
  return pair.orientations[1];
}

const CoronalPoint& minimum_coronal_kappa(
    const CoronalPair& pair) {
  return compare_kappa(
             pair.orientations[0], pair.orientations[1]) <= 0
      ? pair.orientations[0]
      : pair.orientations[1];
}

CoronalPair evaluate_coronal_pair(
    const CoreMatrix& core, std::int64_t determinant,
    HopperStatistics& statistics) {
  CoronalPair result = coronal_pair(core, determinant);
  ++statistics.coronal_pair_evaluations;
  statistics.coronal_formula_checks += 2;
  return result;
}

std::uint64_t stable_sketch_hash(const GramSketch& sketch) {
  std::uint64_t hash = UINT64_C(0x6a09e667f3bcc909);
  for (const std::uint16_t label : sketch.labels) {
    hash = splitmix64(
        hash ^ static_cast<std::uint64_t>(label));
  }
  return hash;
}

std::string hexadecimal(std::uint64_t value) {
  std::ostringstream output;
  output << std::hex << std::setw(16) << std::setfill('0')
         << value;
  return output.str();
}

SignMatrix transpose_sign_matrix(const SignMatrix& matrix) {
  SignMatrix result{};
  for (int row = 0; row < kSignOrder; ++row) {
    for (int column = 0; column < kSignOrder; ++column) {
      result[row][column] = matrix[column][row];
    }
  }
  return result;
}

int hexadecimal_digit(char character) {
  if (character >= '0' && character <= '9') {
    return character - '0';
  }
  if (character >= 'a' && character <= 'f') {
    return 10 + character - 'a';
  }
  throw std::runtime_error(
      "invalid hexadecimal digit in pinned core");
}

CoreMatrix core_from_hex(std::string_view encoded) {
  constexpr std::size_t kEncodedLength =
      static_cast<std::size_t>(kCoreEntries / 4);
  static_assert(kCoreEntries % 4 == 0);
  if (encoded.size() != kEncodedLength) {
    throw std::runtime_error(
        "pinned core hexadecimal length mismatch");
  }
  CoreMatrix core{};
  int flat = 0;
  for (const char character : encoded) {
    const int nibble = hexadecimal_digit(character);
    for (int shift = 3; shift >= 0; --shift) {
      core[flat / kCoreOrder][flat % kCoreOrder] =
          static_cast<std::uint8_t>((nibble >> shift) & 1);
      ++flat;
    }
  }
  if (flat != kCoreEntries) {
    throw std::runtime_error(
        "pinned core hexadecimal decode mismatch");
  }
  return core;
}

void check_pinned_coronal(
    std::string_view encoded_core, std::uint64_t expected_q,
    std::string_view expected_det_m,
    std::string_view expected_kappa_numerator,
    std::string_view expected_kappa_denominator) {
  const CoreMatrix core = core_from_hex(encoded_core);
  const std::int64_t determinant =
      exact_core_determinant(core);
  if (magnitude(determinant) != expected_q) {
    throw std::runtime_error(
        "pinned coronal core quotient mismatch");
  }
  const CoronalPair pair = coronal_pair(core, determinant);
  const CoronalPoint expected{
      decimal_big(expected_det_m),
      decimal_big(expected_kappa_numerator),
      decimal_big(expected_kappa_denominator)};
  for (const CoronalPoint& point : pair.orientations) {
    if (!(point == expected)) {
      throw std::runtime_error(
          "pinned exact det(M)/kappa coronal check failed");
    }
  }
}

void run_pinned_coronal_self_tests() {
  check_pinned_coronal(
      "70c3ffcf300df0333bc0333f7198fec997c2da5f0796c4f16b23ca50"
      "ecea83735543370d0ce3ca56a3eaa64f9a647e5652fa6a7d2a55f899"
      "abca596f1",
      UINT64_C(662671875), "267647730468750000", "6455",
      "2736");
  check_pinned_coronal(
      "0607f7a303d25feb2a27da955179c0d9776e56334d7019eb57563894e"
      "287af9705790db872c6bb338e74a6e4a3a51dad4c9e46b8cfe0b243dc"
      "ce699aa",
      UINT64_C(622339200), "293534443238400000", "2962",
      "1105");

  const CoreMatrix asymmetric_core = core_from_hex(
      "0931b0d6bcd01bd8e3d6c71e63e2f4450a8294704d94216f274ad0d29"
      "2f08c036043853805eb382280c1fd490c6059377435c27dd316ffa80702"
      "087b4");
  const std::int64_t asymmetric_determinant =
      exact_core_determinant(asymmetric_core);
  if (magnitude(asymmetric_determinant) != UINT64_C(29302)) {
    throw std::runtime_error(
        "pinned asymmetric coronal quotient mismatch");
  }
  const CoronalPair asymmetric =
      coronal_pair(asymmetric_core, asymmetric_determinant);
  const CoronalPair asymmetric_expected{{
      CoronalPoint{
          INT64_C(1677254768), INT64_C(208943281),
          INT64_C(59901956)},
      CoronalPoint{
          INT64_C(4473391164), INT64_C(8050547),
          INT64_C(2114079)}}};
  if (!(asymmetric == asymmetric_expected)) {
    throw std::runtime_error(
        "pinned asymmetric row/column coronal pair failed");
  }
  if (!(
          coronal_pair(
              transpose_sign_matrix(
                  core_to_sign(asymmetric_core)),
              magnitude(asymmetric_determinant)) ==
          asymmetric)) {
    throw std::runtime_error(
        "pinned asymmetric coronal transpose check failed");
  }
}

void run_coronal_selection_self_tests();

void run_hopper_logic_self_tests() {
  HopperArguments arguments;
  if (
      archive_gate_applies(
          arguments, ObservationKind::kSeed, true) ||
      archive_gate_applies(
          arguments, ObservationKind::kKick, false) ||
      archive_gate_applies(
          arguments, ObservationKind::kPeak, false) ||
      !archive_gate_applies(
          arguments, ObservationKind::kTransient, false) ||
      !archive_gate_applies(
          arguments, ObservationKind::kPolished, false)) {
    throw std::runtime_error(
        "default archive-gate scope self-test failed");
  }
  arguments.gate_all = true;
  if (
      archive_gate_applies(
          arguments, ObservationKind::kSeed, true) ||
      !archive_gate_applies(
          arguments, ObservationKind::kKick, false) ||
      !archive_gate_applies(
          arguments, ObservationKind::kPeak, false) ||
      !archive_gate_applies(
          arguments, ObservationKind::kTransient, false) ||
      !archive_gate_applies(
          arguments, ObservationKind::kPolished, false)) {
    throw std::runtime_error(
        "gate-all scope self-test failed");
  }
  if (
      search_added(UINT64_C(7), UINT64_C(3), true) !=
          UINT64_C(4) ||
      search_added(UINT64_C(7), UINT64_C(3), false) !=
          UINT64_C(0)) {
    throw std::runtime_error(
        "initialization-baseline self-test failed");
  }
  if (
      !path_is_within_or_equal(
          fs::path("/tmp/hopper/archive/elite"),
          fs::path("/tmp/hopper/archive")) ||
      !path_is_within_or_equal(
          fs::path("/tmp/hopper/archive"),
          fs::path("/tmp/hopper/archive")) ||
      path_is_within_or_equal(
          fs::path("/tmp/hopper/archive-other"),
          fs::path("/tmp/hopper/archive"))) {
    throw std::runtime_error(
        "path-containment self-test failed");
  }
}

void run_descriptor_self_tests(
    int rounds, std::uint64_t seed,
    HopperStatistics& statistics, bool test_coronal) {
  run_hopper_logic_self_tests();
  if (test_coronal) {
    run_pinned_coronal_self_tests();
    run_coronal_selection_self_tests();
  }
  std::mt19937_64 randomizer(seed);
  int coronal_rounds = 0;
  for (int round = 0; round < rounds; ++round) {
    SignMatrix matrix{};
    for (auto& row : matrix) {
      for (int& value : row) {
        value =
            (randomizer() & UINT64_C(1)) == 0U ? -1 : 1;
      }
    }
    const GramSketch expected = gram_sketch(matrix);

    std::array<int, kSignOrder> row_permutation{};
    std::array<int, kSignOrder> column_permutation{};
    std::iota(
        row_permutation.begin(), row_permutation.end(), 0);
    std::iota(
        column_permutation.begin(), column_permutation.end(), 0);
    std::shuffle(
        row_permutation.begin(), row_permutation.end(),
        randomizer);
    std::shuffle(
        column_permutation.begin(), column_permutation.end(),
        randomizer);
    std::array<int, kSignOrder> row_signs{};
    std::array<int, kSignOrder> column_signs{};
    for (int index = 0; index < kSignOrder; ++index) {
      row_signs[index] =
          (randomizer() & UINT64_C(1)) == 0U ? -1 : 1;
      column_signs[index] =
          (randomizer() & UINT64_C(1)) == 0U ? -1 : 1;
    }

    SignMatrix transformed{};
    for (int row = 0; row < kSignOrder; ++row) {
      for (int column = 0; column < kSignOrder; ++column) {
        transformed[row][column] =
            row_signs[row] * column_signs[column] *
            matrix[row_permutation[row]]
                  [column_permutation[column]];
      }
    }
    if (!(gram_sketch(transformed) == expected)) {
      throw std::runtime_error(
          "Gram sketch failed signed-permutation invariance");
    }
    if (!(gram_sketch(transpose_sign_matrix(matrix)) == expected)) {
      throw std::runtime_error(
          "Gram sketch failed transpose invariance");
    }
    if (!(
            gram_sketch(transpose_sign_matrix(transformed)) ==
            expected)) {
      throw std::runtime_error(
          "Gram sketch failed combined invariance");
    }
    statistics.descriptor_invariance_checks += 3;

    if (test_coronal && coronal_rounds < 4) {
      const CoreMatrix core = dephase_to_core(matrix);
      const std::int64_t determinant =
          exact_core_determinant(core);
      if (determinant != 0) {
        const CoronalPair coronal =
            coronal_pair(matrix, magnitude(determinant));
        if (!(
                coronal_pair(
                    transformed, magnitude(determinant)) ==
                coronal)) {
          throw std::runtime_error(
              "coronal pair failed signed-permutation invariance");
        }
        if (!(
                coronal_pair(
                    transpose_sign_matrix(matrix),
                    magnitude(determinant)) == coronal)) {
          throw std::runtime_error(
              "coronal pair failed transpose invariance");
        }
        statistics.coronal_invariance_checks += 2;
        ++coronal_rounds;
      }
    }
  }
  if (test_coronal && rounds > 0 && coronal_rounds == 0) {
    throw std::runtime_error(
        "coronal invariance self-test found no nonsingular matrix");
  }
}

std::string flatten_core(const CoreMatrix& core) {
  std::string result;
  result.reserve(kCoreEntries);
  for (const auto& row : core) {
    for (const std::uint8_t value : row) {
      result.push_back(static_cast<char>('0' + value));
    }
  }
  return result;
}

bool elite_is_better(
    std::int64_t determinant, const CoreMatrix& core,
    const SketchElite& incumbent) {
  const std::uint64_t candidate_magnitude =
      magnitude(determinant);
  const std::uint64_t incumbent_magnitude =
      magnitude(incumbent.determinant);
  if (candidate_magnitude != incumbent_magnitude) {
    return candidate_magnitude > incumbent_magnitude;
  }
  return flatten_core(core) < flatten_core(incumbent.core);
}

const CoronalPair& require_coronal(const SketchElite& elite) {
  if (!elite.coronal.has_value()) {
    throw std::runtime_error(
        "coronal-Pareto archive elite lacks an exact descriptor");
  }
  return *elite.coronal;
}

bool coronal_det_m_is_larger(
    const SketchElite& first, const SketchElite& second) {
  return maximum_coronal_det_m(require_coronal(first))
             .determinant_m >
      maximum_coronal_det_m(require_coronal(second))
          .determinant_m;
}

bool coronal_kappa_is_smaller(
    const SketchElite& first, const SketchElite& second) {
  return compare_kappa(
             minimum_coronal_kappa(require_coronal(first)),
             minimum_coronal_kappa(require_coronal(second))) < 0;
}

bool coronal_dominates(
    const SketchElite& first, const SketchElite& second) {
  const std::uint64_t first_quality =
      magnitude(first.determinant);
  const std::uint64_t second_quality =
      magnitude(second.determinant);
  const BigInt& first_det_m =
      maximum_coronal_det_m(require_coronal(first))
          .determinant_m;
  const BigInt& second_det_m =
      maximum_coronal_det_m(require_coronal(second))
          .determinant_m;
  const int kappa_comparison = compare_kappa(
      minimum_coronal_kappa(require_coronal(first)),
      minimum_coronal_kappa(require_coronal(second)));
  const bool no_worse =
      first_quality >= second_quality &&
      first_det_m >= second_det_m &&
      kappa_comparison <= 0;
  const bool strict =
      first_quality > second_quality ||
      first_det_m > second_det_m ||
      kappa_comparison < 0;
  return no_worse && strict;
}

std::size_t coronal_pareto_front_size(
    const std::map<GramSketch, SketchElite>& archive) {
  std::size_t count = 0;
  for (auto candidate = archive.begin();
       candidate != archive.end(); ++candidate) {
    bool dominated = false;
    for (auto other = archive.begin(); other != archive.end();
         ++other) {
      if (
          candidate != other &&
          coronal_dominates(other->second, candidate->second)) {
        dominated = true;
        break;
      }
    }
    if (!dominated) ++count;
  }
  return count;
}

std::map<GramSketch, SketchElite>::iterator
worst_coronal_elite(
    std::map<GramSketch, SketchElite>& archive) {
  if (archive.size() < 4) {
    throw std::runtime_error(
        "coronal-Pareto eviction needs at least four elites");
  }
  using Iterator =
      std::map<GramSketch, SketchElite>::iterator;
  std::vector<Iterator> entries;
  entries.reserve(archive.size());
  for (auto iterator = archive.begin();
       iterator != archive.end(); ++iterator) {
    entries.push_back(iterator);
  }

  std::vector<std::size_t> dominated_by(entries.size(), 0);
  for (std::size_t first = 0; first < entries.size(); ++first) {
    for (std::size_t second = 0; second < entries.size();
         ++second) {
      if (
          first != second &&
          coronal_dominates(
              entries[second]->second, entries[first]->second)) {
        ++dominated_by[first];
      }
    }
  }

  std::size_t quality_extreme = 0;
  std::size_t det_m_extreme = 0;
  std::size_t kappa_extreme = 0;
  for (std::size_t index = 1; index < entries.size(); ++index) {
    if (
        magnitude(entries[index]->second.determinant) >
        magnitude(entries[quality_extreme]->second.determinant)) {
      quality_extreme = index;
    }
    if (coronal_det_m_is_larger(
            entries[index]->second,
            entries[det_m_extreme]->second)) {
      det_m_extreme = index;
    }
    if (coronal_kappa_is_smaller(
            entries[index]->second,
            entries[kappa_extreme]->second)) {
      kappa_extreme = index;
    }
  }

  auto protected_extreme = [&](std::size_t index) {
    return index == quality_extreme || index == det_m_extreme ||
        index == kappa_extreme;
  };
  auto is_worse = [&](std::size_t candidate,
                      std::size_t incumbent) {
    if (dominated_by[candidate] != dominated_by[incumbent]) {
      return dominated_by[candidate] > dominated_by[incumbent];
    }
    const SketchElite& candidate_elite =
        entries[candidate]->second;
    const SketchElite& incumbent_elite =
        entries[incumbent]->second;
    const std::uint64_t candidate_quality =
        magnitude(candidate_elite.determinant);
    const std::uint64_t incumbent_quality =
        magnitude(incumbent_elite.determinant);
    if (candidate_quality != incumbent_quality) {
      return candidate_quality < incumbent_quality;
    }
    const BigInt& candidate_det_m =
        maximum_coronal_det_m(require_coronal(candidate_elite))
            .determinant_m;
    const BigInt& incumbent_det_m =
        maximum_coronal_det_m(require_coronal(incumbent_elite))
            .determinant_m;
    if (candidate_det_m != incumbent_det_m) {
      return candidate_det_m < incumbent_det_m;
    }
    const int kappa_comparison = compare_kappa(
        minimum_coronal_kappa(require_coronal(candidate_elite)),
        minimum_coronal_kappa(require_coronal(incumbent_elite)));
    if (kappa_comparison != 0) return kappa_comparison > 0;
    if (candidate_elite.selections != incumbent_elite.selections) {
      return candidate_elite.selections >
          incumbent_elite.selections;
    }
    return candidate_elite.last_updated_epoch <
        incumbent_elite.last_updated_epoch;
  };

  std::optional<std::size_t> worst;
  for (std::size_t index = 0; index < entries.size(); ++index) {
    if (protected_extreme(index)) continue;
    if (!worst.has_value() || is_worse(index, *worst)) {
      worst = index;
    }
  }
  if (!worst.has_value()) {
    throw std::runtime_error(
        "coronal-Pareto capacity cannot preserve its extremes");
  }
  return entries[*worst];
}

void run_coronal_selection_self_tests() {
  auto make_elite = [](
                        std::int64_t determinant,
                        std::int64_t determinant_m,
                        std::int64_t kappa_numerator,
                        std::int64_t kappa_denominator,
                        std::string origin) {
    const CoronalPoint point{
        determinant_m, kappa_numerator, kappa_denominator};
    CoronalPair pair{{point, point}};
    SketchElite elite;
    elite.determinant = determinant;
    elite.origin = std::move(origin);
    elite.coronal = std::move(pair);
    return elite;
  };
  std::map<GramSketch, SketchElite> archive;
  archive.emplace(
      GramSketch{{0}}, make_elite(100, 100, 2, 1, "quality"));
  archive.emplace(
      GramSketch{{1}}, make_elite(90, 120, 3, 1, "det_m"));
  archive.emplace(
      GramSketch{{2}}, make_elite(80, 80, 1, 1, "kappa"));
  archive.emplace(
      GramSketch{{3}}, make_elite(70, 70, 3, 1, "dominated"));
  if (coronal_pareto_front_size(archive) != 3) {
    throw std::runtime_error(
        "coronal Pareto-front self-test failed");
  }
  const auto victim = worst_coronal_elite(archive);
  if (victim->second.origin != "dominated") {
    throw std::runtime_error(
        "coronal extreme-preserving eviction self-test failed");
  }
}

std::map<GramSketch, SketchElite>::iterator worst_elite(
    std::map<GramSketch, SketchElite>& archive) {
  if (archive.empty()) {
    throw std::runtime_error("cannot rank an empty archive");
  }
  auto worst = archive.begin();
  for (auto iterator = std::next(archive.begin());
       iterator != archive.end(); ++iterator) {
    const auto worst_key = std::make_tuple(
        magnitude(worst->second.determinant),
        std::numeric_limits<std::uint64_t>::max() -
            worst->second.selections,
        worst->second.last_updated_epoch);
    const auto candidate_key = std::make_tuple(
        magnitude(iterator->second.determinant),
        std::numeric_limits<std::uint64_t>::max() -
            iterator->second.selections,
        iterator->second.last_updated_epoch);
    if (candidate_key < worst_key) worst = iterator;
  }
  return worst;
}

const char* hopper_engine_name(
    const HopperArguments& arguments) {
  return arguments.coronal_pareto
      ? "exact-core-gram-sketch-basin-hopper-coronal-pareto-v1"
      : "exact-core-gram-sketch-basin-hopper-v1";
}

std::string coronal_point_json(const CoronalPoint& point) {
  std::ostringstream output;
  output
      << "{\"det_m\":\""
      << big_to_string(point.determinant_m)
      << "\",\"kappa_numerator\":\""
      << big_to_string(point.kappa_numerator)
      << "\",\"kappa_denominator\":\""
      << big_to_string(point.kappa_denominator) << "\"}";
  return output.str();
}

std::string coronal_pair_json(const CoronalPair& pair) {
  return "[" + coronal_point_json(pair.orientations[0]) +
      "," + coronal_point_json(pair.orientations[1]) + "]";
}

std::string parent_probability_json(
    const HopperArguments& arguments) {
  if (arguments.coronal_pareto) {
    return
        "      \"quality\": 0.30,\n"
        "      \"coronal_det_m\": 0.20,\n"
        "      \"coronal_kappa\": 0.20,\n"
        "      \"coronal_pareto\": 0.20,\n"
        "      \"least_used\": 0.05,\n"
        "      \"uniform\": 0.05\n";
  }
  return
      "      \"quality\": 0.50,\n"
      "      \"coronal_det_m\": 0.00,\n"
      "      \"coronal_kappa\": 0.00,\n"
      "      \"coronal_pareto\": 0.00,\n"
      "      \"least_used\": 0.35,\n"
      "      \"uniform\": 0.15\n";
}

std::string coronal_archive_json(
    const HopperArguments& arguments,
    const HopperStatistics& statistics,
    const std::map<GramSketch, SketchElite>& archive) {
  std::ostringstream output;
  output
      << "  \"coronal_pareto\": {\n"
      << "    \"enabled\": "
      << (arguments.coronal_pareto ? "true" : "false")
      << ",\n"
      << "    \"identity\": "
         "\"W=(G-(24I-J))/4; M=6I+W; "
         "q^2=det(M)*(4-kappa)\",\n"
      << "    \"kappa_derivation\": "
         "\"4-q^2/det(M); no inverse\",\n"
      << "    \"arithmetic\": "
         "\"exact signed-128 Bareiss and rational cross-products\",\n"
      << "    \"selection_objectives\": "
         "[\"larger_core_quotient\",\"larger_det_m\","
         "\"smaller_kappa\"],\n"
      << "    \"exact_pair_evaluations\": "
      << statistics.coronal_pair_evaluations << ",\n"
      << "    \"exact_formula_checks\": "
      << statistics.coronal_formula_checks << ",\n"
      << "    \"invariance_checks\": "
      << statistics.coronal_invariance_checks << ",\n"
      << "    \"capacity_rejections\": "
      << statistics.coronal_capacity_rejections << ",\n"
      << "    \"pareto_evictions\": "
      << statistics.coronal_pareto_evictions;
  if (!arguments.coronal_pareto || archive.empty()) {
    output
        << ",\n"
        << "    \"pareto_front_size\": 0,\n"
        << "    \"orientation_pair\": null\n"
        << "  },\n";
    return output.str();
  }

  auto quality = archive.begin();
  auto det_m = archive.begin();
  auto kappa = archive.begin();
  for (auto iterator = std::next(archive.begin());
       iterator != archive.end(); ++iterator) {
    if (
        magnitude(iterator->second.determinant) >
        magnitude(quality->second.determinant)) {
      quality = iterator;
    }
    if (coronal_det_m_is_larger(
            iterator->second, det_m->second)) {
      det_m = iterator;
    }
    if (coronal_kappa_is_smaller(
            iterator->second, kappa->second)) {
      kappa = iterator;
    }
  }
  output
      << ",\n"
      << "    \"pareto_front_size\": "
      << coronal_pareto_front_size(archive) << ",\n"
      << "    \"orientation_pair\": "
         "\"sorted_row_column_coronal_points\",\n"
      << "    \"quality_extreme\": {\"core_quotient\": "
      << magnitude(quality->second.determinant)
      << ",\"points\":"
      << coronal_pair_json(require_coronal(quality->second))
      << "},\n"
      << "    \"det_m_extreme\": {\"core_quotient\": "
      << magnitude(det_m->second.determinant)
      << ",\"point\":"
      << coronal_point_json(maximum_coronal_det_m(
             require_coronal(det_m->second)))
      << "},\n"
      << "    \"kappa_extreme\": {\"core_quotient\": "
      << magnitude(kappa->second.determinant)
      << ",\"point\":"
      << coronal_point_json(minimum_coronal_kappa(
             require_coronal(kappa->second)))
      << "}\n"
      << "  },\n";
  return output.str();
}

std::string hopper_log_json(
    std::string_view event, double elapsed,
    const HopperArguments& arguments,
    const HopperStatistics& statistics,
    const std::map<GramSketch, SketchElite>& archive,
    std::uint64_t best_magnitude) {
  const std::uint64_t initialization_sketches =
      statistics.initialization_complete
          ? statistics.initialization_sketch_discoveries
          : statistics.sketch_discoveries;
  const std::uint64_t initialization_promotions =
      statistics.initialization_complete
          ? statistics.initialization_promotions
          : statistics.promotions;
  const double directions_per_second =
      elapsed > 0.0
          ? static_cast<double>(
                statistics.kernel.candidate_evaluations) /
                elapsed
          : 0.0;
  std::ostringstream output;
  output
      << "{\"event\":\"" << event << "\""
      << ",\"engine\":\"" << hopper_engine_name(arguments)
      << "\""
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(6) << elapsed
      << ",\"seed\":" << arguments.seed
      << ",\"epochs_started\":" << statistics.epochs_started
      << ",\"epochs_completed\":"
      << statistics.epochs_completed
      << ",\"tabu_moves\":" << statistics.tabu_moves
      << ",\"polish_moves\":" << statistics.polish_moves
      << ",\"candidate_evaluations\":"
      << statistics.kernel.candidate_evaluations
      << ",\"directions_per_second\":"
      << directions_per_second
      << ",\"kick_evaluations\":"
      << statistics.kick_evaluations
      << ",\"observations\":" << statistics.observations
      << ",\"sketch_discoveries\":"
      << statistics.sketch_discoveries
      << ",\"initialization_sketch_discoveries\":"
      << initialization_sketches
      << ",\"search_added_sketch_discoveries\":"
      << search_added(
             statistics.sketch_discoveries,
             initialization_sketches,
             statistics.initialization_complete)
      << ",\"archive_size\":" << archive.size()
      << ",\"archive_insertions\":"
      << statistics.archive_insertions
      << ",\"archive_replacements\":"
      << statistics.archive_replacements
      << ",\"archive_evictions\":"
      << statistics.archive_evictions
      << ",\"coronal_pareto\":"
      << (arguments.coronal_pareto ? "true" : "false")
      << ",\"coronal_pair_evaluations\":"
      << statistics.coronal_pair_evaluations
      << ",\"coronal_formula_checks\":"
      << statistics.coronal_formula_checks
      << ",\"parent_coronal_det_m\":"
      << statistics.parent_coronal_det_m
      << ",\"parent_coronal_kappa\":"
      << statistics.parent_coronal_kappa
      << ",\"parent_coronal_pareto\":"
      << statistics.parent_coronal_pareto
      << ",\"coronal_capacity_rejections\":"
      << statistics.coronal_capacity_rejections
      << ",\"coronal_pareto_evictions\":"
      << statistics.coronal_pareto_evictions
      << ",\"promotions\":" << statistics.promotions
      << ",\"initialization_promotions\":"
      << initialization_promotions
      << ",\"search_added_promotions\":"
      << search_added(
             statistics.promotions, initialization_promotions,
             statistics.initialization_complete)
      << ",\"quotient_gate\":" << arguments.quotient_gate
      << ",\"gate_scope\":\""
      << gate_scope_name(arguments) << "\""
      << ",\"strict_target_states\":"
      << statistics.strict_target_states
      << ",\"best_core_quotient\":" << best_magnitude
      << ",\"strict_target\":" << kStrictCoreTarget
      << ",\"best_absolute_determinant\":\""
      << wide_to_string(
             static_cast<Wide>(best_magnitude) * kScale)
      << "\"}\n";
  return output.str();
}

std::string hopper_summary_json(
    const HopperArguments& arguments,
    const HopperStatistics& statistics,
    const std::map<GramSketch, SketchElite>& archive,
    double elapsed, std::uint64_t best_magnitude,
    bool complete, std::string_view reason,
    const std::vector<fs::path>& archive_outputs) {
  const std::uint64_t initialization_sketches =
      statistics.initialization_complete
          ? statistics.initialization_sketch_discoveries
          : statistics.sketch_discoveries;
  const std::uint64_t initialization_promotions =
      statistics.initialization_complete
          ? statistics.initialization_promotions
          : statistics.promotions;
  const std::uint64_t initialization_archive_size =
      statistics.initialization_complete
          ? statistics.initialization_archive_size
          : static_cast<std::uint64_t>(archive.size());
  const std::uint64_t all_moves =
      statistics.tabu_moves + statistics.polish_moves;
  const double moves_per_second =
      elapsed > 0.0
          ? static_cast<double>(all_moves) / elapsed
          : 0.0;
  const double directions_per_second =
      elapsed > 0.0
          ? static_cast<double>(
                statistics.kernel.candidate_evaluations) /
                elapsed
          : 0.0;
  std::ostringstream output;
  output
      << "{\n"
      << "  \"schema_version\": 1,\n"
      << "  \"engine\": \""
      << hopper_engine_name(arguments) << "\",\n"
      << "  \"complete\": " << (complete ? "true" : "false")
      << ",\n"
      << "  \"reason\": \"" << json_escape(reason) << "\",\n"
      << "  \"elapsed_seconds\": " << std::fixed
      << std::setprecision(6) << elapsed << ",\n"
      << "  \"seed\": " << arguments.seed << ",\n"
      << "  \"known_core_frontier\": " << kKnownCoreFrontier
      << ",\n"
      << "  \"strict_core_target\": " << kStrictCoreTarget
      << ",\n"
      << "  \"best_core_quotient\": " << best_magnitude
      << ",\n"
      << "  \"best_absolute_determinant\": \""
      << wide_to_string(
             static_cast<Wide>(best_magnitude) * kScale)
      << "\",\n"
      << "  \"strict_target_reached_unverified\": "
      << (best_magnitude >= kStrictCoreTarget ? "true" : "false")
      << ",\n"
      << "  \"epochs_started\": " << statistics.epochs_started
      << ",\n"
      << "  \"epochs_completed\": "
      << statistics.epochs_completed << ",\n"
      << "  \"tabu_moves\": " << statistics.tabu_moves
      << ",\n"
      << "  \"polish_moves\": " << statistics.polish_moves
      << ",\n"
      << "  \"moves_per_second\": " << moves_per_second
      << ",\n"
      << "  \"candidate_evaluations\": "
      << statistics.kernel.candidate_evaluations << ",\n"
      << "  \"directions_per_second\": "
      << directions_per_second << ",\n"
      << "  \"kick_evaluations\": "
      << statistics.kick_evaluations << ",\n"
      << "  \"kick_singular\": " << statistics.kick_singular
      << ",\n"
      << "  \"kick_endpoints\": " << statistics.kick_endpoints
      << ",\n"
      << "  \"parent_selection\": {\n"
      << "    \"configured_probabilities\": {\n"
      << parent_probability_json(arguments)
      << "    },\n"
      << "    \"quality\": " << statistics.parent_quality
      << ",\n"
      << "    \"coronal_det_m\": "
      << statistics.parent_coronal_det_m << ",\n"
      << "    \"coronal_kappa\": "
      << statistics.parent_coronal_kappa << ",\n"
      << "    \"coronal_pareto\": "
      << statistics.parent_coronal_pareto << ",\n"
      << "    \"least_used\": "
      << statistics.parent_least_used << ",\n"
      << "    \"uniform\": " << statistics.parent_uniform
      << "\n"
      << "  },\n"
      << "  \"observations\": " << statistics.observations
      << ",\n"
      << "  \"observation_origins\": {\n"
      << "    \"seed\": " << statistics.seed_observations
      << ",\n"
      << "    \"kick_before_repair\": "
      << statistics.kick_observations << ",\n"
      << "    \"transient\": "
      << statistics.transient_observations << ",\n"
      << "    \"epoch_peak\": "
      << statistics.peak_observations << ",\n"
      << "    \"polished\": "
      << statistics.polished_observations << "\n"
      << "  },\n"
      << "  \"quotient_gate\": " << arguments.quotient_gate
      << ",\n"
      << "  \"gate_all\": "
      << (arguments.gate_all ? "true" : "false") << ",\n"
      << "  \"gate_scope\": \""
      << gate_scope_name(arguments) << "\",\n"
      << "  \"gate_rejections\": "
      << statistics.gate_rejections << ",\n"
      << "  \"sketch_evaluations\": "
      << statistics.sketch_evaluations << ",\n"
      << "  \"sketch_discoveries\": "
      << statistics.sketch_discoveries << ",\n"
      << "  \"initialization_baseline\": {\n"
      << "    \"sketch_discoveries\": "
      << initialization_sketches << ",\n"
      << "    \"promotions\": "
      << initialization_promotions << ",\n"
      << "    \"archive_size\": "
      << initialization_archive_size << "\n"
      << "  },\n"
      << "  \"search_added_sketch_discoveries\": "
      << search_added(
             statistics.sketch_discoveries,
             initialization_sketches,
             statistics.initialization_complete)
      << ",\n"
      << "  \"archive_size\": " << archive.size() << ",\n"
      << "  \"archive_capacity\": "
      << arguments.archive_capacity << ",\n"
      << "  \"archive_insertions\": "
      << statistics.archive_insertions << ",\n"
      << "  \"archive_replacements\": "
      << statistics.archive_replacements << ",\n"
      << "  \"archive_rejections\": "
      << statistics.archive_rejections << ",\n"
      << "  \"archive_capacity_rejections\": "
      << statistics.archive_capacity_rejections << ",\n"
      << "  \"archive_evictions\": "
      << statistics.archive_evictions << ",\n"
      << "  \"archive_output_count\": "
      << archive_outputs.size() << ",\n"
      << "  \"promotions\": " << statistics.promotions
      << ",\n"
      << "  \"search_added_promotions\": "
      << search_added(
             statistics.promotions, initialization_promotions,
             statistics.initialization_complete)
      << ",\n"
      << "  \"strict_target_states\": "
      << statistics.strict_target_states << ",\n"
      << "  \"identity_checks\": "
      << statistics.identity_checks << ",\n"
      << "  \"determinant_checks\": "
      << statistics.determinant_checks << ",\n"
      << "  \"descriptor_invariance_checks\": "
      << statistics.descriptor_invariance_checks << ",\n"
      << coronal_archive_json(
             arguments, statistics, archive)
      << "  \"search\": {\n"
      << "    \"epoch_move_range\": ["
      << arguments.epoch_min_moves << ", "
      << arguments.epoch_max_moves << "],\n"
      << "    \"kick_flips\": " << arguments.kick_flips
      << ",\n"
      << "    \"kick_attempts\": " << arguments.kick_attempts
      << ",\n"
      << "    \"transient_stride\": "
      << arguments.transient_stride << ",\n"
      << "    \"polish_move_limit\": "
      << arguments.polish_moves << ",\n"
      << "    \"coronal_pareto_parent_selection\": "
      << (arguments.coronal_pareto ? "true" : "false")
      << ",\n"
      << "    \"gate_scope\": \""
      << gate_scope_name(arguments) << "\"\n"
      << "  },\n"
      << "  \"gram_sketch\": {\n"
      << "    \"canonical\": false,\n"
      << "    \"signed_line_orientation\": "
         "\"odd-order line product\",\n"
      << "    \"mod4_normalized_labels\": true,\n"
      << "    \"global_label_histogram\": true,\n"
      << "    \"sorted_incident_label_histograms\": true,\n"
      << "    \"row_column_pair_unordered\": true,\n"
      << "    \"lifetime_discovery_set_bounded\": false\n"
      << "  },\n"
      << "  \"claim_boundary\": [\n"
      << "    \"Every reported determinant objective is an exact integer.\",\n"
      << "    \"The Gram sketch is invariant but not canonical; collisions are possible.\",\n"
      << "    \"Distinct sketches imply distinct signed Gram orbits, but matching sketches do not prove equivalence.\",\n"
      << "    \"Coronal det(M) and kappa coordinates are exact diversity heuristics and never replace determinant ranking.\",\n"
      << "    \"The coronal descriptor is the sorted row/column pair; a Pareto point is not a proof of a distinct H/HT class.\",\n"
      << "    \"Any strict target must pass the independent arena verifier before publication or a claim of improvement.\",\n"
      << "    \"This stochastic search makes no optimality or literature-novelty claim.\"\n"
      << "  ]\n"
      << "}\n";
  return output.str();
}

const GramSketch& select_parent(
    std::map<GramSketch, SketchElite>& archive,
    const HopperArguments& arguments,
    ParentMode mode, std::mt19937_64& randomizer,
    HopperStatistics& statistics) {
  if (archive.empty()) {
    throw std::runtime_error("Gram-sketch archive is empty");
  }
  auto random_iterator = [&]() {
    auto iterator = archive.begin();
    const std::size_t offset = static_cast<std::size_t>(
        randomizer() % archive.size());
    std::advance(
        iterator, static_cast<std::ptrdiff_t>(offset));
    return iterator;
  };

  auto chosen = random_iterator();
  if (mode == ParentMode::kQuality) {
    for (int trial = 1; trial < arguments.quality_tournament;
         ++trial) {
      auto candidate = random_iterator();
      if (
          magnitude(candidate->second.determinant) >
          magnitude(chosen->second.determinant)) {
        chosen = candidate;
      }
    }
    ++statistics.parent_quality;
  } else if (mode == ParentMode::kCoronalDetM) {
    if (!arguments.coronal_pareto) {
      throw std::runtime_error(
          "coronal det(M) parent mode is not enabled");
    }
    for (int trial = 1; trial < arguments.quality_tournament;
         ++trial) {
      auto candidate = random_iterator();
      if (coronal_det_m_is_larger(
              candidate->second, chosen->second)) {
        chosen = candidate;
      }
    }
    ++statistics.parent_coronal_det_m;
  } else if (mode == ParentMode::kCoronalKappa) {
    if (!arguments.coronal_pareto) {
      throw std::runtime_error(
          "coronal kappa parent mode is not enabled");
    }
    for (int trial = 1; trial < arguments.quality_tournament;
         ++trial) {
      auto candidate = random_iterator();
      if (coronal_kappa_is_smaller(
              candidate->second, chosen->second)) {
        chosen = candidate;
      }
    }
    ++statistics.parent_coronal_kappa;
  } else if (mode == ParentMode::kCoronalPareto) {
    if (!arguments.coronal_pareto) {
      throw std::runtime_error(
          "coronal Pareto parent mode is not enabled");
    }
    std::vector<decltype(archive.begin())> pareto;
    pareto.reserve(archive.size());
    for (auto candidate = archive.begin();
         candidate != archive.end(); ++candidate) {
      bool dominated = false;
      for (auto other = archive.begin();
           other != archive.end(); ++other) {
        if (
            candidate != other &&
            coronal_dominates(
                other->second, candidate->second)) {
          dominated = true;
          break;
        }
      }
      if (!dominated) pareto.push_back(candidate);
    }
    if (pareto.empty()) {
      throw std::runtime_error(
          "coronal Pareto front is unexpectedly empty");
    }
    chosen = pareto[static_cast<std::size_t>(
        randomizer() % pareto.size())];
    ++statistics.parent_coronal_pareto;
  } else if (mode == ParentMode::kLeastUsed) {
    for (auto iterator = archive.begin();
         iterator != archive.end(); ++iterator) {
      if (
          iterator->second.selections <
              chosen->second.selections ||
          (iterator->second.selections ==
               chosen->second.selections &&
           magnitude(iterator->second.determinant) >
               magnitude(chosen->second.determinant))) {
        chosen = iterator;
      }
    }
    ++statistics.parent_least_used;
  } else {
    ++statistics.parent_uniform;
  }
  ++chosen->second.selections;
  return chosen->first;
}

template <typename StopPredicate>
std::optional<std::pair<CoreMatrix, std::int64_t>> best_exact_kick(
    const CoreMatrix& parent, const HopperArguments& arguments,
    std::mt19937_64& randomizer,
    HopperStatistics& statistics,
    StopPredicate&& stop_or_deadline) {
  std::array<int, kCoreEntries> coordinates{};
  std::iota(coordinates.begin(), coordinates.end(), 0);
  std::optional<std::pair<CoreMatrix, std::int64_t>> best;
  for (int attempt = 0; attempt < arguments.kick_attempts;
       ++attempt) {
    if (stop_or_deadline()) break;
    std::shuffle(
        coordinates.begin(), coordinates.end(), randomizer);
    CoreMatrix candidate = parent;
    for (int index = 0; index < arguments.kick_flips; ++index) {
      const int flat =
          coordinates[static_cast<std::size_t>(index)];
      candidate[flat / kCoreOrder][flat % kCoreOrder] ^= 1U;
    }
    const std::int64_t determinant =
        exact_core_determinant(candidate);
    ++statistics.kick_evaluations;
    if (determinant == 0) {
      ++statistics.kick_singular;
      continue;
    }
    if (
        !best.has_value() ||
        magnitude(determinant) >
            magnitude(best->second)) {
      best = std::make_pair(candidate, determinant);
    }
  }
  return best;
}

std::optional<Move> best_improving_move(
    const State& state, HopperStatistics& statistics) {
  std::optional<Move> best;
  std::uint64_t best_magnitude = magnitude(state.determinant);
  auto consider = [&](Move candidate) {
    ++statistics.kernel.candidate_evaluations;
    if (candidate.determinant == 0) {
      ++statistics.kernel.singular_candidates;
      return;
    }
    const std::uint64_t candidate_magnitude =
        magnitude(candidate.determinant);
    if (candidate_magnitude > best_magnitude) {
      best = candidate;
      best_magnitude = candidate_magnitude;
    }
  };
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      consider(Move{
          MoveKind::kBit, row, column,
          row * kCoreOrder + column,
          bit_candidate_determinant(state, row, column), false});
    }
  }
  for (int row = 0; row < kCoreOrder; ++row) {
    consider(Move{
        MoveKind::kRowComplement, row, -1,
        kCoreEntries + row,
        row_complement_determinant(state, row), false});
  }
  for (int column = 0; column < kCoreOrder; ++column) {
    consider(Move{
        MoveKind::kColumnComplement, column, -1,
        kCoreEntries + kCoreOrder + column,
        column_complement_determinant(state, column), false});
  }
  consider(Move{
      MoveKind::kWholeComplement, -1, -1, kMoveCount - 1,
      whole_complement_determinant(state), false});
  return best;
}

std::vector<fs::path> export_archive(
    const HopperArguments& arguments,
    const std::map<GramSketch, SketchElite>& archive,
    std::uint64_t& nonce, HopperStatistics& statistics) {
  using Entry =
      std::pair<const GramSketch*, const SketchElite*>;
  std::vector<Entry> ranked;
  ranked.reserve(archive.size());
  for (const auto& item : archive) {
    ranked.emplace_back(&item.first, &item.second);
  }
  std::sort(
      ranked.begin(), ranked.end(),
      [](const Entry& first, const Entry& second) {
        const std::uint64_t first_magnitude =
            magnitude(first.second->determinant);
        const std::uint64_t second_magnitude =
            magnitude(second.second->determinant);
        if (first_magnitude != second_magnitude) {
          return first_magnitude > second_magnitude;
        }
        return first.second->sketch_hash <
            second.second->sketch_hash;
      });

  const std::size_t count = std::min(
      ranked.size(),
      static_cast<std::size_t>(arguments.archive_outputs));
  std::vector<fs::path> outputs;
  outputs.reserve(count);
  std::ostringstream manifest;
  manifest
      << "{\n"
      << "  \"schema_version\": 1,\n"
      << "  \"descriptor\": \"invariant-noncanonical-gram-sketch-v1\",\n"
      << "  \"coronal_pareto\": "
      << (arguments.coronal_pareto ? "true" : "false")
      << ",\n"
      << "  \"elites\": [\n";
  for (std::size_t index = 0; index < count; ++index) {
    const GramSketch& sketch = *ranked[index].first;
    const SketchElite& elite = *ranked[index].second;
    const std::int64_t checked =
        exact_core_determinant(elite.core);
    if (checked != elite.determinant) {
      throw std::runtime_error(
          "archive elite determinant invariant failed");
    }
    ++statistics.determinant_checks;
    if (!(gram_sketch(elite.core) == sketch)) {
      throw std::runtime_error(
          "archive elite sketch invariant failed");
    }
    ++statistics.sketch_evaluations;
    if (arguments.coronal_pareto) {
      const CoronalPair checked_coronal =
          evaluate_coronal_pair(
              elite.core, elite.determinant, statistics);
      if (!(checked_coronal == require_coronal(elite))) {
        throw std::runtime_error(
            "archive elite coronal invariant failed");
      }
    }
    std::ostringstream filename;
    filename
        << "elite-" << std::setw(3) << std::setfill('0')
        << index << "-q" << magnitude(elite.determinant)
        << "-s" << hexadecimal(elite.sketch_hash) << ".txt";
    const fs::path path =
        arguments.archive_directory / filename.str();
    atomic_write(
        path,
        sign_matrix_text(core_to_sign(elite.core)), nonce++);
    outputs.push_back(path);
    if (index != 0) manifest << ",\n";
    manifest
        << "    {\"rank\": " << index
        << ", \"path\": \"" << json_escape(path.string())
        << "\", \"core_quotient\": "
        << magnitude(elite.determinant)
        << ", \"absolute_determinant\": \""
        << wide_to_string(
               static_cast<Wide>(
                   magnitude(elite.determinant)) *
               kScale)
        << "\", \"sketch_hash\": \""
        << hexadecimal(elite.sketch_hash)
        << "\", \"origin\": \""
        << json_escape(elite.origin) << "\"";
    if (arguments.coronal_pareto) {
      manifest
          << ", \"coronal_points\": "
          << coronal_pair_json(require_coronal(elite));
    }
    manifest << "}";
  }
  manifest << "\n  ]\n}\n";
  atomic_write(
      arguments.archive_directory / "manifest.json",
      manifest.str(), nonce++);
  return outputs;
}

int run_hopper(const HopperArguments& arguments) {
  HopperStatistics statistics;
  run_descriptor_self_tests(
      64, arguments.seed ^ UINT64_C(0x4752414d534b4554),
      statistics, arguments.coronal_pareto);

  for (const fs::path& path :
       {arguments.output, arguments.log, arguments.summary}) {
    if (!path.parent_path().empty()) {
      fs::create_directories(path.parent_path());
    }
  }
  if (!arguments.archive_directory.parent_path().empty()) {
    fs::create_directories(
        arguments.archive_directory.parent_path());
  }
  if (!fs::create_directory(arguments.archive_directory)) {
    throw std::runtime_error(
        "archive directory appeared after path validation");
  }
  std::ofstream log(
      arguments.log, std::ios::out | std::ios::trunc);
  if (!log) {
    throw std::runtime_error("cannot create basin-hopper log");
  }

  std::map<GramSketch, SketchElite> archive;
  std::unordered_set<GramSketch, SketchHash> seen_sketches;
  std::unordered_set<std::string> strict_cores;
  std::mt19937_64 randomizer(arguments.seed);
  const auto zobrist = make_zobrist(
      arguments.seed ^ UINT64_C(0x424153494e484f50));
  CoreMatrix best_core{};
  std::uint64_t best_magnitude = 0;
  std::uint64_t nonce = 0;

  const auto started = Clock::now();
  const auto deadline =
      started + std::chrono::duration<double>(arguments.seconds);
  auto next_heartbeat =
      started + std::chrono::duration<double>(
                    arguments.heartbeat_seconds);
  bool deadline_reached = false;
  std::vector<fs::path> no_archive_outputs;

  auto elapsed_seconds = [&]() {
    return std::chrono::duration<double>(
               Clock::now() - started)
        .count();
  };
  auto stop_or_deadline = [&]() {
    if (Clock::now() >= deadline) deadline_reached = true;
    return stop_requested || deadline_reached;
  };
  auto write_summary = [&](bool complete,
                           std::string_view reason) {
    atomic_write(
        arguments.summary,
        hopper_summary_json(
            arguments, statistics, archive, elapsed_seconds(),
            best_magnitude, complete, reason,
            no_archive_outputs),
        nonce++);
  };
  auto write_log = [&](std::string_view event) {
    log << hopper_log_json(
        event, elapsed_seconds(), arguments, statistics,
        archive, best_magnitude);
    log.flush();
    if (!log) {
      throw std::runtime_error(
          "cannot write basin-hopper log");
    }
  };
  auto write_seed_coronal = [&](
                                 std::string_view source,
                                 std::uint64_t quotient,
                                 const CoronalPair& pair) {
    log
        << "{\"event\":\"seed_coronal\""
        << ",\"engine\":\"" << hopper_engine_name(arguments)
        << "\",\"elapsed_seconds\":" << std::fixed
        << std::setprecision(6) << elapsed_seconds()
        << ",\"source\":\"" << json_escape(source) << "\""
        << ",\"core_quotient\":" << quotient
        << ",\"absolute_determinant\":\""
        << wide_to_string(
               static_cast<Wide>(quotient) * kScale)
        << "\",\"points\":" << coronal_pair_json(pair)
        << ",\"ranking_role\":\"diversity_only\"}\n";
    log.flush();
    if (!log) {
      throw std::runtime_error(
          "cannot write seed-coronal log");
    }
    std::cout
        << "coronal seed=" << source
        << " q=" << quotient;
    for (const CoronalPoint& point : pair.orientations) {
      std::cout
          << " [detM=" << big_to_string(point.determinant_m)
          << " kappa="
          << big_to_string(point.kappa_numerator) << '/'
          << big_to_string(point.kappa_denominator) << ']';
    }
    std::cout << '\n' << std::flush;
  };
  std::string active_seed_source;
  auto promote = [&](const CoreMatrix& core,
                     std::int64_t determinant,
                     ObservationKind kind) {
    const std::uint64_t candidate_magnitude =
        magnitude(determinant);
    if (candidate_magnitude <= best_magnitude) return;
    best_magnitude = candidate_magnitude;
    best_core = core;
    ++statistics.promotions;
    atomic_write(
        arguments.output,
        sign_matrix_text(core_to_sign(best_core)), nonce++);
    write_log("new_best");
    std::cout
        << "new best core quotient=" << best_magnitude
        << " absolute="
        << wide_to_string(
               static_cast<Wide>(best_magnitude) * kScale)
        << " origin=" << observation_name(kind) << '\n'
        << std::flush;
  };
  auto observe = [&](
                     const CoreMatrix& core,
                     std::int64_t determinant,
                     ObservationKind kind, bool force_seed) {
    ++statistics.observations;
    switch (kind) {
      case ObservationKind::kSeed:
        ++statistics.seed_observations;
        break;
      case ObservationKind::kKick:
        ++statistics.kick_observations;
        break;
      case ObservationKind::kTransient:
        ++statistics.transient_observations;
        break;
      case ObservationKind::kPeak:
        ++statistics.peak_observations;
        break;
      case ObservationKind::kPolished:
        ++statistics.polished_observations;
        break;
    }
    promote(core, determinant, kind);
    const std::uint64_t candidate_magnitude =
        magnitude(determinant);
    if (candidate_magnitude >= kStrictCoreTarget) {
      if (strict_cores.insert(flatten_core(core)).second) {
        ++statistics.strict_target_states;
        write_log("strict_target_unverified");
      }
    }
    if (
        archive_gate_applies(arguments, kind, force_seed) &&
        candidate_magnitude < arguments.quotient_gate &&
        candidate_magnitude < kStrictCoreTarget) {
      ++statistics.gate_rejections;
      return false;
    }

    std::optional<CoronalPair> eager_coronal;
    if (
        arguments.coronal_pareto &&
        kind == ObservationKind::kSeed) {
      eager_coronal = evaluate_coronal_pair(
          core, determinant, statistics);
      write_seed_coronal(
          active_seed_source, candidate_magnitude,
          *eager_coronal);
    }

    GramSketch sketch = gram_sketch(core);
    ++statistics.sketch_evaluations;
    const bool new_discovery =
        seen_sketches.insert(sketch).second;
    if (new_discovery) ++statistics.sketch_discoveries;
    const std::uint64_t sketch_hash =
        stable_sketch_hash(sketch);
    auto found = archive.find(sketch);
    if (found != archive.end()) {
      if (elite_is_better(
              determinant, core, found->second)) {
        const std::uint64_t selections =
            found->second.selections;
        std::optional<CoronalPair> coronal;
        if (arguments.coronal_pareto) {
          coronal = eager_coronal.has_value()
              ? std::move(eager_coronal)
              : std::optional<CoronalPair>(
                    evaluate_coronal_pair(
                        core, determinant, statistics));
        }
        found->second = SketchElite{
            core,
            determinant,
            selections,
            found->second.inserted_at_epoch,
            statistics.epochs_started,
            sketch_hash,
            observation_name(kind),
            std::move(coronal)};
        ++statistics.archive_replacements;
        return true;
      }
      ++statistics.archive_rejections;
      return false;
    }

    std::optional<CoronalPair> coronal;
    if (arguments.coronal_pareto) {
      coronal = eager_coronal.has_value()
          ? std::move(eager_coronal)
          : std::optional<CoronalPair>(
                evaluate_coronal_pair(
                    core, determinant, statistics));
    }
    SketchElite candidate{
        core,
        determinant,
        0,
        statistics.epochs_started,
        statistics.epochs_started,
        sketch_hash,
        observation_name(kind),
        std::move(coronal)};

    if (
        archive.size() >=
        static_cast<std::size_t>(
            arguments.archive_capacity)) {
      if (arguments.coronal_pareto) {
        auto inserted = archive.emplace(
            std::move(sketch), std::move(candidate));
        if (!inserted.second) {
          throw std::runtime_error(
              "new Gram sketch collided during coronal insertion");
        }
        auto victim = worst_coronal_elite(archive);
        if (victim == inserted.first) {
          archive.erase(victim);
          ++statistics.archive_capacity_rejections;
          ++statistics.coronal_capacity_rejections;
          return false;
        }
        archive.erase(victim);
        ++statistics.archive_evictions;
        ++statistics.coronal_pareto_evictions;
        ++statistics.archive_insertions;
        return true;
      }
      auto worst = worst_elite(archive);
      if (!elite_is_better(
              determinant, core, worst->second)) {
        ++statistics.archive_capacity_rejections;
        return false;
      }
      archive.erase(worst);
      ++statistics.archive_evictions;
    }
    archive.emplace(
        std::move(sketch), std::move(candidate));
    ++statistics.archive_insertions;
    return true;
  };

  for (const fs::path& seed_path : arguments.seeds) {
    State seed = read_exact_state(seed_path);
    ++statistics.identity_checks;
    ++statistics.determinant_checks;
    active_seed_source = seed_path.string();
    static_cast<void>(observe(
        seed.core, seed.determinant, ObservationKind::kSeed,
        true));
    active_seed_source.clear();
  }
  if (best_magnitude == 0 || archive.empty()) {
    throw std::runtime_error(
        "no nonsingular seed entered the archive");
  }
  statistics.initialization_sketch_discoveries =
      statistics.sketch_discoveries;
  statistics.initialization_promotions = statistics.promotions;
  statistics.initialization_archive_size =
      static_cast<std::uint64_t>(archive.size());
  statistics.initialization_complete = true;
  atomic_write(
      arguments.output,
      sign_matrix_text(core_to_sign(best_core)), nonce++);
  write_log("start");
  write_summary(false, "running");
  std::signal(SIGINT, request_stop);
  std::signal(SIGTERM, request_stop);

  while (!stop_requested) {
    if (Clock::now() >= deadline) {
      deadline_reached = true;
      break;
    }
    if (
        arguments.maximum_epochs != 0 &&
        statistics.epochs_started >=
            arguments.maximum_epochs) {
      break;
    }
    ++statistics.epochs_started;
    const std::uint64_t parent_roll =
        randomizer() % UINT64_C(100);
    ParentMode mode = ParentMode::kUniform;
    if (arguments.coronal_pareto) {
      mode =
          parent_roll < UINT64_C(30)
          ? ParentMode::kQuality
          : (parent_roll < UINT64_C(50)
                 ? ParentMode::kCoronalDetM
                 : (parent_roll < UINT64_C(70)
                        ? ParentMode::kCoronalKappa
                        : (parent_roll < UINT64_C(90)
                               ? ParentMode::kCoronalPareto
                               : (parent_roll < UINT64_C(95)
                                      ? ParentMode::kLeastUsed
                                      : ParentMode::kUniform))));
    } else {
      mode =
          parent_roll < UINT64_C(50)
          ? ParentMode::kQuality
          : (parent_roll < UINT64_C(85)
                 ? ParentMode::kLeastUsed
                 : ParentMode::kUniform);
    }
    const GramSketch parent_key = select_parent(
        archive, arguments, mode, randomizer, statistics);
    const CoreMatrix parent_core =
        archive.at(parent_key).core;

    const auto kicked = best_exact_kick(
        parent_core, arguments, randomizer, statistics,
        stop_or_deadline);
    if (!kicked.has_value()) {
      if (stop_or_deadline()) break;
      write_log("all_kicks_singular");
      continue;
    }
    ++statistics.kick_endpoints;
    static_cast<void>(observe(
        kicked->first, kicked->second, ObservationKind::kKick,
        false));
    if (stop_or_deadline()) break;

    State state = exact_state_from_core(kicked->first);
    if (state.determinant != kicked->second) {
      throw std::runtime_error(
          "kick determinant changed during state rebuild");
    }
    ++statistics.identity_checks;
    ++statistics.determinant_checks;
    std::uint64_t state_hash = core_hash(state.core, zobrist);
    std::array<std::uint64_t, kMoveCount> tabu_until{};
    std::vector<Visit> visits(kVisitTableSize);
    visits[state_hash & (kVisitTableSize - 1U)] =
        Visit{state_hash, 0, true};
    const int baseline_tenure =
        kMinimumTenure +
        static_cast<int>(randomizer() % UINT64_C(5));
    int tenure = baseline_tenure;
    std::uint64_t last_cycle_iteration = 0;

    State epoch_peak = state;
    std::uint64_t epoch_peak_magnitude =
        magnitude(state.determinant);
    std::uniform_int_distribution<int> epoch_length_distribution(
        arguments.epoch_min_moves,
        arguments.epoch_max_moves);
    const int epoch_move_limit =
        epoch_length_distribution(randomizer);
    int epoch_moves = 0;
    for (int local_move = 1; local_move <= epoch_move_limit;
         ++local_move) {
      if (
          (local_move & 255) == 0 &&
          stop_or_deadline()) {
        break;
      }
      const std::uint64_t iteration =
          static_cast<std::uint64_t>(local_move);
      const Move move = choose_move(
          state, tabu_until, iteration, epoch_peak_magnitude,
          randomizer, statistics.kernel);
      apply_move(
          state, move, state_hash, zobrist, statistics.kernel);
      ++statistics.tabu_moves;
      ++epoch_moves;

      const int jitter =
          static_cast<int>(randomizer() % UINT64_C(5));
      tabu_until[move.id] =
          iteration +
          static_cast<std::uint64_t>(tenure + jitter + 1);
      Visit& visit =
          visits[state_hash & (kVisitTableSize - 1U)];
      if (
          visit.occupied && visit.hash == state_hash &&
          iteration > visit.iteration) {
        const std::uint64_t cycle_length =
            iteration - visit.iteration;
        if (
            cycle_length <=
            static_cast<std::uint64_t>(
                4 * kMaximumTenure)) {
          ++statistics.kernel.cycles;
          last_cycle_iteration = iteration;
          tenure = std::min(
              kMaximumTenure,
              tenure + 2 +
                  static_cast<int>(
                      std::min<std::uint64_t>(
                          cycle_length / UINT64_C(8),
                          UINT64_C(8))));
        }
      }
      visit = Visit{state_hash, iteration, true};
      if (
          iteration - last_cycle_iteration >= UINT64_C(512) &&
          (iteration & UINT64_C(127)) == 0 &&
          tenure > baseline_tenure) {
        --tenure;
      }

      const std::uint64_t current_magnitude =
          magnitude(state.determinant);
      if (current_magnitude > epoch_peak_magnitude) {
        epoch_peak = state;
        epoch_peak_magnitude = current_magnitude;
      }
      const bool scheduled_transient =
          local_move % arguments.transient_stride == 0 &&
          current_magnitude >= arguments.quotient_gate;
      if (
          scheduled_transient ||
          current_magnitude >= kStrictCoreTarget) {
        static_cast<void>(observe(
            state.core, state.determinant,
            ObservationKind::kTransient, false));
      }
      if (
          iteration % kEpochIdentityCheckInterval == 0) {
        check_adjugate_identity(state);
        ++statistics.identity_checks;
      }
      if (
          iteration % kEpochDeterminantCheckInterval == 0) {
        if (
            exact_core_determinant(state.core) !=
            state.determinant) {
          throw std::runtime_error(
              "tabu incremental determinant invariant failed");
        }
        ++statistics.determinant_checks;
      }

      const auto now = Clock::now();
      if (
          arguments.heartbeat_seconds > 0.0 &&
          now >= next_heartbeat) {
        write_log("heartbeat");
        write_summary(false, "running");
        next_heartbeat =
            now + std::chrono::duration<double>(
                      arguments.heartbeat_seconds);
      }
    }
    check_adjugate_identity(state);
    ++statistics.identity_checks;
    if (
        exact_core_determinant(state.core) !=
        state.determinant) {
      throw std::runtime_error(
          "epoch-final determinant invariant failed");
    }
    ++statistics.determinant_checks;
    static_cast<void>(observe(
        epoch_peak.core, epoch_peak.determinant,
        ObservationKind::kPeak, false));

    State polished = epoch_peak;
    std::uint64_t polished_hash =
        core_hash(polished.core, zobrist);
    for (int move_index = 0;
         move_index < arguments.polish_moves; ++move_index) {
      if (stop_or_deadline()) break;
      const std::optional<Move> move =
          best_improving_move(polished, statistics);
      if (!move.has_value()) break;
      apply_move(
          polished, *move, polished_hash, zobrist,
          statistics.kernel);
      ++statistics.polish_moves;
    }
    check_adjugate_identity(polished);
    ++statistics.identity_checks;
    if (
        exact_core_determinant(polished.core) !=
        polished.determinant) {
      throw std::runtime_error(
          "polished determinant invariant failed");
    }
    ++statistics.determinant_checks;
    if (
        magnitude(polished.determinant) >=
            arguments.quotient_gate ||
        magnitude(polished.determinant) >=
            kStrictCoreTarget) {
      static_cast<void>(observe(
          polished.core, polished.determinant,
          ObservationKind::kPolished, false));
    }
    if (epoch_moves == epoch_move_limit) {
      ++statistics.epochs_completed;
    }
    write_log("epoch");
    if (stop_requested || deadline_reached) break;
  }

  const std::int64_t checked_best_determinant =
      exact_core_determinant(best_core);
  ++statistics.determinant_checks;
  if (magnitude(checked_best_determinant) != best_magnitude) {
    throw std::runtime_error(
        "final best-core determinant magnitude invariant failed");
  }
  const std::vector<fs::path> archive_outputs =
      export_archive(
          arguments, archive, nonce, statistics);
  no_archive_outputs = archive_outputs;
  atomic_write(
      arguments.output,
      sign_matrix_text(core_to_sign(best_core)), nonce++);
  const std::string reason = stop_requested
      ? "signal"
      : (deadline_reached ? "time_limit" : "epoch_limit");
  atomic_write(
      arguments.summary,
      hopper_summary_json(
          arguments, statistics, archive, elapsed_seconds(),
          best_magnitude, !stop_requested, reason,
          archive_outputs),
      nonce++);
  write_log(stop_requested ? "stopped" : "finished");
  std::cout
      << "finished best_core_quotient=" << best_magnitude
      << " strict_target=" << kStrictCoreTarget
      << " archive_size=" << archive.size()
      << " sketch_discoveries="
      << statistics.sketch_discoveries
      << " tabu_moves=" << statistics.tabu_moves
      << " candidate_evaluations="
      << statistics.kernel.candidate_evaluations << '\n';
  if (best_magnitude >= kStrictCoreTarget) {
    std::cout
        << "WARNING: strict target reached but not independently "
           "arena-verified; verify "
        << arguments.output.string() << '\n';
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const HopperArguments arguments =
        parse_hopper_arguments(argc, argv);
    if (arguments.self_test_rounds != 0) {
      HopperStatistics statistics;
      run_descriptor_self_tests(
          arguments.self_test_rounds, arguments.seed,
          statistics, true);
      std::cout
          << "Gram-sketch self-test passed rounds="
          << arguments.self_test_rounds
          << " invariance_checks="
          << statistics.descriptor_invariance_checks
          << " coronal_invariance_checks="
          << statistics.coronal_invariance_checks << '\n';
      return 0;
    }
    return run_hopper(arguments);
  } catch (const std::exception& error) {
    std::cerr
        << "core_gram_basin_hopper: " << error.what() << '\n';
    return 2;
  }
}
