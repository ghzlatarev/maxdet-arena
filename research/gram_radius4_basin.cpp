#include <atomic>
#include <gmpxx.h>
#include <mutex>
#include <thread>

// Keep the radius-3 catalog parser, exact CRT implementation, reference-graph
// validation, and automorphism construction in one audited implementation.
// Renaming its CLI entry point lets this translation unit specialize the next
// search without copying that mathematical trusted core.
#define main gram_radius3_embedded_main
#include "gram_radius3_orbits.cpp"
#undef main

namespace radius4_basin {

using Big = mpz_class;

constexpr std::string_view kFrozenCatalogSha256 =
    "b05098cd1a3bec294f67422b9c796885e9aef6d807b450f5ee1a3a38d607d6ed";
constexpr std::uint64_t kFrozenRadius3Representatives = 9'967'496;
constexpr std::uint64_t kExpectedAboveParents = 1'958;
constexpr std::uint64_t kExpectedAboveParentLabeled = 148'736;
constexpr std::uint64_t kOutwardRemovalsPerParent = 42;
constexpr std::uint64_t kOutwardAdditionsPerParent = 205;
constexpr std::uint64_t kTransitionsPerParent =
    kOutwardRemovalsPerParent * kOutwardAdditionsPerParent;
constexpr std::uint64_t kExpectedTransitions =
    kExpectedAboveParents * kTransitionsPerParent;
constexpr std::uint64_t kIndependentPrime = 998'244'353;

enum class UpdateMethod {
  kParentRank4,
  kBaseRank16,
};

struct BasinArguments {
  bool self_test = false;
  bool screen = false;
  std::filesystem::path reference =
      "references/orrick-et-al-2003/matrix.txt";
  std::filesystem::path catalog;
  std::filesystem::path output;
  std::filesystem::path route_snapshot;
  std::filesystem::path parents_snapshot;
  std::string catalog_sha256;
  UpdateMethod method = UpdateMethod::kParentRank4;
  std::uint64_t threads = 1;
  double heartbeat_seconds = 10.0;
};

struct Radius4Key {
  std::array<std::uint8_t, 4> removed{};
  std::array<std::uint8_t, 4> added{};

  auto operator<=>(const Radius4Key&) const = default;
};

struct AboveParent {
  std::uint64_t catalog_index = 0;
  std::uint64_t orbit_size = 0;
  Wide determinant = 0;
  OrbitRecord record;
};

struct RelevantHit {
  Radius4Key key;
  Radius4Key canonical_key;
  Wide determinant = 0;
  Wide root = 0;
  bool root_divisible = false;
  bool positive_definite = false;
};

struct WorkerSummary {
  std::uint64_t parents = 0;
  std::uint64_t transitions = 0;
  std::uint64_t positive_determinant_transitions = 0;
  std::uint64_t above_frontier_transitions = 0;
  std::uint64_t square_transitions = 0;
  std::uint64_t above_frontier_square_transitions = 0;
  std::uint64_t divisible_above_frontier_square_transitions = 0;
  std::uint64_t positive_definite_route_transitions = 0;
  std::uint64_t singular_parent_prime_fallbacks = 0;
  std::uint64_t direct_spot_checks = 0;
  bool has_maximum = false;
  Wide maximum_determinant = 0;
  Radius4Key maximum_key;
  std::map<Radius4Key, RelevantHit> relevant_hits;
};

struct BasinSummary : WorkerSummary {
  std::uint64_t radius3_catalog_records = 0;
  std::uint64_t radius3_catalog_orbit_size_sum = 0;
  std::uint64_t above_parent_labeled = 0;
  std::uint64_t unique_relevant_states = 0;
  std::uint64_t unique_relevant_orbits = 0;
  std::uint64_t unique_route_states = 0;
  std::uint64_t unique_route_orbits = 0;
  bool maximum_positive_definite = false;
  std::vector<RelevantHit> canonical_hits;
  double recovery_seconds = 0.0;
  double expansion_seconds = 0.0;
};

std::string method_name(UpdateMethod method) {
  return method == UpdateMethod::kParentRank4 ? "parent-rank4"
                                              : "base-rank16";
}

std::string sha256_file(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("cannot open " + path.string());
  }
  Sha256 digest;
  std::array<char, 1 << 20> buffer{};
  while (input) {
    input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    const std::streamsize count = input.gcount();
    if (count > 0) {
      digest.update(std::string_view(
          buffer.data(), static_cast<std::size_t>(count)));
    }
  }
  if (!input.eof()) {
    throw std::runtime_error("cannot read " + path.string());
  }
  return digest.finish();
}

std::string next_value(int& index, int argc, char** argv) {
  if (index + 1 >= argc) {
    throw std::runtime_error(
        "missing value after " + std::string(argv[index]));
  }
  return argv[++index];
}

BasinArguments parse_basin_arguments(int argc, char** argv) {
  BasinArguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (option == "--self-test") {
      arguments.self_test = true;
    } else if (option == "--screen") {
      arguments.screen = true;
    } else if (option == "--reference") {
      arguments.reference = next_value(index, argc, argv);
    } else if (option == "--catalog") {
      arguments.catalog = next_value(index, argc, argv);
    } else if (option == "--catalog-sha256") {
      arguments.catalog_sha256 = next_value(index, argc, argv);
    } else if (option == "--output") {
      arguments.output = next_value(index, argc, argv);
    } else if (option == "--route-snapshot") {
      arguments.route_snapshot = next_value(index, argc, argv);
    } else if (option == "--parents-snapshot") {
      arguments.parents_snapshot = next_value(index, argc, argv);
    } else if (option == "--method") {
      const std::string method = next_value(index, argc, argv);
      if (method == "parent-rank4") {
        arguments.method = UpdateMethod::kParentRank4;
      } else if (method == "base-rank16") {
        arguments.method = UpdateMethod::kBaseRank16;
      } else {
        throw std::runtime_error(
            "--method must be parent-rank4 or base-rank16");
      }
    } else if (option == "--threads") {
      arguments.threads =
          parse_unsigned(next_value(index, argc, argv), option);
    } else if (option == "--heartbeat-seconds") {
      arguments.heartbeat_seconds = parse_nonnegative_double(
          next_value(index, argc, argv), option);
    } else if (option == "--help") {
      std::cout
          << "Usage:\n"
          << "  gram_radius4_basin --self-test [--reference FILE]\n"
          << "  gram_radius4_basin --screen --catalog FILE "
             "--catalog-sha256 HEX --output FILE "
             "--parents-snapshot FILE --route-snapshot FILE "
             "[--method parent-rank4|base-rank16] [--threads N]\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (arguments.self_test == arguments.screen) {
    throw std::runtime_error("select exactly one of --self-test or --screen");
  }
  if (arguments.threads == 0 || arguments.threads > 64) {
    throw std::runtime_error("--threads must be in [1,64]");
  }
  if (arguments.screen &&
      (arguments.catalog.empty() || arguments.output.empty() ||
       arguments.route_snapshot.empty() ||
       arguments.parents_snapshot.empty() ||
       arguments.catalog_sha256.empty())) {
    throw std::runtime_error(
        "--screen requires catalog, catalog SHA, output, parent snapshot, "
        "and route snapshot");
  }
  if (arguments.screen &&
      arguments.catalog_sha256 != kFrozenCatalogSha256) {
    throw std::runtime_error(
        "--catalog-sha256 is not the frozen radius-3 catalog hash");
  }
  const std::array<std::filesystem::path, 3> outputs{
      arguments.output, arguments.route_snapshot,
      arguments.parents_snapshot};
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

std::uint64_t wide_mod(Wide value, std::uint64_t modulus) {
  Wide residue = value % static_cast<Wide>(modulus);
  if (residue < 0) residue += modulus;
  return static_cast<std::uint64_t>(residue);
}

Big bareiss_determinant(const Gram& matrix, int order = kOrder) {
  std::array<std::array<Big, kOrder>, kOrder> work{};
  for (int row = 0; row < order; ++row) {
    for (int column = 0; column < order; ++column) {
      work[row][column] = matrix[row][column];
    }
  }
  Big sign = 1;
  Big previous = 1;
  for (int pivot = 0; pivot < order - 1; ++pivot) {
    int pivot_row = pivot;
    while (pivot_row < order && work[pivot_row][pivot] == 0) {
      ++pivot_row;
    }
    if (pivot_row == order) return 0;
    if (pivot_row != pivot) {
      std::swap(work[pivot_row], work[pivot]);
      sign = -sign;
    }
    const Big pivot_value = work[pivot][pivot];
    for (int row = pivot + 1; row < order; ++row) {
      for (int column = pivot + 1; column < order; ++column) {
        const Big numerator =
            work[row][column] * pivot_value -
            work[row][pivot] * work[pivot][column];
        if (pivot > 0 && numerator % previous != 0) {
          throw std::runtime_error("Bareiss exact division failed");
        }
        work[row][column] =
            pivot == 0 ? numerator : numerator / previous;
      }
      work[row][pivot] = 0;
    }
    previous = pivot_value;
  }
  return sign * work[order - 1][order - 1];
}

void require_bareiss_match(const Gram& gram, Wide determinant) {
  const Big independent = bareiss_determinant(gram);
  if (independent.get_str() != decimal(determinant)) {
    throw std::runtime_error(
        "independent multiprecision Bareiss determinant disagrees");
  }
}

Radius4Key add_outward_edges(
    const OrbitRecord& parent, std::uint8_t removed,
    std::uint8_t added) {
  Radius4Key key;
  for (int index = 0; index < 3; ++index) {
    key.removed[index] = parent.removed_global.values[index];
    key.added[index] = parent.added_global.values[index];
  }
  key.removed[3] = removed;
  key.added[3] = added;
  std::sort(key.removed.begin(), key.removed.end());
  std::sort(key.added.begin(), key.added.end());
  if (std::adjacent_find(key.removed.begin(), key.removed.end()) !=
          key.removed.end() ||
      std::adjacent_find(key.added.begin(), key.added.end()) !=
          key.added.end()) {
    throw std::runtime_error("outward radius-4 update contains duplicates");
  }
  return key;
}

Gram radius4_gram(const Gram& base, const Radius4Key& key) {
  Gram result = base;
  const auto& edges = global_edges();
  for (const std::uint8_t edge_number : key.removed) {
    const Edge edge = edges[edge_number];
    if (result[edge.first][edge.second] != 3) {
      throw std::runtime_error("radius-4 removal is not a base edge");
    }
    result[edge.first][edge.second] = -1;
    result[edge.second][edge.first] = -1;
  }
  for (const std::uint8_t edge_number : key.added) {
    const Edge edge = edges[edge_number];
    if (result[edge.first][edge.second] != -1) {
      throw std::runtime_error("radius-4 addition is not a base nonedge");
    }
    result[edge.first][edge.second] = 3;
    result[edge.second][edge.first] = 3;
  }
  return result;
}

template <std::size_t ChangeCount>
std::array<std::uint64_t, kCrtPrimes.size()> update_residues(
    const std::array<std::uint8_t, ChangeCount>& edge_numbers,
    const std::array<int, ChangeCount>& coefficients,
    const ModularBases& bases) {
  std::array<bool, kOrder> used{};
  const auto& edges = global_edges();
  for (const std::uint8_t edge_number : edge_numbers) {
    const Edge edge = edges[edge_number];
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
  if (dimension <= 0 ||
      dimension > static_cast<int>(2 * ChangeCount)) {
    throw std::runtime_error("invalid rank-update endpoint dimension");
  }
  std::array<std::array<int, kOrder>, kOrder> delta{};
  for (std::size_t index = 0; index < ChangeCount; ++index) {
    const Edge edge = edges[edge_numbers[index]];
    const int first = position[edge.first];
    const int second = position[edge.second];
    delta[first][second] += coefficients[index];
    delta[second][first] += coefficients[index];
  }
  std::array<std::uint64_t, kCrtPrimes.size()> residues{};
  for (std::size_t prime_index = 0;
       prime_index < bases.size(); ++prime_index) {
    const ModularBase& modular = bases[prime_index];
    std::array<std::array<std::uint64_t, kOrder>, kOrder> lemma{};
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
                  : modular.prime -
                        static_cast<std::uint64_t>(-coefficient);
          value += coefficient_mod *
                   modular.inverse[vertices[inner]][vertices[column]] %
                   modular.prime;
          if (value >= modular.prime) value -= modular.prime;
        }
        lemma[row][column] = value;
      }
    }
    const std::uint64_t correction =
        determinant_modulo_residue_matrix(
            lemma, dimension, modular.prime);
    residues[prime_index] =
        modular.determinant * correction % modular.prime;
  }
  return residues;
}

std::array<std::uint64_t, kCrtPrimes.size()> base_update_residues(
    const Radius4Key& key, const ModularBases& bases) {
  std::array<std::uint8_t, 8> changes{};
  std::array<int, 8> coefficients{};
  for (int index = 0; index < 4; ++index) {
    changes[index] = key.removed[index];
    coefficients[index] = -4;
    changes[4 + index] = key.added[index];
    coefficients[4 + index] = 4;
  }
  return update_residues(changes, coefficients, bases);
}

std::array<std::uint64_t, kCrtPrimes.size()> parent_update_residues(
    std::uint8_t removed, std::uint8_t added,
    const ModularBases& parent_bases) {
  const std::array<std::uint8_t, 2> changes{removed, added};
  const std::array<int, 2> coefficients{-4, 4};
  return update_residues(changes, coefficients, parent_bases);
}

struct RecoveryResult {
  std::vector<AboveParent> parents;
  std::uint64_t record_count = 0;
  std::uint64_t orbit_size_sum = 0;
  std::uint64_t above_labeled = 0;
  std::uint64_t direct_checks = 0;
  double elapsed_seconds = 0.0;
};

std::string catalog_header(
    std::ifstream& input, std::string_view expected_key) {
  std::string line;
  if (!std::getline(input, line)) {
    throw std::runtime_error("catalog is truncated");
  }
  const std::string prefix = std::string(expected_key) + " ";
  if (!line.starts_with(prefix)) {
    throw std::runtime_error(
        "catalog expected header " + std::string(expected_key));
  }
  return line.substr(prefix.size());
}

RecoveryResult recover_above_parents(
    const std::filesystem::path& catalog_path,
    std::string_view reference_sha256,
    std::string_view reference_gram_sha256, const Gram& base,
    const ModularBases& base_bases, double heartbeat_seconds) {
  const Clock::time_point started = Clock::now();
  Clock::time_point next_heartbeat =
      started + std::chrono::duration_cast<Clock::duration>(
                    std::chrono::duration<double>(heartbeat_seconds));
  std::ifstream input(catalog_path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("cannot open " + catalog_path.string());
  }
  std::string line;
  if (!std::getline(input, line) ||
      line != "GRAM_RADIUS3_ORBIT_CATALOG_V1") {
    throw std::runtime_error("unsupported catalog schema");
  }
  if (catalog_header(input, "reference_raw_sha256") != reference_sha256 ||
      catalog_header(input, "reference_gram_sha256") !=
          reference_gram_sha256) {
    throw std::runtime_error("catalog reference hash mismatch");
  }
  if (parse_unsigned(
          catalog_header(input, "group_order"),
          "catalog group_order") != kAutomorphismGroupOrder) {
    throw std::runtime_error("catalog group order mismatch");
  }
  if (parse_unsigned(catalog_header(input, "radius"), "catalog radius") !=
      3) {
    throw std::runtime_error("catalog radius is not three");
  }
  if (parse_unsigned(
          catalog_header(input, "labeled_count"),
          "catalog labeled_count") != kRadius3Labeled) {
    throw std::runtime_error("catalog labeled count mismatch");
  }
  if (parse_unsigned(
          catalog_header(input, "representative_count"),
          "catalog representative_count") !=
      kFrozenRadius3Representatives) {
    throw std::runtime_error("catalog representative count mismatch");
  }
  if (catalog_header(input, "edge_numbering") !=
      "lexicographic-zero-based-pairs-on-23-vertices") {
    throw std::runtime_error("catalog edge numbering mismatch");
  }
  if (!std::getline(input, line) ||
      line != "records index orbit_size removed_edge_ids added_edge_ids") {
    throw std::runtime_error("catalog records header mismatch");
  }

  const std::vector<Edge> reference_edges = expected_reference_edges();
  const std::set<Edge> present(
      reference_edges.begin(), reference_edges.end());
  const Wide frontier_squared =
      static_cast<Wide>(kFrontierRoot) * kFrontierRoot;
  RecoveryResult result;
  std::pair<Subset, Subset> previous_key;
  bool have_previous = false;

  while (std::getline(input, line)) {
    if (line.empty()) {
      throw std::runtime_error("catalog contains a blank record");
    }
    std::istringstream fields(line);
    std::string index_text;
    std::string orbit_size_text;
    std::string removed_text;
    std::string added_text;
    std::string trailing;
    if (!(fields >> index_text >> orbit_size_text >>
          removed_text >> added_text) ||
        (fields >> trailing)) {
      throw std::runtime_error("malformed catalog record");
    }
    OrbitRecord record;
    record.index = parse_unsigned(index_text, "catalog record index");
    record.orbit_size =
        parse_unsigned(orbit_size_text, "catalog record orbit size");
    record.removed_global = parse_subset_csv(removed_text, 3);
    record.added_global = parse_subset_csv(added_text, 3);
    if (record.index != result.record_count ||
        record.orbit_size == 0 ||
        kAutomorphismGroupOrder % record.orbit_size != 0) {
      throw std::runtime_error("invalid catalog record index/orbit size");
    }
    const std::pair<Subset, Subset> key{
        record.removed_global, record.added_global};
    if (have_previous && !(previous_key < key)) {
      throw std::runtime_error(
          "catalog records are not in strict canonical order");
    }
    previous_key = key;
    have_previous = true;
    for (int index = 0; index < 3; ++index) {
      if (!present.contains(
              global_edges()[record.removed_global.values[index]]) ||
          present.contains(
              global_edges()[record.added_global.values[index]])) {
        throw std::runtime_error(
            "catalog record violates reference edge colors");
      }
    }

    const auto residues =
        fast_determinant_residues(record, base_bases);
    const Wide determinant = reconstruct_crt(residues);
    const bool periodic_check =
        (record.index & ((std::uint64_t{1} << 17) - 1)) == 0;
    if (periodic_check || determinant > frontier_squared) {
      const Gram candidate = candidate_gram(base, record);
      for (std::size_t prime_index = 0;
           prime_index < kCrtPrimes.size(); ++prime_index) {
        if (residues[prime_index] != determinant_modulo(
                candidate, kOrder, kCrtPrimes[prime_index])) {
          throw std::runtime_error(
              "parent fast/full determinant residue mismatch");
        }
      }
      if (exact_determinant(candidate) != determinant ||
          determinant_modulo(candidate, kOrder, kIndependentPrime) !=
              wide_mod(determinant, kIndependentPrime)) {
        throw std::runtime_error(
            "parent exact/independent-prime determinant mismatch");
      }
      ++result.direct_checks;
    }
    if (determinant > frontier_squared) {
      AboveParent parent;
      parent.catalog_index = record.index;
      parent.orbit_size = record.orbit_size;
      parent.determinant = determinant;
      parent.record = record;
      result.parents.push_back(parent);
      result.above_labeled += record.orbit_size;
    }
    ++result.record_count;
    result.orbit_size_sum += record.orbit_size;

    const Clock::time_point now = Clock::now();
    if (heartbeat_seconds > 0.0 && now >= next_heartbeat) {
      const double elapsed =
          std::chrono::duration<double>(now - started).count();
      std::cerr
          << "{\"above_parents\":" << result.parents.size()
          << ",\"elapsed_seconds\":" << std::fixed
          << std::setprecision(3) << elapsed
          << ",\"event\":\"radius4-basin-parent-recovery-heartbeat\""
          << ",\"records\":" << result.record_count << "}\n";
      next_heartbeat =
          now + std::chrono::duration_cast<Clock::duration>(
                    std::chrono::duration<double>(heartbeat_seconds));
    }
  }
  if (!input.eof()) {
    throw std::runtime_error("catalog read failed");
  }
  if (result.record_count != kFrozenRadius3Representatives ||
      result.orbit_size_sum != kRadius3Labeled ||
      result.parents.size() != kExpectedAboveParents ||
      result.above_labeled != kExpectedAboveParentLabeled) {
    throw std::runtime_error(
        "recovered parent family disagrees with frozen radius-3 screen");
  }
  result.elapsed_seconds =
      std::chrono::duration<double>(Clock::now() - started).count();
  return result;
}

std::string parents_snapshot_bytes(
    const RecoveryResult& recovery, std::string_view catalog_sha256) {
  std::ostringstream output;
  output << "GRAM_RADIUS4_BASIN_PARENTS_V1\n";
  output << "catalog_sha256 " << catalog_sha256 << '\n';
  output << "frontier_squared "
         << decimal(static_cast<Wide>(kFrontierRoot) * kFrontierRoot)
         << '\n';
  output << "selection determinant_strictly_greater_than_frontier_squared\n";
  output << "representative_count " << recovery.parents.size() << '\n';
  output << "labeled_count " << recovery.above_labeled << '\n';
  output << "records catalog_index orbit_size determinant "
            "removed_edge_ids added_edge_ids\n";
  for (const AboveParent& parent : recovery.parents) {
    output << parent.catalog_index << ' ' << parent.orbit_size << ' '
           << decimal(parent.determinant) << ' '
           << subset_csv(parent.record.removed_global) << ' '
           << subset_csv(parent.record.added_global) << '\n';
  }
  return output.str();
}

bool contains_edge(
    const Subset& subset, std::uint8_t edge_number) {
  return std::find(
             subset.values.begin(),
             subset.values.begin() + subset.size,
             edge_number) != subset.values.begin() + subset.size;
}

std::vector<std::uint8_t> global_edge_ids(
    const std::vector<Edge>& edges) {
  const auto numbers = global_edge_numbers();
  std::vector<std::uint8_t> result;
  result.reserve(edges.size());
  for (const Edge edge : edges) {
    result.push_back(numbers[edge.first][edge.second]);
  }
  return result;
}

std::vector<Edge> absent_reference_edges(
    const std::vector<Edge>& present) {
  const std::set<Edge> present_set(present.begin(), present.end());
  std::vector<Edge> absent;
  for (const Edge edge : global_edges()) {
    if (!present_set.contains(edge)) absent.push_back(edge);
  }
  if (present.size() != 45 || absent.size() != 208) {
    throw std::runtime_error("reference edge partition mismatch");
  }
  return absent;
}

Radius4Key canonical_key(
    const Radius4Key& key,
    const std::vector<Permutation>& group) {
  const auto edge_numbers = global_edge_numbers();
  const auto& edges = global_edges();
  Radius4Key best = key;
  for (const Permutation& permutation : group) {
    Radius4Key image;
    for (int index = 0; index < 4; ++index) {
      const Edge removed_image =
          permute_edge(edges[key.removed[index]], permutation);
      const Edge added_image =
          permute_edge(edges[key.added[index]], permutation);
      image.removed[index] =
          edge_numbers[removed_image.first][removed_image.second];
      image.added[index] =
          edge_numbers[added_image.first][added_image.second];
    }
    std::sort(image.removed.begin(), image.removed.end());
    std::sort(image.added.begin(), image.added.end());
    if (image < best) best = image;
  }
  return best;
}

std::vector<Edge> defect_edges(
    const std::vector<Edge>& reference_present,
    const Radius4Key& key) {
  std::set<Edge> defects(
      reference_present.begin(), reference_present.end());
  const auto& edges = global_edges();
  for (const std::uint8_t edge_number : key.removed) {
    if (defects.erase(edges[edge_number]) != 1) {
      throw std::runtime_error("cannot remove radius-4 defect edge");
    }
  }
  for (const std::uint8_t edge_number : key.added) {
    if (!defects.insert(edges[edge_number]).second) {
      throw std::runtime_error("cannot add radius-4 defect edge");
    }
  }
  if (defects.size() != 45) {
    throw std::runtime_error("radius-4 defect edge count is not 45");
  }
  return std::vector<Edge>(defects.begin(), defects.end());
}

void merge_worker_summary(
    BasinSummary& target, WorkerSummary&& source) {
  target.parents += source.parents;
  target.transitions += source.transitions;
  target.positive_determinant_transitions +=
      source.positive_determinant_transitions;
  target.above_frontier_transitions +=
      source.above_frontier_transitions;
  target.square_transitions += source.square_transitions;
  target.above_frontier_square_transitions +=
      source.above_frontier_square_transitions;
  target.divisible_above_frontier_square_transitions +=
      source.divisible_above_frontier_square_transitions;
  target.positive_definite_route_transitions +=
      source.positive_definite_route_transitions;
  target.singular_parent_prime_fallbacks +=
      source.singular_parent_prime_fallbacks;
  target.direct_spot_checks += source.direct_spot_checks;
  if (source.has_maximum &&
      (!target.has_maximum ||
       source.maximum_determinant > target.maximum_determinant ||
       (source.maximum_determinant == target.maximum_determinant &&
        source.maximum_key < target.maximum_key))) {
    target.has_maximum = true;
    target.maximum_determinant = source.maximum_determinant;
    target.maximum_key = source.maximum_key;
  }
  for (auto& [key, hit] : source.relevant_hits) {
    const auto [iterator, inserted] =
        target.relevant_hits.emplace(key, std::move(hit));
    if (!inserted &&
        (iterator->second.determinant != hit.determinant ||
         iterator->second.root != hit.root ||
         iterator->second.root_divisible != hit.root_divisible ||
         iterator->second.positive_definite != hit.positive_definite)) {
      throw std::runtime_error(
          "duplicate relevant state has inconsistent invariants");
    }
  }
}

WorkerSummary screen_parent_partition(
    std::size_t thread_index, std::size_t thread_count,
    const std::vector<AboveParent>& parents, const Gram& base,
    const ModularBases& base_bases,
    const std::vector<std::uint8_t>& present_ids,
    const std::vector<std::uint8_t>& absent_ids,
    UpdateMethod method, std::atomic<std::uint64_t>& shared_transitions,
    const Clock::time_point& expansion_started,
    double heartbeat_seconds, std::mutex& heartbeat_mutex,
    Clock::time_point& next_heartbeat) {
  WorkerSummary summary;
  const Wide frontier_squared =
      static_cast<Wide>(kFrontierRoot) * kFrontierRoot;

  for (std::size_t parent_position = thread_index;
       parent_position < parents.size();
       parent_position += thread_count) {
    const AboveParent& parent = parents[parent_position];
    const Gram parent_gram = candidate_gram(base, parent.record);
    bool parent_update_available =
        method == UpdateMethod::kParentRank4;
    for (const std::uint64_t prime : kCrtPrimes) {
      if (wide_mod(parent.determinant, prime) == 0) {
        parent_update_available = false;
      }
    }
    ModularBases parent_bases{};
    if (parent_update_available) {
      for (std::size_t prime_index = 0;
           prime_index < kCrtPrimes.size(); ++prime_index) {
        parent_bases[prime_index] =
            modular_base(parent_gram, kCrtPrimes[prime_index]);
        if (parent_bases[prime_index].determinant !=
            wide_mod(parent.determinant, kCrtPrimes[prime_index])) {
          throw std::runtime_error(
              "parent modular-base determinant mismatch");
        }
      }
    } else if (method == UpdateMethod::kParentRank4) {
      ++summary.singular_parent_prime_fallbacks;
    }

    std::uint64_t local_transition = 0;
    std::uint64_t removal_count = 0;
    for (const std::uint8_t removed : present_ids) {
      if (contains_edge(parent.record.removed_global, removed)) continue;
      ++removal_count;
      std::uint64_t addition_count = 0;
      for (const std::uint8_t added : absent_ids) {
        if (contains_edge(parent.record.added_global, added)) continue;
        ++addition_count;
        const Radius4Key key =
            add_outward_edges(parent.record, removed, added);
        const auto residues =
            parent_update_available
                ? parent_update_residues(removed, added, parent_bases)
                : base_update_residues(key, base_bases);
        const Wide determinant = reconstruct_crt(residues);

        const std::uint64_t global_transition =
            parent_position * kTransitionsPerParent + local_transition;
        if (global_transition % 131'071 == 0) {
          const Gram candidate = radius4_gram(base, key);
          for (std::size_t prime_index = 0;
               prime_index < kCrtPrimes.size(); ++prime_index) {
            if (residues[prime_index] != determinant_modulo(
                    candidate, kOrder, kCrtPrimes[prime_index])) {
              throw std::runtime_error(
                  "radius-4 update/full determinant residue mismatch");
            }
          }
          if (exact_determinant(candidate) != determinant ||
              determinant_modulo(
                  candidate, kOrder, kIndependentPrime) !=
                  wide_mod(determinant, kIndependentPrime)) {
            throw std::runtime_error(
                "radius-4 exact/independent-prime spot-check failed");
          }
          ++summary.direct_spot_checks;
        }

        ++summary.transitions;
        ++local_transition;
        if (!summary.has_maximum ||
            determinant > summary.maximum_determinant ||
            (determinant == summary.maximum_determinant &&
             key < summary.maximum_key)) {
          summary.has_maximum = true;
          summary.maximum_determinant = determinant;
          summary.maximum_key = key;
        }
        if (determinant > 0) {
          ++summary.positive_determinant_transitions;
        }
        if (determinant > frontier_squared) {
          ++summary.above_frontier_transitions;
        }
        if (determinant > 0) {
          const Wide root = integer_square_root(determinant);
          if (root * root == determinant) {
            ++summary.square_transitions;
            if (root > static_cast<Wide>(kFrontierRoot)) {
              ++summary.above_frontier_square_transitions;
              const Gram candidate = radius4_gram(base, key);
              if (exact_determinant(candidate) != determinant ||
                  determinant_modulo(
                      candidate, kOrder, kIndependentPrime) !=
                      wide_mod(determinant, kIndependentPrime)) {
                throw std::runtime_error(
                    "relevant square determinant cross-check failed");
              }
              require_bareiss_match(candidate, determinant);
              const bool divisible =
                  root % kRequiredDivisor == 0;
              const bool positive_definite =
                  exact_positive_definite(candidate);
              if (divisible) {
                ++summary
                     .divisible_above_frontier_square_transitions;
              }
              if (divisible && positive_definite) {
                ++summary.positive_definite_route_transitions;
              }
              RelevantHit hit;
              hit.key = key;
              hit.determinant = determinant;
              hit.root = root;
              hit.root_divisible = divisible;
              hit.positive_definite = positive_definite;
              const auto [iterator, inserted] =
                  summary.relevant_hits.emplace(key, hit);
              if (!inserted &&
                  (iterator->second.determinant != determinant ||
                   iterator->second.root != root ||
                   iterator->second.root_divisible != divisible ||
                   iterator->second.positive_definite !=
                       positive_definite)) {
                throw std::runtime_error(
                    "duplicate relevant transition is inconsistent");
              }
            }
          }
        }
      }
      if (addition_count != kOutwardAdditionsPerParent) {
        throw std::runtime_error(
            "parent does not have exactly 205 outward additions");
      }
    }
    if (removal_count != kOutwardRemovalsPerParent ||
        local_transition != kTransitionsPerParent) {
      throw std::runtime_error(
          "parent does not have exactly 8,610 outward transitions");
    }
    ++summary.parents;
    const std::uint64_t completed =
        shared_transitions.fetch_add(
            kTransitionsPerParent, std::memory_order_relaxed) +
        kTransitionsPerParent;
    if (heartbeat_seconds > 0.0) {
      const Clock::time_point now = Clock::now();
      if (now >= next_heartbeat) {
        std::lock_guard<std::mutex> lock(heartbeat_mutex);
        if (now >= next_heartbeat) {
          const double elapsed =
              std::chrono::duration<double>(
                  now - expansion_started).count();
          const double rate =
              elapsed > 0.0 ? completed / elapsed : 0.0;
          const double eta =
              rate > 0.0
                  ? (kExpectedTransitions - completed) / rate
                  : 0.0;
          std::cerr
              << "{\"elapsed_seconds\":" << std::fixed
              << std::setprecision(3) << elapsed
              << ",\"eta_seconds\":" << eta
              << ",\"event\":\"radius4-basin-expansion-heartbeat\""
              << ",\"method\":" << json_escape(method_name(method))
              << ",\"rate\":" << rate
              << ",\"transitions\":" << completed << "}\n";
          next_heartbeat =
              now + std::chrono::duration_cast<Clock::duration>(
                        std::chrono::duration<double>(
                            heartbeat_seconds));
        }
      }
    }
  }
  return summary;
}

BasinSummary screen_basin(
    const RecoveryResult& recovery, const Gram& base,
    const ModularBases& base_bases,
    const std::vector<Edge>& present, UpdateMethod method,
    std::uint64_t requested_threads, double heartbeat_seconds) {
  BasinSummary summary;
  summary.radius3_catalog_records = recovery.record_count;
  summary.radius3_catalog_orbit_size_sum = recovery.orbit_size_sum;
  summary.above_parent_labeled = recovery.above_labeled;
  summary.recovery_seconds = recovery.elapsed_seconds;

  const std::vector<Edge> absent = absent_reference_edges(present);
  const std::vector<std::uint8_t> present_ids =
      global_edge_ids(present);
  const std::vector<std::uint8_t> absent_ids =
      global_edge_ids(absent);
  const std::size_t thread_count =
      std::min<std::size_t>(
          static_cast<std::size_t>(requested_threads),
          recovery.parents.size());
  const Clock::time_point started = Clock::now();
  Clock::time_point next_heartbeat =
      started + std::chrono::duration_cast<Clock::duration>(
                    std::chrono::duration<double>(heartbeat_seconds));
  std::atomic<std::uint64_t> shared_transitions{0};
  std::mutex heartbeat_mutex;
  std::vector<WorkerSummary> worker_summaries(thread_count);
  std::vector<std::exception_ptr> errors(thread_count);
  std::vector<std::thread> workers;
  workers.reserve(thread_count);
  for (std::size_t thread_index = 0;
       thread_index < thread_count; ++thread_index) {
    workers.emplace_back([&, thread_index] {
      try {
        worker_summaries[thread_index] =
            screen_parent_partition(
                thread_index, thread_count, recovery.parents, base,
                base_bases, present_ids, absent_ids, method,
                shared_transitions, started, heartbeat_seconds,
                heartbeat_mutex, next_heartbeat);
      } catch (...) {
        errors[thread_index] = std::current_exception();
      }
    });
  }
  for (std::thread& worker : workers) worker.join();
  for (const std::exception_ptr& error : errors) {
    if (error) std::rethrow_exception(error);
  }
  for (WorkerSummary& worker_summary : worker_summaries) {
    merge_worker_summary(summary, std::move(worker_summary));
  }
  summary.expansion_seconds =
      std::chrono::duration<double>(Clock::now() - started).count();
  if (summary.parents != kExpectedAboveParents ||
      summary.transitions != kExpectedTransitions ||
      shared_transitions != kExpectedTransitions) {
    throw std::runtime_error(
        "radius-4 transition coverage count mismatch");
  }
  if (!summary.has_maximum) {
    throw std::runtime_error("radius-4 screen has no maximum");
  }

  const Gram maximum = radius4_gram(base, summary.maximum_key);
  if (exact_determinant(maximum) != summary.maximum_determinant ||
      determinant_modulo(maximum, kOrder, kIndependentPrime) !=
          wide_mod(summary.maximum_determinant, kIndependentPrime)) {
    throw std::runtime_error("maximum determinant replay mismatch");
  }
  require_bareiss_match(maximum, summary.maximum_determinant);
  summary.maximum_positive_definite =
      exact_positive_definite(maximum);

  if (!summary.relevant_hits.empty()) {
    const std::vector<Permutation> group =
        generate_automorphism_group();
    validate_automorphism_group(group, present);
    std::set<Radius4Key> relevant_orbits;
    std::set<Radius4Key> route_states;
    std::set<Radius4Key> route_orbits;
    for (auto& [key, hit] : summary.relevant_hits) {
      hit.canonical_key = canonical_key(key, group);
      relevant_orbits.insert(hit.canonical_key);
      if (hit.root_divisible && hit.positive_definite) {
        route_states.insert(key);
        route_orbits.insert(hit.canonical_key);
      }
    }
    summary.unique_relevant_states = summary.relevant_hits.size();
    summary.unique_relevant_orbits = relevant_orbits.size();
    summary.unique_route_states = route_states.size();
    summary.unique_route_orbits = route_orbits.size();
    std::set<Radius4Key> emitted;
    for (const auto& [key, hit] : summary.relevant_hits) {
      (void)key;
      if (emitted.insert(hit.canonical_key).second) {
        summary.canonical_hits.push_back(hit);
      }
    }
  }
  return summary;
}

void append_edge_list(
    std::ostream& output, const std::vector<Edge>& edges) {
  output << '[';
  for (std::size_t index = 0; index < edges.size(); ++index) {
    if (index != 0) output << ',';
    output << '['
           << static_cast<unsigned>(edges[index].first) + 1U << ','
           << static_cast<unsigned>(edges[index].second) + 1U << ']';
  }
  output << ']';
}

void append_edge_id_array(
    std::ostream& output,
    const std::array<std::uint8_t, 4>& ids) {
  output << '[';
  for (int index = 0; index < 4; ++index) {
    if (index != 0) output << ',';
    output << static_cast<unsigned>(ids[index]);
  }
  output << ']';
}

void append_key_json(std::ostream& output, const Radius4Key& key) {
  output << "{\"added_edge_ids\":";
  append_edge_id_array(output, key.added);
  output << ",\"removed_edge_ids\":";
  append_edge_id_array(output, key.removed);
  output << '}';
}

std::string report_json(
    const BasinSummary& summary, const BasinArguments& arguments,
    std::string_view reference_sha256,
    std::string_view reference_gram_sha256,
    std::string_view parents_sha256) {
  const Wide frontier_squared =
      static_cast<Wide>(kFrontierRoot) * kFrontierRoot;
  const Wide maximum_root =
      summary.maximum_determinant > 0
          ? integer_square_root(summary.maximum_determinant)
          : 0;
  std::ostringstream output;
  output << std::setprecision(17);
  output << "{\"catalog\":{";
  output << "\"path\":" << json_escape(arguments.catalog.string());
  output << ",\"sha256\":" << json_escape(arguments.catalog_sha256);
  output << ",\"frozen\":true}";
  output << ",\"challenge_id\":\"maxdet-23-v1\"";
  output << ",\"claim_boundary\":"
         << json_escape(
                "Exact only for radius-4 Aut(B0)-orbits having at least "
                "one radius-3 parent with det(G)>frontier^2. The "
                "16,858,380 count is parent-child transitions, not unique "
                "states or orbits. A routed Gram remains only a necessary "
                "condition for a {-1,+1} factor.");
  output << ",\"complete\":true";
  output << ",\"coverage\":{";
  output << "\"above_parent_labeled\":"
         << summary.above_parent_labeled;
  output << ",\"above_parent_representatives\":"
         << kExpectedAboveParents;
  output << ",\"all_parent_radius4_children_enumerated\":true";
  output << ",\"boundary\":"
         << json_escape(
                "Every radius-4 state with at least one strict-above-"
                "frontier-squared radius-3 parent, modulo Aut(B0)");
  output << ",\"children_per_parent\":"
         << kTransitionsPerParent;
  output << ",\"outward_additions_per_parent\":"
         << kOutwardAdditionsPerParent;
  output << ",\"outward_removals_per_parent\":"
         << kOutwardRemovalsPerParent;
  output << ",\"radius3_parents_per_radius4_state\":16";
  output << ",\"transition_count\":" << summary.transitions;
  output << ",\"transition_count_is_unique_states\":false";
  output << ",\"transition_count_is_unique_orbits\":false}";
  output << ",\"elapsed\":{";
  output << "\"expansion_seconds\":" << summary.expansion_seconds;
  output << ",\"parent_recovery_seconds\":"
         << summary.recovery_seconds;
  output << ",\"total_seconds\":"
         << summary.expansion_seconds + summary.recovery_seconds << '}';
  output << ",\"engine\":\"gram-radius4-frontier-basin\"";
  output << ",\"frontier_root\":"
         << json_escape(std::to_string(kFrontierRoot));
  output << ",\"frontier_squared\":"
         << json_escape(decimal(frontier_squared));
  output << ",\"maximum\":{";
  output << "\"determinant\":"
         << json_escape(decimal(summary.maximum_determinant));
  output << ",\"is_square\":"
         << (maximum_root * maximum_root ==
                     summary.maximum_determinant
                 ? "true"
                 : "false");
  output << ",\"key\":";
  append_key_json(output, summary.maximum_key);
  output << ",\"positive_definite\":"
         << (summary.maximum_positive_definite ? "true" : "false");
  output << ",\"square_root_floor\":"
         << json_escape(decimal(maximum_root)) << '}';
  output << ",\"method\":{";
  output << "\"crt_primes\":[";
  for (std::size_t index = 0; index < kCrtPrimes.size(); ++index) {
    if (index != 0) output << ',';
    output << kCrtPrimes[index];
  }
  output << ']';
  output << ",\"independent_prime\":" << kIndependentPrime;
  output << ",\"name\":"
         << json_escape(method_name(arguments.method));
  output << ",\"threads\":" << arguments.threads;
  output << ",\"update_identity\":"
         << json_escape(
                "det(B+P Delta P^T)=det(B)det(I+Delta P^T B^-1 P)");
  output << ",\"uniqueness_bound\":"
         << json_escape(
                "727^(23/2)<2^110; four-prime CRT modulus exceeds 2^123")
         << '}';
  output << ",\"normalization\":\"G=24I-J+4A\"";
  output << ",\"parents_snapshot\":{";
  output << "\"path\":"
         << json_escape(arguments.parents_snapshot.string());
  output << ",\"sha256\":" << json_escape(parents_sha256) << '}';
  output << ",\"reference\":{";
  output << "\"gram_sha256\":"
         << json_escape(reference_gram_sha256);
  output << ",\"path\":" << json_escape(arguments.reference.string());
  output << ",\"raw_sha256\":"
         << json_escape(reference_sha256) << '}';
  output << ",\"relevant_orbit_hits\":[";
  for (std::size_t index = 0;
       index < summary.canonical_hits.size(); ++index) {
    if (index != 0) output << ',';
    const RelevantHit& hit = summary.canonical_hits[index];
    output << "{\"canonical_key\":";
    append_key_json(output, hit.canonical_key);
    output << ",\"determinant\":"
           << json_escape(decimal(hit.determinant));
    output << ",\"positive_definite\":"
           << (hit.positive_definite ? "true" : "false");
    output << ",\"root_divisible_by_2_22\":"
           << (hit.root_divisible ? "true" : "false");
    output << ",\"square_root\":"
           << json_escape(decimal(hit.root)) << '}';
  }
  output << ']';
  output << ",\"schema_version\":1";
  output << ",\"statistics\":{";
  output << "\"above_frontier_square_transitions\":"
         << summary.above_frontier_square_transitions;
  output << ",\"above_frontier_transitions\":"
         << summary.above_frontier_transitions;
  output << ",\"direct_spot_checks\":"
         << summary.direct_spot_checks;
  output << ",\"divisible_above_frontier_square_transitions\":"
         << summary.divisible_above_frontier_square_transitions;
  output << ",\"positive_definite_route_transitions\":"
         << summary.positive_definite_route_transitions;
  output << ",\"positive_determinant_transitions\":"
         << summary.positive_determinant_transitions;
  output << ",\"radius3_catalog_orbit_size_sum\":"
         << summary.radius3_catalog_orbit_size_sum;
  output << ",\"radius3_catalog_records\":"
         << summary.radius3_catalog_records;
  output << ",\"singular_parent_prime_fallbacks\":"
         << summary.singular_parent_prime_fallbacks;
  output << ",\"square_transitions\":"
         << summary.square_transitions;
  output << ",\"unique_above_frontier_square_orbits\":"
         << summary.unique_relevant_orbits;
  output << ",\"unique_above_frontier_square_states\":"
         << summary.unique_relevant_states;
  output << ",\"unique_route_orbits\":"
         << summary.unique_route_orbits;
  output << ",\"unique_route_states\":"
         << summary.unique_route_states << '}';
  output << ",\"termination\":\"completed\"}\n";
  return output.str();
}

std::string route_json(
    const BasinSummary& summary, const BasinArguments& arguments,
    const std::vector<Edge>& present,
    std::string_view parents_sha256) {
  std::ostringstream output;
  output << std::setprecision(17);
  output << "{\"challenge_id\":\"maxdet-23-v1\"";
  output << ",\"claim_boundary\":"
         << json_escape(
                "Exact radius-4 frontier-basin Gram screen. Stored hits "
                "are one exact representative per Aut(B0) orbit and still "
                "require Hasse, shell, and {-1,+1} factor checks.");
  output << ",\"complete\":true";
  // The established Hasse and shell consumers accept this binary-Gram
  // compatibility engine. The mode below preserves the true producer.
  output << ",\"engine\":\"gram-tabu\"";
  output << ",\"frontier_root\":"
         << json_escape(std::to_string(kFrontierRoot));
  output << ",\"frontier_squared\":"
         << json_escape(decimal(
                static_cast<Wide>(kFrontierRoot) * kFrontierRoot));
  output << ",\"hits\":[";
  bool first_hit = true;
  for (const RelevantHit& hit : summary.canonical_hits) {
    if (!hit.root_divisible || !hit.positive_definite) continue;
    if (!first_hit) output << ',';
    first_hit = false;
    const std::vector<Edge> edges =
        defect_edges(present, hit.canonical_key);
    output << "{\"determinant\":"
           << json_escape(decimal(hit.determinant));
    output << ",\"divisible_by_2_22\":"
           << (hit.root_divisible ? "true" : "false");
    output << ",\"edge_count\":" << edges.size();
    output << ",\"edges\":";
    append_edge_list(output, edges);
    output << ",\"positive_definite\":"
           << (hit.positive_definite ? "true" : "false");
    output << ",\"qualified\":"
           << (hit.root_divisible && hit.positive_definite
                   ? "true"
                   : "false");
    output << ",\"square_root\":"
           << json_escape(decimal(hit.root)) << '}';
  }
  output << ']';
  output << ",\"mode\":\"reference-radius4-frontier-basin-exact\"";
  output << ",\"normalization\":\"G=24I-J+4A\"";
  output << ",\"parameters\":{";
  output << "\"catalog_sha256\":"
         << json_escape(arguments.catalog_sha256);
  output << ",\"max_stored_hits\":"
         << summary.unique_route_orbits;
  output << ",\"parents_snapshot_sha256\":"
         << json_escape(parents_sha256) << '}';
  output << ",\"schema_version\":1";
  output << ",\"statistics\":{";
  output << "\"exact_squares\":"
         << summary.above_frontier_square_transitions;
  output << ",\"qualified_survivors\":"
         << summary.positive_definite_route_transitions;
  output << ",\"unrecorded_square_observations\":"
         << (summary.unique_relevant_states -
             summary.unique_route_orbits)
         << '}';
  output << ",\"termination\":\"completed\"}\n";
  return output.str();
}

void run_basin_self_test(const Gram& base, const ModularBases& base_bases) {
  if (sha256("abc") !=
      "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb"
      "410ff61f20015ad") {
    throw std::runtime_error("SHA-256 self-test failed");
  }
  const Wide frontier_squared =
      static_cast<Wide>(kFrontierRoot) * kFrontierRoot;
  if (exact_determinant(base) != frontier_squared ||
      bareiss_determinant(base).get_str() !=
          decimal(frontier_squared)) {
    throw std::runtime_error("reference determinant self-test failed");
  }
  const std::vector<Edge> present = validate_reference_gram(base);
  const std::vector<Edge> absent = absent_reference_edges(present);
  const std::vector<std::uint8_t> present_ids =
      global_edge_ids(present);
  const std::vector<std::uint8_t> absent_ids =
      global_edge_ids(absent);
  OrbitRecord parent;
  parent.removed_global.size = 3;
  parent.added_global.size = 3;
  for (int index = 0; index < 3; ++index) {
    parent.removed_global.values[index] = present_ids[index];
    parent.added_global.values[index] = absent_ids[index];
  }
  const Gram parent_gram = candidate_gram(base, parent);
  const Wide parent_determinant = exact_determinant(parent_gram);
  bool parent_invertible = true;
  for (const std::uint64_t prime : kCrtPrimes) {
    if (wide_mod(parent_determinant, prime) == 0) {
      parent_invertible = false;
    }
  }
  ModularBases parent_bases{};
  if (parent_invertible) {
    for (std::size_t index = 0; index < kCrtPrimes.size(); ++index) {
      parent_bases[index] =
          modular_base(parent_gram, kCrtPrimes[index]);
    }
  }
  const std::uint8_t removed = present_ids[3];
  const std::uint8_t added = absent_ids[3];
  const Radius4Key key =
      add_outward_edges(parent, removed, added);
  const Gram child = radius4_gram(base, key);
  const Wide direct = exact_determinant(child);
  if (reconstruct_crt(base_update_residues(key, base_bases)) != direct) {
    throw std::runtime_error("base-rank16 update self-test failed");
  }
  if (parent_invertible &&
      reconstruct_crt(
          parent_update_residues(removed, added, parent_bases)) !=
          direct) {
    throw std::runtime_error("parent-rank4 update self-test failed");
  }
  if (determinant_modulo(child, kOrder, kIndependentPrime) !=
      wide_mod(direct, kIndependentPrime)) {
    throw std::runtime_error("independent-prime self-test failed");
  }
  require_bareiss_match(child, direct);
  if (kExpectedTransitions != 16'858'380 ||
      kTransitionsPerParent != 8'610) {
    throw std::runtime_error("transition-count arithmetic self-test failed");
  }
  std::cout
      << "SELF-TEST PASS"
      << " base_rank16=true"
      << " parent_rank4=" << (parent_invertible ? "true" : "skipped")
      << " transition_count=" << kExpectedTransitions << '\n';
}

int run_main(int argc, char** argv) {
  const BasinArguments arguments =
      parse_basin_arguments(argc, argv);
  const std::string reference_bytes =
      read_file_bytes(arguments.reference);
  const std::string reference_sha256 = sha256(reference_bytes);
  const SignMatrix reference = parse_sign_matrix(reference_bytes);
  const Gram base = gram_of(reference);
  const std::string reference_gram_sha256 =
      sha256(gram_bytes(base));
  const std::vector<Edge> present =
      validate_reference_gram(base);
  const ModularBases base_bases = modular_bases(base);

  if (arguments.self_test) {
    run_basin_self_test(base, base_bases);
    return 0;
  }
  const std::string actual_catalog_sha256 =
      sha256_file(arguments.catalog);
  if (actual_catalog_sha256 != kFrozenCatalogSha256 ||
      actual_catalog_sha256 != arguments.catalog_sha256) {
    throw std::runtime_error(
        "catalog bytes do not match the frozen SHA-256");
  }
  const RecoveryResult recovery = recover_above_parents(
      arguments.catalog, reference_sha256,
      reference_gram_sha256, base, base_bases,
      arguments.heartbeat_seconds);
  const std::string parents_bytes =
      parents_snapshot_bytes(recovery, actual_catalog_sha256);
  const std::string parents_sha256 = sha256(parents_bytes);
  atomic_write(arguments.parents_snapshot, parents_bytes);

  const BasinSummary summary = screen_basin(
      recovery, base, base_bases, present, arguments.method,
      arguments.threads, arguments.heartbeat_seconds);
  atomic_write(
      arguments.output,
      report_json(
          summary, arguments, reference_sha256,
          reference_gram_sha256, parents_sha256));
  atomic_write(
      arguments.route_snapshot,
      route_json(summary, arguments, present, parents_sha256));
  std::cout
      << "SCREEN COMPLETE"
      << " method=" << method_name(arguments.method)
      << " parents=" << summary.parents
      << " transitions=" << summary.transitions
      << " above_frontier_squares="
      << summary.above_frontier_square_transitions
      << " route_orbits=" << summary.unique_route_orbits
      << " maximum=" << decimal(summary.maximum_determinant)
      << " recovery_seconds=" << std::fixed
      << std::setprecision(6) << summary.recovery_seconds
      << " expansion_seconds=" << summary.expansion_seconds << '\n';
  return 0;
}

}  // namespace radius4_basin

int main(int argc, char** argv) {
  try {
    return radius4_basin::run_main(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "gram_radius4_basin: " << error.what() << '\n';
    return 1;
  }
}
