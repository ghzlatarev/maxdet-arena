// Reuse the independently validated determinant, canonical-connector, and
// JSON primitives without changing the completed distinct-pair engine.
#define main gram_connector12_distinct_pair_main
#include "gram_connector12.cpp"
#undef main

#include <set>

namespace {

constexpr std::uint64_t kExpectedDoublePairGraphs = 765'720;
constexpr std::uint64_t kExpectedTriplePairGraphs = 360;
constexpr std::uint64_t kExpectedUniqueLabeledGraphs =
    kExpectedLabeledConfigurations +
    kExpectedDoublePairGraphs +
    kExpectedTriplePairGraphs;

struct ReuseArguments {
  std::filesystem::path output;
  std::filesystem::path route_snapshot;
  double heartbeat_seconds = 10.0;
};

struct PairUnion2 {
  std::uint16_t edge_mask = 0;
  std::vector<std::array<Connector, 2>> decompositions;
};

struct PairUnion3 {
  std::uint16_t edge_mask = 0;
  std::vector<Configuration> decompositions;
};

ReuseArguments parse_reuse_arguments(int argc, char** argv) {
  ReuseArguments arguments;
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
          << "Usage: gram_connector12_reuse --output FILE "
             "[--route-snapshot FILE] [--heartbeat-seconds S]\n\n"
          << "Exactly enumerate unions of three edge-disjoint K2,2 "
             "connectors, allowing repeated base-block pairs, modulo "
             "S3 x (S4 wr S5).\n";
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

std::uint16_t relation_mask(const Connector& connector) {
  const int second_size = kBlockSizes[connector.second.block];
  std::uint16_t result = 0;
  for (int first = 0; first < kBlockSizes[connector.first.block]; ++first) {
    if (!((connector.first.mask >> first) & 1)) continue;
    for (int second = 0; second < second_size; ++second) {
      if (!((connector.second.mask >> second) & 1)) continue;
      result |= static_cast<std::uint16_t>(
          1U << (first * second_size + second));
    }
  }
  if (std::popcount(result) != 4) {
    throw std::runtime_error("connector relation mask is not K2,2");
  }
  return result;
}

std::vector<PairUnion2> pair_unions2(
    const std::vector<Connector>& options) {
  std::map<std::uint16_t, std::vector<std::array<Connector, 2>>> unions;
  for (std::size_t first = 0; first < options.size(); ++first) {
    const std::uint16_t first_mask = relation_mask(options[first]);
    for (std::size_t second = first + 1; second < options.size(); ++second) {
      const std::uint16_t second_mask = relation_mask(options[second]);
      if ((first_mask & second_mask) != 0) continue;
      unions[static_cast<std::uint16_t>(first_mask | second_mask)]
          .push_back({options[first], options[second]});
    }
  }
  std::vector<PairUnion2> result;
  result.reserve(unions.size());
  for (auto& [mask, decompositions] : unions) {
    if (std::popcount(mask) != 8) {
      throw std::runtime_error("two-connector union does not have 8 edges");
    }
    result.push_back(PairUnion2{mask, std::move(decompositions)});
  }
  return result;
}

std::vector<PairUnion3> pair_unions3(
    const std::vector<Connector>& options) {
  std::map<std::uint16_t, std::vector<Configuration>> unions;
  for (std::size_t first = 0; first < options.size(); ++first) {
    const std::uint16_t first_mask = relation_mask(options[first]);
    for (std::size_t second = first + 1; second < options.size(); ++second) {
      const std::uint16_t second_mask = relation_mask(options[second]);
      if ((first_mask & second_mask) != 0) continue;
      for (std::size_t third = second + 1; third < options.size(); ++third) {
        const std::uint16_t third_mask = relation_mask(options[third]);
        if (((first_mask | second_mask) & third_mask) != 0) continue;
        unions[static_cast<std::uint16_t>(
                   first_mask | second_mask | third_mask)]
            .push_back(Configuration{{
                options[first], options[second], options[third]}});
      }
    }
  }
  std::vector<PairUnion3> result;
  result.reserve(unions.size());
  for (auto& [mask, decompositions] : unions) {
    if (std::popcount(mask) != 12) {
      throw std::runtime_error("three-connector union does not have 12 edges");
    }
    result.push_back(PairUnion3{mask, std::move(decompositions)});
  }
  return result;
}

std::string reuse_route_snapshot_json(
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
                "Compatibility snapshot produced by "
                "gram-connector12-reuse; survivors still require exact "
                "sign-matrix factorization.");
  output << ",\"complete\":true";
  output << ",\"elapsed_seconds\":" << elapsed_seconds;
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
  output << ",\"mode\":\"connector12-reuse-exact-orbit-enumeration\"";
  output << ",\"normalization\":\"G=24I-J+4A\"";
  output << ",\"parameters\":{";
  output << "\"max_stored_hits\":" << route_hits.size();
  output << ",\"orbit_count\":" << orbit_count;
  output << ",\"producer_engine\":\"gram-connector12-reuse\"}";
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
    const ReuseArguments arguments =
        parse_reuse_arguments(argc, argv);
    const Clock::time_point started = Clock::now();
    Clock::time_point next_heartbeat =
        started + std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(
                          arguments.heartbeat_seconds));

    const SequenceCanonicalizer sequence_canonicalizer;
    const std::vector<BlockPair> pairs = block_pairs();
    std::array<std::vector<Connector>, 15> options{};
    std::array<std::vector<PairUnion2>, 15> unions2{};
    std::array<std::vector<PairUnion3>, 15> unions3{};
    for (std::size_t index = 0; index < pairs.size(); ++index) {
      options[index] = connector_options(pairs[index]);
      unions2[index] = pair_unions2(options[index]);
      unions3[index] = pair_unions3(options[index]);
    }

    // These local counts independently characterize all repeated-pair
    // deduplication before global symmetry reduction.
    for (std::size_t index = 0; index < pairs.size(); ++index) {
      const bool k3_k4 = pairs[index].first == 0;
      const std::size_t expected2 = k3_k4 ? 21U : 174U;
      const std::size_t expected3 = k3_k4 ? 0U : 36U;
      if (unions2[index].size() != expected2 ||
          unions3[index].size() != expected3) {
        throw std::runtime_error("local repeated-pair union count mismatch");
      }
    }

    std::unordered_map<std::uint64_t, OrbitData> orbits;
    std::uint64_t distinct_pair_graphs = 0;
    std::uint64_t double_pair_graphs = 0;
    std::uint64_t triple_pair_graphs = 0;

    // Type 4+4+4: each occupied block pair is distinct, so the connector
    // decomposition is unique.
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
                ++distinct_pair_graphs;
              }
            }
          }
        }
      }
      const Clock::time_point now = Clock::now();
      if (now >= next_heartbeat) {
        std::cerr
            << "{\"elapsed_seconds\":"
            << std::chrono::duration<double>(now - started).count()
            << ",\"event\":\"distinct-pair-heartbeat\""
            << ",\"labeled_unique_graphs\":" << distinct_pair_graphs
            << ",\"orbit_count\":" << orbits.size() << "}\n";
        next_heartbeat =
            now + std::chrono::duration_cast<Clock::duration>(
                      std::chrono::duration<double>(
                          arguments.heartbeat_seconds));
      }
    }

    // Type 8+4: the 8-edge relation identifies the repeated pair. Minimize
    // over every K2,2 decomposition of that relation before recording the
    // graph's orbit key.
    for (int repeated_pair = 0; repeated_pair < 15; ++repeated_pair) {
      for (int single_pair = 0; single_pair < 15; ++single_pair) {
        if (single_pair == repeated_pair) continue;
        for (const PairUnion2& local_union : unions2[repeated_pair]) {
          for (const Connector& single : options[single_pair]) {
            std::uint64_t best =
                std::numeric_limits<std::uint64_t>::max();
            for (const auto& decomposition :
                 local_union.decompositions) {
              const Configuration configuration{{
                  decomposition[0], decomposition[1], single}};
              best = std::min(
                  best,
                  canonical_key(
                      configuration, sequence_canonicalizer));
            }
            if (best == std::numeric_limits<std::uint64_t>::max()) {
              throw std::runtime_error("empty two-connector decomposition");
            }
            ++orbits[best].labeled_count;
            ++double_pair_graphs;
          }
        }
      }
    }

    // Type 12: all three connectors occupy one block pair. A K3--K4 pair
    // has no such edge-disjoint union; each K4--K4 union has five
    // decompositions, all minimized here.
    for (int repeated_pair = 0; repeated_pair < 15; ++repeated_pair) {
      for (const PairUnion3& local_union : unions3[repeated_pair]) {
        std::uint64_t best =
            std::numeric_limits<std::uint64_t>::max();
        for (const Configuration& decomposition :
             local_union.decompositions) {
          best = std::min(
              best,
              canonical_key(decomposition, sequence_canonicalizer));
        }
        if (best == std::numeric_limits<std::uint64_t>::max()) {
          throw std::runtime_error("empty three-connector decomposition");
        }
        ++orbits[best].labeled_count;
        ++triple_pair_graphs;
      }
    }

    if (distinct_pair_graphs != kExpectedLabeledConfigurations ||
        double_pair_graphs != kExpectedDoublePairGraphs ||
        triple_pair_graphs != kExpectedTriplePairGraphs) {
      throw std::runtime_error("labeled unique graph count mismatch");
    }
    const std::uint64_t labeled_unique_graphs =
        distinct_pair_graphs +
        double_pair_graphs +
        triple_pair_graphs;
    if (labeled_unique_graphs != kExpectedUniqueLabeledGraphs) {
      throw std::runtime_error("total labeled unique graph count mismatch");
    }

    const std::uint64_t orbit_multiplicity_sum =
        std::accumulate(
            orbits.begin(), orbits.end(), std::uint64_t{0},
            [](std::uint64_t sum, const auto& entry) {
              return sum + entry.second.labeled_count;
            });
    if (orbit_multiplicity_sum != labeled_unique_graphs) {
      throw std::runtime_error(
          "orbit multiplicities do not sum to unique labeled graphs");
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

    const Wide frontier_squared =
        static_cast<Wide>(kFrontierRoot) *
        static_cast<Wide>(kFrontierRoot);
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

    for (const std::uint64_t key : keys) {
      const Configuration configuration = decode_key(key);
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
    }

    const double elapsed_seconds =
        std::chrono::duration<double>(Clock::now() - started).count();
    const Configuration maximum_configuration = decode_key(maximum_key);

    std::ostringstream report;
    report << std::setprecision(17);
    report << "{\"challenge_id\":\"maxdet-23-v1\"";
    report << ",\"claim_boundary\":"
           << json_escape(
                  "Exact Gram-orbit screening only; a qualified square "
                  "still requires an exact factor in {-1,+1}^{23x23}.");
    report << ",\"complete\":true";
    report << ",\"elapsed_seconds\":" << elapsed_seconds;
    report << ",\"engine\":\"gram-connector12-reuse\"";
    report << ",\"family\":{";
    report << "\"base\":\"K3 disjoint-union 5K4\"";
    report << ",\"connector\":\"K2,2\"";
    report << ",\"connector_count\":3";
    report << ",\"connector_edges_pairwise_disjoint\":true";
    report << ",\"distinct_pair_graphs\":" << distinct_pair_graphs;
    report << ",\"double_pair_graphs\":" << double_pair_graphs;
    report << ",\"triple_pair_graphs\":" << triple_pair_graphs;
    report << ",\"unique_labeled_graph_count\":"
           << labeled_unique_graphs;
    report << ",\"symmetry_group\":"
           << json_escape("S3 x (S4 wr S5)");
    report << ",\"symmetry_group_order\":"
           << kBaseAutomorphismGroupOrder;
    report << ",\"symmetry_reduced_orbit_count\":" << orbits.size();
    report << ",\"orbit_size_histogram\":{";
    bool first_histogram_entry = true;
    for (const auto& [orbit_size, count] : orbit_size_histogram) {
      if (!first_histogram_entry) report << ',';
      first_histogram_entry = false;
      report << json_escape(std::to_string(orbit_size))
             << ':' << count;
    }
    report << '}';
    report << ",\"local_repeated_pair_unions\":{";
    report << "\"K3-K4\":{\"two_connectors\":21,"
              "\"three_connectors\":0}";
    report << ",\"K4-K4\":{\"two_connectors\":174,"
              "\"three_connectors\":36}}}";
    report << ",\"frontier_root\":"
           << json_escape(std::to_string(kFrontierRoot));
    report << ",\"frontier_squared\":"
           << json_escape(decimal(frontier_squared));
    report << ",\"maximum\":{";
    report << "\"added_edges\":";
    append_edge_array(
        report, added_edges(maximum_configuration));
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
    const Configuration published = published_configuration();
    const std::uint64_t published_key =
        canonical_key(published, sequence_canonicalizer);
    report << ",\"published_reproduction\":{";
    report << "\"canonical_key\":"
           << json_escape(std::to_string(published_key));
    report << ",\"determinant\":"
           << json_escape(decimal(exact_determinant(gram(published))));
    report << ",\"square_root\":"
           << json_escape(std::to_string(kFrontierRoot));
    report << ",\"verified\":true}";
    report << ",\"normalization\":\"G=24I-J+4A\"";
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
          reuse_route_snapshot_json(
              route_hits, exact_squares, qualified_survivors,
              keys.size(), elapsed_seconds));
    }

    std::cout
        << "complete unique_labeled=" << labeled_unique_graphs
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
    std::cerr
        << "gram_connector12_reuse: " << error.what() << '\n';
    return 2;
  }
}
