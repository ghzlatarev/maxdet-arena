// Differential and invariant audit for core_adjugate_tabu.cpp.
//
// The production source is included into this translation unit with its main
// renamed.  This keeps the audit on exactly the same determinant-direction
// and adjugate-update code without adding test-only branches to the search
// binary.

#define main core_adjugate_tabu_embedded_main
#include "core_adjugate_tabu.cpp"
#undef main

namespace {

struct AuditArguments {
  fs::path start;
  std::uint64_t seed = 32001;
  int random_steps = 1000;
  int exact_adjugate_interval = 100;
};

int strict_int(std::string_view text, std::string_view option) {
  const std::uint64_t value = strict_unsigned(text, option);
  if (value > static_cast<std::uint64_t>(INT_MAX)) {
    throw std::runtime_error(std::string(option) + " is too large");
  }
  return static_cast<int>(value);
}

AuditArguments parse_audit_arguments(int argc, char** argv) {
  AuditArguments arguments;
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
    } else if (option == "--seed") {
      arguments.seed = strict_unsigned(value(), option);
    } else if (option == "--random-steps") {
      arguments.random_steps = strict_int(value(), option);
    } else if (option == "--exact-adjugate-interval") {
      arguments.exact_adjugate_interval =
          strict_int(value(), option);
    } else if (option == "--help") {
      std::cout
          << "usage: core_adjugate_audit --start MATRIX "
             "[--seed N] [--random-steps N] "
             "[--exact-adjugate-interval N]\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (arguments.start.empty()) {
    throw std::runtime_error("--start is required");
  }
  if (arguments.random_steps < 0 ||
      arguments.exact_adjugate_interval < 0) {
    throw std::runtime_error("audit counts must be non-negative");
  }
  return arguments;
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

void check_other_adjugate_identity(const State& state) {
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      Wide value = 0;
      for (int inner = 0; inner < kCoreOrder; ++inner) {
        value +=
            static_cast<Wide>(state.adjugate[row][inner]) *
            state.core[inner][column];
      }
      const Wide expected =
          row == column ? state.determinant : 0;
      if (value != expected) {
        throw std::runtime_error(
            "adj(B)*B invariant failed at (" +
            std::to_string(row) + "," +
            std::to_string(column) + ")");
      }
    }
  }
}

void check_exact_adjugate(const State& state) {
  const Adjugate exact = exact_adjugate(state.core);
  if (exact != state.adjugate) {
    throw std::runtime_error(
        "incremental adjugate differs from cofactor recomputation");
  }
}

Move direction_move(const State& state, int id) {
  if (id < 0 || id >= kMoveCount) {
    throw std::runtime_error("direction id is out of range");
  }
  if (id < kCoreEntries) {
    const int row = id / kCoreOrder;
    const int column = id % kCoreOrder;
    return Move{MoveKind::kBit, row, column, id,
                bit_candidate_determinant(state, row, column),
                false};
  }
  if (id < kCoreEntries + kCoreOrder) {
    const int row = id - kCoreEntries;
    return Move{MoveKind::kRowComplement, row, -1, id,
                row_complement_determinant(state, row), false};
  }
  if (id < kCoreEntries + 2 * kCoreOrder) {
    const int column =
        id - kCoreEntries - kCoreOrder;
    return Move{MoveKind::kColumnComplement, column, -1, id,
                column_complement_determinant(state, column),
                false};
  }
  return Move{MoveKind::kWholeComplement, -1, -1, id,
              whole_complement_determinant(state), false};
}

CoreMatrix directly_apply_direction(const CoreMatrix& source,
                                    int id) {
  CoreMatrix result = source;
  if (id < kCoreEntries) {
    result[id / kCoreOrder][id % kCoreOrder] ^= 1U;
    return result;
  }
  if (id < kCoreEntries + kCoreOrder) {
    const int row = id - kCoreEntries;
    for (int column = 0; column < kCoreOrder; ++column) {
      result[row][column] ^= 1U;
    }
    return result;
  }
  if (id < kCoreEntries + 2 * kCoreOrder) {
    const int column =
        id - kCoreEntries - kCoreOrder;
    for (int row = 0; row < kCoreOrder; ++row) {
      result[row][column] ^= 1U;
    }
    return result;
  }
  for (auto& row : result) {
    for (std::uint8_t& entry : row) {
      entry ^= 1U;
    }
  }
  return result;
}

std::pair<int, int> sign_coordinate_for_direction(int id) {
  if (id < kCoreEntries) {
    return {id / kCoreOrder + 1, id % kCoreOrder + 1};
  }
  if (id < kCoreEntries + kCoreOrder) {
    return {id - kCoreEntries + 1, 0};
  }
  if (id < kCoreEntries + 2 * kCoreOrder) {
    return {0, id - kCoreEntries - kCoreOrder + 1};
  }
  return {0, 0};
}

void check_state(const State& state, std::uint64_t hash,
                 const std::array<std::uint64_t, kCoreEntries>&
                     zobrist,
                 bool recompute_adjugate) {
  const std::int64_t exact_determinant =
      exact_core_determinant(state.core);
  if (exact_determinant != state.determinant) {
    throw std::runtime_error(
        "incremental determinant differs from Bareiss");
  }
  if (hash != core_hash(state.core, zobrist)) {
    throw std::runtime_error("incremental Zobrist hash drifted");
  }
  check_adjugate_identity(state);
  check_other_adjugate_identity(state);
  if (recompute_adjugate) {
    check_exact_adjugate(state);
  }
}

void exhaustive_direction_audit(
    const SignMatrix& input, const State& base,
    const std::array<std::uint64_t, kCoreEntries>& zobrist,
    std::uint64_t& nonsingular_updates,
    std::uint64_t& singular_directions) {
  const std::uint64_t base_hash =
      core_hash(base.core, zobrist);
  for (int id = 0; id < kMoveCount; ++id) {
    SignMatrix flipped = input;
    const auto [sign_row, sign_column] =
        sign_coordinate_for_direction(id);
    flipped[sign_row][sign_column] *= -1;
    const CoreMatrix mapped = dephase_to_core(flipped);
    const CoreMatrix expected =
        directly_apply_direction(base.core, id);
    if (mapped != expected) {
      throw std::runtime_error(
          "sign/core direction mapping failed for id " +
          std::to_string(id));
    }

    const Move move = direction_move(base, id);
    const std::int64_t exact_candidate =
        exact_core_determinant(expected);
    if (move.determinant != exact_candidate) {
      throw std::runtime_error(
          "candidate determinant formula failed for id " +
          std::to_string(id));
    }
    const Wide sign_determinant = bareiss(sign_as_wide(flipped));
    const Wide scaled_core =
        static_cast<Wide>(exact_candidate) * kScale;
    if (sign_determinant != scaled_core &&
        sign_determinant != -scaled_core) {
      throw std::runtime_error(
          "sign/core quotient mapping failed for id " +
          std::to_string(id));
    }

    if (exact_candidate == 0) {
      ++singular_directions;
      continue;
    }
    State updated = base;
    std::uint64_t updated_hash = base_hash;
    Statistics statistics;
    apply_move(updated, move, updated_hash, zobrist, statistics);
    if (updated.core != expected) {
      throw std::runtime_error(
          "updated core failed for direction id " +
          std::to_string(id));
    }
    check_state(updated, updated_hash, zobrist, false);
    ++nonsingular_updates;
  }
}

int random_direction_id(int phase, std::mt19937_64& randomizer) {
  if (phase == 0) {
    return static_cast<int>(
        randomizer() % static_cast<std::uint64_t>(kCoreEntries));
  }
  if (phase == 1) {
    return kCoreEntries +
           static_cast<int>(
               randomizer() %
               static_cast<std::uint64_t>(kCoreOrder));
  }
  if (phase == 2) {
    return kCoreEntries + kCoreOrder +
           static_cast<int>(
               randomizer() %
               static_cast<std::uint64_t>(kCoreOrder));
  }
  return kMoveCount - 1;
}

void random_walk_audit(
    State state,
    const std::array<std::uint64_t, kCoreEntries>& zobrist,
    const AuditArguments& arguments,
    std::uint64_t& applied,
    std::array<std::uint64_t, 4>& kind_counts) {
  std::mt19937_64 randomizer(arguments.seed);
  std::uint64_t hash = core_hash(state.core, zobrist);
  for (int step = 0; step < arguments.random_steps; ++step) {
    const int phase = step % 4;
    const int id = random_direction_id(phase, randomizer);
    const Move move = direction_move(state, id);
    const CoreMatrix expected =
        directly_apply_direction(state.core, id);
    const std::int64_t exact_candidate =
        exact_core_determinant(expected);
    if (move.determinant != exact_candidate) {
      throw std::runtime_error(
          "random-walk candidate determinant mismatch at step " +
          std::to_string(step));
    }
    if (exact_candidate == 0) {
      continue;
    }
    Statistics statistics;
    apply_move(state, move, hash, zobrist, statistics);
    if (state.core != expected) {
      throw std::runtime_error(
          "random-walk core mismatch at step " +
          std::to_string(step));
    }
    const bool recompute_adjugate =
        arguments.exact_adjugate_interval != 0 &&
        (step % arguments.exact_adjugate_interval == 0);
    check_state(state, hash, zobrist, recompute_adjugate);

    if ((step & 31) == 0) {
      const SignMatrix sign = core_to_sign(state.core);
      const Wide sign_determinant =
          bareiss(sign_as_wide(sign));
      if (sign_determinant !=
          static_cast<Wide>(state.determinant) * kScale) {
        throw std::runtime_error(
            "random-walk sign/core determinant mismatch");
      }
    }
    ++applied;
    ++kind_counts[static_cast<std::size_t>(phase)];
  }
  check_state(state, hash, zobrist, true);
}

int run_audit(const AuditArguments& arguments) {
  const auto started = Clock::now();
  const SignMatrix input = read_sign_matrix(arguments.start);
  State base;
  base.core = dephase_to_core(input);
  base.determinant = exact_core_determinant(base.core);
  if (base.determinant == 0) {
    throw std::runtime_error("audit start must be nonsingular");
  }
  base.adjugate = exact_adjugate(base.core);
  const auto zobrist = make_zobrist(arguments.seed);
  check_state(base, core_hash(base.core, zobrist), zobrist, true);

  std::uint64_t nonsingular_updates = 0;
  std::uint64_t singular_directions = 0;
  exhaustive_direction_audit(
      input, base, zobrist, nonsingular_updates,
      singular_directions);

  std::uint64_t random_updates = 0;
  std::array<std::uint64_t, 4> kind_counts{};
  random_walk_audit(base, zobrist, arguments, random_updates,
                    kind_counts);

  const double elapsed =
      std::chrono::duration<double>(Clock::now() - started)
          .count();
  std::cout
      << "{\"audit\":\"core-adjugate-v1\",\"passed\":true"
      << ",\"directions_checked\":" << kMoveCount
      << ",\"nonsingular_direction_updates\":"
      << nonsingular_updates
      << ",\"singular_directions\":" << singular_directions
      << ",\"random_steps_requested\":" << arguments.random_steps
      << ",\"random_updates_applied\":" << random_updates
      << ",\"random_bit_updates\":" << kind_counts[0]
      << ",\"random_row_complements\":" << kind_counts[1]
      << ",\"random_column_complements\":" << kind_counts[2]
      << ",\"random_whole_complements\":" << kind_counts[3]
      << ",\"exact_adjugate_interval\":"
      << arguments.exact_adjugate_interval
      << ",\"elapsed_seconds\":" << std::fixed
      << std::setprecision(6) << elapsed << "}\n";
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    return run_audit(parse_audit_arguments(argc, argv));
  } catch (const std::exception& error) {
    std::cerr << "core_adjugate_audit: " << error.what() << '\n';
    return 1;
  }
}
