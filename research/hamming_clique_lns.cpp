// Exact Hamming-code clique large-neighborhood search for order-23
// {-1,+1} matrices.
//
// The input is switched into a deterministic orientation:
//
//   * row 0 is the orientation anchor;
//   * every off-diagonal row Gram entry is 3 modulo 4; and
//   * column switches make row 0 all +1.
//
// In this orientation, Gram entries 3, -1, and -5 correspond respectively
// to Hamming distances 10, 12, and 14.  Each LNS iteration destroys t
// non-anchor rows, exhaustively filters all 2^23 oriented row words against
// the fixed rows, builds their exact compatibility graph, and searches for a
// t-clique.  Every unpruned completion is ranked with an exact determinant.
//
// Each repair conditions first on its fixed rows.  If D is their Gram
// determinant, N is the exact integer Schur kernel of the candidate pool,
// N=gK, and D=gd, then a partial repair Q with r rows still needed obeys
//
//   det(full Gram)
//       <= g det(K_Q) product(top-r current K_ii) / d^(t-1).
//
// Conditioning can only decrease each remaining Schur diagonal, so this
// strengthened Fischer bound is exact and sound.  A proper greedy coloring
// rejects candidate subgraphs that cannot contain the remaining clique before
// any determinant work.  Both tests reject only branches that cannot strictly
// improve the exact incumbent.

#include <algorithm>
#include <array>
#include <bit>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

namespace {

namespace fs = std::filesystem;
using BigInt = boost::multiprecision::cpp_int;
using Clock = std::chrono::steady_clock;
using Word = std::uint32_t;

constexpr int kOrder = 23;
constexpr int kNonAnchorRows = kOrder - 1;
constexpr Word kWordCount = Word{1} << kOrder;
constexpr std::uint64_t kPoolPollStride = UINT64_C(4096);
constexpr std::uint64_t kGraphPollStride = UINT64_C(16384);
constexpr std::uint64_t kCliquePollStride = UINT64_C(1024);
constexpr std::uint64_t kDefaultPoolLimit = UINT64_C(20000);
constexpr std::size_t kMaximumSchurKernelBytes =
    std::size_t{128} * 1024U * 1024U;

volatile std::sig_atomic_t stop_requested = 0;

void request_stop(int) { stop_requested = 1; }

struct Arguments {
  fs::path start;
  fs::path output;
  fs::path log;
  fs::path summary;
  std::uint64_t seed = UINT64_C(231023);
  double seconds = 60.0;
  double heartbeat_seconds = 5.0;
  std::uint64_t max_iterations = 0;
  int destroy_min = 4;
  int destroy_max = 7;
  std::uint64_t pool_limit = kDefaultPoolLimit;
  bool destroy_without_replacement = false;
  std::uint64_t destroy_shard_count = 1;
  std::uint64_t destroy_shard_index = 0;
  bool allow_distance_14 = false;
  bool transpose_start = false;
  bool self_test = false;
};

struct SearchStatistics {
  std::uint64_t iterations_started = 0;
  std::uint64_t iterations_completed = 0;
  std::uint64_t pool_words_examined = 0;
  std::uint64_t pools_completed = 0;
  std::uint64_t pool_candidates = 0;
  std::uint64_t maximum_pool_size = 0;
  std::uint64_t graph_pair_checks = 0;
  std::uint64_t graph_edges = 0;
  std::uint64_t graphs_completed = 0;
  std::uint64_t clique_nodes = 0;
  std::uint64_t cardinality_prunes = 0;
  std::uint64_t color_checks = 0;
  std::uint64_t color_prunes = 0;
  std::uint64_t bound_checks = 0;
  std::uint64_t bound_prunes = 0;
  std::uint64_t clique_leaves = 0;
  std::uint64_t exact_matrix_scores = 0;
  std::uint64_t improvements = 0;
  std::uint64_t interrupted_iterations = 0;
};

std::string json_escape(std::string_view input) {
  std::ostringstream output;
  for (const char character : input) {
    switch (character) {
      case '"':
        output << "\\\"";
        break;
      case '\\':
        output << "\\\\";
        break;
      case '\b':
        output << "\\b";
        break;
      case '\f':
        output << "\\f";
        break;
      case '\n':
        output << "\\n";
        break;
      case '\r':
        output << "\\r";
        break;
      case '\t':
        output << "\\t";
        break;
      default:
        if (static_cast<unsigned char>(character) < 0x20U) {
          output << "\\u" << std::hex << std::setw(4)
                 << std::setfill('0')
                 << static_cast<unsigned>(
                        static_cast<unsigned char>(character))
                 << std::dec << std::setfill(' ');
        } else {
          output << character;
        }
    }
  }
  return output.str();
}

std::string big_to_string(const BigInt& value) {
  std::ostringstream output;
  output << value;
  return output.str();
}

std::uint64_t strict_unsigned(
    std::string_view text, std::string_view option) {
  if (text.empty() ||
      !std::all_of(
          text.begin(), text.end(),
          [](const unsigned char character) {
            return character >= static_cast<unsigned char>('0') &&
                character <= static_cast<unsigned char>('9');
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

double strict_double(
    std::string_view text, std::string_view option, bool allow_zero) {
  std::size_t consumed = 0;
  const double result = std::stod(std::string(text), &consumed);
  if (
      consumed != text.size() || !std::isfinite(result) ||
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
    } else if (option == "--max-iterations") {
      arguments.max_iterations = strict_unsigned(value(), option);
    } else if (option == "--destroy-min") {
      const std::uint64_t parsed = strict_unsigned(value(), option);
      if (parsed < 1 || parsed > kNonAnchorRows) {
        throw std::runtime_error("--destroy-min must be in [1,22]");
      }
      arguments.destroy_min = static_cast<int>(parsed);
    } else if (option == "--destroy-max") {
      const std::uint64_t parsed = strict_unsigned(value(), option);
      if (parsed < 1 || parsed > kNonAnchorRows) {
        throw std::runtime_error("--destroy-max must be in [1,22]");
      }
      arguments.destroy_max = static_cast<int>(parsed);
    } else if (option == "--pool-limit") {
      const std::uint64_t parsed = strict_unsigned(value(), option);
      if (parsed == 0) {
        throw std::runtime_error("--pool-limit must be positive");
      }
      arguments.pool_limit = parsed;
    } else if (option == "--destroy-without-replacement") {
      arguments.destroy_without_replacement = true;
    } else if (option == "--destroy-shard-count") {
      const std::uint64_t parsed = strict_unsigned(value(), option);
      if (parsed == 0) {
        throw std::runtime_error(
            "--destroy-shard-count must be positive");
      }
      arguments.destroy_shard_count = parsed;
      arguments.destroy_without_replacement = true;
    } else if (option == "--destroy-shard-index") {
      arguments.destroy_shard_index =
          strict_unsigned(value(), option);
      arguments.destroy_without_replacement = true;
    } else if (option == "--allow-distance-14") {
      arguments.allow_distance_14 = true;
    } else if (option == "--transpose-start") {
      arguments.transpose_start = true;
    } else if (option == "--self-test") {
      arguments.self_test = true;
    } else if (option == "--help") {
      std::cout
          << "usage: hamming_clique_lns --start MATRIX --output MATRIX "
             "--log JSONL --summary JSON [--seed N] [--seconds S] "
             "[--heartbeat-seconds S] [--max-iterations N] "
             "[--destroy-min T] [--destroy-max T] [--pool-limit N] "
             "[--destroy-without-replacement] "
             "[--destroy-shard-count N] [--destroy-shard-index I] "
             "[--allow-distance-14] [--transpose-start]\n"
             "       hamming_clique_lns --self-test\n";
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (arguments.destroy_min > arguments.destroy_max) {
    throw std::runtime_error(
        "--destroy-min cannot exceed --destroy-max");
  }
  if (
      arguments.destroy_shard_index >=
      arguments.destroy_shard_count) {
    throw std::runtime_error(
        "--destroy-shard-index must be smaller than "
        "--destroy-shard-count");
  }
  if (arguments.self_test) return arguments;
  if (
      arguments.start.empty() || arguments.output.empty() ||
      arguments.log.empty() || arguments.summary.empty()) {
    throw std::runtime_error(
        "--start, --output, --log, and --summary are required");
  }
  return arguments;
}

using SignMatrix =
    std::array<std::array<int, kOrder>, kOrder>;
using RowWords = std::array<Word, kOrder>;

SignMatrix transpose_sign_matrix(const SignMatrix& matrix) {
  SignMatrix transposed{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      transposed[row][column] = matrix[column][row];
    }
  }
  return transposed;
}

SignMatrix read_sign_matrix(const fs::path& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("cannot open matrix: " + path.string());
  }
  SignMatrix matrix{};
  std::string line;
  for (int row = 0; row < kOrder; ++row) {
    if (!std::getline(input, line)) {
      throw std::runtime_error("matrix has fewer than 23 rows");
    }
    std::istringstream row_input(line);
    for (int column = 0; column < kOrder; ++column) {
      if (
          !(row_input >> matrix[row][column]) ||
          (matrix[row][column] != -1 &&
           matrix[row][column] != 1)) {
        throw std::runtime_error(
            "invalid sign-matrix entry at row " +
            std::to_string(row + 1));
      }
    }
    std::string extra;
    if (row_input >> extra) {
      throw std::runtime_error(
          "matrix row has more than 23 entries");
    }
  }
  while (std::getline(input, line)) {
    if (line.find_first_not_of(" \t\r\n") != std::string::npos) {
      throw std::runtime_error("matrix has more than 23 rows");
    }
  }
  return matrix;
}

std::vector<std::vector<BigInt>> sign_matrix_as_big(
    const std::vector<Word>& rows, int order) {
  std::vector<std::vector<BigInt>> result(
      rows.size(), std::vector<BigInt>(static_cast<std::size_t>(order)));
  for (std::size_t row = 0; row < rows.size(); ++row) {
    for (int column = 0; column < order; ++column) {
      result[row][static_cast<std::size_t>(column)] =
          ((rows[row] >> column) & Word{1}) != 0U ? -1 : 1;
    }
  }
  return result;
}

BigInt bareiss(std::vector<std::vector<BigInt>> work) {
  if (work.empty()) return 1;
  const std::size_t order = work.size();
  for (const auto& row : work) {
    if (row.size() != order) {
      throw std::runtime_error("Bareiss input must be square");
    }
  }
  if (order == 1) return work[0][0];
  BigInt previous = 1;
  int sign = 1;
  for (std::size_t column = 0; column + 1 < order; ++column) {
    std::size_t pivot = column;
    while (pivot < order && work[pivot][column] == 0) ++pivot;
    if (pivot == order) return 0;
    if (pivot != column) {
      std::swap(work[pivot], work[column]);
      sign = -sign;
    }
    const BigInt pivot_value = work[column][column];
    for (std::size_t row = column + 1; row < order; ++row) {
      for (std::size_t inner = column + 1; inner < order; ++inner) {
        const BigInt numerator =
            work[row][inner] * pivot_value -
            work[row][column] * work[column][inner];
        if (numerator % previous != 0) {
          throw std::runtime_error("non-exact Bareiss division");
        }
        work[row][inner] = numerator / previous;
      }
      work[row][column] = 0;
    }
    previous = pivot_value;
  }
  return sign * work.back().back();
}

BigInt absolute(BigInt value) {
  return value < 0 ? -value : value;
}

BigInt exact_gcd(BigInt first, BigInt second) {
  first = absolute(std::move(first));
  second = absolute(std::move(second));
  while (second != 0) {
    BigInt remainder = first % second;
    first = std::move(second);
    second = std::move(remainder);
  }
  return first;
}

BigInt exact_sign_determinant(
    const std::vector<Word>& rows, int order) {
  if (rows.size() != static_cast<std::size_t>(order)) {
    throw std::runtime_error(
        "exact determinant requires a square row set");
  }
  return bareiss(sign_matrix_as_big(rows, order));
}

BigInt exact_sign_determinant(const RowWords& rows) {
  return exact_sign_determinant(
      std::vector<Word>(rows.begin(), rows.end()), kOrder);
}

std::vector<std::vector<BigInt>> partial_gram_matrix(
    const std::vector<Word>& rows, int order) {
  const std::size_t count = rows.size();
  std::vector<std::vector<BigInt>> gram(
      count, std::vector<BigInt>(count));
  for (std::size_t first = 0; first < count; ++first) {
    gram[first][first] = order;
    for (std::size_t second = first + 1; second < count; ++second) {
      const int distance =
          std::popcount(rows[first] ^ rows[second]);
      const int inner_product = order - 2 * distance;
      gram[first][second] = inner_product;
      gram[second][first] = inner_product;
    }
  }
  return gram;
}

BigInt partial_gram_determinant(
    const std::vector<Word>& rows, int order) {
  const BigInt determinant =
      bareiss(partial_gram_matrix(rows, order));
  if (determinant < 0) {
    throw std::runtime_error(
        "exact Gram determinant was negative");
  }
  return determinant;
}

std::vector<std::vector<BigInt>> exact_adjugate(
    const std::vector<std::vector<BigInt>>& matrix) {
  const std::size_t order = matrix.size();
  for (const auto& row : matrix) {
    if (row.size() != order) {
      throw std::runtime_error(
          "adjugate input must be square");
    }
  }
  std::vector<std::vector<BigInt>> adjugate(
      order, std::vector<BigInt>(order));
  if (order == 0) return adjugate;
  if (order == 1) {
    adjugate[0][0] = 1;
    return adjugate;
  }
  for (std::size_t row = 0; row < order; ++row) {
    for (std::size_t column = 0; column < order; ++column) {
      std::vector<std::vector<BigInt>> minor(
          order - 1, std::vector<BigInt>(order - 1));
      std::size_t minor_row = 0;
      for (std::size_t source_row = 0;
           source_row < order; ++source_row) {
        if (source_row == column) continue;
        std::size_t minor_column = 0;
        for (std::size_t source_column = 0;
             source_column < order; ++source_column) {
          if (source_column == row) continue;
          minor[minor_row][minor_column] =
              matrix[source_row][source_column];
          ++minor_column;
        }
        ++minor_row;
      }
      BigInt cofactor = bareiss(std::move(minor));
      if (((row + column) & 1U) != 0U) {
        cofactor = -cofactor;
      }
      adjugate[row][column] = std::move(cofactor);
    }
  }
  return adjugate;
}

// For G' = [[G,v],[v^T,23]], retain D = det(G) and C = adj(G).
// With q = Cv,
//
//   D' = 23D - v^Tq
//   adj(G') = [[(D'C + qq^T)/D, -q], [-q^T, D]].
//
// A child is prepared with q and D' before the Fischer test.  Its adjugate
// is materialized only if that test survives and the DFS must descend.
class BorderedGramState {
 public:
  explicit BorderedGramState(int vector_order)
      : vector_order_(vector_order) {
    if (vector_order_ < 1 || vector_order_ > kOrder) {
      throw std::runtime_error(
          "bordered Gram state order must be in [1,23]");
    }
  }

  const BigInt& push(Word word) {
    static_cast<void>(prepare_push(word));
    materialize_adjugate();
    return determinant();
  }

  const BigInt& prepare_push(Word word) {
    if (size_ >= vector_order_) {
      throw std::runtime_error(
          "bordered Gram state exceeds vector dimension");
    }
    if (size_ != 0) {
      const Frame& current =
          frames_[static_cast<std::size_t>(size_)];
      if (current.determinant == 0) {
        throw std::runtime_error(
            "cannot extend a singular bordered Gram state");
      }
      if (!current.adjugate_ready) {
        throw std::runtime_error(
            "cannot extend an unmaterialized bordered Gram state");
      }
    }
    words_[static_cast<std::size_t>(size_)] = word;
    Frame& next = frames_[static_cast<std::size_t>(size_ + 1)];
    if (size_ == 0) {
      next.determinant = vector_order_;
      next.adjugate[0] = 1;
      next.adjugate_ready = true;
      ++size_;
      return next.determinant;
    }

    const Frame& current =
        frames_[static_cast<std::size_t>(size_)];
    std::array<BigInt, kOrder> border{};
    for (int row = 0; row < size_; ++row) {
      const int distance = std::popcount(
          words_[static_cast<std::size_t>(row)] ^ word);
      border[static_cast<std::size_t>(row)] =
          vector_order_ - 2 * distance;
      next.product[static_cast<std::size_t>(row)] = 0;
    }
    for (int row = 0; row < size_; ++row) {
      for (int column = 0; column < size_; ++column) {
        next.product[static_cast<std::size_t>(row)] +=
            current.at(row, column) *
            border[static_cast<std::size_t>(column)];
      }
    }
    BigInt quadratic = 0;
    for (int row = 0; row < size_; ++row) {
      quadratic +=
          border[static_cast<std::size_t>(row)] *
          next.product[static_cast<std::size_t>(row)];
    }
    next.determinant =
        BigInt(vector_order_) * current.determinant - quadratic;
    if (next.determinant < 0) {
      throw std::runtime_error(
          "bordered Gram determinant was negative");
    }
    next.adjugate_ready = false;
    ++size_;
    return next.determinant;
  }

  void materialize_adjugate() {
    if (size_ == 0) {
      throw std::runtime_error(
          "cannot materialize an empty bordered Gram state");
    }
    Frame& current =
        frames_[static_cast<std::size_t>(size_)];
    if (current.adjugate_ready) return;
    const Frame& parent =
        frames_[static_cast<std::size_t>(size_ - 1)];
    const int parent_size = size_ - 1;
    for (int row = 0; row < parent_size; ++row) {
      for (int column = row;
           column < parent_size; ++column) {
        const BigInt numerator =
            current.determinant * parent.at(row, column) +
            current.product[static_cast<std::size_t>(row)] *
                current.product[static_cast<std::size_t>(column)];
        if (numerator % parent.determinant != 0) {
          throw std::runtime_error(
              "non-exact bordered adjugate division");
        }
        const BigInt value =
            numerator / parent.determinant;
        current.at(row, column) = value;
        current.at(column, row) = value;
      }
      current.at(row, parent_size) =
          -current.product[static_cast<std::size_t>(row)];
      current.at(parent_size, row) =
          -current.product[static_cast<std::size_t>(row)];
    }
    current.at(parent_size, parent_size) =
        parent.determinant;
    current.adjugate_ready = true;
  }

  void pop() {
    if (size_ == 0) {
      throw std::runtime_error(
          "cannot pop an empty bordered Gram state");
    }
    --size_;
  }

  int size() const { return size_; }

  const BigInt& determinant() const {
    if (size_ == 0) return empty_determinant_;
    return frames_[static_cast<std::size_t>(size_)].determinant;
  }

  const BigInt& adjugate(int row, int column) const {
    if (
        row < 0 || column < 0 || row >= size_ ||
        column >= size_) {
      throw std::runtime_error(
          "bordered Gram adjugate index out of range");
    }
    const Frame& current =
        frames_[static_cast<std::size_t>(size_)];
    if (!current.adjugate_ready) {
      throw std::runtime_error(
          "bordered Gram adjugate is not materialized");
    }
    return current.at(row, column);
  }

 private:
  struct Frame {
    BigInt determinant = 1;
    std::array<BigInt, kOrder * kOrder> adjugate{};
    std::array<BigInt, kOrder> product{};
    bool adjugate_ready = true;

    BigInt& at(int row, int column) {
      return adjugate[
          static_cast<std::size_t>(row * kOrder + column)];
    }

    const BigInt& at(int row, int column) const {
      return adjugate[
          static_cast<std::size_t>(row * kOrder + column)];
    }
  };

  int vector_order_ = 0;
  int size_ = 0;
  std::array<Word, kOrder> words_{};
  std::array<Frame, kOrder + 1> frames_{};
  const BigInt empty_determinant_ = 1;
};

// Let A be the Gram matrix of the fixed rows, D = det(A), and
// J = adj(A).  For candidate rows x_i with fixed-row correlation vectors
// b_i, the integer Schur kernel is
//
//   N_ij = D <x_i,x_j> - b_i^T J b_j.
//
// A clique branch can use only diagonal entries and compatibility-graph
// edges.  If g is the gcd of D and exactly those retained N entries, K=N/g
// and d=D/g keep all later determinant comparisons integral while
// substantially reducing operands.  Nonedges are deliberately absent.
class FixedSchurKernel {
 public:
  template <typename Compatible>
  FixedSchurKernel(
      int vector_order, const std::vector<Word>& fixed,
      const std::vector<Word>& candidates,
      std::uint64_t expected_edges, Compatible compatible)
      : vector_order_(vector_order),
        vertices_(candidates.size()) {
    if (
        vector_order_ < 1 || vector_order_ > kOrder ||
        fixed.empty() ||
        fixed.size() >= static_cast<std::size_t>(vector_order_)) {
      throw std::runtime_error(
          "invalid fixed-row dimensions for Schur kernel");
    }
    if (
        vertices_ >
        static_cast<std::size_t>(
            std::numeric_limits<std::uint32_t>::max())) {
      throw std::runtime_error(
          "Schur-kernel vertex index exceeds uint32_t");
    }
    const std::uint64_t possible_edges =
        checked_possible_edges(vertices_);
    if (expected_edges > possible_edges) {
      throw std::runtime_error(
          "Schur-kernel edge count exceeds the possible pairs");
    }
    if (
        expected_edges >
        static_cast<std::uint64_t>(
            std::numeric_limits<std::size_t>::max())) {
      throw std::runtime_error(
          "Schur-kernel edge count exceeds size_t");
    }
    edge_count_ = static_cast<std::size_t>(expected_edges);
    blocks_ = checked_blocks(vertices_);
    const std::size_t bitmap_entries =
        checked_product(
            vertices_, blocks_,
            "Schur-kernel edge bitmap");
    const std::size_t offset_entries =
        checked_add(
            vertices_, std::size_t{1},
            "Schur-kernel row offsets");
    const std::size_t retained_entries =
        checked_add(
            vertices_, edge_count_,
            "Schur-kernel retained entries");
    std::size_t planned_bytes = checked_product(
        retained_entries, sizeof(BigInt),
        "Schur-kernel integer storage");
    planned_bytes = checked_add(
        planned_bytes,
        checked_product(
            bitmap_entries, sizeof(std::uint64_t),
            "Schur-kernel edge bitmap bytes"),
        "Schur-kernel storage");
    planned_bytes = checked_add(
        planned_bytes,
        checked_product(
            bitmap_entries, sizeof(std::uint32_t),
            "Schur-kernel edge-rank bytes"),
        "Schur-kernel storage");
    planned_bytes = checked_add(
        planned_bytes,
        checked_product(
            offset_entries, sizeof(std::size_t),
            "Schur-kernel row-offset bytes"),
        "Schur-kernel storage");
    const std::size_t workspace_entries =
        checked_product(
            vertices_, fixed.size(),
            "Schur-kernel projection workspace");
    planned_bytes = checked_add(
        planned_bytes,
        checked_product(
            workspace_entries,
            sizeof(int) + sizeof(BigInt),
            "Schur-kernel projection-workspace bytes"),
        "Schur-kernel storage");
    if (planned_bytes > kMaximumSchurKernelBytes) {
      throw std::runtime_error(
          "exact sparse Schur kernel would exceed its checked "
          "128 MiB memory safety limit; no candidate pool or "
          "compatibility edge was truncated");
    }
    planned_storage_bytes_ = planned_bytes;
    diagonal_.resize(vertices_);
    edges_.reserve(edge_count_);
    upper_edges_.assign(bitmap_entries, UINT64_C(0));
    rank_before_block_.resize(bitmap_entries);
    edge_offsets_.resize(offset_entries);

    BorderedGramState fixed_state(vector_order_);
    for (const Word word : fixed) {
      if (fixed_state.push(word) == 0) {
        throw std::runtime_error(
            "fixed rows have singular Gram matrix");
      }
    }
    determinant_ = fixed_state.determinant();
    if (determinant_ <= 0) {
      throw std::runtime_error(
          "fixed Gram determinant must be positive");
    }

    const std::size_t fixed_count = fixed.size();
    std::vector<std::vector<int>> correlations(
        vertices_, std::vector<int>(fixed_count));
    std::vector<std::vector<BigInt>> adjugate_products(
        vertices_, std::vector<BigInt>(fixed_count));
    for (std::size_t vertex = 0; vertex < vertices_; ++vertex) {
      for (std::size_t row = 0; row < fixed_count; ++row) {
        correlations[vertex][row] =
            vector_order_ -
            2 * std::popcount(candidates[vertex] ^ fixed[row]);
      }
      for (std::size_t row = 0; row < fixed_count; ++row) {
        BigInt product = 0;
        for (std::size_t column = 0;
             column < fixed_count; ++column) {
          product +=
              fixed_state.adjugate(
                  static_cast<int>(row),
                  static_cast<int>(column)) *
              correlations[vertex][column];
        }
        adjugate_products[vertex][row] = std::move(product);
      }
    }

    gcd_ = determinant_;
    for (std::size_t first = 0; first < vertices_; ++first) {
      BigInt diagonal = determinant_ * vector_order_;
      for (std::size_t row = 0; row < fixed_count; ++row) {
        diagonal -=
            correlations[first][row] *
            adjugate_products[first][row];
      }
      if (diagonal < 0) {
        throw std::runtime_error(
            "Schur-kernel diagonal was negative");
      }
      diagonal_[first] = diagonal;
      gcd_ = exact_gcd(gcd_, diagonal);
      edge_offsets_[first] = edges_.size();
      for (std::size_t second = first + 1;
           second < vertices_; ++second) {
        if (!compatible(first, second)) continue;
        const int correlation =
            vector_order_ -
            2 * std::popcount(
                    candidates[first] ^ candidates[second]);
        BigInt value = determinant_ * correlation;
        for (std::size_t row = 0; row < fixed_count; ++row) {
          value -=
              correlations[first][row] *
              adjugate_products[second][row];
        }
        const std::size_t bitmap_offset =
            first * blocks_ + second / 64U;
        const std::uint64_t bit =
            UINT64_C(1) << (second % 64U);
        if ((upper_edges_[bitmap_offset] & bit) != 0U) {
          throw std::runtime_error(
              "duplicate Schur-kernel compatibility edge");
        }
        if (edges_.size() >= edge_count_) {
          throw std::runtime_error(
              "Schur-kernel predicate produced more edges than "
              "the exact graph declared");
        }
        upper_edges_[bitmap_offset] |= bit;
        edges_.push_back(std::move(value));
        gcd_ = exact_gcd(gcd_, edges_.back());
      }
    }
    edge_offsets_[vertices_] = edges_.size();
    if (edges_.size() != edge_count_) {
      throw std::runtime_error(
          "Schur-kernel compatibility-edge count disagrees "
          "with the exact graph");
    }
    for (std::size_t vertex = 0;
         vertex < vertices_; ++vertex) {
      std::uint64_t rank = 0;
      for (std::size_t block = 0; block < blocks_; ++block) {
        if (
            rank >
            static_cast<std::uint64_t>(
                std::numeric_limits<std::uint32_t>::max())) {
          throw std::runtime_error(
              "Schur-kernel edge rank exceeds uint32_t");
        }
        const std::size_t bitmap_offset =
            vertex * blocks_ + block;
        rank_before_block_[bitmap_offset] =
            static_cast<std::uint32_t>(rank);
        rank += static_cast<std::uint64_t>(
            std::popcount(upper_edges_[bitmap_offset]));
      }
      const std::size_t retained_in_row =
          edge_offsets_[vertex + 1] -
          edge_offsets_[vertex];
      if (rank != retained_in_row) {
        throw std::runtime_error(
            "Schur-kernel edge-rank index is inconsistent");
      }
    }
    if (gcd_ <= 0 || determinant_ % gcd_ != 0) {
      throw std::runtime_error(
          "invalid Schur-kernel gcd");
    }
    reduced_fixed_determinant_ = determinant_ / gcd_;
    for (BigInt& value : diagonal_) {
      if (value % gcd_ != 0) {
        throw std::runtime_error(
            "Schur-kernel gcd did not divide a diagonal");
      }
      value /= gcd_;
    }
    for (BigInt& value : edges_) {
      if (value % gcd_ != 0) {
        throw std::runtime_error(
            "Schur-kernel gcd did not divide an edge");
      }
      value /= gcd_;
    }
  }

  std::size_t size() const { return vertices_; }
  std::size_t stored_edges() const { return edge_count_; }
  std::size_t planned_storage_bytes() const {
    return planned_storage_bytes_;
  }

  const BigInt& at(
      std::size_t first, std::size_t second) const {
    if (first >= vertices_ || second >= vertices_) {
      throw std::runtime_error(
          "Schur-kernel vertex out of range");
    }
    if (first == second) return diagonal_[first];
    if (first > second) std::swap(first, second);
    const std::size_t block = second / 64U;
    const std::size_t bitmap_offset =
        first * blocks_ + block;
    const std::uint64_t word = upper_edges_[bitmap_offset];
    const std::uint64_t edge_bit =
        UINT64_C(1) << (second % 64U);
    if ((word & edge_bit) == 0U) {
      throw std::runtime_error(
          "Schur-kernel off-diagonal access is not a "
          "compatibility edge");
    }
    const unsigned bit_index =
        static_cast<unsigned>(second % 64U);
    const std::uint64_t lower_mask =
        bit_index == 0U
        ? UINT64_C(0)
        : (UINT64_C(1) << bit_index) - 1U;
    const std::size_t rank =
        static_cast<std::size_t>(
            rank_before_block_[bitmap_offset]) +
        static_cast<std::size_t>(
            std::popcount(word & lower_mask));
    const std::size_t edge_index =
        edge_offsets_[first] + rank;
    if (edge_index >= edge_offsets_[first + 1]) {
      throw std::runtime_error(
          "Schur-kernel edge rank escaped its row");
    }
    return edges_[edge_index];
  }

  const BigInt& determinant() const { return determinant_; }
  const BigInt& gcd() const { return gcd_; }
  const BigInt& reduced_fixed_determinant() const {
    return reduced_fixed_determinant_;
  }

 private:
  static std::size_t checked_add(
      std::size_t first, std::size_t second,
      std::string_view context) {
    if (
        first >
        std::numeric_limits<std::size_t>::max() - second) {
      throw std::runtime_error(
          std::string(context) + " size overflow");
    }
    return first + second;
  }

  static std::size_t checked_product(
      std::size_t first, std::size_t second,
      std::string_view context) {
    if (
        second != 0 &&
        first >
            std::numeric_limits<std::size_t>::max() / second) {
      throw std::runtime_error(
          std::string(context) + " size overflow");
    }
    return first * second;
  }

  static std::size_t checked_blocks(std::size_t vertices) {
    return checked_add(
               vertices, std::size_t{63},
               "Schur-kernel bitmap block count") /
        64U;
  }

  static std::uint64_t checked_possible_edges(
      std::size_t vertices) {
    if (
        vertices >
        static_cast<std::size_t>(
            std::numeric_limits<std::uint64_t>::max())) {
      throw std::runtime_error(
          "Schur-kernel vertex count exceeds uint64_t");
    }
    const std::uint64_t count =
        static_cast<std::uint64_t>(vertices);
    if (
        count != 0 &&
        count - 1 >
            std::numeric_limits<std::uint64_t>::max() / count) {
      throw std::runtime_error(
          "Schur-kernel possible-edge count overflow");
    }
    return count * (count - 1) / 2;
  }

  int vector_order_ = 0;
  std::size_t vertices_ = 0;
  std::size_t blocks_ = 0;
  std::size_t edge_count_ = 0;
  std::size_t planned_storage_bytes_ = 0;
  BigInt determinant_ = 0;
  BigInt gcd_ = 0;
  BigInt reduced_fixed_determinant_ = 0;
  std::vector<BigInt> diagonal_;
  std::vector<BigInt> edges_;
  std::vector<std::uint64_t> upper_edges_;
  std::vector<std::uint32_t> rank_before_block_;
  std::vector<std::size_t> edge_offsets_;
};

// Lazy exact determinant/adjugate updates for principal submatrices of a
// fixed symmetric Schur kernel.  The determinant is prepared before a bound
// check; the adjugate is materialized only when the DFS actually descends.
class BorderedKernelState {
 public:
  explicit BorderedKernelState(const FixedSchurKernel& kernel)
      : kernel_(kernel) {}

  const BigInt& push(std::size_t vertex) {
    static_cast<void>(prepare_push(vertex));
    materialize_adjugate();
    return determinant();
  }

  const BigInt& prepare_push(std::size_t vertex) {
    if (
        vertex >= kernel_.size() ||
        size_ >= kOrder) {
      throw std::runtime_error(
          "invalid bordered-kernel push");
    }
    if (size_ != 0) {
      const Frame& current =
          frames_[static_cast<std::size_t>(size_)];
      if (current.determinant == 0) {
        throw std::runtime_error(
            "cannot extend a singular bordered kernel");
      }
      if (!current.adjugate_ready) {
        throw std::runtime_error(
            "cannot extend an unmaterialized bordered kernel");
      }
    }
    vertices_[static_cast<std::size_t>(size_)] = vertex;
    Frame& next = frames_[static_cast<std::size_t>(size_ + 1)];
    if (size_ == 0) {
      next.determinant = kernel_.at(vertex, vertex);
      if (next.determinant < 0) {
        throw std::runtime_error(
            "bordered-kernel determinant was negative");
      }
      next.adjugate[0] = 1;
      next.adjugate_ready = true;
      ++size_;
      return next.determinant;
    }

    const Frame& current =
        frames_[static_cast<std::size_t>(size_)];
    for (int column = 0; column < size_; ++column) {
      next.border[static_cast<std::size_t>(column)] =
          kernel_.at(
              vertices_[static_cast<std::size_t>(column)],
              vertex);
    }
    for (int row = 0; row < size_; ++row) {
      next.product[static_cast<std::size_t>(row)] = 0;
    }
    for (int row = 0; row < size_; ++row) {
      for (int column = 0; column < size_; ++column) {
        next.product[static_cast<std::size_t>(row)] +=
            current.at(row, column) *
            next.border[static_cast<std::size_t>(column)];
      }
    }
    BigInt quadratic = 0;
    for (int row = 0; row < size_; ++row) {
      quadratic +=
          next.border[static_cast<std::size_t>(row)] *
          next.product[static_cast<std::size_t>(row)];
    }
    next.determinant =
        kernel_.at(vertex, vertex) * current.determinant -
        quadratic;
    if (next.determinant < 0) {
      throw std::runtime_error(
          "bordered-kernel determinant was negative");
    }
    next.adjugate_ready = false;
    ++size_;
    return next.determinant;
  }

  void materialize_adjugate() {
    if (size_ == 0) {
      throw std::runtime_error(
          "cannot materialize an empty bordered kernel");
    }
    Frame& current =
        frames_[static_cast<std::size_t>(size_)];
    if (current.adjugate_ready) return;
    const Frame& parent =
        frames_[static_cast<std::size_t>(size_ - 1)];
    const int parent_size = size_ - 1;
    for (int row = 0; row < parent_size; ++row) {
      for (int column = row;
           column < parent_size; ++column) {
        const BigInt numerator =
            current.determinant * parent.at(row, column) +
            current.product[static_cast<std::size_t>(row)] *
                current.product[static_cast<std::size_t>(column)];
        if (numerator % parent.determinant != 0) {
          throw std::runtime_error(
              "non-exact bordered-kernel adjugate division");
        }
        const BigInt value =
            numerator / parent.determinant;
        current.at(row, column) = value;
        current.at(column, row) = value;
      }
      current.at(row, parent_size) =
          -current.product[static_cast<std::size_t>(row)];
      current.at(parent_size, row) =
          -current.product[static_cast<std::size_t>(row)];
    }
    current.at(parent_size, parent_size) =
        parent.determinant;
    current.adjugate_ready = true;
  }

  void pop() {
    if (size_ == 0) {
      throw std::runtime_error(
          "cannot pop an empty bordered kernel");
    }
    --size_;
  }

  int size() const { return size_; }

  const BigInt& determinant() const {
    if (size_ == 0) return empty_determinant_;
    return frames_[static_cast<std::size_t>(size_)].determinant;
  }

  const BigInt& adjugate(int row, int column) const {
    if (
        row < 0 || column < 0 || row >= size_ ||
        column >= size_) {
      throw std::runtime_error(
          "bordered-kernel adjugate index out of range");
    }
    const Frame& current =
        frames_[static_cast<std::size_t>(size_)];
    if (!current.adjugate_ready) {
      throw std::runtime_error(
          "bordered-kernel adjugate is not materialized");
    }
    return current.at(row, column);
  }

 private:
  struct Frame {
    BigInt determinant = 1;
    std::array<BigInt, kOrder * kOrder> adjugate{};
    std::array<BigInt, kOrder> border{};
    std::array<BigInt, kOrder> product{};
    bool adjugate_ready = true;

    BigInt& at(int row, int column) {
      return adjugate[
          static_cast<std::size_t>(row * kOrder + column)];
    }

    const BigInt& at(int row, int column) const {
      return adjugate[
          static_cast<std::size_t>(row * kOrder + column)];
    }
  };

  const FixedSchurKernel& kernel_;
  int size_ = 0;
  std::array<std::size_t, kOrder> vertices_{};
  std::array<Frame, kOrder + 1> frames_{};
  const BigInt empty_determinant_ = 1;
};

BigInt principal_kernel_determinant(
    const FixedSchurKernel& kernel,
    const std::vector<std::size_t>& vertices) {
  std::vector<std::vector<BigInt>> principal(
      vertices.size(),
      std::vector<BigInt>(vertices.size()));
  for (std::size_t row = 0; row < vertices.size(); ++row) {
    for (std::size_t column = 0;
         column < vertices.size(); ++column) {
      principal[row][column] =
          kernel.at(vertices[row], vertices[column]);
    }
  }
  const BigInt determinant = bareiss(std::move(principal));
  if (determinant < 0) {
    throw std::runtime_error(
        "Schur-kernel principal determinant was negative");
  }
  return determinant;
}

BigInt largest_schur_diagonal_product(
    const FixedSchurKernel& kernel,
    const std::vector<std::size_t>& candidates,
    int count) {
  if (
      count < 0 ||
      candidates.size() < static_cast<std::size_t>(count)) {
      throw std::runtime_error(
        "invalid Schur diagonal-product count");
  }
  if (count == 0) return 1;
  std::vector<const BigInt*> largest;
  largest.reserve(static_cast<std::size_t>(count));
  for (const std::size_t vertex : candidates) {
    const BigInt& diagonal = kernel.at(vertex, vertex);
    if (diagonal < 0) {
      throw std::runtime_error(
          "negative Schur-kernel diagonal");
    }
    const auto position = std::lower_bound(
        largest.begin(), largest.end(), &diagonal,
        [](const BigInt* first, const BigInt* second) {
          return *first > *second;
        });
    largest.insert(position, &diagonal);
    if (largest.size() > static_cast<std::size_t>(count)) {
      largest.pop_back();
    }
  }
  BigInt result = 1;
  for (const BigInt* diagonal : largest) {
    result *= *diagonal;
  }
  return result;
}

BigInt integer_power(int base, int exponent) {
  BigInt result = 1;
  BigInt factor = base;
  while (exponent > 0) {
    if ((exponent & 1) != 0) result *= factor;
    exponent >>= 1;
    if (exponent != 0) factor *= factor;
  }
  return result;
}

BigInt big_power(BigInt base, int exponent) {
  BigInt result = 1;
  while (exponent > 0) {
    if ((exponent & 1) != 0) result *= base;
    exponent >>= 1;
    if (exponent != 0) base *= base;
  }
  return result;
}

int positive_mod4(int value) {
  const int remainder = value % 4;
  return remainder < 0 ? remainder + 4 : remainder;
}

RowWords matrix_to_words(const SignMatrix& matrix) {
  RowWords rows{};
  for (int row = 0; row < kOrder; ++row) {
    Word word = 0;
    for (int column = 0; column < kOrder; ++column) {
      if (matrix[row][column] == -1) {
        word |= Word{1} << column;
      }
    }
    rows[static_cast<std::size_t>(row)] = word;
  }
  return rows;
}

RowWords normalize_hamming_orientation(const SignMatrix& input) {
  SignMatrix matrix = input;
  for (int row = 1; row < kOrder; ++row) {
    int inner_product = 0;
    for (int column = 0; column < kOrder; ++column) {
      inner_product +=
          matrix[row][column] * matrix[0][column];
    }
    const int residue = positive_mod4(inner_product);
    if (residue == 1) {
      for (int column = 0; column < kOrder; ++column) {
        matrix[row][column] = -matrix[row][column];
      }
    } else if (residue != 3) {
      throw std::runtime_error(
          "odd-order row inner product was not odd");
    }
  }
  for (int column = 0; column < kOrder; ++column) {
    const int column_switch = matrix[0][column];
    for (int row = 0; row < kOrder; ++row) {
      matrix[row][column] *= column_switch;
    }
  }
  for (int column = 0; column < kOrder; ++column) {
    if (matrix[0][column] != 1) {
      throw std::runtime_error(
          "normalization failed to anchor row 1");
    }
  }
  for (int first = 0; first < kOrder; ++first) {
    for (int second = first + 1; second < kOrder; ++second) {
      int inner_product = 0;
      for (int column = 0; column < kOrder; ++column) {
        inner_product +=
            matrix[first][column] * matrix[second][column];
      }
      if (positive_mod4(inner_product) != 3) {
        throw std::runtime_error(
            "normalization failed the Gram residue invariant");
      }
    }
  }
  return matrix_to_words(matrix);
}

std::string row_words_text(const RowWords& rows) {
  std::ostringstream output;
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      if (column != 0) output << ' ';
      output
          << ((((rows[static_cast<std::size_t>(row)] >> column) &
                 Word{1}) != 0U)
                  ? -1
                  : 1);
    }
    output << '\n';
  }
  return output.str();
}

bool distance_allowed(int distance, bool allow_distance_14) {
  return distance == 10 || distance == 12 ||
      (allow_distance_14 && distance == 14);
}

void validate_hamming_code(
    const RowWords& rows, bool allow_distance_14) {
  if (rows[0] != 0U) {
    throw std::runtime_error(
        "normalized anchor row is not all +1");
  }
  for (int first = 0; first < kOrder; ++first) {
    for (int second = first + 1; second < kOrder; ++second) {
      const int distance = std::popcount(
          rows[static_cast<std::size_t>(first)] ^
          rows[static_cast<std::size_t>(second)]);
      if (!distance_allowed(distance, allow_distance_14)) {
        throw std::runtime_error(
            "normalized seed is outside the selected Hamming code: "
            "rows " +
            std::to_string(first + 1) + " and " +
            std::to_string(second + 1) + " have distance " +
            std::to_string(distance) +
            (allow_distance_14
                 ? " (allowed: 10,12,14)"
                 : " (allowed: 10,12; try --allow-distance-14 only "
                   "for the separate {-5 Gram} arm)"));
      }
    }
  }
}

fs::path normalized_absolute_path(const fs::path& path) {
  return fs::weakly_canonical(fs::absolute(path));
}

bool path_entry_exists(const fs::path& path) {
  std::error_code error;
  const fs::file_status status = fs::symlink_status(path, error);
  if (error == std::errc::no_such_file_or_directory) return false;
  if (error) {
    throw std::runtime_error(
        "cannot inspect path " + path.string() + ": " +
        error.message());
  }
  return status.type() != fs::file_type::not_found;
}

void sync_directory(const fs::path& directory) {
  const int descriptor =
      ::open(directory.c_str(), O_RDONLY | O_DIRECTORY);
  if (descriptor < 0) {
    throw std::runtime_error(
        "cannot open output directory for sync: " +
        std::string(std::strerror(errno)));
  }
  if (::fsync(descriptor) != 0) {
    const int saved = errno;
    ::close(descriptor);
    throw std::runtime_error(
        "cannot sync output directory: " +
        std::string(std::strerror(saved)));
  }
  if (::close(descriptor) != 0) {
    throw std::runtime_error(
        "cannot close output directory: " +
        std::string(std::strerror(errno)));
  }
}

class AtomicArtifact {
 public:
  explicit AtomicArtifact(fs::path path) : path_(std::move(path)) {}

  AtomicArtifact(const AtomicArtifact&) = delete;
  AtomicArtifact& operator=(const AtomicArtifact&) = delete;

  void write(std::string_view content) {
    const fs::path directory =
        path_.parent_path().empty() ? fs::path(".")
                                    : path_.parent_path();
    fs::create_directories(directory);
    const fs::path temporary =
        directory /
        (path_.filename().string() + ".tmp." +
         std::to_string(static_cast<unsigned long>(::getpid())) +
         "." + std::to_string(++nonce_));
    int descriptor = -1;
    bool installed_or_renamed = false;
    try {
      descriptor = ::open(
          temporary.c_str(), O_WRONLY | O_CREAT | O_EXCL,
          S_IRUSR | S_IWUSR);
      if (descriptor < 0) {
        throw std::runtime_error(
            "cannot create atomic temporary file: " +
            std::string(std::strerror(errno)));
      }
      std::size_t offset = 0;
      while (offset < content.size()) {
        const ssize_t written = ::write(
            descriptor, content.data() + offset,
            content.size() - offset);
        if (written < 0) {
          if (errno == EINTR) continue;
          throw std::runtime_error(
              "cannot write atomic temporary file: " +
              std::string(std::strerror(errno)));
        }
        offset += static_cast<std::size_t>(written);
      }
      if (::fsync(descriptor) != 0) {
        throw std::runtime_error(
            "cannot sync atomic temporary file: " +
            std::string(std::strerror(errno)));
      }
      if (::close(descriptor) != 0) {
        descriptor = -1;
        throw std::runtime_error(
            "cannot close atomic temporary file: " +
            std::string(std::strerror(errno)));
      }
      descriptor = -1;
      if (!installed_) {
        if (::link(temporary.c_str(), path_.c_str()) != 0) {
          throw std::runtime_error(
              "refusing to overwrite artifact " + path_.string() +
              ": " + std::string(std::strerror(errno)));
        }
        if (::unlink(temporary.c_str()) != 0) {
          throw std::runtime_error(
              "cannot remove linked temporary file: " +
              std::string(std::strerror(errno)));
        }
        installed_ = true;
      } else {
        if (::rename(temporary.c_str(), path_.c_str()) != 0) {
          throw std::runtime_error(
              "cannot atomically replace owned artifact: " +
              std::string(std::strerror(errno)));
        }
      }
      installed_or_renamed = true;
      sync_directory(directory);
    } catch (...) {
      if (descriptor >= 0) ::close(descriptor);
      if (!installed_or_renamed) {
        std::error_code ignored;
        fs::remove(temporary, ignored);
      }
      throw;
    }
  }

  const fs::path& path() const { return path_; }

 private:
  fs::path path_;
  std::uint64_t nonce_ = 0;
  bool installed_ = false;
};

void preflight_paths(const Arguments& arguments) {
  const std::array<fs::path, 4> targets{
      arguments.output, arguments.log, arguments.summary,
      arguments.start};
  std::array<fs::path, 4> normalized{};
  for (std::size_t index = 0; index < targets.size(); ++index) {
    normalized[index] = normalized_absolute_path(targets[index]);
  }
  for (std::size_t first = 0; first < normalized.size(); ++first) {
    for (std::size_t second = first + 1;
         second < normalized.size(); ++second) {
      if (normalized[first] == normalized[second]) {
        throw std::runtime_error(
            "input and artifact paths must be distinct");
      }
    }
  }
  for (std::size_t index = 0; index < 3; ++index) {
    if (path_entry_exists(targets[index])) {
      throw std::runtime_error(
          "refusing to overwrite existing artifact: " +
          targets[index].string());
    }
  }
}

struct PoolResult {
  std::vector<Word> words;
  std::uint64_t examined = 0;
  bool complete = true;
};

PoolResult enumerate_compatible_words(
    int order, const std::vector<Word>& fixed,
    const std::function<bool(int)>& allowed,
    std::uint64_t pool_limit,
    const std::function<bool()>& should_stop) {
  if (order < 1 || order > 30) {
    throw std::runtime_error(
        "word enumeration supports orders in [1,30]");
  }
  const Word word_count = Word{1} << order;
  PoolResult result;
  result.words.reserve(
      static_cast<std::size_t>(
          std::min<std::uint64_t>(pool_limit, UINT64_C(1024))));
  for (Word word = 0; word < word_count; ++word) {
    ++result.examined;
    bool compatible = true;
    for (const Word fixed_word : fixed) {
      if (!allowed(std::popcount(word ^ fixed_word))) {
        compatible = false;
        break;
      }
    }
    if (compatible) {
      result.words.push_back(word);
      if (
          result.words.size() >
          static_cast<std::size_t>(pool_limit)) {
        throw std::runtime_error(
            "exact candidate pool exceeded --pool-limit " +
            std::to_string(pool_limit) +
            "; no truncated graph was searched");
      }
    }
    if (
        (static_cast<std::uint64_t>(word) &
         (kPoolPollStride - 1U)) == 0U &&
        should_stop()) {
      result.complete = false;
      return result;
    }
  }
  return result;
}

class BitGraph {
 public:
  explicit BitGraph(std::size_t vertices)
      : vertices_(vertices),
        blocks_((vertices + 63U) / 64U),
        adjacency_(checked_word_count(vertices, blocks_), 0),
        degrees_(vertices, 0) {}

  void add_edge(std::size_t first, std::size_t second) {
    if (first == second || first >= vertices_ || second >= vertices_) {
      throw std::runtime_error("invalid bit-graph edge");
    }
    const std::size_t first_offset = first * blocks_;
    const std::size_t second_offset = second * blocks_;
    adjacency_[first_offset + second / 64U] |=
        UINT64_C(1) << (second % 64U);
    adjacency_[second_offset + first / 64U] |=
        UINT64_C(1) << (first % 64U);
    ++degrees_[first];
    ++degrees_[second];
  }

  bool adjacent(std::size_t first, std::size_t second) const {
    if (first >= vertices_ || second >= vertices_) {
      throw std::runtime_error("bit-graph vertex out of range");
    }
    return (
        adjacency_[first * blocks_ + second / 64U] >>
        (second % 64U)) &
        UINT64_C(1);
  }

  std::size_t size() const { return vertices_; }
  std::uint32_t degree(std::size_t vertex) const {
    return degrees_.at(vertex);
  }

 private:
  static std::size_t checked_word_count(
      std::size_t vertices, std::size_t blocks) {
    if (
        blocks != 0 &&
        vertices >
            std::numeric_limits<std::size_t>::max() / blocks) {
      throw std::runtime_error(
          "bit-graph allocation size overflow");
    }
    return vertices * blocks;
  }

  std::size_t vertices_ = 0;
  std::size_t blocks_ = 0;
  std::vector<std::uint64_t> adjacency_;
  std::vector<std::uint32_t> degrees_;
};

struct GraphResult {
  BitGraph graph;
  std::uint64_t pair_checks = 0;
  std::uint64_t edges = 0;
  bool complete = true;

  explicit GraphResult(std::size_t vertices) : graph(vertices) {}
};

GraphResult build_compatibility_graph(
    const std::vector<Word>& words,
    const std::function<bool(int)>& allowed,
    const std::function<bool()>& should_stop) {
  GraphResult result(words.size());
  for (std::size_t first = 0; first < words.size(); ++first) {
    for (std::size_t second = first + 1;
         second < words.size(); ++second) {
      ++result.pair_checks;
      if (allowed(std::popcount(words[first] ^ words[second]))) {
        result.graph.add_edge(first, second);
        ++result.edges;
      }
      if (
          (result.pair_checks & (kGraphPollStride - 1U)) == 0U &&
          should_stop()) {
        result.complete = false;
        return result;
      }
    }
  }
  return result;
}

std::vector<std::vector<std::size_t>>
greedy_independent_color_classes(
    const BitGraph& graph,
    const std::vector<std::size_t>& vertices) {
  std::vector<std::vector<std::size_t>> classes;
  for (const std::size_t vertex : vertices) {
    if (vertex >= graph.size()) {
      throw std::runtime_error(
          "greedy coloring vertex out of range");
    }
    bool placed = false;
    for (std::vector<std::size_t>& color_class : classes) {
      bool independent = true;
      for (const std::size_t member : color_class) {
        if (graph.adjacent(vertex, member)) {
          independent = false;
          break;
        }
      }
      if (independent) {
        color_class.push_back(vertex);
        placed = true;
        break;
      }
    }
    if (!placed) {
      classes.push_back(
          std::vector<std::size_t>{vertex});
    }
  }
  return classes;
}

struct CliqueWalkStatistics {
  std::uint64_t nodes = 0;
  std::uint64_t cardinality_prunes = 0;
  std::uint64_t color_checks = 0;
  std::uint64_t color_prunes = 0;
  std::uint64_t leaves = 0;
  bool complete = true;
};

class CliqueWalker {
 public:
  using Prune =
      std::function<bool(
          const std::vector<std::size_t>&,
          const std::vector<std::size_t>&)>;
  using Visit =
      std::function<void(const std::vector<std::size_t>&)>;
  using Stop = std::function<bool()>;
  using Push = std::function<void(std::size_t)>;
  using Commit = std::function<void()>;
  using Pop = std::function<void()>;

  CliqueWalker(
      const BitGraph& graph, int target, Prune prune,
      Visit visit, Stop stop, Push push, Commit commit, Pop pop)
      : graph_(graph),
        target_(target),
        prune_(std::move(prune)),
        visit_(std::move(visit)),
        stop_(std::move(stop)),
        push_(std::move(push)),
        commit_(std::move(commit)),
        pop_(std::move(pop)) {}

  CliqueWalkStatistics run(std::vector<std::size_t> candidates) {
    if (target_ < 0) {
      throw std::runtime_error("negative clique target");
    }
    chosen_.clear();
    depth_first(candidates, false);
    return statistics_;
  }

 private:
  void depth_first(
      const std::vector<std::size_t>& candidates,
      bool push_chosen_vertex) {
    if (!statistics_.complete) return;
    ++statistics_.nodes;
    if (
        (statistics_.nodes & (kCliquePollStride - 1U)) == 0U &&
        stop_()) {
      statistics_.complete = false;
      return;
    }
    if (
        chosen_.size() + candidates.size() <
        static_cast<std::size_t>(target_)) {
      ++statistics_.cardinality_prunes;
      return;
    }
    const int remaining =
        target_ - static_cast<int>(chosen_.size());
    if (remaining > 0) {
      ++statistics_.color_checks;
      const std::size_t colors =
          greedy_independent_color_classes(
              graph_, candidates)
              .size();
      if (colors < static_cast<std::size_t>(remaining)) {
        ++statistics_.color_prunes;
        return;
      }
    }
    if (push_chosen_vertex) {
      push_(chosen_.back());
    }
    if (prune_(chosen_, candidates)) {
      if (push_chosen_vertex) pop_();
      return;
    }
    if (chosen_.size() == static_cast<std::size_t>(target_)) {
      ++statistics_.leaves;
      visit_(chosen_);
      if (push_chosen_vertex) pop_();
      return;
    }
    if (push_chosen_vertex) commit_();
    for (std::size_t position = 0; position < candidates.size();
         ++position) {
      if (
          chosen_.size() + (candidates.size() - position) <
          static_cast<std::size_t>(target_)) {
        ++statistics_.cardinality_prunes;
        break;
      }
      const std::size_t vertex = candidates[position];
      std::vector<std::size_t> next;
      next.reserve(candidates.size() - position - 1);
      for (std::size_t inner = position + 1;
           inner < candidates.size(); ++inner) {
        if (graph_.adjacent(vertex, candidates[inner])) {
          next.push_back(candidates[inner]);
        }
      }
      chosen_.push_back(vertex);
      depth_first(next, true);
      chosen_.pop_back();
      if (!statistics_.complete) break;
    }
    if (push_chosen_vertex) pop_();
  }

  const BitGraph& graph_;
  int target_;
  Prune prune_;
  Visit visit_;
  Stop stop_;
  Push push_;
  Commit commit_;
  Pop pop_;
  std::vector<std::size_t> chosen_;
  CliqueWalkStatistics statistics_{};
};

std::string integer_array_json(const std::vector<int>& values) {
  std::ostringstream output;
  output << '[';
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index != 0) output << ',';
    output << values[index];
  }
  output << ']';
  return output.str();
}

std::uint64_t binomial_coefficient(int n, int k) {
  if (n < 0 || k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  std::uint64_t result = 1;
  for (int index = 1; index <= k; ++index) {
    const std::uint64_t numerator =
        static_cast<std::uint64_t>(n - k + index);
    if (
        result >
        std::numeric_limits<std::uint64_t>::max() /
            numerator) {
      throw std::runtime_error(
          "binomial coefficient overflow");
    }
    result *= numerator;
    result /= static_cast<std::uint64_t>(index);
  }
  return result;
}

std::uint64_t destroy_mask_rank(
    Word mask, int row_count, int destroy_count) {
  if (
      row_count < 1 || row_count > 30 ||
      destroy_count < 0 || destroy_count > row_count ||
      (mask >> row_count) != 0U ||
      std::popcount(mask) != destroy_count) {
    throw std::runtime_error(
        "invalid destroy mask for ranking");
  }
  std::uint64_t rank = 0;
  int selected = 0;
  for (int bit = 0; bit < row_count; ++bit) {
    if (((mask >> bit) & Word{1}) == 0U) continue;
    ++selected;
    rank += binomial_coefficient(bit, selected);
  }
  return rank;
}

Word destroy_mask_from_rank(
    std::uint64_t rank, int row_count, int destroy_count) {
  const std::uint64_t total =
      binomial_coefficient(row_count, destroy_count);
  if (
      row_count < 1 || row_count > 30 ||
      destroy_count < 0 || destroy_count > row_count ||
      rank >= total) {
    throw std::runtime_error(
        "invalid destroy-mask rank");
  }
  Word mask = 0;
  int maximum_bit = row_count - 1;
  for (int selected = destroy_count;
       selected >= 1; --selected) {
    int bit = maximum_bit;
    while (
        bit >= selected - 1 &&
        binomial_coefficient(bit, selected) > rank) {
      --bit;
    }
    if (bit < selected - 1) {
      throw std::runtime_error(
          "destroy-mask unranking failed");
    }
    mask |= Word{1} << bit;
    rank -= binomial_coefficient(bit, selected);
    maximum_bit = bit - 1;
  }
  if (rank != 0) {
    throw std::runtime_error(
        "destroy-mask unranking left a remainder");
  }
  return mask;
}

std::uint64_t splitmix64(std::uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value =
      (value ^ (value >> 30U)) *
      UINT64_C(0xbf58476d1ce4e5b9);
  value =
      (value ^ (value >> 27U)) *
      UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31U);
}

std::uint64_t deterministic_uniform_index(
    std::mt19937_64& generator, std::uint64_t bound) {
  if (bound == 0) {
    throw std::runtime_error(
        "uniform-index bound must be positive");
  }
  const std::uint64_t threshold =
      (UINT64_C(0) - bound) % bound;
  for (;;) {
    const std::uint64_t value = generator();
    if (value >= threshold) return value % bound;
  }
}

std::uint64_t destroy_mask_order_key(
    std::uint64_t seed, std::uint64_t parent_generation,
    int destroy_count, Word mask) {
  std::uint64_t key = splitmix64(seed);
  key = splitmix64(
      key ^ splitmix64(parent_generation));
  key = splitmix64(
      key ^ splitmix64(
                static_cast<std::uint64_t>(destroy_count)));
  return splitmix64(key ^ static_cast<std::uint64_t>(mask));
}

struct DestroyMaskChoice {
  int destroy_count = 0;
  Word mask = 0;
  std::uint64_t rank = 0;
  std::uint64_t parent_generation = 0;
};

class DestroyMaskSchedule {
 public:
  DestroyMaskSchedule(
      int row_count, int destroy_min, int destroy_max,
      std::uint64_t seed, std::uint64_t shard_count,
      std::uint64_t shard_index)
      : row_count_(row_count),
        destroy_min_(destroy_min),
        destroy_max_(destroy_max),
        seed_(seed),
        shard_count_(shard_count),
        shard_index_(shard_index) {
    if (
        row_count_ < 1 || row_count_ > 30 ||
        destroy_min_ < 1 || destroy_min_ > destroy_max_ ||
        destroy_max_ > row_count_ || shard_count_ == 0 ||
        shard_index_ >= shard_count_) {
      throw std::runtime_error(
          "invalid no-replacement destroy schedule");
    }
    reset(0);
  }

  void reset(std::uint64_t parent_generation) {
    parent_generation_ = parent_generation;
    sizes_.clear();
    sizes_.reserve(static_cast<std::size_t>(
        destroy_max_ - destroy_min_ + 1));
    for (int destroy_count = destroy_min_;
         destroy_count <= destroy_max_; ++destroy_count) {
      SizeSchedule schedule;
      schedule.destroy_count = destroy_count;
      const std::uint64_t total =
          binomial_coefficient(row_count_, destroy_count);
      const std::uint64_t shard_items =
          shard_index_ >= total
          ? 0
          : 1 + (total - 1 - shard_index_) /
                    shard_count_;
      if (
          shard_items >
          std::numeric_limits<std::size_t>::max()) {
        throw std::runtime_error(
            "destroy schedule allocation size overflow");
      }
      schedule.items.reserve(
          static_cast<std::size_t>(shard_items));
      for (std::uint64_t rank = 0; rank < total; ++rank) {
        if (rank % shard_count_ != shard_index_) continue;
        const Word mask = destroy_mask_from_rank(
            rank, row_count_, destroy_count);
        schedule.items.push_back(
            Item{
                mask, rank,
                destroy_mask_order_key(
                    seed_, parent_generation_,
                    destroy_count, mask)});
      }
      std::sort(
          schedule.items.begin(), schedule.items.end(),
          [](const Item& first, const Item& second) {
            if (first.order_key != second.order_key) {
              return first.order_key < second.order_key;
            }
            return first.mask < second.mask;
          });
      sizes_.push_back(std::move(schedule));
    }
  }

  std::optional<DestroyMaskChoice> next(
      std::mt19937_64& generator) {
    std::vector<std::size_t> nonexhausted;
    nonexhausted.reserve(sizes_.size());
    for (std::size_t index = 0; index < sizes_.size(); ++index) {
      if (sizes_[index].cursor < sizes_[index].items.size()) {
        nonexhausted.push_back(index);
      }
    }
    if (nonexhausted.empty()) return std::nullopt;
    const std::uint64_t selected =
        deterministic_uniform_index(
            generator,
            static_cast<std::uint64_t>(
                nonexhausted.size()));
    SizeSchedule& schedule =
        sizes_[nonexhausted[
            static_cast<std::size_t>(selected)]];
    const Item& item =
        schedule.items[schedule.cursor++];
    return DestroyMaskChoice{
        schedule.destroy_count, item.mask, item.rank,
        parent_generation_};
  }

  std::uint64_t remaining() const {
    std::uint64_t result = 0;
    for (const SizeSchedule& schedule : sizes_) {
      result += static_cast<std::uint64_t>(
          schedule.items.size() - schedule.cursor);
    }
    return result;
  }

 private:
  struct Item {
    Word mask = 0;
    std::uint64_t rank = 0;
    std::uint64_t order_key = 0;
  };

  struct SizeSchedule {
    int destroy_count = 0;
    std::vector<Item> items;
    std::size_t cursor = 0;
  };

  int row_count_ = 0;
  int destroy_min_ = 0;
  int destroy_max_ = 0;
  std::uint64_t seed_ = 0;
  std::uint64_t shard_count_ = 1;
  std::uint64_t shard_index_ = 0;
  std::uint64_t parent_generation_ = 0;
  std::vector<SizeSchedule> sizes_;
};

std::vector<int> destroyed_rows_from_mask(
    Word mask, int row_count) {
  if (
      row_count < 1 || row_count > 30 ||
      (mask >> row_count) != 0U) {
    throw std::runtime_error(
        "invalid destroy mask dimensions");
  }
  std::vector<int> rows;
  for (int bit = 0; bit < row_count; ++bit) {
    if (((mask >> bit) & Word{1}) != 0U) {
      rows.push_back(bit + 1);
    }
  }
  return rows;
}

class SearchRunner {
 public:
  explicit SearchRunner(Arguments arguments)
      : arguments_(std::move(arguments)),
        output_(arguments_.output),
        log_(arguments_.log),
        summary_(arguments_.summary),
        generator_(arguments_.seed) {
    if (arguments_.destroy_without_replacement) {
      destroy_schedule_.emplace(
          kNonAnchorRows, arguments_.destroy_min,
          arguments_.destroy_max, arguments_.seed,
          arguments_.destroy_shard_count,
          arguments_.destroy_shard_index);
    }
  }

  int run() {
    preflight_paths(arguments_);
    SignMatrix input = read_sign_matrix(arguments_.start);
    if (arguments_.transpose_start) {
      input = transpose_sign_matrix(input);
    }
    const RowWords input_words = matrix_to_words(input);
    const BigInt input_determinant =
        absolute(exact_sign_determinant(input_words));
    best_rows_ = normalize_hamming_orientation(input);
    best_determinant_ =
        absolute(exact_sign_determinant(best_rows_));
    if (best_determinant_ != input_determinant) {
      throw std::runtime_error(
          "switch normalization changed the exact determinant");
    }
    validate_hamming_code(
        best_rows_, arguments_.allow_distance_14);
    start_determinant_ = best_determinant_;

    started_at_ = Clock::now();
    deadline_ = started_at_ + std::chrono::duration_cast<Clock::duration>(
        std::chrono::duration<double>(arguments_.seconds));
    next_heartbeat_ =
        arguments_.heartbeat_seconds == 0.0
        ? Clock::time_point::max()
        : started_at_ +
              std::chrono::duration_cast<Clock::duration>(
                  std::chrono::duration<double>(
                      arguments_.heartbeat_seconds));

    output_.write(row_words_text(best_rows_));
    record_start_event();
    flush_log();
    write_summary("running", "");

    int result = 0;
    std::string final_status;
    std::string final_reason;
    try {
      search_loop();
      if (stop_reason_ == "signal") {
        final_status = "interrupted";
        final_reason = "SIGINT_or_SIGTERM";
        result = 130;
      } else if (stop_reason_ == "deadline") {
        final_status = "deadline";
        final_reason = "wall_clock_deadline";
      } else if (
          stop_reason_ == "destroy_schedule_exhausted") {
        final_status = "completed";
        final_reason = "destroy_schedule_exhausted";
      } else {
        final_status = "completed";
        final_reason = "iteration_limit";
      }
    } catch (const std::exception& error) {
      final_status = "aborted";
      final_reason = error.what();
      result = 2;
    }

    verify_owned_output();
    record_finish_event(final_status, final_reason);
    flush_log();
    write_summary(final_status, final_reason);
    if (result == 2) {
      std::cerr << "error: " << final_reason << '\n';
    }
    return result;
  }

 private:
  double elapsed_seconds() const {
    return std::chrono::duration<double>(
               Clock::now() - started_at_)
        .count();
  }

  std::string allowed_distances_json() const {
    return arguments_.allow_distance_14
        ? "[10,12,14]" : "[10,12]";
  }

  void record_start_event() {
    std::ostringstream event;
    event << "{\"schema_version\":1,\"event\":\"start\","
          << "\"elapsed_seconds\":" << std::fixed
          << std::setprecision(6) << elapsed_seconds() << ','
          << "\"start\":\""
          << json_escape(arguments_.start.string()) << "\","
          << "\"output\":\""
          << json_escape(arguments_.output.string()) << "\","
          << "\"seed\":" << arguments_.seed << ','
          << "\"seconds\":" << arguments_.seconds << ','
          << "\"destroy_min\":" << arguments_.destroy_min << ','
          << "\"destroy_max\":" << arguments_.destroy_max << ','
          << "\"destroy_without_replacement\":"
          << (arguments_.destroy_without_replacement
                  ? "true" : "false")
          << ','
          << "\"destroy_shard_count\":"
          << arguments_.destroy_shard_count << ','
          << "\"destroy_shard_index\":"
          << arguments_.destroy_shard_index << ','
          << "\"pool_limit\":" << arguments_.pool_limit << ','
          << "\"allowed_distances\":"
          << allowed_distances_json() << ','
          << "\"start_absolute_determinant\":\""
          << big_to_string(start_determinant_) << "\"}";
    events_.push_back(event.str());
  }

  void record_heartbeat() {
    std::ostringstream event;
    event << "{\"schema_version\":1,\"event\":\"heartbeat\","
          << "\"elapsed_seconds\":" << std::fixed
          << std::setprecision(6) << elapsed_seconds() << ','
          << "\"phase\":\"" << json_escape(phase_) << "\","
          << "\"iterations_started\":"
          << statistics_.iterations_started << ','
          << "\"iterations_completed\":"
          << statistics_.iterations_completed << ','
          << "\"pool_words_examined\":"
          << statistics_.pool_words_examined << ','
          << "\"clique_nodes\":" << statistics_.clique_nodes << ','
          << "\"color_checks\":" << statistics_.color_checks << ','
          << "\"color_prunes\":" << statistics_.color_prunes << ','
          << "\"parent_generation\":"
          << parent_generation_ << ','
          << "\"best_absolute_determinant\":\""
          << big_to_string(best_determinant_) << "\"}";
    events_.push_back(event.str());
  }

  void record_iteration_start(
      std::uint64_t iteration, const std::vector<int>& destroyed,
      const DestroyMaskChoice& choice) {
    std::vector<int> one_based = destroyed;
    for (int& row : one_based) ++row;
    std::ostringstream event;
    event << "{\"schema_version\":1,"
          << "\"event\":\"iteration_start\","
          << "\"elapsed_seconds\":" << std::fixed
          << std::setprecision(6) << elapsed_seconds() << ','
          << "\"iteration\":" << iteration << ','
          << "\"destroy_mask\":"
          << static_cast<std::uint64_t>(choice.mask) << ','
          << "\"destroy_rank\":" << choice.rank << ','
          << "\"parent_generation\":"
          << choice.parent_generation << ','
          << "\"destroy_without_replacement\":"
          << (arguments_.destroy_without_replacement
                  ? "true" : "false")
          << ','
          << "\"destroy_shard_count\":"
          << arguments_.destroy_shard_count << ','
          << "\"destroy_shard_index\":"
          << arguments_.destroy_shard_index << ','
          << "\"destroyed_rows\":"
          << integer_array_json(one_based) << '}';
    events_.push_back(event.str());
  }

  void record_pool_event(
      std::uint64_t iteration, std::size_t pool_size,
      std::uint64_t examined) {
    std::ostringstream event;
    event << "{\"schema_version\":1,\"event\":\"pool\","
          << "\"elapsed_seconds\":" << std::fixed
          << std::setprecision(6) << elapsed_seconds() << ','
          << "\"iteration\":" << iteration << ','
          << "\"words_examined\":" << examined << ','
          << "\"candidate_words\":" << pool_size << '}';
    events_.push_back(event.str());
  }

  void record_graph_event(
      std::uint64_t iteration, const GraphResult& graph) {
    std::ostringstream event;
    event << "{\"schema_version\":1,\"event\":\"graph\","
          << "\"elapsed_seconds\":" << std::fixed
          << std::setprecision(6) << elapsed_seconds() << ','
          << "\"iteration\":" << iteration << ','
          << "\"vertices\":" << graph.graph.size() << ','
          << "\"pair_checks\":" << graph.pair_checks << ','
          << "\"edges\":" << graph.edges << '}';
    events_.push_back(event.str());
  }

  void record_improvement(
      std::uint64_t iteration, const BigInt& old_determinant) {
    std::ostringstream event;
    event << "{\"schema_version\":1,\"event\":\"improvement\","
          << "\"elapsed_seconds\":" << std::fixed
          << std::setprecision(6) << elapsed_seconds() << ','
          << "\"iteration\":" << iteration << ','
          << "\"parent_generation\":"
          << parent_generation_ << ','
          << "\"previous_absolute_determinant\":\""
          << big_to_string(old_determinant) << "\","
          << "\"absolute_determinant\":\""
          << big_to_string(best_determinant_) << "\"}";
    events_.push_back(event.str());
  }

  void record_iteration_finish(
      std::uint64_t iteration, const CliqueWalkStatistics& walk,
      std::uint64_t scored_before, std::uint64_t bounds_before,
      std::uint64_t prunes_before,
      std::uint64_t colors_before,
      std::uint64_t color_prunes_before) {
    std::ostringstream event;
    event << "{\"schema_version\":1,"
          << "\"event\":\"iteration_finish\","
          << "\"elapsed_seconds\":" << std::fixed
          << std::setprecision(6) << elapsed_seconds() << ','
          << "\"iteration\":" << iteration << ','
          << "\"clique_nodes\":" << walk.nodes << ','
          << "\"clique_leaves\":" << walk.leaves << ','
          << "\"color_checks\":"
          << (statistics_.color_checks - colors_before) << ','
          << "\"color_prunes\":"
          << (statistics_.color_prunes -
              color_prunes_before)
          << ','
          << "\"bound_checks\":"
          << (statistics_.bound_checks - bounds_before) << ','
          << "\"bound_prunes\":"
          << (statistics_.bound_prunes - prunes_before) << ','
          << "\"exact_matrix_scores\":"
          << (statistics_.exact_matrix_scores - scored_before)
          << ",\"best_absolute_determinant\":\""
          << big_to_string(best_determinant_) << "\"}";
    events_.push_back(event.str());
  }

  void record_finish_event(
      std::string_view status, std::string_view reason) {
    std::ostringstream event;
    event << "{\"schema_version\":1,\"event\":\"finish\","
          << "\"elapsed_seconds\":" << std::fixed
          << std::setprecision(6) << elapsed_seconds() << ','
          << "\"status\":\"" << json_escape(status) << "\","
          << "\"reason\":\"" << json_escape(reason) << "\","
          << "\"iterations_completed\":"
          << statistics_.iterations_completed << ','
          << "\"best_absolute_determinant\":\""
          << big_to_string(best_determinant_) << "\"}";
    events_.push_back(event.str());
  }

  std::string render_log() const {
    std::ostringstream output;
    for (const std::string& event : events_) {
      output << event << '\n';
    }
    return output.str();
  }

  std::string render_summary(
      std::string_view status, std::string_view reason) const {
    std::ostringstream output;
    output << "{\n"
           << "  \"schema_version\": 1,\n"
           << "  \"algorithm\": "
              "\"exact-hamming-code-clique-lns-v2\",\n"
           << "  \"status\": \"" << json_escape(status) << "\",\n"
           << "  \"reason\": \"" << json_escape(reason) << "\",\n"
           << "  \"elapsed_seconds\": " << std::fixed
           << std::setprecision(6) << elapsed_seconds() << ",\n"
           << "  \"start\": \""
           << json_escape(arguments_.start.string()) << "\",\n"
           << "  \"output\": \""
           << json_escape(arguments_.output.string()) << "\",\n"
           << "  \"log\": \""
           << json_escape(arguments_.log.string()) << "\",\n"
           << "  \"seed\": " << arguments_.seed << ",\n"
           << "  \"transpose_start\": "
           << (arguments_.transpose_start ? "true" : "false") << ",\n"
           << "  \"seconds\": " << arguments_.seconds << ",\n"
           << "  \"destroy_min\": " << arguments_.destroy_min
           << ",\n"
           << "  \"destroy_max\": " << arguments_.destroy_max
           << ",\n"
           << "  \"destroy_without_replacement\": "
           << (arguments_.destroy_without_replacement
                   ? "true" : "false")
           << ",\n"
           << "  \"destroy_shard_count\": "
           << arguments_.destroy_shard_count << ",\n"
           << "  \"destroy_shard_index\": "
           << arguments_.destroy_shard_index << ",\n"
           << "  \"parent_generation\": "
           << parent_generation_ << ",\n"
           << "  \"pool_limit\": " << arguments_.pool_limit
           << ",\n"
           << "  \"allowed_distances\": "
           << allowed_distances_json() << ",\n"
           << "  \"candidate_universe_per_completed_pool\": "
           << static_cast<std::uint64_t>(kWordCount) << ",\n"
           << "  \"start_absolute_determinant\": \""
           << big_to_string(start_determinant_) << "\",\n"
           << "  \"best_absolute_determinant\": \""
           << big_to_string(best_determinant_) << "\",\n"
           << "  \"exactness\": {\n"
           << "    \"pool\": \"all oriented 23-bit words; never "
              "truncated\",\n"
           << "    \"repair\": \"exact clique branch-and-bound until "
              "the recorded stop condition\",\n"
           << "    \"ranking\": \"cpp_int Bareiss determinant\",\n"
           << "    \"pruning\": \"exact Hadamard-Fischer upper "
              "bounds with a fixed-Gram Schur residual bound and "
              "proper greedy-color feasibility\",\n"
           << "    \"claim_boundary\": \"a completed iteration is "
              "exhaustive only for its recorded destroyed rows and "
              "distance alphabet\"\n"
           << "  },\n"
           << "  \"statistics\": {\n"
           << "    \"iterations_started\": "
           << statistics_.iterations_started << ",\n"
           << "    \"iterations_completed\": "
           << statistics_.iterations_completed << ",\n"
           << "    \"interrupted_iterations\": "
           << statistics_.interrupted_iterations << ",\n"
           << "    \"pool_words_examined\": "
           << statistics_.pool_words_examined << ",\n"
           << "    \"pools_completed\": "
           << statistics_.pools_completed << ",\n"
           << "    \"pool_candidates\": "
           << statistics_.pool_candidates << ",\n"
           << "    \"maximum_pool_size\": "
           << statistics_.maximum_pool_size << ",\n"
           << "    \"graph_pair_checks\": "
           << statistics_.graph_pair_checks << ",\n"
           << "    \"graph_edges\": "
           << statistics_.graph_edges << ",\n"
           << "    \"graphs_completed\": "
           << statistics_.graphs_completed << ",\n"
           << "    \"clique_nodes\": "
           << statistics_.clique_nodes << ",\n"
           << "    \"cardinality_prunes\": "
           << statistics_.cardinality_prunes << ",\n"
           << "    \"color_checks\": "
           << statistics_.color_checks << ",\n"
           << "    \"color_prunes\": "
           << statistics_.color_prunes << ",\n"
           << "    \"bound_checks\": "
           << statistics_.bound_checks << ",\n"
           << "    \"bound_prunes\": "
           << statistics_.bound_prunes << ",\n"
           << "    \"clique_leaves\": "
           << statistics_.clique_leaves << ",\n"
           << "    \"exact_matrix_scores\": "
           << statistics_.exact_matrix_scores << ",\n"
           << "    \"improvements\": "
           << statistics_.improvements << "\n"
           << "  }\n"
           << "}\n";
    return output.str();
  }

  void flush_log() { log_.write(render_log()); }

  void write_summary(
      std::string_view status, std::string_view reason) {
    summary_.write(render_summary(status, reason));
  }

  bool poll() {
    if (stop_requested != 0) {
      stop_reason_ = "signal";
      return true;
    }
    const Clock::time_point now = Clock::now();
    if (now >= deadline_) {
      stop_reason_ = "deadline";
      return true;
    }
    if (now >= next_heartbeat_) {
      record_heartbeat();
      flush_log();
      write_summary("running", "");
      next_heartbeat_ =
          now + std::chrono::duration_cast<Clock::duration>(
                    std::chrono::duration<double>(
                        arguments_.heartbeat_seconds));
    }
    return false;
  }

  std::vector<int> choose_destroyed_rows(int destroy_count) {
    std::vector<int> rows(kNonAnchorRows);
    std::iota(rows.begin(), rows.end(), 1);
    std::shuffle(rows.begin(), rows.end(), generator_);
    rows.resize(static_cast<std::size_t>(destroy_count));
    std::sort(rows.begin(), rows.end());
    return rows;
  }

  std::optional<DestroyMaskChoice> choose_destroy_mask(
      std::uniform_int_distribution<int>& destroy_distribution) {
    if (destroy_schedule_.has_value()) {
      return destroy_schedule_->next(generator_);
    }
    const int destroy_count = destroy_distribution(generator_);
    const std::vector<int> destroyed =
        choose_destroyed_rows(destroy_count);
    Word mask = 0;
    for (const int row : destroyed) {
      mask |= Word{1} << (row - 1);
    }
    return DestroyMaskChoice{
        destroy_count, mask,
        destroy_mask_rank(
            mask, kNonAnchorRows, destroy_count),
        parent_generation_};
  }

  void audit_seed_clique(
      const RowWords& parent, const std::vector<int>& destroyed,
      const std::vector<Word>& pool, const BitGraph& graph) const {
    std::vector<std::size_t> locations;
    locations.reserve(destroyed.size());
    for (const int row : destroyed) {
      const Word word = parent[static_cast<std::size_t>(row)];
      const auto found = std::lower_bound(
          pool.begin(), pool.end(), word);
      if (found == pool.end() || *found != word) {
        throw std::runtime_error(
            "exact pool omitted a destroyed seed row");
      }
      locations.push_back(
          static_cast<std::size_t>(found - pool.begin()));
    }
    for (std::size_t first = 0; first < locations.size(); ++first) {
      for (std::size_t second = first + 1;
           second < locations.size(); ++second) {
        if (!graph.adjacent(locations[first], locations[second])) {
          throw std::runtime_error(
              "destroyed seed rows do not form a repair clique");
        }
      }
    }
  }

  void search_loop() {
    std::uniform_int_distribution<int> destroy_distribution(
        arguments_.destroy_min, arguments_.destroy_max);
    while (
        arguments_.max_iterations == 0 ||
        statistics_.iterations_started <
            arguments_.max_iterations) {
      phase_ = "between_iterations";
      if (poll()) return;

      const std::uint64_t iteration =
          statistics_.iterations_started + 1;
      const std::optional<DestroyMaskChoice> choice =
          choose_destroy_mask(destroy_distribution);
      if (!choice.has_value()) {
        stop_reason_ = "destroy_schedule_exhausted";
        return;
      }
      ++statistics_.iterations_started;
      const int destroy_count = choice->destroy_count;
      const std::vector<int> destroyed =
          destroyed_rows_from_mask(
              choice->mask, kNonAnchorRows);
      if (
          destroyed.size() !=
          static_cast<std::size_t>(destroy_count)) {
        throw std::runtime_error(
            "destroy schedule emitted a mask of wrong weight");
      }
      const RowWords parent = best_rows_;
      std::array<bool, kOrder> is_destroyed{};
      for (const int row : destroyed) {
        is_destroyed[static_cast<std::size_t>(row)] = true;
      }
      std::vector<Word> fixed;
      fixed.reserve(
          static_cast<std::size_t>(kOrder - destroy_count));
      for (int row = 0; row < kOrder; ++row) {
        if (!is_destroyed[static_cast<std::size_t>(row)]) {
          fixed.push_back(parent[static_cast<std::size_t>(row)]);
        }
      }
      record_iteration_start(iteration, destroyed, *choice);
      flush_log();

      const auto allowed = [this](const int distance) {
        return distance_allowed(
            distance, arguments_.allow_distance_14);
      };

      phase_ = "enumerating_pool";
      PoolResult pool = enumerate_compatible_words(
          kOrder, fixed, allowed, arguments_.pool_limit,
          [this]() { return poll(); });
      statistics_.pool_words_examined += pool.examined;
      if (!pool.complete) {
        ++statistics_.interrupted_iterations;
        return;
      }
      if (pool.examined != static_cast<std::uint64_t>(kWordCount)) {
        throw std::runtime_error(
            "completed pool did not examine all 2^23 words");
      }
      ++statistics_.pools_completed;
      statistics_.pool_candidates += pool.words.size();
      statistics_.maximum_pool_size = std::max(
          statistics_.maximum_pool_size,
          static_cast<std::uint64_t>(pool.words.size()));
      if (
          pool.words.size() <
          static_cast<std::size_t>(destroy_count)) {
        throw std::runtime_error(
            "candidate pool is too small to restore the seed");
      }
      record_pool_event(iteration, pool.words.size(), pool.examined);
      flush_log();

      phase_ = "building_graph";
      GraphResult graph = build_compatibility_graph(
          pool.words, allowed, [this]() { return poll(); });
      statistics_.graph_pair_checks += graph.pair_checks;
      statistics_.graph_edges += graph.edges;
      if (!graph.complete) {
        ++statistics_.interrupted_iterations;
        return;
      }
      ++statistics_.graphs_completed;
      audit_seed_clique(parent, destroyed, pool.words, graph.graph);
      record_graph_event(iteration, graph);
      flush_log();

      std::vector<std::size_t> root(pool.words.size());
      std::iota(root.begin(), root.end(), std::size_t{0});
      std::sort(
          root.begin(), root.end(),
          [&graph, &pool](
              const std::size_t first,
              const std::size_t second) {
            const std::uint32_t first_degree =
                graph.graph.degree(first);
            const std::uint32_t second_degree =
                graph.graph.degree(second);
            if (first_degree != second_degree) {
              return first_degree > second_degree;
            }
            return pool.words[first] < pool.words[second];
          });

      const std::uint64_t scored_before =
          statistics_.exact_matrix_scores;
      const std::uint64_t bounds_before =
          statistics_.bound_checks;
      const std::uint64_t prunes_before =
          statistics_.bound_prunes;
      const std::uint64_t colors_before =
          statistics_.color_checks;
      const std::uint64_t color_prunes_before =
          statistics_.color_prunes;
      phase_ = "clique_repair";

      const FixedSchurKernel schur_kernel(
          kOrder, fixed, pool.words, graph.edges,
          [&graph](
              const std::size_t first,
              const std::size_t second) {
            return graph.graph.adjacent(first, second);
          });
      BorderedKernelState kernel_state(schur_kernel);
      const BigInt reduced_fixed_power =
          big_power(
              schur_kernel.reduced_fixed_determinant(),
              destroy_count - 1);

      const auto prune =
          [this, &kernel_state, &schur_kernel,
           &reduced_fixed_power, destroy_count](
              const std::vector<std::size_t>& chosen,
              const std::vector<std::size_t>& candidates) {
            if (
                kernel_state.size() !=
                    static_cast<int>(chosen.size())) {
              throw std::runtime_error(
                  "clique/Schur-kernel depth mismatch");
            }
            ++statistics_.bound_checks;
            const int remaining =
                destroy_count -
                static_cast<int>(chosen.size());
            if (
                remaining < 0 ||
                candidates.size() <
                    static_cast<std::size_t>(remaining)) {
              throw std::runtime_error(
                  "invalid Schur-bound remaining depth");
            }
            const BigInt residual_product =
                largest_schur_diagonal_product(
                    schur_kernel, candidates, remaining);
            const BigInt left =
                schur_kernel.gcd() *
                kernel_state.determinant() *
                residual_product;
            const BigInt incumbent_squared =
                best_determinant_ * best_determinant_;
            const BigInt right =
                incumbent_squared * reduced_fixed_power;
            if (left <= right) {
              ++statistics_.bound_prunes;
              return true;
            }
            return false;
          };

      const auto visit =
          [this, &parent, &destroyed, &pool, iteration](
              const std::vector<std::size_t>& chosen) {
            ++statistics_.exact_matrix_scores;
            RowWords completion = parent;
            for (std::size_t index = 0;
                 index < destroyed.size(); ++index) {
              completion[static_cast<std::size_t>(
                  destroyed[index])] =
                  pool.words[chosen[index]];
            }
            const BigInt determinant =
                absolute(exact_sign_determinant(completion));
            if (determinant > best_determinant_) {
              const BigInt gram_determinant =
                  partial_gram_determinant(
                      std::vector<Word>(
                          completion.begin(), completion.end()),
                      kOrder);
              if (
                  gram_determinant !=
                  determinant * determinant) {
                throw std::runtime_error(
                    "matrix/Gram determinant identity failed");
              }
              const BigInt old_determinant = best_determinant_;
              best_rows_ = completion;
              best_determinant_ = determinant;
              ++statistics_.improvements;
              ++parent_generation_;
              if (destroy_schedule_.has_value()) {
                destroy_schedule_->reset(parent_generation_);
              }
              output_.write(row_words_text(best_rows_));
              record_improvement(iteration, old_determinant);
              flush_log();
              write_summary("running", "");
            }
          };

      CliqueWalker walker(
          graph.graph, destroy_count, prune, visit,
          [this]() { return poll(); },
          [&kernel_state](const std::size_t vertex) {
            static_cast<void>(
                kernel_state.prepare_push(vertex));
          },
          [&kernel_state]() {
            kernel_state.materialize_adjugate();
          },
          [&kernel_state]() { kernel_state.pop(); });
      const CliqueWalkStatistics walk =
          walker.run(std::move(root));
      statistics_.clique_nodes += walk.nodes;
      statistics_.cardinality_prunes +=
          walk.cardinality_prunes;
      statistics_.color_checks += walk.color_checks;
      statistics_.color_prunes += walk.color_prunes;
      statistics_.clique_leaves += walk.leaves;
      if (!walk.complete) {
        ++statistics_.interrupted_iterations;
        return;
      }
      ++statistics_.iterations_completed;
      record_iteration_finish(
          iteration, walk, scored_before, bounds_before,
          prunes_before, colors_before,
          color_prunes_before);
      flush_log();
      write_summary("running", "");
    }
  }

  void verify_owned_output() const {
    const SignMatrix written = read_sign_matrix(output_.path());
    const RowWords written_words = matrix_to_words(written);
    if (written_words != best_rows_) {
      throw std::runtime_error(
          "atomically written output differs from incumbent");
    }
    validate_hamming_code(
        written_words, arguments_.allow_distance_14);
    const BigInt determinant =
        absolute(exact_sign_determinant(written_words));
    if (determinant != best_determinant_) {
      throw std::runtime_error(
          "final exact output determinant disagrees with incumbent");
    }
  }

  Arguments arguments_;
  AtomicArtifact output_;
  AtomicArtifact log_;
  AtomicArtifact summary_;
  std::mt19937_64 generator_;
  std::optional<DestroyMaskSchedule> destroy_schedule_;
  Clock::time_point started_at_{};
  Clock::time_point deadline_{};
  Clock::time_point next_heartbeat_{};
  RowWords best_rows_{};
  BigInt start_determinant_ = 0;
  BigInt best_determinant_ = 0;
  std::uint64_t parent_generation_ = 0;
  SearchStatistics statistics_{};
  std::vector<std::string> events_;
  std::string phase_ = "initializing";
  std::string stop_reason_;
};

void require_test(bool condition, std::string_view message) {
  if (!condition) {
    throw std::runtime_error(
        "self-test failed: " + std::string(message));
  }
}

std::size_t brute_clique_number(const BitGraph& graph) {
  if (graph.size() > 20) {
    throw std::runtime_error(
        "brute clique self-test graph is too large");
  }
  const std::uint64_t subsets =
      UINT64_C(1) << graph.size();
  std::size_t maximum = 0;
  for (std::uint64_t subset = 0; subset < subsets; ++subset) {
    const std::size_t cardinality =
        static_cast<std::size_t>(std::popcount(subset));
    if (cardinality <= maximum) continue;
    bool clique = true;
    for (std::size_t first = 0;
         first < graph.size() && clique; ++first) {
      if (
          ((subset >> first) & UINT64_C(1)) == 0U) {
        continue;
      }
      for (std::size_t second = first + 1;
           second < graph.size(); ++second) {
        if (
            ((subset >> second) & UINT64_C(1)) != 0U &&
            !graph.adjacent(first, second)) {
          clique = false;
          break;
        }
      }
    }
    if (clique) maximum = cardinality;
  }
  return maximum;
}

void self_test_normalization_and_determinant() {
  SignMatrix matrix{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      matrix[row][column] = row == column ? -1 : 1;
    }
  }
  const BigInt expected =
      BigInt(21) * (BigInt(1) << 22);
  require_test(
      absolute(exact_sign_determinant(matrix_to_words(matrix))) ==
          expected,
      "det(J-2I) has the known exact value");
  const RowWords base_normalized =
      normalize_hamming_orientation(matrix);

  for (int row = 0; row < kOrder; ++row) {
    if ((row % 3) == 1) {
      for (int column = 0; column < kOrder; ++column) {
        matrix[row][column] = -matrix[row][column];
      }
    }
  }
  for (int column = 0; column < kOrder; ++column) {
    if ((column % 4) == 2) {
      for (int row = 0; row < kOrder; ++row) {
        matrix[row][column] = -matrix[row][column];
      }
    }
  }
  const BigInt switched_determinant =
      absolute(exact_sign_determinant(matrix_to_words(matrix)));
  require_test(
      absolute(exact_sign_determinant(
          matrix_to_words(transpose_sign_matrix(matrix)))) ==
          switched_determinant,
      "explicit input transposition preserves absolute determinant");
  const RowWords normalized =
      normalize_hamming_orientation(matrix);
  require_test(
      normalized == base_normalized,
      "normalization removes arbitrary row and column switches");
  require_test(
      normalized[0] == 0U,
      "normalization makes the anchor row all +1");
  require_test(
      absolute(exact_sign_determinant(normalized)) ==
          switched_determinant,
      "normalization preserves absolute determinant");
  for (int first = 0; first < kOrder; ++first) {
    for (int second = first + 1; second < kOrder; ++second) {
      const int distance = std::popcount(
          normalized[static_cast<std::size_t>(first)] ^
          normalized[static_cast<std::size_t>(second)]);
      require_test(
          positive_mod4(kOrder - 2 * distance) == 3,
          "normalized Gram entries are 3 modulo 4");
    }
  }
}

void self_test_pool_enumeration() {
  constexpr int order = 7;
  const std::vector<Word> fixed{
      0U, UINT32_C(0b0000011), UINT32_C(0b0011100)};
  const auto allowed = [](const int distance) {
    return distance == 2 || distance == 4;
  };
  const PoolResult pool = enumerate_compatible_words(
      order, fixed, allowed, UINT64_C(1024),
      []() { return false; });
  std::vector<Word> brute;
  for (Word word = 0; word < (Word{1} << order); ++word) {
    const int anchor_distance = std::popcount(word);
    const int second_distance =
        std::popcount(word ^ fixed[1]);
    const int third_distance =
        std::popcount(word ^ fixed[2]);
    if (
        (anchor_distance == 2 || anchor_distance == 4) &&
        (second_distance == 2 || second_distance == 4) &&
        (third_distance == 2 || third_distance == 4)) {
      brute.push_back(word);
    }
  }
  require_test(pool.complete, "small pool completes");
  require_test(
      pool.examined == (UINT64_C(1) << order),
      "small pool examines its full oriented universe");
  require_test(
      pool.words == brute,
      "pool enumeration agrees with independent brute predicate");

  bool refused_truncation = false;
  try {
    static_cast<void>(enumerate_compatible_words(
        order, std::vector<Word>{0U},
        [](const int distance) {
          return distance == 2 || distance == 4;
        },
        UINT64_C(1), []() { return false; }));
  } catch (const std::runtime_error& error) {
    refused_truncation =
        std::string(error.what()).find("no truncated graph") !=
        std::string::npos;
  }
  require_test(
      refused_truncation,
      "pool limit aborts rather than truncating");
}

void self_test_graph_and_cliques() {
  constexpr int order = 6;
  std::vector<Word> words;
  for (Word word = 0; word < (Word{1} << order); ++word) {
    if (std::popcount(word) == 2) words.push_back(word);
  }
  const auto allowed = [](const int distance) {
    return distance == 2 || distance == 4;
  };
  const GraphResult graph = build_compatibility_graph(
      words, allowed, []() { return false; });
  require_test(graph.complete, "small graph completes");
  for (std::size_t first = 0; first < words.size(); ++first) {
    for (std::size_t second = 0; second < words.size(); ++second) {
      const bool expected =
          first != second &&
          allowed(std::popcount(words[first] ^ words[second]));
      require_test(
          graph.graph.adjacent(first, second) == expected,
          "bit graph has exactly the compatibility edges");
    }
  }

  std::vector<std::vector<Word>> walked;
  std::vector<std::size_t> root(words.size());
  std::iota(root.begin(), root.end(), std::size_t{0});
  CliqueWalker walker(
      graph.graph, 3,
      [](
          const std::vector<std::size_t>&,
          const std::vector<std::size_t>&) {
        return false;
      },
      [&walked, &words](
          const std::vector<std::size_t>& clique) {
        std::vector<Word> selected;
        for (const std::size_t vertex : clique) {
          selected.push_back(words[vertex]);
        }
        std::sort(selected.begin(), selected.end());
        walked.push_back(std::move(selected));
      },
      []() { return false; },
      [](std::size_t) {},
      []() {},
      []() {});
  const CliqueWalkStatistics walk = walker.run(std::move(root));
  require_test(walk.complete, "small clique walk completes");

  std::vector<std::vector<Word>> brute;
  for (std::size_t first = 0; first < words.size(); ++first) {
    for (std::size_t second = first + 1;
         second < words.size(); ++second) {
      for (std::size_t third = second + 1;
           third < words.size(); ++third) {
        if (
            graph.graph.adjacent(first, second) &&
            graph.graph.adjacent(first, third) &&
            graph.graph.adjacent(second, third)) {
          brute.push_back(
              std::vector<Word>{
                  words[first], words[second], words[third]});
        }
      }
    }
  }
  std::sort(walked.begin(), walked.end());
  std::sort(brute.begin(), brute.end());
  require_test(
      walked == brute,
      "clique branch-and-bound agrees with brute combinations");
  require_test(
      walk.leaves == brute.size(),
      "clique leaf count is exact");

  BitGraph color_graph(9);
  for (std::size_t first = 0; first < color_graph.size(); ++first) {
    for (std::size_t second = first + 1;
         second < color_graph.size(); ++second) {
      if (
          ((first * 7U + second * 11U +
            first * second) % 5U) < 2U) {
        color_graph.add_edge(first, second);
      }
    }
  }
  const std::size_t omega = brute_clique_number(color_graph);
  std::vector<std::size_t> color_order(color_graph.size());
  std::iota(
      color_order.begin(), color_order.end(), std::size_t{0});
  for (int permutation = 0; permutation < 4; ++permutation) {
    if (permutation != 0) {
      std::rotate(
          color_order.begin(), color_order.begin() + 2,
          color_order.end());
    }
    const auto classes =
        greedy_independent_color_classes(
            color_graph, color_order);
    for (const auto& color_class : classes) {
      for (std::size_t first = 0;
           first < color_class.size(); ++first) {
        for (std::size_t second = first + 1;
             second < color_class.size(); ++second) {
          require_test(
              !color_graph.adjacent(
                  color_class[first], color_class[second]),
              "every greedy color class is independent");
        }
      }
    }
    require_test(
        omega <= classes.size(),
        "proper greedy color count upper-bounds omega");
  }
}

void self_test_hadamard_fischer() {
  std::mt19937_64 generator(UINT64_C(230014));
  for (int order : {5, 7}) {
    for (int trial = 0; trial < 16; ++trial) {
      std::vector<Word> rows(static_cast<std::size_t>(order));
      const Word mask = (Word{1} << order) - 1U;
      for (Word& row : rows) {
        row = static_cast<Word>(generator()) & mask;
      }
      const BigInt determinant =
          absolute(exact_sign_determinant(rows, order));
      const BigInt determinant_squared =
          determinant * determinant;
      for (int count = 1; count <= order; ++count) {
        const std::vector<Word> partial(
            rows.begin(), rows.begin() + count);
        const BigInt upper_squared =
            partial_gram_determinant(partial, order) *
            integer_power(order, order - count);
        require_test(
            determinant_squared <= upper_squared,
            "Hadamard-Fischer bound dominates exact determinant");
        if (count == order) {
          require_test(
              determinant_squared == upper_squared,
              "full Gram determinant equals det(A)^2");
        }
      }
    }
  }
}

void validate_bordered_gram_state(
    const std::vector<Word>& rows, int order,
    const BorderedGramState& state,
    std::string_view context) {
  const auto message =
      [context](std::string_view detail) {
        return std::string(context) + ": " + std::string(detail);
      };
  require_test(
      state.size() == static_cast<int>(rows.size()),
      message("prefix size agrees"));
  const std::vector<std::vector<BigInt>> gram =
      partial_gram_matrix(rows, order);
  const BigInt determinant = bareiss(gram);
  require_test(
      determinant >= 0,
      message("independent Gram determinant is non-negative"));
  require_test(
      state.determinant() == determinant,
      message("bordered determinant agrees with Bareiss"));
  const std::vector<std::vector<BigInt>> adjugate =
      exact_adjugate(gram);
  for (std::size_t row = 0; row < rows.size(); ++row) {
    for (std::size_t column = 0;
         column < rows.size(); ++column) {
      require_test(
          state.adjugate(
              static_cast<int>(row),
              static_cast<int>(column)) ==
              adjugate[row][column],
          message("bordered adjugate agrees with cofactors"));
      require_test(
          state.adjugate(
              static_cast<int>(row),
              static_cast<int>(column)) ==
              state.adjugate(
                  static_cast<int>(column),
                  static_cast<int>(row)),
          message("bordered adjugate is symmetric"));
      BigInt product = 0;
      for (std::size_t inner = 0;
           inner < rows.size(); ++inner) {
        product +=
            gram[row][inner] *
            state.adjugate(
                static_cast<int>(inner),
                static_cast<int>(column));
      }
      const BigInt expected =
          row == column ? determinant : BigInt(0);
      require_test(
          product == expected,
          message("G times adj(G) equals det(G) times I"));
    }
  }
}

void validate_pending_bordered_determinant(
    const std::vector<Word>& rows, int order,
    const BorderedGramState& state,
    std::string_view context) {
  require_test(
      state.size() == static_cast<int>(rows.size()),
      std::string(context) + ": pending prefix size agrees");
  require_test(
      state.determinant() ==
          bareiss(partial_gram_matrix(rows, order)),
      std::string(context) +
          ": pending bordered determinant agrees with Bareiss");
}

void self_test_bordered_gram_state() {
  std::mt19937_64 generator(UINT64_C(230023));
  for (const int order : {5, 7, 9}) {
    const Word mask = (Word{1} << order) - 1U;
    for (int trial = 0; trial < 6; ++trial) {
      BorderedGramState state(order);
      std::vector<Word> rows;
      for (int count = 0; count < order; ++count) {
        if (!rows.empty() && count < order) {
          rows.push_back(rows.front());
          static_cast<void>(
              state.prepare_push(rows.back()));
          validate_pending_bordered_determinant(
              rows, order, state,
              "pending duplicate singular extension");
          state.materialize_adjugate();
          validate_bordered_gram_state(
              rows, order, state,
              "randomized duplicate singular extension");
          require_test(
              state.determinant() == 0,
              "duplicate row makes bordered Gram singular");
          state.pop();
          rows.pop_back();
          validate_bordered_gram_state(
              rows, order, state,
              "state restored after singular extension");
        }

        bool accepted = false;
        for (int attempt = 0; attempt < 4096; ++attempt) {
          const Word word =
              static_cast<Word>(generator()) & mask;
          rows.push_back(word);
          static_cast<void>(state.prepare_push(word));
          validate_pending_bordered_determinant(
              rows, order, state,
              "pending randomized prefix");
          state.materialize_adjugate();
          validate_bordered_gram_state(
              rows, order, state,
              "randomized nonsingular-prefix trial");
          if (state.determinant() != 0) {
            accepted = true;
            break;
          }
          state.pop();
          rows.pop_back();
        }
        require_test(
            accepted,
            "randomized test found a nonsingular next row");
      }

      while (!rows.empty()) {
        state.pop();
        rows.pop_back();
        if (rows.empty()) {
          require_test(
              state.size() == 0 && state.determinant() == 1,
              "empty bordered Gram state has determinant one");
        } else {
          validate_bordered_gram_state(
              rows, order, state,
              "randomized pop restores exact prefix");
        }
      }
    }
  }

  const std::vector<Word> dependent_unique{
      UINT32_C(0b0000), UINT32_C(0b0011),
      UINT32_C(0b0100), UINT32_C(0b1000)};
  BorderedGramState singular_state(4);
  std::vector<Word> singular_prefix;
  for (const Word word : dependent_unique) {
    singular_prefix.push_back(word);
    static_cast<void>(singular_state.prepare_push(word));
    validate_pending_bordered_determinant(
        singular_prefix, 4, singular_state,
        "pending unique-row singular-prefix fixture");
    singular_state.materialize_adjugate();
    validate_bordered_gram_state(
        singular_prefix, 4, singular_state,
        "unique-row singular-prefix fixture");
  }
  require_test(
      singular_state.determinant() == 0,
      "unique rows with repeated sign columns are singular");
}

void validate_bordered_kernel_state(
    const FixedSchurKernel& kernel,
    const std::vector<std::size_t>& vertices,
    const BorderedKernelState& state,
    std::string_view context) {
  require_test(
      state.size() == static_cast<int>(vertices.size()),
      std::string(context) +
          ": bordered-kernel prefix size agrees");
  std::vector<std::vector<BigInt>> principal(
      vertices.size(),
      std::vector<BigInt>(vertices.size()));
  for (std::size_t row = 0; row < vertices.size(); ++row) {
    for (std::size_t column = 0;
         column < vertices.size(); ++column) {
      principal[row][column] =
          kernel.at(vertices[row], vertices[column]);
    }
  }
  const BigInt determinant = bareiss(principal);
  require_test(
      determinant >= 0 &&
          state.determinant() == determinant &&
          principal_kernel_determinant(kernel, vertices) ==
              determinant,
      std::string(context) +
          ": bordered-kernel determinant agrees with Bareiss");
  const auto adjugate = exact_adjugate(principal);
  for (std::size_t row = 0; row < vertices.size(); ++row) {
    for (std::size_t column = 0;
         column < vertices.size(); ++column) {
      require_test(
          state.adjugate(
              static_cast<int>(row),
              static_cast<int>(column)) ==
              adjugate[row][column],
          std::string(context) +
              ": bordered-kernel adjugate agrees with cofactors");
    }
  }
}

void self_test_schur_residual_bound() {
  std::mt19937_64 generator(UINT64_C(230031));
  for (const int order : {5, 7}) {
    const int repair_count = 3;
    const int fixed_count = order - repair_count;
    const Word mask = (Word{1} << order) - 1U;
    for (int trial = 0; trial < 8; ++trial) {
      std::vector<Word> full_rows(
          static_cast<std::size_t>(order));
      bool nonsingular = false;
      for (int attempt = 0; attempt < 4096; ++attempt) {
        for (Word& word : full_rows) {
          word = static_cast<Word>(generator()) & mask;
        }
        if (exact_sign_determinant(full_rows, order) != 0) {
          nonsingular = true;
          break;
        }
      }
      require_test(
          nonsingular,
          "Schur-bound fixture found a nonsingular sign matrix");

      const std::vector<Word> fixed(
          full_rows.begin(),
          full_rows.begin() + fixed_count);
      std::vector<Word> candidates(
          full_rows.begin() + fixed_count, full_rows.end());
      candidates.push_back(
          static_cast<Word>(generator()) & mask);
      candidates.push_back(
          static_cast<Word>(generator()) & mask);
      candidates.push_back(
          static_cast<Word>(generator()) & mask);
      const std::uint64_t complete_edges =
          static_cast<std::uint64_t>(candidates.size()) *
          static_cast<std::uint64_t>(candidates.size() - 1) /
          2;
      const FixedSchurKernel kernel(
          order, fixed, candidates, complete_edges,
          [](std::size_t first, std::size_t second) {
            return first != second;
          });
      require_test(
          kernel.determinant() ==
              kernel.gcd() *
                  kernel.reduced_fixed_determinant(),
          "Schur gcd exactly factors the fixed determinant");

      BorderedKernelState state(kernel);
      std::vector<std::size_t> state_vertices;
      for (int index = 0; index < repair_count; ++index) {
        state_vertices.push_back(
            static_cast<std::size_t>(index));
        static_cast<void>(state.prepare_push(
            static_cast<std::size_t>(index)));
        require_test(
            state.determinant() ==
                principal_kernel_determinant(
                    kernel, state_vertices),
            "pending bordered-kernel determinant agrees with "
            "Bareiss");
        state.materialize_adjugate();
        validate_bordered_kernel_state(
            kernel, state_vertices, state,
            "randomized Schur-kernel prefix");
        require_test(
            state.determinant() > 0,
            "nonsingular full-row prefix has positive "
            "Schur-kernel determinant");
      }
      while (!state_vertices.empty()) {
        state.pop();
        state_vertices.pop_back();
        if (!state_vertices.empty()) {
          validate_bordered_kernel_state(
              kernel, state_vertices, state,
              "bordered-kernel pop");
        }
      }

      const BigInt full_gram =
          partial_gram_determinant(full_rows, order);
      for (int chosen_count = 0;
           chosen_count <= repair_count; ++chosen_count) {
        std::vector<std::size_t> chosen;
        std::vector<Word> partial = fixed;
        for (int index = 0; index < chosen_count; ++index) {
          chosen.push_back(static_cast<std::size_t>(index));
          partial.push_back(
              candidates[static_cast<std::size_t>(index)]);
        }
        const BigInt kernel_determinant =
            principal_kernel_determinant(kernel, chosen);
        const BigInt partial_determinant =
            partial_gram_determinant(partial, order);
        if (chosen_count == 0) {
          require_test(
              partial_determinant ==
                  kernel.gcd() *
                      kernel.reduced_fixed_determinant(),
              "empty Schur selection recovers fixed determinant");
        } else {
          require_test(
              partial_determinant *
                      big_power(
                          kernel.reduced_fixed_determinant(),
                          chosen_count - 1) ==
                  kernel.gcd() * kernel_determinant,
              "Schur principal determinant identity");
        }

        std::vector<std::size_t> current;
        for (std::size_t index =
                 static_cast<std::size_t>(chosen_count);
             index < candidates.size(); ++index) {
          current.push_back(index);
        }
        const int remaining =
            repair_count - chosen_count;
        const BigInt strong_numerator =
            kernel.gcd() * kernel_determinant *
            largest_schur_diagonal_product(
                kernel, current, remaining);
        const BigInt denominator =
            big_power(
                kernel.reduced_fixed_determinant(),
                repair_count - 1);
        require_test(
            full_gram * denominator <= strong_numerator,
            "strengthened Schur-Fischer bound dominates the "
            "known completion");
        require_test(
            full_gram <=
                partial_determinant *
                    integer_power(order, remaining),
            "old Hadamard-Fischer bound dominates the known "
            "completion");
        if (remaining == 0) {
          require_test(
              full_gram * denominator == strong_numerator,
              "full Schur-Fischer endpoint is exact");
        }

        BigInt brute_maximum = 0;
        std::vector<std::size_t> extension;
        const std::function<void(std::size_t, int)> enumerate =
            [&](const std::size_t begin, const int needed) {
              if (needed == 0) {
                std::vector<Word> completion = partial;
                for (const std::size_t vertex : extension) {
                  completion.push_back(candidates[vertex]);
                }
                const BigInt determinant =
                    partial_gram_determinant(
                        completion, order);
                brute_maximum =
                    std::max(brute_maximum, determinant);
                require_test(
                    determinant <=
                        partial_determinant *
                            integer_power(order, remaining),
                    "old Fischer bound dominates a brute "
                    "completion");
                require_test(
                    determinant * denominator <=
                        strong_numerator,
                    "Schur-Fischer bound dominates a brute "
                    "completion");
                return;
              }
              for (std::size_t vertex = begin;
                   vertex + static_cast<std::size_t>(needed) <=
                       candidates.size();
                   ++vertex) {
                extension.push_back(vertex);
                enumerate(vertex + 1, needed - 1);
                extension.pop_back();
              }
            };
        enumerate(
            static_cast<std::size_t>(chosen_count), remaining);
        require_test(
            brute_maximum * denominator <= strong_numerator,
            "Schur-Fischer bound dominates brute maximum");
      }
    }
  }

  const std::vector<Word> singular_fixed{
      UINT32_C(0b00000), UINT32_C(0b00011)};
  const std::vector<Word> singular_candidates{
      UINT32_C(0b00000), UINT32_C(0b00101),
      UINT32_C(0b01001)};
  const FixedSchurKernel singular_kernel(
      5, singular_fixed, singular_candidates, 3,
      [](std::size_t first, std::size_t second) {
        return first != second;
      });
  require_test(
      singular_kernel.at(0, 0) == 0 &&
          principal_kernel_determinant(
              singular_kernel, {0}) == 0,
      "candidate in the fixed span has zero Schur residual");
  const BigInt zero_bound =
      singular_kernel.gcd() *
      principal_kernel_determinant(
          singular_kernel, {0}) *
      largest_schur_diagonal_product(
          singular_kernel, {1, 2}, 2);
  require_test(
      zero_bound == 0,
      "singular Schur prefix has a zero strengthened bound");
}

void self_test_sparse_schur_storage() {
  constexpr int order = 7;
  const std::vector<Word> fixed{
      UINT32_C(0b0000000), UINT32_C(0b0000011)};
  std::vector<Word> candidates;
  candidates.reserve(70);
  for (std::size_t vertex = 0; vertex < 70; ++vertex) {
    candidates.push_back(
        static_cast<Word>((vertex * 37U + 19U) & 127U));
  }

  BitGraph graph(candidates.size());
  std::uint64_t edge_count = 0;
  for (std::size_t first = 0; first < candidates.size(); ++first) {
    for (std::size_t second = first + 1;
         second < candidates.size(); ++second) {
      const bool patterned =
          ((first * 13U + second * 17U +
            first * second) % 11U) < 3U;
      const bool block_boundary =
          (first == 0 && second == 64) ||
          (first == 2 && second == 69);
      if (patterned || block_boundary) {
        graph.add_edge(first, second);
        ++edge_count;
      }
    }
  }
  const FixedSchurKernel kernel(
      order, fixed, candidates, edge_count,
      [&graph](
          const std::size_t first,
          const std::size_t second) {
        return graph.adjacent(first, second);
      });
  require_test(
      kernel.stored_edges() == edge_count,
      "sparse Schur kernel retains every declared edge");
  require_test(
      kernel.planned_storage_bytes() <
          candidates.size() * candidates.size() *
              sizeof(BigInt),
      "sparse Schur fixture plans less storage than dense");

  BorderedGramState fixed_state(order);
  for (const Word word : fixed) {
    static_cast<void>(fixed_state.push(word));
  }
  require_test(
      fixed_state.determinant() > 0,
      "sparse Schur fixture has nonsingular fixed rows");
  std::vector<std::array<int, 2>> correlations(
      candidates.size());
  std::vector<std::array<BigInt, 2>> products(
      candidates.size());
  for (std::size_t vertex = 0;
       vertex < candidates.size(); ++vertex) {
    for (std::size_t row = 0; row < fixed.size(); ++row) {
      correlations[vertex][row] =
          order -
          2 * std::popcount(
                  candidates[vertex] ^ fixed[row]);
    }
    for (std::size_t row = 0; row < fixed.size(); ++row) {
      for (std::size_t column = 0;
           column < fixed.size(); ++column) {
        products[vertex][row] +=
            fixed_state.adjugate(
                static_cast<int>(row),
                static_cast<int>(column)) *
            correlations[vertex][column];
      }
    }
  }
  const auto raw_entry =
      [&fixed_state, &correlations, &products, &candidates](
          const std::size_t first,
          const std::size_t second) {
        BigInt value =
            fixed_state.determinant() *
            (order -
             2 * std::popcount(
                     candidates[first] ^ candidates[second]));
        for (std::size_t row = 0; row < 2; ++row) {
          value -=
              correlations[first][row] *
              products[second][row];
        }
        return value;
      };

  BigInt expected_gcd = fixed_state.determinant();
  for (std::size_t vertex = 0;
       vertex < candidates.size(); ++vertex) {
    const BigInt raw = raw_entry(vertex, vertex);
    expected_gcd = exact_gcd(expected_gcd, raw);
  }
  for (std::size_t first = 0; first < candidates.size(); ++first) {
    for (std::size_t second = first + 1;
         second < candidates.size(); ++second) {
      if (graph.adjacent(first, second)) {
        expected_gcd =
            exact_gcd(expected_gcd, raw_entry(first, second));
      }
    }
  }
  require_test(
      kernel.gcd() == expected_gcd,
      "sparse Schur gcd covers D, diagonals, and every edge");
  for (std::size_t vertex = 0;
       vertex < candidates.size(); ++vertex) {
    require_test(
        kernel.gcd() * kernel.at(vertex, vertex) ==
            raw_entry(vertex, vertex),
        "sparse Schur diagonal matches independent formula");
  }

  bool found_nonedge = false;
  std::vector<std::size_t> triangle;
  for (std::size_t first = 0; first < candidates.size(); ++first) {
    for (std::size_t second = first + 1;
         second < candidates.size(); ++second) {
      if (graph.adjacent(first, second)) {
        require_test(
            kernel.at(first, second) ==
                    kernel.at(second, first) &&
                kernel.gcd() * kernel.at(first, second) ==
                    raw_entry(first, second),
            "sparse Schur edge rank maps to its exact entry");
      } else if (!found_nonedge) {
        bool rejected = false;
        try {
          static_cast<void>(kernel.at(first, second));
        } catch (const std::runtime_error& error) {
          rejected =
              std::string(error.what()).find(
                  "not a compatibility edge") !=
              std::string::npos;
        }
        require_test(
            rejected,
            "sparse Schur storage rejects a nonedge access");
        found_nonedge = true;
      }
      if (!triangle.empty() || !graph.adjacent(first, second)) {
        continue;
      }
      for (std::size_t third = second + 1;
           third < candidates.size(); ++third) {
        if (
            graph.adjacent(first, third) &&
            graph.adjacent(second, third)) {
          triangle = {first, second, third};
          break;
        }
      }
    }
  }
  require_test(
      found_nonedge,
      "sparse Schur fixture contains a guarded nonedge");
  require_test(
      kernel.gcd() * kernel.at(0, 64) ==
          raw_entry(0, 64),
      "sparse edge rank crosses a 64-bit block boundary");
  require_test(
      !triangle.empty(),
      "sparse Schur fixture contains a test triangle");
  std::vector<Word> triangle_rows = fixed;
  for (const std::size_t vertex : triangle) {
    triangle_rows.push_back(candidates[vertex]);
  }
  const BigInt triangle_gram =
      partial_gram_determinant(triangle_rows, order);
  require_test(
      triangle_gram *
              big_power(
                  kernel.reduced_fixed_determinant(), 2) ==
          kernel.gcd() *
              principal_kernel_determinant(kernel, triangle),
      "sparse clique preserves the reduced Schur identity");
  BorderedKernelState triangle_state(kernel);
  std::vector<std::size_t> triangle_prefix;
  for (const std::size_t vertex : triangle) {
    triangle_prefix.push_back(vertex);
    static_cast<void>(triangle_state.prepare_push(vertex));
    require_test(
        triangle_state.determinant() ==
            principal_kernel_determinant(
                kernel, triangle_prefix),
        "sparse bordered-kernel determinant agrees with Bareiss");
    triangle_state.materialize_adjugate();
  }

  const std::vector<Word> small_candidates{
      UINT32_C(0b0000101), UINT32_C(0b0001001),
      UINT32_C(0b0010001)};
  const FixedSchurKernel empty_kernel(
      order, fixed, small_candidates, 0,
      [](std::size_t, std::size_t) { return false; });
  require_test(
      empty_kernel.stored_edges() == 0,
      "sparse Schur storage accepts an empty graph");
  bool empty_nonedge_rejected = false;
  try {
    static_cast<void>(empty_kernel.at(0, 1));
  } catch (const std::runtime_error&) {
    empty_nonedge_rejected = true;
  }
  require_test(
      empty_nonedge_rejected,
      "empty sparse Schur graph rejects off-diagonal access");

  bool count_mismatch_rejected = false;
  try {
    static_cast<void>(FixedSchurKernel(
        order, fixed, small_candidates, 0,
        [](std::size_t first, std::size_t second) {
          return first == 0 && second == 1;
        }));
  } catch (const std::runtime_error& error) {
    count_mismatch_rejected =
        std::string(error.what()).find(
            "more edges") != std::string::npos;
  }
  require_test(
      count_mismatch_rejected,
      "sparse Schur storage rejects excess predicate edges");

  bool impossible_count_rejected = false;
  try {
    static_cast<void>(FixedSchurKernel(
        order, fixed, small_candidates, 4,
        [](std::size_t, std::size_t) { return true; }));
  } catch (const std::runtime_error& error) {
    impossible_count_rejected =
        std::string(error.what()).find(
            "possible pairs") != std::string::npos;
  }
  require_test(
      impossible_count_rejected,
      "sparse Schur storage rejects an impossible edge count");

  std::vector<Word> oversized_candidates(3000, 0U);
  bool memory_guard_rejected = false;
  bool predicate_called = false;
  try {
    const std::uint64_t complete_edges =
        static_cast<std::uint64_t>(
            oversized_candidates.size()) *
        static_cast<std::uint64_t>(
            oversized_candidates.size() - 1) /
        2;
    static_cast<void>(FixedSchurKernel(
        order, fixed, oversized_candidates, complete_edges,
        [&predicate_called](std::size_t, std::size_t) {
          predicate_called = true;
          return true;
        }));
  } catch (const std::runtime_error& error) {
    memory_guard_rejected =
        std::string(error.what()).find("128 MiB") !=
        std::string::npos;
  }
  require_test(
      memory_guard_rejected && !predicate_called,
      "sparse Schur byte guard fires before edge construction");
}

void self_test_destroy_schedules() {
  for (int row_count = 1; row_count <= 8; ++row_count) {
    for (int destroy_count = 1;
         destroy_count <= row_count; ++destroy_count) {
      const std::uint64_t total =
          binomial_coefficient(row_count, destroy_count);
      Word previous = 0;
      for (std::uint64_t rank = 0; rank < total; ++rank) {
        const Word mask = destroy_mask_from_rank(
            rank, row_count, destroy_count);
        require_test(
            destroy_mask_rank(
                mask, row_count, destroy_count) == rank,
            "destroy-mask rank/unrank roundtrip");
        if (rank != 0) {
          require_test(
              previous < mask,
              "destroy-mask ranks increase with numeric mask");
        }
        previous = mask;
      }
    }
  }

  using MaskIdentity = std::pair<int, Word>;
  std::set<MaskIdentity> shard_union;
  std::array<std::set<MaskIdentity>, 3> memberships;
  for (std::uint64_t shard = 0; shard < 3; ++shard) {
    DestroyMaskSchedule schedule(
        6, 2, 4, UINT64_C(117), 3, shard);
    std::mt19937_64 generator(UINT64_C(8000) + shard);
    while (const auto choice = schedule.next(generator)) {
      const MaskIdentity identity{
          choice->destroy_count, choice->mask};
      require_test(
          memberships[static_cast<std::size_t>(shard)]
              .insert(identity)
              .second,
          "a no-replacement shard has no repeated masks");
      require_test(
          choice->rank % 3 == shard,
          "destroy-mask rank selects its shard");
      require_test(
          destroy_mask_from_rank(
              choice->rank, 6, choice->destroy_count) ==
              choice->mask,
          "scheduled destroy-mask rank roundtrips");
      shard_union.insert(identity);
    }
    require_test(
        schedule.remaining() == 0,
        "drained destroy schedule reports no remaining masks");

    DestroyMaskSchedule different_seed(
        6, 2, 4, UINT64_C(999), 3, shard);
    std::mt19937_64 other_generator(UINT64_C(9100) + shard);
    std::set<MaskIdentity> other_membership;
    while (const auto choice =
               different_seed.next(other_generator)) {
      other_membership.emplace(
          choice->destroy_count, choice->mask);
    }
    require_test(
        other_membership ==
            memberships[static_cast<std::size_t>(shard)],
        "shard membership is independent of seed");
  }
  for (std::size_t first = 0; first < memberships.size(); ++first) {
    for (std::size_t second = first + 1;
         second < memberships.size(); ++second) {
      std::vector<MaskIdentity> overlap;
      std::set_intersection(
          memberships[first].begin(), memberships[first].end(),
          memberships[second].begin(), memberships[second].end(),
          std::back_inserter(overlap));
      require_test(
          overlap.empty(),
          "different destroy shards are disjoint");
    }
  }
  std::size_t expected_union = 0;
  for (int destroy_count = 2; destroy_count <= 4;
       ++destroy_count) {
    expected_union += static_cast<std::size_t>(
        binomial_coefficient(6, destroy_count));
  }
  require_test(
      shard_union.size() == expected_union,
      "destroy shards unite to the complete mask family");

  DestroyMaskSchedule first(
      7, 2, 5, UINT64_C(81234), 2, 1);
  DestroyMaskSchedule second(
      7, 2, 5, UINT64_C(81234), 2, 1);
  std::mt19937_64 first_generator(UINT64_C(551));
  std::mt19937_64 second_generator(UINT64_C(551));
  for (;;) {
    const auto first_choice = first.next(first_generator);
    const auto second_choice = second.next(second_generator);
    require_test(
        first_choice.has_value() ==
            second_choice.has_value(),
        "destroy schedule exhaustion is reproducible");
    if (!first_choice.has_value()) break;
    require_test(
        first_choice->destroy_count ==
                second_choice->destroy_count &&
            first_choice->mask == second_choice->mask &&
            first_choice->rank == second_choice->rank &&
            first_choice->parent_generation ==
                second_choice->parent_generation,
        "destroy schedule sequence is reproducible");
  }

  first.reset(4);
  second.reset(4);
  first_generator.seed(UINT64_C(992));
  second_generator.seed(UINT64_C(992));
  for (int index = 0; index < 12; ++index) {
    const auto first_choice = first.next(first_generator);
    const auto second_choice = second.next(second_generator);
    require_test(
        first_choice.has_value() &&
            second_choice.has_value() &&
            first_choice->mask == second_choice->mask &&
            first_choice->destroy_count ==
                second_choice->destroy_count &&
            first_choice->parent_generation == 4 &&
            second_choice->parent_generation == 4,
        "parent-generation reset is deterministic");
  }

  bool invalid_schedule = false;
  try {
    static_cast<void>(DestroyMaskSchedule(
        6, 2, 4, 1, 0, 0));
  } catch (const std::runtime_error&) {
    invalid_schedule = true;
  }
  require_test(
      invalid_schedule,
      "zero-count destroy shard is rejected");
  invalid_schedule = false;
  try {
    static_cast<void>(DestroyMaskSchedule(
        6, 2, 4, 1, 2, 2));
  } catch (const std::runtime_error&) {
    invalid_schedule = true;
  }
  require_test(
      invalid_schedule,
      "out-of-range destroy shard is rejected");

  const auto parse_test_arguments =
      [](std::vector<std::string> values) {
        std::vector<char*> pointers;
        pointers.reserve(values.size());
        for (std::string& value : values) {
          pointers.push_back(value.data());
        }
        return parse_arguments(
            static_cast<int>(pointers.size()),
            pointers.data());
      };
  const Arguments implicit_sharding = parse_test_arguments(
      {"test", "--self-test", "--destroy-shard-count", "2"});
  require_test(
      implicit_sharding.destroy_without_replacement &&
          implicit_sharding.destroy_shard_count == 2 &&
          implicit_sharding.destroy_shard_index == 0,
      "shard count enables no-replacement scheduling");
  const Arguments implicit_index = parse_test_arguments(
      {"test", "--self-test", "--destroy-shard-index", "0"});
  require_test(
      implicit_index.destroy_without_replacement,
      "shard index enables no-replacement scheduling");
  bool invalid_arguments = false;
  try {
    static_cast<void>(parse_test_arguments(
        {"test", "--self-test", "--destroy-shard-index", "1"}));
  } catch (const std::runtime_error&) {
    invalid_arguments = true;
  }
  require_test(
      invalid_arguments,
      "CLI rejects a shard index outside the default count");
  invalid_arguments = false;
  try {
    static_cast<void>(parse_test_arguments(
        {"test", "--self-test", "--destroy-shard-count", "2",
         "--destroy-shard-index", "2"}));
  } catch (const std::runtime_error&) {
    invalid_arguments = true;
  }
  require_test(
      invalid_arguments,
      "CLI rejects a shard index equal to its count");
}

void run_self_tests() {
  self_test_normalization_and_determinant();
  self_test_pool_enumeration();
  self_test_graph_and_cliques();
  self_test_hadamard_fischer();
  self_test_bordered_gram_state();
  self_test_schur_residual_bound();
  self_test_sparse_schur_storage();
  self_test_destroy_schedules();
  std::cout
      << "self-test passed: normalization, exact determinant, pool, "
         "bit graph, clique enumeration/coloring, old and "
         "Schur-residual Fischer bounds, sparse edge-only "
         "Schur storage, bordered Gram/kernel determinants, and "
         "deterministic sharded destroy schedules\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    if (arguments.self_test) {
      run_self_tests();
      return 0;
    }
    stop_requested = 0;
    std::signal(SIGINT, request_stop);
    std::signal(SIGTERM, request_stop);
    SearchRunner runner(arguments);
    return runner.run();
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
