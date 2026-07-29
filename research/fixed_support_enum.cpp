// Exhaustive fixed-support switch enumerator for order-23 MaxDet matrices.
//
// The search fixes nine row labels and nine column labels, with labeled
// degrees (2,2,2,1,1,1,1,1,1) on both sides, and enumerates every simple
// bipartite mask with those degrees.  Every selected entry is sign-flipped.
//
// Scores are exact without a floating-point gate.  For every +/-1 matrix of
// order 23, det(A) is divisible by 2^22.  Modulo p = 2^32 - 5, the quotient
// det(A)/2^22 is unique inside the order-23 Hadamard bound because
//
//   2 * floor(sqrt(23^23) / 2^22) < p.
//
// Thus one modular determinant recovers the signed integer determinant.  A
// rank-nine determinant lemma and a 3+6 exterior-product meet-in-the-middle
// make the complete 6,508,620-mask search practical.

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <chrono>
#include <csignal>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr int kOrder = 23;
constexpr int kSupport = 9;
constexpr int kDegreeTwo = 3;
constexpr int kDegreeOne = 6;
constexpr int kFlipCount = 12;
constexpr std::uint64_t kTwo22 = UINT64_C(1) << 22;
constexpr std::uint64_t kPrime = UINT64_C(4294967291);  // 2^32 - 5, prime.

std::atomic<bool> g_stop{false};

void on_signal(int) {
  g_stop.store(true, std::memory_order_relaxed);
}

std::uint64_t add_mod(std::uint64_t a, std::uint64_t b) {
  const std::uint64_t s = a + b;
  return s >= kPrime ? s - kPrime : s;
}

std::uint64_t sub_mod(std::uint64_t a, std::uint64_t b) {
  return a >= b ? a - b : a + kPrime - b;
}

// Fast reduction modulo 2^32 - 5.  Inputs to mul_mod are below p, so their
// product fits in uint64_t.  Two pseudo-Mersenne folds suffice.
std::uint64_t reduce_mod(std::uint64_t x) {
  std::uint64_t r = static_cast<std::uint32_t>(x) + 5 * (x >> 32);
  r = static_cast<std::uint32_t>(r) + 5 * (r >> 32);
  if (r >= kPrime) {
    r -= kPrime;
  }
  return r;
}

std::uint64_t mul_mod(std::uint64_t a, std::uint64_t b) {
  return reduce_mod(a * b);
}

std::uint64_t pow_mod(std::uint64_t a, std::uint64_t e) {
  std::uint64_t result = 1;
  while (e != 0) {
    if ((e & 1U) != 0) {
      result = mul_mod(result, a);
    }
    a = mul_mod(a, a);
    e >>= 1;
  }
  return result;
}

bool is_prime_u32(std::uint64_t n) {
  if (n < 2) {
    return false;
  }
  if ((n & 1U) == 0) {
    return n == 2;
  }
  for (std::uint64_t d = 3; d * d <= n; d += 2) {
    if (n % d == 0) {
      return false;
    }
  }
  return true;
}

using U128 = unsigned __int128;

U128 pow_u128(std::uint64_t base, int exponent) {
  U128 result = 1;
  for (int i = 0; i < exponent; ++i) {
    result *= base;
  }
  return result;
}

std::uint64_t isqrt_u128(U128 value) {
  std::uint64_t lo = 0;
  std::uint64_t hi = UINT64_C(1) << 63;
  while (lo + 1 < hi) {
    const std::uint64_t mid = lo + (hi - lo) / 2;
    if (static_cast<U128>(mid) * mid <= value) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return lo;
}

std::string json_escape(std::string_view input) {
  std::ostringstream out;
  for (const unsigned char ch : input) {
    switch (ch) {
      case '"':
        out << "\\\"";
        break;
      case '\\':
        out << "\\\\";
        break;
      case '\b':
        out << "\\b";
        break;
      case '\f':
        out << "\\f";
        break;
      case '\n':
        out << "\\n";
        break;
      case '\r':
        out << "\\r";
        break;
      case '\t':
        out << "\\t";
        break;
      default:
        if (ch < 0x20) {
          out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
              << static_cast<int>(ch) << std::dec;
        } else {
          out << static_cast<char>(ch);
        }
    }
  }
  return out.str();
}

std::string u128_string(U128 value) {
  if (value == 0) {
    return "0";
  }
  std::string result;
  while (value != 0) {
    result.push_back(static_cast<char>('0' + value % 10));
    value /= 10;
  }
  std::reverse(result.begin(), result.end());
  return result;
}

using Matrix = std::array<std::array<int, kOrder>, kOrder>;
using ModMatrix =
    std::array<std::array<std::uint64_t, kOrder>, kOrder>;
using Vec9 = std::array<std::uint64_t, kSupport>;

Matrix read_matrix(const fs::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("cannot open matrix: " + path.string());
  }
  Matrix matrix{};
  std::string line;
  for (int row = 0; row < kOrder; ++row) {
    if (!std::getline(in, line)) {
      throw std::runtime_error("matrix has fewer than 23 rows: " +
                               path.string());
    }
    std::istringstream row_in(line);
    for (int col = 0; col < kOrder; ++col) {
      if (!(row_in >> matrix[row][col]) ||
          (matrix[row][col] != -1 && matrix[row][col] != 1)) {
        throw std::runtime_error("invalid matrix entry: " + path.string());
      }
    }
    std::string extra;
    if (row_in >> extra) {
      throw std::runtime_error("matrix row has more than 23 entries: " +
                               path.string());
    }
  }
  while (std::getline(in, line)) {
    if (line.find_first_not_of(" \t\r\n") != std::string::npos) {
      throw std::runtime_error("matrix has more than 23 rows: " +
                               path.string());
    }
  }
  return matrix;
}

void write_atomic(const fs::path& destination, const std::string& content) {
  const fs::path temporary = destination.string() + ".tmp";
  {
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    if (!out) {
      throw std::runtime_error("cannot create temporary output: " +
                               temporary.string());
    }
    out << content;
    out.flush();
    if (!out) {
      throw std::runtime_error("failed writing temporary output: " +
                               temporary.string());
    }
  }
  fs::rename(temporary, destination);
}

std::string matrix_text(const Matrix& base,
                        const std::array<std::array<bool, kSupport>, kSupport>&
                            mask,
                        const std::array<int, kSupport>& rows,
                        const std::array<int, kSupport>& cols) {
  Matrix result = base;
  for (int i = 0; i < kSupport; ++i) {
    for (int j = 0; j < kSupport; ++j) {
      if (mask[i][j]) {
        result[rows[i]][cols[j]] *= -1;
      }
    }
  }
  std::ostringstream out;
  for (int i = 0; i < kOrder; ++i) {
    for (int j = 0; j < kOrder; ++j) {
      if (j != 0) {
        out << ' ';
      }
      out << result[i][j];
    }
    out << '\n';
  }
  return out.str();
}

std::string mask_text(
    const std::array<std::array<bool, kSupport>, kSupport>& mask,
    const std::array<int, kSupport>& rows,
    const std::array<int, kSupport>& cols) {
  std::vector<std::pair<int, int>> coordinates;
  for (int i = 0; i < kSupport; ++i) {
    for (int j = 0; j < kSupport; ++j) {
      if (mask[i][j]) {
        coordinates.emplace_back(rows[i] + 1, cols[j] + 1);
      }
    }
  }
  std::sort(coordinates.begin(), coordinates.end());
  std::ostringstream out;
  for (const auto& [row, col] : coordinates) {
    out << row << ' ' << col << '\n';
  }
  return out.str();
}

struct InverseResult {
  std::uint64_t determinant = 0;
  ModMatrix inverse{};
};

InverseResult invert_matrix(const Matrix& matrix) {
  std::array<std::array<std::uint64_t, 2 * kOrder>, kOrder> work{};
  for (int i = 0; i < kOrder; ++i) {
    for (int j = 0; j < kOrder; ++j) {
      work[i][j] =
          matrix[i][j] == 1 ? 1 : static_cast<std::uint64_t>(kPrime - 1);
    }
    work[i][kOrder + i] = 1;
  }

  std::uint64_t determinant = 1;
  bool negative = false;
  for (int col = 0; col < kOrder; ++col) {
    int pivot = col;
    while (pivot < kOrder && work[pivot][col] == 0) {
      ++pivot;
    }
    if (pivot == kOrder) {
      throw std::runtime_error("base matrix is singular modulo exact prime");
    }
    if (pivot != col) {
      std::swap(work[pivot], work[col]);
      negative = !negative;
    }
    const std::uint64_t pivot_value = work[col][col];
    determinant = mul_mod(determinant, pivot_value);
    const std::uint64_t inverse_pivot = pow_mod(pivot_value, kPrime - 2);
    for (int j = 0; j < 2 * kOrder; ++j) {
      work[col][j] = mul_mod(work[col][j], inverse_pivot);
    }
    for (int row = 0; row < kOrder; ++row) {
      if (row == col || work[row][col] == 0) {
        continue;
      }
      const std::uint64_t factor = work[row][col];
      for (int j = 0; j < 2 * kOrder; ++j) {
        work[row][j] =
            sub_mod(work[row][j], mul_mod(factor, work[col][j]));
      }
    }
  }
  if (negative && determinant != 0) {
    determinant = kPrime - determinant;
  }

  InverseResult result;
  result.determinant = determinant;
  for (int i = 0; i < kOrder; ++i) {
    for (int j = 0; j < kOrder; ++j) {
      result.inverse[i][j] = work[i][kOrder + j];
    }
  }
  return result;
}

struct ExteriorTables {
  std::array<std::vector<std::uint16_t>, kSupport + 1> subsets;
  std::array<std::array<int, 1 << kSupport>, kSupport + 1> index{};

  ExteriorTables() {
    for (auto& by_degree : index) {
      by_degree.fill(-1);
    }
    for (std::uint16_t mask = 0; mask < (1U << kSupport); ++mask) {
      const int degree = __builtin_popcount(mask);
      index[degree][mask] = static_cast<int>(subsets[degree].size());
      subsets[degree].push_back(mask);
    }
  }
};

using Exterior = std::array<std::uint64_t, 126>;

Exterior wedge_append(const Exterior& previous,
                      int previous_degree,
                      const Vec9& row,
                      const ExteriorTables& tables) {
  Exterior next{};
  const int next_degree = previous_degree + 1;
  for (std::size_t out_index = 0;
       out_index < tables.subsets[next_degree].size(); ++out_index) {
    const std::uint16_t subset = tables.subsets[next_degree][out_index];
    std::uint64_t value = 0;
    int position = 0;
    for (int col = 0; col < kSupport; ++col) {
      if ((subset & (1U << col)) == 0) {
        continue;
      }
      const std::uint16_t smaller =
          static_cast<std::uint16_t>(subset ^ (1U << col));
      const int input_index = tables.index[previous_degree][smaller];
      const std::uint64_t term =
          mul_mod(previous[input_index], row[col]);
      if (((previous_degree + position) & 1) == 0) {
        value = add_mod(value, term);
      } else {
        value = sub_mod(value, term);
      }
      ++position;
    }
    next[out_index] = value;
  }
  return next;
}

struct SupportSpec {
  std::array<int, kSupport> labels{};
  std::array<int, kSupport> degrees{};
};

SupportSpec parse_support(std::string_view text, std::string_view kind) {
  std::vector<std::pair<int, int>> entries;
  std::size_t start = 0;
  while (start < text.size()) {
    const std::size_t comma = text.find(',', start);
    const std::string_view token =
        text.substr(start, comma == std::string_view::npos
                               ? text.size() - start
                               : comma - start);
    const std::size_t colon = token.find(':');
    if (colon == std::string_view::npos) {
      throw std::runtime_error(std::string(kind) +
                               " support token must be LABEL:DEGREE");
    }
    const int label = std::stoi(std::string(token.substr(0, colon)));
    const int degree = std::stoi(std::string(token.substr(colon + 1)));
    if (label < 1 || label > kOrder || (degree != 1 && degree != 2)) {
      throw std::runtime_error("invalid " + std::string(kind) +
                               " support entry");
    }
    entries.emplace_back(label - 1, degree);
    if (comma == std::string_view::npos) {
      break;
    }
    start = comma + 1;
  }
  if (entries.size() != kSupport) {
    throw std::runtime_error(std::string(kind) +
                             " support must contain exactly nine entries");
  }
  std::sort(entries.begin(), entries.end(),
            [](const auto& left, const auto& right) {
              if (left.second != right.second) {
                return left.second > right.second;
              }
              return left.first < right.first;
            });
  if (std::adjacent_find(entries.begin(), entries.end(),
                         [](const auto& left, const auto& right) {
                           return left.first == right.first;
                         }) != entries.end()) {
    throw std::runtime_error(std::string(kind) +
                             " support labels must be unique");
  }
  // Sorting by degree can separate equal labels, so also check via a bitmap.
  std::array<bool, kOrder> seen{};
  int degree_two = 0;
  int degree_sum = 0;
  SupportSpec result;
  for (int i = 0; i < kSupport; ++i) {
    const auto [label, degree] = entries[i];
    if (seen[label]) {
      throw std::runtime_error(std::string(kind) +
                               " support labels must be unique");
    }
    seen[label] = true;
    result.labels[i] = label;
    result.degrees[i] = degree;
    degree_two += degree == 2;
    degree_sum += degree;
  }
  if (degree_two != kDegreeTwo || degree_sum != kFlipCount) {
    throw std::runtime_error(std::string(kind) +
                             " degrees must be exactly 2,2,2,1,1,1,1,1,1");
  }
  return result;
}

struct BaseSpec {
  std::string label;
  fs::path path;
};

struct Options {
  SupportSpec rows;
  SupportSpec cols;
  bool have_rows = false;
  bool have_cols = false;
  std::vector<BaseSpec> bases;
  fs::path output_dir;
  bool have_output = false;
  std::uint64_t frontier = UINT64_C(2779447296000000);
  std::uint64_t shard_count = 1;
  std::uint64_t shard_index = 0;
  double heartbeat_seconds = 10.0;
  std::optional<std::uint64_t> limit;
};

void usage(const char* program) {
  std::cerr
      << "usage: " << program << " --rows R:D,... --cols C:D,...\\n"
      << "  --base LABEL=PATH [--base LABEL=PATH ...] --out-dir DIR\\n"
      << "  [--frontier N] [--shard-count N --shard-index I]\\n"
      << "  [--heartbeat-seconds N] [--limit N]\\n";
}

Options parse_options(int argc, char** argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto value = [&](std::string_view name) -> std::string {
      if (i + 1 >= argc) {
        throw std::runtime_error("missing value after " + std::string(name));
      }
      return argv[++i];
    };
    if (arg == "--rows") {
      options.rows = parse_support(value(arg), "row");
      options.have_rows = true;
    } else if (arg == "--cols") {
      options.cols = parse_support(value(arg), "column");
      options.have_cols = true;
    } else if (arg == "--base") {
      const std::string spec = value(arg);
      const std::size_t equals = spec.find('=');
      if (equals == std::string::npos || equals == 0 ||
          equals + 1 == spec.size()) {
        throw std::runtime_error("--base must be LABEL=PATH");
      }
      options.bases.push_back(
          {spec.substr(0, equals), fs::path(spec.substr(equals + 1))});
    } else if (arg == "--out-dir") {
      options.output_dir = value(arg);
      options.have_output = true;
    } else if (arg == "--frontier") {
      options.frontier = std::stoull(value(arg));
    } else if (arg == "--shard-count") {
      options.shard_count = std::stoull(value(arg));
    } else if (arg == "--shard-index") {
      options.shard_index = std::stoull(value(arg));
    } else if (arg == "--heartbeat-seconds") {
      options.heartbeat_seconds = std::stod(value(arg));
    } else if (arg == "--limit") {
      options.limit = std::stoull(value(arg));
    } else if (arg == "--help" || arg == "-h") {
      usage(argv[0]);
      std::exit(0);
    } else {
      throw std::runtime_error("unknown option: " + arg);
    }
  }
  if (!options.have_rows || !options.have_cols || options.bases.empty() ||
      !options.have_output) {
    throw std::runtime_error(
        "--rows, --cols, at least one --base, and --out-dir are required");
  }
  if (options.shard_count == 0 ||
      options.shard_index >= options.shard_count) {
    throw std::runtime_error("invalid shard index/count");
  }
  if (!(options.heartbeat_seconds > 0.0)) {
    throw std::runtime_error("--heartbeat-seconds must be positive");
  }
  std::set<std::string> labels;
  for (const auto& base : options.bases) {
    if (base.label.empty() || base.label.size() > 255) {
      throw std::runtime_error("base label must contain 1..255 characters");
    }
    for (const unsigned char ch : base.label) {
      if (!(std::isalnum(ch) || ch == '-' || ch == '_')) {
        throw std::runtime_error(
          "base label may contain only letters, digits, '-' and '_'");
      }
    }
    if (!labels.insert(base.label).second) {
      throw std::runtime_error("base labels must be unique");
    }
  }
  return options;
}

std::uint64_t encode_counts(
    const std::array<std::uint8_t, kSupport>& counts,
    const std::array<int, kSupport>& capacities) {
  std::uint64_t key = 0;
  std::uint64_t multiplier = 1;
  for (int i = 0; i < kSupport; ++i) {
    key += counts[i] * multiplier;
    multiplier *= static_cast<std::uint64_t>(capacities[i] + 1);
  }
  return key;
}

struct CountMemoKey {
  int row = 0;
  std::uint64_t remaining = 0;

  bool operator==(const CountMemoKey& other) const {
    return row == other.row && remaining == other.remaining;
  }
};

struct CountMemoHash {
  std::size_t operator()(const CountMemoKey& key) const {
    return std::hash<std::uint64_t>{}(
        key.remaining ^ (static_cast<std::uint64_t>(key.row) << 48));
  }
};

std::uint64_t count_masks_recursive(
    int row,
    const std::array<int, kSupport>& row_degrees,
    std::array<std::uint8_t, kSupport>& remaining,
    const std::array<int, kSupport>& capacities,
    std::unordered_map<CountMemoKey, std::uint64_t, CountMemoHash>& memo) {
  if (row == kSupport) {
    return std::all_of(remaining.begin(), remaining.end(),
                       [](std::uint8_t value) { return value == 0; })
               ? 1
               : 0;
  }
  const CountMemoKey memo_key{row, encode_counts(remaining, capacities)};
  if (const auto it = memo.find(memo_key); it != memo.end()) {
    return it->second;
  }
  std::uint64_t count = 0;
  if (row_degrees[row] == 1) {
    for (int col = 0; col < kSupport; ++col) {
      if (remaining[col] == 0) {
        continue;
      }
      --remaining[col];
      count += count_masks_recursive(row + 1, row_degrees, remaining,
                                     capacities, memo);
      ++remaining[col];
    }
  } else {
    for (int first = 0; first < kSupport; ++first) {
      if (remaining[first] == 0) {
        continue;
      }
      --remaining[first];
      for (int second = first + 1; second < kSupport; ++second) {
        if (remaining[second] == 0) {
          continue;
        }
        --remaining[second];
        count += count_masks_recursive(row + 1, row_degrees, remaining,
                                       capacities, memo);
        ++remaining[second];
      }
      ++remaining[first];
    }
  }
  memo.emplace(memo_key, count);
  return count;
}

std::uint64_t count_masks(const SupportSpec& rows, const SupportSpec& cols) {
  std::array<std::uint8_t, kSupport> remaining{};
  for (int i = 0; i < kSupport; ++i) {
    remaining[i] = static_cast<std::uint8_t>(cols.degrees[i]);
  }
  std::unordered_map<CountMemoKey, std::uint64_t, CountMemoHash> memo;
  return count_masks_recursive(0, rows.degrees, remaining, cols.degrees, memo);
}

struct SuffixRecord {
  std::array<std::uint32_t, 84> dual{};
  std::array<std::uint8_t, kDegreeOne> choices{};
};

using SuffixGroups = std::array<std::vector<SuffixRecord>, 1728>;

std::array<std::array<Vec9, 36>, kDegreeTwo> build_pair_rows(
    const Matrix& base,
    const InverseResult& inverse,
    const SupportSpec& rows,
    const SupportSpec& cols) {
  std::array<std::array<Vec9, 36>, kDegreeTwo> result{};
  for (int row_index = 0; row_index < kDegreeTwo; ++row_index) {
    int option = 0;
    for (int first = 0; first < kSupport; ++first) {
      for (int second = first + 1; second < kSupport; ++second) {
        Vec9 vector{};
        vector[row_index] = 1;
        for (const int col_index : {first, second}) {
          const int value =
              base[rows.labels[row_index]][cols.labels[col_index]];
          const std::uint64_t delta = value == 1 ? kPrime - 2 : 2;
          for (int k = 0; k < kSupport; ++k) {
            const std::uint64_t inverse_entry =
                inverse.inverse[cols.labels[col_index]][rows.labels[k]];
            vector[k] =
                add_mod(vector[k], mul_mod(delta, inverse_entry));
          }
        }
        result[row_index][option++] = vector;
      }
    }
  }
  return result;
}

std::array<std::array<Vec9, kSupport>, kDegreeOne> build_single_rows(
    const Matrix& base,
    const InverseResult& inverse,
    const SupportSpec& rows,
    const SupportSpec& cols) {
  std::array<std::array<Vec9, kSupport>, kDegreeOne> result{};
  for (int local_row = 0; local_row < kDegreeOne; ++local_row) {
    const int row_index = kDegreeTwo + local_row;
    for (int col_index = 0; col_index < kSupport; ++col_index) {
      Vec9 vector{};
      vector[row_index] = 1;
      const int value =
          base[rows.labels[row_index]][cols.labels[col_index]];
      const std::uint64_t delta = value == 1 ? kPrime - 2 : 2;
      for (int k = 0; k < kSupport; ++k) {
        const std::uint64_t inverse_entry =
            inverse.inverse[cols.labels[col_index]][rows.labels[k]];
        vector[k] = add_mod(vector[k], mul_mod(delta, inverse_entry));
      }
      result[local_row][col_index] = vector;
    }
  }
  return result;
}

struct SuffixBuilder {
  const std::array<std::array<Vec9, kSupport>, kDegreeOne>& row_vectors;
  const SupportSpec& cols;
  const ExteriorTables& tables;
  SuffixGroups& groups;
  std::array<std::uint8_t, kSupport> counts{};
  std::array<std::uint8_t, kDegreeOne> choices{};
  std::uint64_t records = 0;

  void recurse(int depth, const Exterior& exterior) {
    if (depth == kDegreeOne) {
      SuffixRecord record;
      record.choices = choices;
      constexpr std::uint16_t all = (1U << kSupport) - 1;
      for (std::size_t prefix_index = 0;
           prefix_index < tables.subsets[3].size(); ++prefix_index) {
        const std::uint16_t prefix_subset =
            tables.subsets[3][prefix_index];
        const std::uint16_t complement =
            static_cast<std::uint16_t>(all ^ prefix_subset);
        const int suffix_index = tables.index[6][complement];
        std::uint64_t value = exterior[suffix_index];
        int inversions = 0;
        int rank = 0;
        for (int col = 0; col < kSupport; ++col) {
          if ((prefix_subset & (1U << col)) != 0) {
            inversions += col - rank;
            ++rank;
          }
        }
        if ((inversions & 1) != 0 && value != 0) {
          value = kPrime - value;
        }
        record.dual[prefix_index] = static_cast<std::uint32_t>(value);
      }
      groups[encode_counts(counts, cols.degrees)].push_back(record);
      ++records;
      return;
    }
    for (int col = 0; col < kSupport; ++col) {
      if (counts[col] >= cols.degrees[col]) {
        continue;
      }
      ++counts[col];
      choices[depth] = static_cast<std::uint8_t>(col);
      const Exterior next =
          wedge_append(exterior, depth, row_vectors[depth][col], tables);
      recurse(depth + 1, next);
      --counts[col];
    }
  }
};

std::pair<int, int> pair_from_option(int option) {
  int current = 0;
  for (int first = 0; first < kSupport; ++first) {
    for (int second = first + 1; second < kSupport; ++second) {
      if (current == option) {
        return {first, second};
      }
      ++current;
    }
  }
  throw std::runtime_error("internal pair-option index error");
}

std::array<std::array<bool, kSupport>, kSupport> make_mask(
    const std::array<std::uint8_t, kDegreeTwo>& prefix_options,
    const std::array<std::uint8_t, kDegreeOne>& suffix_choices) {
  std::array<std::array<bool, kSupport>, kSupport> mask{};
  for (int row = 0; row < kDegreeTwo; ++row) {
    const auto [first, second] = pair_from_option(prefix_options[row]);
    mask[row][first] = true;
    mask[row][second] = true;
  }
  for (int local_row = 0; local_row < kDegreeOne; ++local_row) {
    mask[kDegreeTwo + local_row][suffix_choices[local_row]] = true;
  }
  return mask;
}

std::uint64_t absolute_u64(std::int64_t value) {
  return value < 0 ? static_cast<std::uint64_t>(-value)
                   : static_cast<std::uint64_t>(value);
}

struct Logger {
  std::ofstream stream;

  explicit Logger(const fs::path& path) : stream(path) {
    if (!stream) {
      throw std::runtime_error("cannot create log: " + path.string());
    }
  }

  void line(const std::string& text) {
    stream << text << '\n';
    stream.flush();
    if (!stream) {
      throw std::runtime_error("failed writing JSONL log");
    }
  }
};

std::string support_json(const SupportSpec& spec) {
  std::ostringstream out;
  out << '[';
  for (int i = 0; i < kSupport; ++i) {
    if (i != 0) {
      out << ',';
    }
    out << "{\"label\":" << spec.labels[i] + 1
        << ",\"degree\":" << spec.degrees[i] << '}';
  }
  out << ']';
  return out.str();
}

struct EnumerationResult {
  std::uint64_t examined = 0;
  std::uint64_t prefixes = 0;
  std::uint64_t ties = 0;
  std::uint64_t wins = 0;
  std::uint64_t best_score = 0;
  std::int64_t best_determinant = 0;
  std::array<std::uint8_t, kDegreeTwo> best_prefix{};
  std::array<std::uint8_t, kDegreeOne> best_suffix{};
  bool have_best = false;
  bool complete = true;
};

EnumerationResult enumerate_base(
    const Options& options,
    const BaseSpec& base_spec,
    Logger& logger,
    const ExteriorTables& tables,
    std::uint64_t inverse_two22,
    std::uint64_t quotient_bound,
    std::uint64_t expected_space) {
  const Matrix base = read_matrix(base_spec.path);
  const InverseResult inverse = invert_matrix(base);
  const std::uint64_t base_q_mod =
      mul_mod(inverse.determinant, inverse_two22);
  const std::int64_t base_q =
      base_q_mod <= kPrime / 2
          ? static_cast<std::int64_t>(base_q_mod)
          : static_cast<std::int64_t>(base_q_mod) -
                static_cast<std::int64_t>(kPrime);
  if (absolute_u64(base_q) > quotient_bound) {
    throw std::runtime_error(
        "base determinant residue violates the Hadamard quotient bound");
  }
  const std::int64_t base_determinant =
      base_q * static_cast<std::int64_t>(kTwo22);
  const std::uint64_t base_score = absolute_u64(base_determinant);

  {
    std::ostringstream event;
    event << "{\"event\":\"base_started\",\"label\":\""
          << json_escape(base_spec.label) << "\",\"path\":\""
          << json_escape(base_spec.path.string())
          << "\",\"base_determinant\":\"" << base_determinant
          << "\",\"base_score\":\"" << base_score
          << "\",\"determinant_mod_prime\":" << inverse.determinant << '}';
    logger.line(event.str());
  }

  const auto pair_rows =
      build_pair_rows(base, inverse, options.rows, options.cols);
  const auto single_rows =
      build_single_rows(base, inverse, options.rows, options.cols);

  SuffixGroups suffix_groups;
  Exterior identity{};
  identity[0] = 1;
  SuffixBuilder suffix_builder{single_rows, options.cols, tables,
                               suffix_groups};
  suffix_builder.recurse(0, identity);

  {
    std::uint64_t nonempty_groups = 0;
    std::uint64_t largest_group = 0;
    for (const auto& group : suffix_groups) {
      if (!group.empty()) {
        ++nonempty_groups;
        largest_group = std::max<std::uint64_t>(largest_group, group.size());
      }
    }
    std::ostringstream event;
    event << "{\"event\":\"suffix_table_built\",\"label\":\""
          << json_escape(base_spec.label) << "\",\"records\":"
          << suffix_builder.records << ",\"nonempty_groups\":"
          << nonempty_groups << ",\"largest_group\":" << largest_group
          << '}';
    logger.line(event.str());
  }

  const fs::path base_dir = options.output_dir / base_spec.label;
  fs::create_directory(base_dir);
  const fs::path tie_dir = base_dir / "ties";
  const fs::path win_dir = base_dir / "wins";
  fs::create_directory(tie_dir);
  fs::create_directory(win_dir);

  EnumerationResult result;
  const auto started = std::chrono::steady_clock::now();
  auto next_heartbeat =
      started + std::chrono::duration<double>(options.heartbeat_seconds);
  std::uint64_t valid_prefix_ordinal = 0;
  bool hit_limit = false;

  for (int first_option = 0; first_option < 36 && !g_stop.load(); ++first_option) {
    const auto first_pair = pair_from_option(first_option);
    std::array<std::uint8_t, kSupport> counts{};
    ++counts[first_pair.first];
    ++counts[first_pair.second];
    Exterior exterior1{};
    for (int col = 0; col < kSupport; ++col) {
      exterior1[tables.index[1][1U << col]] =
          pair_rows[0][first_option][col];
    }
    for (int second_option = 0;
         second_option < 36 && !g_stop.load(); ++second_option) {
      const auto second_pair = pair_from_option(second_option);
      if (counts[second_pair.first] >= options.cols.degrees[second_pair.first] ||
          counts[second_pair.second] >=
              options.cols.degrees[second_pair.second]) {
        continue;
      }
      ++counts[second_pair.first];
      ++counts[second_pair.second];
      const Exterior exterior2 =
          wedge_append(exterior1, 1, pair_rows[1][second_option], tables);
      for (int third_option = 0;
           third_option < 36 && !g_stop.load(); ++third_option) {
        const auto third_pair = pair_from_option(third_option);
        if (counts[third_pair.first] >= options.cols.degrees[third_pair.first] ||
            counts[third_pair.second] >=
                options.cols.degrees[third_pair.second]) {
          continue;
        }
        ++counts[third_pair.first];
        ++counts[third_pair.second];
        const std::uint64_t this_prefix = valid_prefix_ordinal++;
        if (this_prefix % options.shard_count != options.shard_index) {
          --counts[third_pair.first];
          --counts[third_pair.second];
          continue;
        }
        ++result.prefixes;
        std::array<std::uint8_t, kSupport> remaining{};
        for (int col = 0; col < kSupport; ++col) {
          remaining[col] = static_cast<std::uint8_t>(
              options.cols.degrees[col] - counts[col]);
        }
        const std::uint64_t key =
            encode_counts(remaining, options.cols.degrees);
        const Exterior exterior3 =
            wedge_append(exterior2, 2, pair_rows[2][third_option], tables);
        const std::array<std::uint8_t, kDegreeTwo> prefix_options{
            static_cast<std::uint8_t>(first_option),
            static_cast<std::uint8_t>(second_option),
            static_cast<std::uint8_t>(third_option)};

        for (const SuffixRecord& suffix : suffix_groups[key]) {
          std::uint64_t small_determinant = 0;
          for (int index = 0; index < 84; ++index) {
            small_determinant =
                add_mod(small_determinant,
                        mul_mod(exterior3[index], suffix.dual[index]));
          }
          const std::uint64_t full_mod =
              mul_mod(inverse.determinant, small_determinant);
          const std::uint64_t q_mod = mul_mod(full_mod, inverse_two22);
          const std::int64_t quotient =
              q_mod <= kPrime / 2
                  ? static_cast<std::int64_t>(q_mod)
                  : static_cast<std::int64_t>(q_mod) -
                        static_cast<std::int64_t>(kPrime);
          if (absolute_u64(quotient) > quotient_bound) {
            throw std::runtime_error(
                "candidate determinant residue violates Hadamard bound");
          }
          const std::int64_t determinant =
              quotient * static_cast<std::int64_t>(kTwo22);
          const std::uint64_t score = absolute_u64(determinant);
          ++result.examined;

          if (!result.have_best || score > result.best_score) {
            result.have_best = true;
            result.best_score = score;
            result.best_determinant = determinant;
            result.best_prefix = prefix_options;
            result.best_suffix = suffix.choices;
          }

          const bool tie = score == base_score;
          const bool win = score > options.frontier;
          if (tie || win) {
            std::uint64_t ordinal = 0;
            fs::path directory;
            std::string kind;
            if (tie) {
              ordinal = ++result.ties;
              directory = tie_dir;
              kind = "tie";
            }
            if (win) {
              ordinal = ++result.wins;
              directory = win_dir;
              kind = "win";
            }
            const auto mask = make_mask(prefix_options, suffix.choices);
            std::ostringstream stem;
            stem << kind << '-' << std::setw(6) << std::setfill('0')
                 << ordinal;
            const fs::path matrix_path =
                directory / (stem.str() + ".matrix.txt");
            const fs::path mask_path =
                directory / (stem.str() + ".mask.txt");
            write_atomic(matrix_path,
                         matrix_text(base, mask, options.rows.labels,
                                     options.cols.labels));
            write_atomic(mask_path,
                         mask_text(mask, options.rows.labels,
                                   options.cols.labels));
            std::ostringstream event;
            event << "{\"event\":\"" << kind << "\",\"label\":\""
                  << json_escape(base_spec.label)
                  << "\",\"candidate_index\":" << result.examined
                  << ",\"determinant\":\"" << determinant
                  << "\",\"score\":\"" << score
                  << "\",\"matrix\":\""
                  << json_escape(matrix_path.string()) << "\",\"mask\":\""
                  << json_escape(mask_path.string()) << "\"}";
            logger.line(event.str());
          }

          if (options.limit && result.examined >= *options.limit) {
            hit_limit = true;
            break;
          }
          const auto now = std::chrono::steady_clock::now();
          if (now >= next_heartbeat) {
            const double elapsed =
                std::chrono::duration<double>(now - started).count();
            std::ostringstream event;
            event << "{\"event\":\"heartbeat\",\"label\":\""
                  << json_escape(base_spec.label)
                  << "\",\"examined\":" << result.examined
                  << ",\"prefixes\":" << result.prefixes
                  << ",\"ties\":" << result.ties << ",\"wins\":"
                  << result.wins << ",\"best_score\":\""
                  << result.best_score << "\",\"elapsed_seconds\":"
                  << std::fixed << std::setprecision(6) << elapsed << '}';
            logger.line(event.str());
            next_heartbeat =
                now + std::chrono::duration<double>(
                          options.heartbeat_seconds);
          }
        }
        --counts[third_pair.first];
        --counts[third_pair.second];
        if (hit_limit) {
          break;
        }
      }
      --counts[second_pair.first];
      --counts[second_pair.second];
      if (hit_limit) {
        break;
      }
    }
    if (hit_limit) {
      break;
    }
  }

  result.complete = !g_stop.load() && !hit_limit;
  if (result.have_best) {
    const auto best_mask =
        make_mask(result.best_prefix, result.best_suffix);
    write_atomic(base_dir / "best.matrix.txt",
                 matrix_text(base, best_mask, options.rows.labels,
                             options.cols.labels));
    write_atomic(base_dir / "best.mask.txt",
                 mask_text(best_mask, options.rows.labels,
                           options.cols.labels));
  }

  const auto stopped = std::chrono::steady_clock::now();
  const double elapsed =
      std::chrono::duration<double>(stopped - started).count();
  std::ostringstream event;
  event << "{\"event\":\"base_finished\",\"label\":\""
        << json_escape(base_spec.label) << "\",\"complete\":"
        << (result.complete ? "true" : "false")
        << ",\"examined\":" << result.examined
        << ",\"expected_full_space\":" << expected_space
        << ",\"prefixes\":" << result.prefixes << ",\"ties\":"
        << result.ties << ",\"wins\":" << result.wins
        << ",\"best_determinant\":\"" << result.best_determinant
        << "\",\"best_score\":\"" << result.best_score
        << "\",\"elapsed_seconds\":" << std::fixed << std::setprecision(6)
        << elapsed << '}';
  logger.line(event.str());
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse_options(argc, argv);
    if (fs::exists(options.output_dir)) {
      throw std::runtime_error("output directory already exists: " +
                               options.output_dir.string());
    }
    fs::create_directories(options.output_dir);
    Logger logger(options.output_dir / "run.jsonl");
    std::signal(SIGINT, on_signal);
    std::signal(SIGTERM, on_signal);

    if (!is_prime_u32(kPrime)) {
      throw std::runtime_error("configured exact modulus is not prime");
    }
    const U128 hadamard_squared = pow_u128(kOrder, kOrder);
    const std::uint64_t hadamard_bound = isqrt_u128(hadamard_squared);
    const std::uint64_t quotient_bound = hadamard_bound / kTwo22;
    if (2 * quotient_bound >= kPrime) {
      throw std::runtime_error(
          "exact residue recovery is ambiguous inside Hadamard bound");
    }
    const std::uint64_t inverse_two22 =
        pow_mod(kTwo22 % kPrime, kPrime - 2);
    const std::uint64_t exact_space = count_masks(options.rows, options.cols);

    {
      std::ostringstream event;
      event << "{\"event\":\"started\",\"algorithm\":"
            << "\"fixed-support-exterior-mitm-v1\""
            << ",\"prime\":" << kPrime
            << ",\"prime_verified\":true,\"divisor\":" << kTwo22
            << ",\"hadamard_squared\":\""
            << u128_string(hadamard_squared)
            << "\",\"hadamard_floor\":\"" << hadamard_bound
            << "\",\"quotient_bound\":" << quotient_bound
            << ",\"twice_quotient_bound\":" << 2 * quotient_bound
            << ",\"exact_recovery\":"
            << "\"centered det/2^22 residue is unique under Hadamard bound\""
            << ",\"row_support\":" << support_json(options.rows)
            << ",\"column_support\":" << support_json(options.cols)
            << ",\"exact_full_space\":" << exact_space
            << ",\"shard_count\":" << options.shard_count
            << ",\"shard_index\":" << options.shard_index
            << ",\"frontier\":\"" << options.frontier << "\"";
      if (options.limit) {
        event << ",\"limit\":" << *options.limit;
      }
      event << '}';
      logger.line(event.str());
    }

    ExteriorTables tables;
    std::uint64_t total_examined = 0;
    std::uint64_t total_ties = 0;
    std::uint64_t total_wins = 0;
    std::uint64_t bases_completed = 0;
    bool complete = true;
    for (const BaseSpec& base : options.bases) {
      if (g_stop.load()) {
        complete = false;
        break;
      }
      const EnumerationResult result =
          enumerate_base(options, base, logger, tables, inverse_two22,
                         quotient_bound, exact_space);
      total_examined += result.examined;
      total_ties += result.ties;
      total_wins += result.wins;
      ++bases_completed;
      complete = complete && result.complete;
    }
    std::ostringstream event;
    event << "{\"event\":\"finished\",\"complete\":"
          << (complete ? "true" : "false")
          << ",\"bases_completed\":" << bases_completed
          << ",\"total_examined\":" << total_examined
          << ",\"total_ties\":" << total_ties
          << ",\"total_wins\":" << total_wins << '}';
    logger.line(event.str());
    return complete ? 0 : 2;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    usage(argv[0]);
    return 1;
  }
}
