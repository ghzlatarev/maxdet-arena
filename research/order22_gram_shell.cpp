#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

// Enumerate the exact normalized sign-column shell of an order-22 factor.
//
// If G = R R^T and R is nonsingular, every sign column s of R satisfies
//
//                  s^T adj(G) s = det(G).
//
// We fix s_0 = +1, enumerate the remaining 2^21 sign vectors in Gray-code
// order, and export every exact shell mask. Determinants are reconstructed
// modulo four 31-bit primes; the combined modulus exceeds 2e37. The known
// D-optimal Gram determinant and cofactors are below 1e30, and the exact
// G*adj(G)=det(G)I check independently catches a bad reconstruction.

namespace {

namespace fs = std::filesystem;

constexpr int kOrder = 22;
constexpr std::uint64_t kAssignments = std::uint64_t{1} << (kOrder - 1);
using Wide = __int128_t;
using UnsignedWide = __uint128_t;
using SignMatrix = std::array<std::array<int, kOrder>, kOrder>;

template <std::size_t Order>
using ExactMatrix = std::array<std::array<Wide, Order>, Order>;

struct Options {
  fs::path matrix;
  fs::path output;
};

std::string wide_to_string(Wide value) {
  if (value == 0) return "0";
  const bool negative = value < 0;
  if (negative) value = -value;
  std::string result;
  while (value != 0) {
    result.push_back(
        static_cast<char>('0' + static_cast<int>(value % 10)));
    value /= 10;
  }
  if (negative) result.push_back('-');
  std::reverse(result.begin(), result.end());
  return result;
}

std::string json_escape(const std::string& value) {
  std::string result;
  for (const char character : value) {
    switch (character) {
      case '\\':
        result += "\\\\";
        break;
      case '"':
        result += "\\\"";
        break;
      case '\n':
        result += "\\n";
        break;
      case '\r':
        result += "\\r";
        break;
      case '\t':
        result += "\\t";
        break;
      default:
        if (static_cast<unsigned char>(character) < 0x20U) {
          throw std::runtime_error("control character in JSON string");
        }
        result.push_back(character);
    }
  }
  return result;
}

std::uint64_t modular_power(
    std::uint64_t base,
    std::uint64_t exponent,
    std::uint64_t modulus) {
  std::uint64_t result = 1;
  while (exponent != 0) {
    if ((exponent & 1U) != 0U) {
      result = (result * base) % modulus;
    }
    base = (base * base) % modulus;
    exponent >>= 1U;
  }
  return result;
}

template <std::size_t Order>
std::uint64_t determinant_modulo(
    const ExactMatrix<Order>& source,
    std::uint64_t prime) {
  std::array<std::array<std::uint64_t, Order>, Order> work{};
  for (std::size_t row = 0; row < Order; ++row) {
    for (std::size_t column = 0; column < Order; ++column) {
      Wide residue = source[row][column] % static_cast<Wide>(prime);
      if (residue < 0) residue += prime;
      work[row][column] = static_cast<std::uint64_t>(residue);
    }
  }

  std::uint64_t determinant = 1;
  for (std::size_t column = 0; column < Order; ++column) {
    std::size_t pivot = column;
    while (pivot < Order && work[pivot][column] == 0) ++pivot;
    if (pivot == Order) return 0;
    if (pivot != column) {
      std::swap(work[pivot], work[column]);
      determinant = prime - determinant;
    }
    const std::uint64_t pivot_value = work[column][column];
    determinant = (determinant * pivot_value) % prime;
    const std::uint64_t inverse =
        modular_power(pivot_value, prime - 2U, prime);
    for (std::size_t row = column + 1; row < Order; ++row) {
      if (work[row][column] == 0) continue;
      const std::uint64_t factor =
          (work[row][column] * inverse) % prime;
      for (std::size_t inner = column + 1; inner < Order; ++inner) {
        const std::uint64_t product =
            (factor * work[column][inner]) % prime;
        work[row][inner] =
            work[row][inner] >= product
                ? work[row][inner] - product
                : work[row][inner] + prime - product;
      }
      work[row][column] = 0;
    }
  }
  return determinant;
}

template <std::size_t Order>
Wide exact_determinant(const ExactMatrix<Order>& matrix) {
  static constexpr std::array<std::uint64_t, 4> primes{
      2'147'483'647ULL,
      2'147'483'629ULL,
      2'147'483'587ULL,
      2'147'483'579ULL};
  UnsignedWide reconstructed = 0;
  UnsignedWide modulus = 1;
  for (const std::uint64_t prime : primes) {
    const std::uint64_t residue =
        determinant_modulo(matrix, prime);
    const std::uint64_t current =
        static_cast<std::uint64_t>(reconstructed % prime);
    const std::uint64_t difference =
        residue >= current ? residue - current
                           : residue + prime - current;
    const std::uint64_t inverse = modular_power(
        static_cast<std::uint64_t>(modulus % prime),
        prime - 2U, prime);
    const std::uint64_t multiplier =
        (difference * inverse) % prime;
    reconstructed += modulus * multiplier;
    modulus *= prime;
  }
  if (reconstructed > modulus / 2U) {
    return static_cast<Wide>(reconstructed) -
           static_cast<Wide>(modulus);
  }
  return static_cast<Wide>(reconstructed);
}

SignMatrix read_sign_matrix(const fs::path& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("cannot open matrix: " + path.string());
  }
  SignMatrix matrix{};
  for (auto& row : matrix) {
    for (int& value : row) {
      if (!(input >> value) || (value != -1 && value != 1)) {
        throw std::runtime_error(
            "matrix must contain exactly 22x22 signs");
      }
    }
  }
  std::string extra;
  if (input >> extra) {
    throw std::runtime_error("matrix contains extra data");
  }
  return matrix;
}

ExactMatrix<kOrder> exact_sign_matrix(const SignMatrix& matrix) {
  ExactMatrix<kOrder> result{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      result[row][column] = matrix[row][column];
    }
  }
  return result;
}

ExactMatrix<kOrder> gram(const SignMatrix& factor) {
  ExactMatrix<kOrder> result{};
  for (int row = 0; row < kOrder; ++row) {
    for (int other = row; other < kOrder; ++other) {
      Wide value = 0;
      for (int column = 0; column < kOrder; ++column) {
        value += static_cast<Wide>(factor[row][column]) *
                 factor[other][column];
      }
      result[row][other] = value;
      result[other][row] = value;
    }
  }
  return result;
}

ExactMatrix<kOrder> exact_adjugate(
    const ExactMatrix<kOrder>& matrix) {
  ExactMatrix<kOrder> adjugate{};
  for (int removed_row = 0; removed_row < kOrder; ++removed_row) {
    for (int removed_column = 0; removed_column < kOrder;
         ++removed_column) {
      ExactMatrix<kOrder - 1> minor{};
      int target_row = 0;
      for (int row = 0; row < kOrder; ++row) {
        if (row == removed_row) continue;
        int target_column = 0;
        for (int column = 0; column < kOrder; ++column) {
          if (column == removed_column) continue;
          minor[target_row][target_column++] = matrix[row][column];
        }
        ++target_row;
      }
      Wide cofactor = exact_determinant(minor);
      if ((removed_row + removed_column) % 2 != 0) {
        cofactor = -cofactor;
      }
      adjugate[removed_column][removed_row] = cofactor;
    }
  }
  return adjugate;
}

void validate_adjugate(
    const ExactMatrix<kOrder>& matrix,
    const ExactMatrix<kOrder>& adjugate,
    Wide determinant) {
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      Wide value = 0;
      for (int inner = 0; inner < kOrder; ++inner) {
        value += matrix[row][inner] * adjugate[inner][column];
      }
      const Wide expected = row == column ? determinant : 0;
      if (value != expected) {
        throw std::runtime_error("Gram-adjugate identity failed");
      }
    }
  }
}

std::vector<std::uint32_t> enumerate_shell(
    const ExactMatrix<kOrder>& adjugate,
    Wide determinant) {
  std::array<int, kOrder> signs{};
  signs.fill(1);
  std::array<Wide, kOrder> weighted{};
  for (int row = 0; row < kOrder; ++row) {
    for (int column = 0; column < kOrder; ++column) {
      weighted[row] += adjugate[row][column];
    }
  }
  Wide quadratic = 0;
  for (int row = 0; row < kOrder; ++row) {
    quadratic += signs[row] * weighted[row];
  }

  std::uint32_t mask = (std::uint32_t{1} << kOrder) - 1U;
  std::vector<std::uint32_t> shell;
  for (std::uint64_t assignment = 0; assignment < kAssignments;
       ++assignment) {
    if (quadratic == determinant) shell.push_back(mask);
    if (assignment + 1 == kAssignments) break;

    const unsigned bit = std::countr_zero(assignment + 1U) + 1U;
    const Wide delta = -2 * signs[bit];
    quadratic +=
        2 * delta * weighted[bit] +
        delta * delta * adjugate[bit][bit];
    for (int row = 0; row < kOrder; ++row) {
      weighted[row] += delta * adjugate[row][bit];
    }
    signs[bit] = -signs[bit];
    mask ^= std::uint32_t{1} << bit;
  }
  return shell;
}

std::vector<std::uint32_t> normalized_factor_columns(
    const SignMatrix& factor) {
  std::vector<std::uint32_t> masks;
  for (int column = 0; column < kOrder; ++column) {
    const int switch_sign = factor[0][column];
    std::uint32_t mask = 0;
    for (int row = 0; row < kOrder; ++row) {
      if (factor[row][column] * switch_sign == 1) {
        mask |= std::uint32_t{1} << row;
      }
    }
    masks.push_back(mask);
  }
  return masks;
}

void atomic_write(const fs::path& path, const std::string& contents) {
  if (!path.parent_path().empty()) {
    fs::create_directories(path.parent_path());
  }
  if (fs::exists(path) && !fs::is_regular_file(path)) {
    throw std::runtime_error("output path is not a regular file");
  }
  const fs::path temporary = path.string() + ".tmp";
  {
    std::ofstream output(temporary, std::ios::trunc);
    if (!output) throw std::runtime_error("cannot create output");
    output << contents;
    output.flush();
    if (!output) throw std::runtime_error("cannot flush output");
  }
  std::error_code error;
  fs::rename(temporary, path, error);
  if (error) {
    fs::remove(temporary);
    throw std::runtime_error("cannot install output");
  }
}

Options parse_options(int argc, char** argv) {
  Options options;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (option == "--help") {
      std::cout
          << "usage: order22_gram_shell --matrix MATRIX --output JSON\n";
      std::exit(0);
    }
    if (index + 1 >= argc) {
      throw std::runtime_error("missing value for " + option);
    }
    const fs::path value = argv[++index];
    if (option == "--matrix") {
      options.matrix = value;
    } else if (option == "--output") {
      options.output = value;
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (options.matrix.empty() || options.output.empty()) {
    throw std::runtime_error("--matrix and --output are required");
  }
  options.matrix = fs::absolute(options.matrix).lexically_normal();
  options.output = fs::absolute(options.output).lexically_normal();
  if (options.matrix == options.output) {
    throw std::runtime_error("matrix and output paths must differ");
  }
  return options;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse_options(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    const SignMatrix factor = read_sign_matrix(options.matrix);
    const Wide factor_determinant =
        exact_determinant(exact_sign_matrix(factor));
    if (factor_determinant == 0) {
      throw std::runtime_error("factor is singular");
    }
    const ExactMatrix<kOrder> factor_gram = gram(factor);
    const Wide gram_determinant = exact_determinant(factor_gram);
    if (gram_determinant != factor_determinant * factor_determinant) {
      throw std::runtime_error("det(R R^T) != det(R)^2");
    }
    const ExactMatrix<kOrder> adjugate =
        exact_adjugate(factor_gram);
    validate_adjugate(factor_gram, adjugate, gram_determinant);
    const std::vector<std::uint32_t> shell =
        enumerate_shell(adjugate, gram_determinant);

    const std::set<std::uint32_t> shell_set(
        shell.begin(), shell.end());
    const std::vector<std::uint32_t> columns =
        normalized_factor_columns(factor);
    if (std::set<std::uint32_t>(columns.begin(), columns.end()).size() !=
        kOrder) {
      throw std::runtime_error("known nonsingular factor has duplicate columns");
    }
    if (!std::all_of(
            columns.begin(), columns.end(),
            [&](std::uint32_t mask) {
              return shell_set.contains(mask);
            })) {
      throw std::runtime_error("known factor column is absent from shell");
    }

    const double elapsed = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - started)
                               .count();
    std::ostringstream output;
    output << "{\"assignments_completed\":" << kAssignments
           << ",\"complete\":true,\"elapsed_seconds\":" << elapsed
           << ",\"engine\":\"order22-gram-shell-v1\""
           << ",\"factor_determinant\":\""
           << json_escape(wide_to_string(factor_determinant)) << "\""
           << ",\"factor_columns_in_shell\":true"
           << ",\"gram_determinant\":\""
           << json_escape(wide_to_string(gram_determinant)) << "\""
           << ",\"matrix\":\"" << json_escape(options.matrix.string()) << "\""
           << ",\"order\":" << kOrder
           << ",\"shell_masks\":[";
    for (std::size_t index = 0; index < shell.size(); ++index) {
      if (index != 0) output << ',';
      output << shell[index];
    }
    output << "],\"shell_size\":" << shell.size()
           << ",\"known_factor_column_masks\":[";
    for (std::size_t index = 0; index < columns.size(); ++index) {
      if (index != 0) output << ',';
      output << columns[index];
    }
    output << "]}\n";
    atomic_write(options.output, output.str());
    std::cout << "completed " << kAssignments
              << " exact sign vectors; shell=" << shell.size()
              << " elapsed_seconds=" << elapsed << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "order22_gram_shell: " << error.what() << '\n';
    return 1;
  }
}
