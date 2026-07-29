#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace {

constexpr int kCoreOrder = 22;
constexpr int kOutputOrder = 23;
constexpr std::uint64_t kAssignments =
    std::uint64_t{1} << (kCoreOrder - 1);
using Wide = __int128_t;

template <std::size_t Order>
using Matrix = std::array<std::array<int, Order>, Order>;

using Core = Matrix<kCoreOrder>;
using Candidate = Matrix<kOutputOrder>;
using Adjugate =
    std::array<std::array<Wide, kCoreOrder>, kCoreOrder>;

struct Cell {
  int row = 0;
  int column = 0;

  auto key() const { return std::pair{row, column}; }
};

struct FlipSpec {
  std::vector<Cell> cells;
  std::size_t line_number = 0;
};

struct Options {
  std::filesystem::path start;
  std::filesystem::path flips;
  std::filesystem::path output;
  std::filesystem::path core_output;
  std::filesystem::path log;
  int radius = 0;
  std::uint64_t progress_every = 100;
};

struct BorderResult {
  Candidate matrix{};
  Wide score = 0;
  std::uint64_t best_assignment = 0;
  std::uint64_t best_multiplicity = 0;
  std::uint64_t improvements = 0;
};

std::string wide_to_string(Wide value) {
  if (value == 0) return "0";
  const bool negative = value < 0;
  if (negative) value = -value;
  std::string result;
  while (value != 0) {
    result.push_back(static_cast<char>('0' + value % 10));
    value /= 10;
  }
  if (negative) result.push_back('-');
  std::reverse(result.begin(), result.end());
  return result;
}

Wide absolute(Wide value) { return value < 0 ? -value : value; }

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

std::uint64_t parse_unsigned(
    const std::string& text,
    const std::string& option) {
  if (text.empty() ||
      !std::all_of(text.begin(), text.end(), [](unsigned char character) {
        return character >= '0' && character <= '9';
      })) {
    throw std::runtime_error("invalid value for " + option);
  }
  std::size_t consumed = 0;
  const unsigned long long parsed = std::stoull(text, &consumed);
  if (consumed != text.size()) {
    throw std::runtime_error("invalid value for " + option);
  }
  return static_cast<std::uint64_t>(parsed);
}

Options parse_options(int argc, char** argv) {
  Options options;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--help") {
      std::cout
          << "Usage: order22_perturb_border --start CORE --flips TSV "
             "--radius 1|2\n"
          << "       --output MATRIX --core-output CORE_MATRIX --log JSONL "
             "[--progress-every CORES]\n"
          << "TSV records are unique, 1-based tab-separated coordinates:\n"
          << "  radius 1: row column\n"
          << "  radius 2: row1 column1 row2 column2\n";
      std::exit(0);
    }
    if (index + 1 >= argc) {
      throw std::runtime_error("missing value for " + argument);
    }
    const std::string value = argv[++index];
    if (argument == "--start") {
      options.start = value;
    } else if (argument == "--flips") {
      options.flips = value;
    } else if (argument == "--output") {
      options.output = value;
    } else if (argument == "--core-output") {
      options.core_output = value;
    } else if (argument == "--log") {
      options.log = value;
    } else if (argument == "--radius") {
      const std::uint64_t parsed = parse_unsigned(value, argument);
      if (parsed != 1 && parsed != 2) {
        throw std::runtime_error("--radius must be 1 or 2");
      }
      options.radius = static_cast<int>(parsed);
    } else if (argument == "--progress-every") {
      options.progress_every = parse_unsigned(value, argument);
    } else {
      throw std::runtime_error("unknown argument: " + argument);
    }
  }

  if (options.start.empty() || options.flips.empty() ||
      options.output.empty() || options.core_output.empty() ||
      options.log.empty() || options.radius == 0) {
    throw std::runtime_error(
        "--start, --flips, --radius, --output, --core-output, and --log "
        "are required");
  }

  std::array<std::filesystem::path*, 5> paths = {
      &options.start,
      &options.flips,
      &options.output,
      &options.core_output,
      &options.log,
  };
  for (std::filesystem::path* path : paths) {
    *path = std::filesystem::absolute(*path).lexically_normal();
  }
  for (std::size_t first = 0; first < paths.size(); ++first) {
    for (std::size_t second = first + 1; second < paths.size(); ++second) {
      if (*paths[first] == *paths[second]) {
        throw std::runtime_error("all input and output paths must differ");
      }
    }
  }
  return options;
}

template <std::size_t Order>
Wide exact_determinant(Matrix<Order> matrix) {
  constexpr int size = static_cast<int>(Order);
  std::array<std::array<Wide, Order>, Order> work{};
  for (int row = 0; row < size; ++row) {
    for (int column = 0; column < size; ++column) {
      work[row][column] = matrix[row][column];
    }
  }

  Wide previous_pivot = 1;
  int sign = 1;
  for (int column = 0; column < size - 1; ++column) {
    int pivot_row = column;
    while (pivot_row < size && work[pivot_row][column] == 0) {
      ++pivot_row;
    }
    if (pivot_row == size) return 0;
    if (pivot_row != column) {
      std::swap(work[pivot_row], work[column]);
      sign = -sign;
    }

    const Wide pivot = work[column][column];
    for (int row = column + 1; row < size; ++row) {
      for (int inner = column + 1; inner < size; ++inner) {
        const Wide numerator =
            work[row][inner] * pivot -
            work[row][column] * work[column][inner];
        if (column != 0 && numerator % previous_pivot != 0) {
          throw std::runtime_error("exact Bareiss division failed");
        }
        work[row][inner] =
            column == 0 ? numerator : numerator / previous_pivot;
      }
      work[row][column] = 0;
    }
    previous_pivot = pivot;
  }
  return sign * work[size - 1][size - 1];
}

Core read_core(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open core matrix");
  Core matrix{};
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      if (!(input >> matrix[row][column]) ||
          (matrix[row][column] != -1 && matrix[row][column] != 1)) {
        throw std::runtime_error("invalid 22 by 22 sign core");
      }
    }
  }
  std::string extra;
  if (input >> extra) throw std::runtime_error("extra core matrix data");
  return matrix;
}

std::vector<FlipSpec> read_flip_specs(
    const std::filesystem::path& path,
    int radius) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open flip TSV");

  std::vector<FlipSpec> result;
  std::set<std::vector<int>> seen_keys;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    std::vector<std::string> fields;
    std::size_t field_start = 0;
    while (true) {
      const std::size_t tab = line.find('\t', field_start);
      fields.push_back(line.substr(field_start, tab - field_start));
      if (tab == std::string::npos) break;
      field_start = tab + 1;
    }
    if (fields.size() != static_cast<std::size_t>(2 * radius)) {
      throw std::runtime_error(
          "flip TSV line " + std::to_string(line_number) +
          " must contain exactly " + std::to_string(2 * radius) +
          " tab-separated integer fields");
    }

    FlipSpec spec;
    spec.line_number = line_number;
    std::vector<int> key;
    for (int index = 0; index < radius; ++index) {
      const std::uint64_t row = parse_unsigned(
          fields[2 * index],
          "row field in flip TSV line " + std::to_string(line_number));
      const std::uint64_t column = parse_unsigned(
          fields[2 * index + 1],
          "column field in flip TSV line " + std::to_string(line_number));
      if (row < 1 || row > kCoreOrder ||
          column < 1 || column > kCoreOrder) {
        throw std::runtime_error(
            "coordinate outside 1..22 in flip TSV line " +
            std::to_string(line_number));
      }
      const Cell cell{
          static_cast<int>(row - 1),
          static_cast<int>(column - 1),
      };
      spec.cells.push_back(cell);
      key.push_back(cell.row);
      key.push_back(cell.column);
    }
    if (radius == 2 && !(spec.cells[0].key() < spec.cells[1].key())) {
      throw std::runtime_error(
          "radius-2 cells must be distinct and sorted in flip TSV line " +
          std::to_string(line_number));
    }
    if (!seen_keys.insert(key).second) {
      throw std::runtime_error(
          "duplicate flip TSV record at line " +
          std::to_string(line_number));
    }
    result.push_back(std::move(spec));
  }
  if (!input.eof()) throw std::runtime_error("cannot read flip TSV");
  if (result.empty()) throw std::runtime_error("flip TSV is empty");
  return result;
}

Adjugate exact_adjugate(const Core& core) {
  Adjugate adjugate{};
  for (int excluded_row = 0; excluded_row < kCoreOrder; ++excluded_row) {
    for (int excluded_column = 0; excluded_column < kCoreOrder;
         ++excluded_column) {
      Matrix<kCoreOrder - 1> minor{};
      int minor_row = 0;
      for (int row = 0; row < kCoreOrder; ++row) {
        if (row == excluded_row) continue;
        int minor_column = 0;
        for (int column = 0; column < kCoreOrder; ++column) {
          if (column == excluded_column) continue;
          minor[minor_row][minor_column] = core[row][column];
          ++minor_column;
        }
        ++minor_row;
      }
      Wide cofactor = exact_determinant(minor);
      if ((excluded_row + excluded_column) % 2 != 0) {
        cofactor = -cofactor;
      }
      adjugate[excluded_column][excluded_row] = cofactor;
    }
  }
  return adjugate;
}

void validate_adjugate(
    const Core& core,
    const Adjugate& adjugate,
    Wide determinant) {
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      Wide product = 0;
      for (int inner = 0; inner < kCoreOrder; ++inner) {
        product +=
            static_cast<Wide>(core[row][inner]) *
            adjugate[inner][column];
      }
      const Wide expected = row == column ? determinant : 0;
      if (product != expected) {
        throw std::runtime_error("exact adjugate identity check failed");
      }
    }
  }
}

Candidate bordered_candidate(
    const Core& core,
    const Adjugate& adjugate,
    const std::array<int, kCoreOrder>& border_column,
    Wide core_determinant) {
  Candidate candidate{};
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      candidate[row][column] = core[row][column];
    }
    candidate[row][kCoreOrder] = border_column[row];
  }

  for (int row = 0; row < kCoreOrder; ++row) {
    Wide value = 0;
    for (int column = 0; column < kCoreOrder; ++column) {
      value += adjugate[row][column] * border_column[column];
    }
    candidate[kCoreOrder][row] = value > 0 ? -1 : 1;
  }
  candidate[kCoreOrder][kCoreOrder] = core_determinant > 0 ? 1 : -1;
  return candidate;
}

template <std::size_t Order>
void write_matrix_atomic(
    const std::filesystem::path& path,
    const Matrix<Order>& matrix) {
  if (!path.parent_path().empty()) {
    std::filesystem::create_directories(path.parent_path());
  }
  if (std::filesystem::exists(path) &&
      !std::filesystem::is_regular_file(path)) {
    throw std::runtime_error("output path is not a regular file");
  }
  const std::filesystem::path temporary = path.string() + ".tmp";
  {
    std::ofstream output(temporary, std::ios::trunc);
    if (!output) throw std::runtime_error("cannot create output matrix");
    for (const auto& row : matrix) {
      for (std::size_t column = 0; column < Order; ++column) {
        if (column != 0) output << ' ';
        output << row[column];
      }
      output << '\n';
    }
    output.flush();
    if (!output) throw std::runtime_error("cannot flush output matrix");
  }
  std::error_code error;
  std::filesystem::rename(temporary, path, error);
  if (error) {
    std::filesystem::remove(temporary);
    throw std::runtime_error("cannot replace output matrix");
  }
}

BorderResult exhaust_borders(
    const Core& core,
    const Adjugate& adjugate,
    Wide core_determinant) {
  std::array<int, kCoreOrder> border_column{};
  border_column.fill(1);
  std::array<Wide, kCoreOrder> products{};
  for (int row = 0; row < kCoreOrder; ++row) {
    for (int column = 0; column < kCoreOrder; ++column) {
      products[row] += adjugate[row][column];
    }
  }

  BorderResult result;
  result.score = absolute(core_determinant);
  for (const Wide value : products) result.score += absolute(value);
  result.matrix =
      bordered_candidate(core, adjugate, border_column, core_determinant);
  if (absolute(exact_determinant(result.matrix)) != result.score) {
    throw std::runtime_error("initial border identity check failed");
  }
  result.best_multiplicity = 1;

  for (std::uint64_t assignment = 1; assignment < kAssignments;
       ++assignment) {
    const unsigned bit = std::countr_zero(assignment) + 1U;
    const int old_sign = border_column[bit];
    border_column[bit] = -old_sign;
    for (int row = 0; row < kCoreOrder; ++row) {
      products[row] -=
          static_cast<Wide>(2 * old_sign) * adjugate[row][bit];
    }

    Wide score = absolute(core_determinant);
    for (const Wide value : products) score += absolute(value);
    if (score < result.score) continue;
    if (score == result.score) {
      ++result.best_multiplicity;
      continue;
    }

    const Candidate candidate = bordered_candidate(
        core, adjugate, border_column, core_determinant);
    if (absolute(exact_determinant(candidate)) != score) {
      throw std::runtime_error("border determinant identity check failed");
    }
    result.score = score;
    result.matrix = candidate;
    result.best_assignment = assignment;
    result.best_multiplicity = 1;
    ++result.improvements;
  }

  if (absolute(exact_determinant(result.matrix)) != result.score) {
    throw std::runtime_error("final exact determinant check failed");
  }
  return result;
}

std::string flips_json(const FlipSpec& spec) {
  std::string result = "[";
  for (std::size_t index = 0; index < spec.cells.size(); ++index) {
    if (index != 0) result += ',';
    result +=
        "[" + std::to_string(spec.cells[index].row + 1) + "," +
        std::to_string(spec.cells[index].column + 1) + "]";
  }
  result += ']';
  return result;
}

void refuse_existing_output(const std::filesystem::path& path) {
  if (std::filesystem::exists(path)) {
    throw std::runtime_error(
        "refusing to overwrite existing campaign output: " +
        path.string());
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse_options(argc, argv);
    refuse_existing_output(options.output);
    refuse_existing_output(options.core_output);
    refuse_existing_output(options.log);

    const Core base = read_core(options.start);
    const Wide base_determinant = exact_determinant(base);
    const std::vector<FlipSpec> specs =
        read_flip_specs(options.flips, options.radius);
    if (specs.size() >
        std::numeric_limits<std::uint64_t>::max() / kAssignments) {
      throw std::runtime_error("campaign assignment count overflows uint64");
    }
    const std::uint64_t total_assignments =
        static_cast<std::uint64_t>(specs.size()) * kAssignments;

    if (!options.log.parent_path().empty()) {
      std::filesystem::create_directories(options.log.parent_path());
    }
    std::ofstream log(options.log, std::ios::trunc);
    if (!log) throw std::runtime_error("cannot create JSONL log");

    const auto start_time = std::chrono::steady_clock::now();
    log << "{\"assignments_per_core\":" << kAssignments
        << ",\"assignments_total\":" << total_assignments
        << ",\"base_absolute_determinant\":\""
        << wide_to_string(absolute(base_determinant))
        << "\",\"core_count\":" << specs.size()
        << ",\"event\":\"campaign_start\",\"flips\":\""
        << json_escape(options.flips.string())
        << "\",\"radius\":" << options.radius
        << ",\"start\":\"" << json_escape(options.start.string()) << "\"}\n";
    log.flush();

    Wide global_best = -1;
    std::uint64_t global_best_core_count = 0;
    std::size_t global_best_core_index = 0;
    std::uint64_t completed_assignments = 0;

    for (std::size_t index = 0; index < specs.size(); ++index) {
      const FlipSpec& spec = specs[index];
      Core core = base;
      for (const Cell cell : spec.cells) {
        core[cell.row][cell.column] = -core[cell.row][cell.column];
      }
      const auto core_start = std::chrono::steady_clock::now();
      const Wide core_determinant = exact_determinant(core);
      log << "{\"core_index\":" << (index + 1)
          << ",\"event\":\"core_start\",\"flips_1_based\":"
          << flips_json(spec)
          << ",\"source_line\":" << spec.line_number << "}\n";

      const Adjugate adjugate = exact_adjugate(core);
      validate_adjugate(core, adjugate, core_determinant);
      const BorderResult result =
          exhaust_borders(core, adjugate, core_determinant);
      completed_assignments += kAssignments;

      if (result.score > global_best) {
        global_best = result.score;
        global_best_core_count = 1;
        global_best_core_index = index + 1;
        write_matrix_atomic(options.core_output, core);
        write_matrix_atomic(options.output, result.matrix);
        log << "{\"absolute_determinant\":\""
            << wide_to_string(result.score)
            << "\",\"core_index\":" << (index + 1)
            << ",\"event\":\"new_global_best\",\"flips_1_based\":"
            << flips_json(spec) << "}\n";
      } else if (result.score == global_best) {
        ++global_best_core_count;
      }

      const double core_elapsed =
          std::chrono::duration<double>(
              std::chrono::steady_clock::now() - core_start)
              .count();
      log << "{\"absolute_core_determinant\":\""
          << wide_to_string(absolute(core_determinant))
          << "\",\"absolute_determinant\":\""
          << wide_to_string(result.score)
          << "\",\"assignments_completed\":" << kAssignments
          << ",\"best_assignment\":" << result.best_assignment
          << ",\"best_border_columns_up_to_global_sign\":"
          << result.best_multiplicity
          << ",\"core_index\":" << (index + 1)
          << ",\"elapsed_seconds\":" << core_elapsed
          << ",\"event\":\"core_finished\",\"flips_1_based\":"
          << flips_json(spec)
          << ",\"improvements\":" << result.improvements
          << ",\"source_line\":" << spec.line_number << "}\n";
      log.flush();
      if (!log) throw std::runtime_error("cannot flush JSONL log");

      const std::uint64_t completed = index + 1;
      if (options.progress_every != 0 &&
          (completed % options.progress_every == 0 ||
           completed == specs.size())) {
        const double elapsed =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start_time)
                .count();
        log << "{\"assignments_completed\":" << completed_assignments
            << ",\"cores_completed\":" << completed
            << ",\"cores_total\":" << specs.size()
            << ",\"elapsed_seconds\":" << elapsed
            << ",\"event\":\"progress\",\"global_best\":\""
            << wide_to_string(global_best) << "\"}\n";
        log.flush();
        std::cout << "completed=" << completed << '/' << specs.size()
                  << " assignments=" << completed_assignments
                  << " global_best=" << wide_to_string(global_best)
                  << " elapsed_seconds=" << elapsed << '\n';
      }
    }

    const double elapsed =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time)
            .count();
    log << "{\"assignments_completed\":" << completed_assignments
        << ",\"complete\":true,\"cores_completed\":" << specs.size()
        << ",\"elapsed_seconds\":" << elapsed
        << ",\"event\":\"finished\",\"global_best_absolute_determinant\":\""
        << wide_to_string(global_best)
        << "\",\"global_best_core_count\":" << global_best_core_count
        << ",\"global_best_core_index\":" << global_best_core_index
        << "}\n";
    log.flush();
    if (!log) throw std::runtime_error("cannot flush final JSONL event");

    std::cout << "complete cores=" << specs.size()
              << " assignments=" << completed_assignments
              << " global_best=" << wide_to_string(global_best)
              << " elapsed_seconds=" << elapsed << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "order22_perturb_border: " << error.what() << '\n';
    return 1;
  }
}
