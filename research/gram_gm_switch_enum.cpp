#include <algorithm>
#include <bit>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

struct Arguments {
  int order = 0;
  std::string edge_mask;
  bool help = false;
};

std::uint64_t parse_hex(std::string_view text) {
  if (text.empty() || text.size() > 16) {
    throw std::runtime_error("hex word must contain 1..16 digits");
  }
  std::uint64_t value = 0;
  for (const char character : text) {
    unsigned digit = 0;
    if (character >= '0' && character <= '9') {
      digit = static_cast<unsigned>(character - '0');
    } else if (character >= 'a' && character <= 'f') {
      digit = static_cast<unsigned>(character - 'a' + 10);
    } else if (character >= 'A' && character <= 'F') {
      digit = static_cast<unsigned>(character - 'A' + 10);
    } else {
      throw std::runtime_error("edge mask contains a non-hex digit");
    }
    value = (value << 4U) | digit;
  }
  return value;
}

std::vector<std::uint64_t> parse_edge_mask(
    std::string_view text, int edge_count) {
  if (text.empty()) {
    throw std::runtime_error("--edge-mask must not be empty");
  }
  if (text.starts_with("0x") || text.starts_with("0X")) {
    text.remove_prefix(2);
  }
  if (text.empty()) {
    throw std::runtime_error("--edge-mask must contain hex digits");
  }
  const std::size_t words =
      static_cast<std::size_t>((edge_count + 63) / 64);
  if (text.size() > words * 16) {
    throw std::runtime_error("--edge-mask has too many hex digits");
  }
  std::vector<std::uint64_t> result(words, 0);
  std::size_t end = text.size();
  for (std::size_t word = 0; end != 0; ++word) {
    const std::size_t begin = end > 16 ? end - 16 : 0;
    result[word] = parse_hex(text.substr(begin, end - begin));
    end = begin;
  }
  if (edge_count % 64 != 0) {
    const std::uint64_t allowed =
        (std::uint64_t{1} << (edge_count % 64)) - 1;
    if ((result.back() & ~allowed) != 0) {
      throw std::runtime_error(
          "--edge-mask sets a bit outside the graph");
    }
  }
  return result;
}

std::string format_edge_mask(
    const std::vector<std::uint64_t>& words, int edge_count) {
  const int digits = (edge_count + 3) / 4;
  std::ostringstream output;
  output << std::hex << std::nouppercase << std::setfill('0');
  for (std::size_t word = words.size(); word-- > 0;) {
    if (word == words.size() - 1) {
      const int leading_digits = digits - 16 * static_cast<int>(word);
      output << std::setw(leading_digits) << words[word];
    } else {
      output << std::setw(16) << words[word];
    }
  }
  return output.str();
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    auto value = [&]() -> std::string {
      if (++index >= argc) {
        throw std::runtime_error("missing value for " + option);
      }
      return argv[index];
    };
    if (option == "--order") {
      const std::string text = value();
      std::size_t consumed = 0;
      arguments.order = std::stoi(text, &consumed);
      if (consumed != text.size()) {
        throw std::runtime_error("--order must be an integer");
      }
    } else if (option == "--edge-mask") {
      arguments.edge_mask = value();
    } else if (option == "--help" || option == "-h") {
      arguments.help = true;
    } else {
      throw std::runtime_error("unknown option: " + option);
    }
  }
  if (arguments.help) return arguments;
  if (arguments.order < 2 || arguments.order > 62) {
    throw std::runtime_error("--order must be in 2..62");
  }
  if (arguments.edge_mask.empty()) {
    throw std::runtime_error("--edge-mask is required");
  }
  return arguments;
}

void print_usage() {
  std::cout
      << "Usage: gram_gm_switch_enum --order N --edge-mask HEX\n"
      << "Exhaust every nontrivial even Godsil--McKay switching set.\n"
      << "Output rows are: switching-set-hex<TAB>mate-edge-mask-hex.\n"
      << "A final #summary row gives masks_examined, valid_sets, "
         "and labeled_mates.\n";
}

int edge_index(int order, int first, int second) {
  if (first > second) std::swap(first, second);
  // Number of upper-triangle entries before row first.
  return first * (2 * order - first - 1) / 2 + second - first - 1;
}

bool edge_bit(
    const std::vector<std::uint64_t>& words, int index) {
  return ((words[static_cast<std::size_t>(index / 64)] >>
           static_cast<unsigned>(index % 64)) &
          1U) != 0;
}

void toggle_edge(
    std::vector<std::uint64_t>& words, int index) {
  words[static_cast<std::size_t>(index / 64)] ^=
      std::uint64_t{1} << static_cast<unsigned>(index % 64);
}

int run(const Arguments& arguments) {
  const int order = arguments.order;
  const int edge_count = order * (order - 1) / 2;
  const std::vector<std::uint64_t> edge_words =
      parse_edge_mask(arguments.edge_mask, edge_count);
  std::vector<std::uint64_t> adjacency(
      static_cast<std::size_t>(order), 0);
  for (int first = 0; first < order; ++first) {
    for (int second = first + 1; second < order; ++second) {
      if (!edge_bit(
              edge_words, edge_index(order, first, second))) {
        continue;
      }
      adjacency[static_cast<std::size_t>(first)] |=
          std::uint64_t{1} << static_cast<unsigned>(second);
      adjacency[static_cast<std::size_t>(second)] |=
          std::uint64_t{1} << static_cast<unsigned>(first);
    }
  }

  const std::uint64_t limit =
      std::uint64_t{1} << static_cast<unsigned>(order);
  std::uint64_t masks_examined = 0;
  std::uint64_t valid_sets = 0;
  std::uint64_t labeled_mates = 0;
  for (std::uint64_t set = 1; set < limit; ++set) {
    const int size = std::popcount(set);
    // A nontrivial switch needs an outside vertex with |X|/2
    // neighbors, hence |X| is even and X is a proper subset.
    if (size < 2 || (size & 1) != 0 || size == order) continue;
    ++masks_examined;

    const int first_vertex = std::countr_zero(set);
    const int induced_degree = std::popcount(
        adjacency[static_cast<std::size_t>(first_vertex)] & set);
    bool regular = true;
    for (std::uint64_t remaining =
             set & ~(std::uint64_t{1}
                     << static_cast<unsigned>(first_vertex));
         remaining != 0;
         remaining &= remaining - 1) {
      const int vertex = std::countr_zero(remaining);
      if (std::popcount(
              adjacency[static_cast<std::size_t>(vertex)] & set) !=
          induced_degree) {
        regular = false;
        break;
      }
    }
    if (!regular) continue;

    const int half = size / 2;
    std::uint64_t half_vertices = 0;
    bool admissible = true;
    const std::uint64_t outside = (limit - 1) ^ set;
    for (std::uint64_t remaining = outside;
         remaining != 0;
         remaining &= remaining - 1) {
      const int vertex = std::countr_zero(remaining);
      const int neighbors = std::popcount(
          adjacency[static_cast<std::size_t>(vertex)] & set);
      if (neighbors == half) {
        half_vertices |=
            std::uint64_t{1} << static_cast<unsigned>(vertex);
      } else if (neighbors != 0 && neighbors != size) {
        admissible = false;
        break;
      }
    }
    if (!admissible || half_vertices == 0) continue;
    ++valid_sets;

    std::vector<std::uint64_t> mate = edge_words;
    for (std::uint64_t outer = half_vertices;
         outer != 0;
         outer &= outer - 1) {
      const int outside_vertex = std::countr_zero(outer);
      for (std::uint64_t inner = set;
           inner != 0;
           inner &= inner - 1) {
        const int inside_vertex = std::countr_zero(inner);
        toggle_edge(
            mate,
            edge_index(
                order, inside_vertex, outside_vertex));
      }
    }
    ++labeled_mates;
    std::cout << std::hex << set << '\t'
              << format_edge_mask(mate, edge_count) << '\n';
  }
  std::cout << std::dec
            << "#summary\tmasks_examined=" << masks_examined
            << "\tvalid_sets=" << valid_sets
            << "\tlabeled_mates=" << labeled_mates << '\n';
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    if (arguments.help) {
      print_usage();
      return 0;
    }
    return run(arguments);
  } catch (const std::exception& error) {
    std::cerr << "gram_gm_switch_enum: " << error.what() << '\n';
    return 2;
  }
}
