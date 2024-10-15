// TODO: put this into a unit test instead

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include "elecstruct/atoms.hpp"
#include "elecstruct/elecstruct.hpp"
#include "elecstruct/input_file_parser.hpp"

auto main() -> int
{
    namespace fs = std::filesystem;

    const auto path = fs::path {"/home/a68ibrah/research/elecstruct/playground/example_input_file.toml"};
    auto toml_stream = std::ifstream {path};
    auto parser = elec::InputFileParser {toml_stream};

    if (!parser.is_valid()) {
        std::cout << "Unable to parse the input file\n";
        std::cout << parser.error_message() << '\n';
        std::exit(EXIT_FAILURE);
    }

    for (const auto& atom_info : parser.atom_information) {
        const auto atom_name = elec::atom_name_from_label(atom_info.label);
        std::cout << atom_name << " : ";
        std::cout << atom_info.x << ", " << atom_info.y << ", " << atom_info.z << '\n';
    }

    return 0;
}
