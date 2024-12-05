#include <filesystem>
#include <iostream>

#include <elecstruct/elecstruct.hpp>

// TODO: move this into another directory

auto main(int argc, const char** argv) -> int
{
    if (argc != 2) {
        std::cerr << "./a.out path/to/input.toml\n";
        std::exit(EXIT_FAILURE);
    }

    const auto path = std::filesystem::path {argv[1]};
    auto toml_stream = std::ifstream {path};
    if (!toml_stream.is_open()) {
        std::cerr << "ERROR: failed to open the input toml file.\n";
        std::exit(EXIT_FAILURE);
    }

    auto parser = elec::InputFileParser {toml_stream};
    parser.parse_all();

    const auto& info = parser.parsed_information();
    
    const auto atoms = info.atom_information();
    const auto inital_fock = info.initial_fock_guess();
    const auto max_iterations = info.max_hartree_fock_iterations();
    const auto tol_change_density_matrix = info.tol_change_density_matrix();
    const auto tol_change_hartree_fock_energy = info.tol_change_hartree_fock_energy();
    const auto n_electrons = info.n_electrons();
    const auto verbose = info.verbose();

    const auto basis = elec::create_atomic_orbitals_sto3g(atoms);

    return 0;
}
