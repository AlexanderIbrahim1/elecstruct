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
    
    auto atoms = info.atom_information();
    elec::fill_atomic_orbitals_sto3g(atoms);

    const auto basis = elec::create_atomic_orbitals_sto3g(atoms);
    const auto initial_fock = info.initial_fock_guess();
    const auto n_electrons = info.n_electrons();
    const auto max_iter = info.max_hartree_fock_iterations();
    const auto tolerance = info.tol_change_density_matrix();
    const auto verbose = info.verbose();

    elec::perform_restricted_hartree_fock(atoms, basis, initial_fock, n_electrons, max_iter, tolerance, verbose);

    return 0;
}
