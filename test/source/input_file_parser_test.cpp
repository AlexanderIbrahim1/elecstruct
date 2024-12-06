#include <iomanip>
#include <sstream>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "elecstruct/geometry.hpp"
#include "elecstruct/input_file_parser/input_file_parser.hpp"

TEST_CASE("parse ATOM_INFORMATION")
{
    constexpr auto ABS_TOL = double {1.0e-8};

    SECTION("one atom at origin")
    {
        auto input_stream = std::stringstream {};
        input_stream << R"(
        positions = [
            ["H", 0.0, 0.0, 0.0],
        ]
        )";

        auto parser = elec::InputFileParser {input_stream};
        parser.parse(elec::InputFileKey::ATOM_INFORMATION);

        const auto& info = parser.parsed_information();
        const auto atom_information = info.atom_information();

        REQUIRE(atom_information.size() == 1);
        REQUIRE(atom_information[0].label == elec::AtomLabel::H);
        REQUIRE(coord::almost_equals(atom_information[0].position, {0.0, 0.0, 0.0}, ABS_TOL));
    }

    SECTION("one atom not at origin")
    {
        auto input_stream = std::stringstream {};
        input_stream << R"(
        positions = [
            ["He", 1.2, 3.4, 5.6],
        ]
        )";

        auto parser = elec::InputFileParser {input_stream};
        parser.parse(elec::InputFileKey::ATOM_INFORMATION);

        const auto& info = parser.parsed_information();
        const auto atom_information = info.atom_information();

        REQUIRE(atom_information.size() == 1);
        REQUIRE(atom_information[0].label == elec::AtomLabel::He);
        REQUIRE(coord::almost_equals(atom_information[0].position, {1.2, 3.4, 5.6}, ABS_TOL));
    }

    SECTION("two atoms, identical")
    {
        auto input_stream = std::stringstream {};
        input_stream << R"(
        positions = [
            ["H", 1.1, 2.2, 3.3],
            ["H", 4.4, 5.5, 6.6],
        ]
        )";

        auto parser = elec::InputFileParser {input_stream};
        parser.parse(elec::InputFileKey::ATOM_INFORMATION);

        const auto& info = parser.parsed_information();
        const auto atom_information = info.atom_information();

        REQUIRE(atom_information.size() == 2);

        REQUIRE(atom_information[0].label == elec::AtomLabel::H);
        REQUIRE(coord::almost_equals(atom_information[0].position, {1.1, 2.2, 3.3}, ABS_TOL));

        REQUIRE(atom_information[1].label == elec::AtomLabel::H);
        REQUIRE(coord::almost_equals(atom_information[1].position, {4.4, 5.5, 6.6}, ABS_TOL));
    }

    SECTION("three atoms, different")
    {
        auto input_stream = std::stringstream {};
        input_stream << R"(
        positions = [
            ["H", 1.1, 2.2, 3.3],
            ["Li", 4.4, 5.5, 6.6],
            ["C", 7.7, 8.8, 9.9],
        ]
        )";

        auto parser = elec::InputFileParser {input_stream};
        parser.parse(elec::InputFileKey::ATOM_INFORMATION);

        const auto& info = parser.parsed_information();
        const auto atom_information = info.atom_information();

        REQUIRE(atom_information.size() == 3);

        REQUIRE(atom_information[0].label == elec::AtomLabel::H);
        REQUIRE(coord::almost_equals(atom_information[0].position, {1.1, 2.2, 3.3}, ABS_TOL));

        REQUIRE(atom_information[1].label == elec::AtomLabel::Li);
        REQUIRE(coord::almost_equals(atom_information[1].position, {4.4, 5.5, 6.6}, ABS_TOL));

        REQUIRE(atom_information[2].label == elec::AtomLabel::C);
        REQUIRE(coord::almost_equals(atom_information[2].position, {7.7, 8.8, 9.9}, ABS_TOL));
    }
}

TEST_CASE("parse INITIAL_FOCK_GUESS")
{
    using IFG = elec::InitialFockGuess;

    SECTION("valid input")
    {
        struct TestPair
        {
            std::string input;
            IFG expected;
        };

        const auto pair = GENERATE(
            TestPair {"zero", IFG::ZERO_MATRIX},
            TestPair {"extended_huckel", IFG::EXTENDED_HUCKEL_MATRIX},
            TestPair {"core_hamiltonian", IFG::CORE_HAMILTONIAN_MATRIX}
        );

        auto input_stream = std::stringstream {};
        input_stream << "initial_fock_guess = "
                     << "\"" << pair.input << "\"" << '\n';

        auto parser = elec::InputFileParser {input_stream};
        parser.parse(elec::InputFileKey::INITIAL_FOCK_GUESS);

        const auto& info = parser.parsed_information();
        const auto guess = info.initial_fock_guess();

        REQUIRE(guess == pair.expected);
    }

    SECTION("invalid throws")
    {
        auto input_stream = std::stringstream {};
        input_stream << R"(initial_fock_guess = "invalid_type"\n)";

        auto parser = elec::InputFileParser {input_stream};

        REQUIRE_THROWS_AS(parser.parse(elec::InputFileKey::INITIAL_FOCK_GUESS), std::runtime_error);
    }
}

TEST_CASE("parse MAX_HARTREE_FOCK_ITERATIONS")
{
    using IFG = elec::InputFileKey;

    SECTION("valid example")
    {
        auto input_stream = std::stringstream {};
        input_stream << "max_hartree_fock_iterations = 5\n";

        auto parser = elec::InputFileParser {input_stream};
        parser.parse(IFG::MAX_HARTREE_FOCK_ITERATIONS);

        const auto& info = parser.parsed_information();
        const auto max_iter = info.max_hartree_fock_iterations();

        REQUIRE(max_iter == 5);
    }

    SECTION("not there")
    {
        auto input_stream = std::stringstream {};
        input_stream << "\n";

        auto parser = elec::InputFileParser {input_stream};
        REQUIRE_THROWS_AS(parser.parse(IFG::MAX_HARTREE_FOCK_ITERATIONS), std::runtime_error);
    }

    SECTION("negative argument")
    {
        auto input_stream = std::stringstream {};
        input_stream << "max_hartree_fock_iterations = -10\n";

        auto parser = elec::InputFileParser {input_stream};
        REQUIRE_THROWS_AS(parser.parse(IFG::MAX_HARTREE_FOCK_ITERATIONS), std::runtime_error);
    }

    SECTION("zero argument")
    {
        auto input_stream = std::stringstream {};
        input_stream << "max_hartree_fock_iterations = 0\n";

        auto parser = elec::InputFileParser {input_stream};
        REQUIRE_THROWS_AS(parser.parse(IFG::MAX_HARTREE_FOCK_ITERATIONS), std::runtime_error);
    }
}

TEST_CASE("parse TOL_CHANGE_DENSITY_MATRIX")
{
    using IFG = elec::InputFileKey;

    SECTION("valid example")
    {
        auto input_stream = std::stringstream {};
        input_stream << "tol_change_density_matrix = 1.0e-5\n";

        auto parser = elec::InputFileParser {input_stream};
        parser.parse(IFG::TOL_CHANGE_DENSITY_MATRIX);

        const auto& info = parser.parsed_information();
        const auto tol_change = info.tol_change_density_matrix();

        REQUIRE_THAT(tol_change, Catch::Matchers::WithinRel(1.0e-5));
    }

    SECTION("not there")
    {
        auto input_stream = std::stringstream {};
        input_stream << "\n";

        auto parser = elec::InputFileParser {input_stream};
        REQUIRE_THROWS_AS(parser.parse(IFG::TOL_CHANGE_DENSITY_MATRIX), std::runtime_error);
    }

    SECTION("negative argument")
    {
        auto input_stream = std::stringstream {};
        input_stream << "tol_change_density_matrix = -4.5e-3\n";

        auto parser = elec::InputFileParser {input_stream};
        REQUIRE_THROWS_AS(parser.parse(IFG::TOL_CHANGE_DENSITY_MATRIX), std::runtime_error);
    }

    SECTION("zero argument")
    {
        auto input_stream = std::stringstream {};
        input_stream << "tol_change_density_matrix = 0.0\n";

        auto parser = elec::InputFileParser {input_stream};
        REQUIRE_THROWS_AS(parser.parse(IFG::TOL_CHANGE_DENSITY_MATRIX), std::runtime_error);
    }
}

TEST_CASE("parse TOL_CHANGE_HARTREE_FOCK_ENERGY")
{
    using IFG = elec::InputFileKey;

    SECTION("valid example")
    {
        auto input_stream = std::stringstream {};
        input_stream << "tol_change_hartree_fock_energy = 1.0e-5\n";

        auto parser = elec::InputFileParser {input_stream};
        parser.parse(IFG::TOL_CHANGE_HARTREE_FOCK_ENERGY);

        const auto& info = parser.parsed_information();
        const auto tol_change = info.tol_change_hartree_fock_energy();

        REQUIRE_THAT(tol_change, Catch::Matchers::WithinRel(1.0e-5));
    }

    SECTION("not there")
    {
        auto input_stream = std::stringstream {};
        input_stream << "\n";

        auto parser = elec::InputFileParser {input_stream};
        REQUIRE_THROWS_AS(parser.parse(IFG::TOL_CHANGE_HARTREE_FOCK_ENERGY), std::runtime_error);
    }

    SECTION("negative argument")
    {
        auto input_stream = std::stringstream {};
        input_stream << "tol_change_hartree_fock_energy = -4.5e-3\n";

        auto parser = elec::InputFileParser {input_stream};
        REQUIRE_THROWS_AS(parser.parse(IFG::TOL_CHANGE_HARTREE_FOCK_ENERGY), std::runtime_error);
    }

    SECTION("zero argument")
    {
        auto input_stream = std::stringstream {};
        input_stream << "tol_change_hartree_fock_energy = 0.0\n";

        auto parser = elec::InputFileParser {input_stream};
        REQUIRE_THROWS_AS(parser.parse(IFG::TOL_CHANGE_HARTREE_FOCK_ENERGY), std::runtime_error);
    }
}

TEST_CASE("parse N_ELECTRONS")
{
    using IFG = elec::InputFileKey;

    SECTION("valid example")
    {
        const auto input = static_cast<std::size_t>(GENERATE(0, 5, 10));

        auto input_stream = std::stringstream {};
        input_stream << "n_electrons = " << input << '\n';

        auto parser = elec::InputFileParser {input_stream};
        parser.parse(IFG::N_ELECTRONS);

        const auto& info = parser.parsed_information();
        const auto n_electrons = info.n_electrons();

        REQUIRE(n_electrons == input);
    }

    SECTION("not there")
    {
        auto input_stream = std::stringstream {};
        input_stream << "\n";

        auto parser = elec::InputFileParser {input_stream};
        REQUIRE_THROWS_AS(parser.parse(IFG::N_ELECTRONS), std::runtime_error);
    }

    SECTION("negative argument")
    {
        auto input_stream = std::stringstream {};
        input_stream << "n_electrons = -10\n";

        auto parser = elec::InputFileParser {input_stream};
        REQUIRE_THROWS_AS(parser.parse(IFG::N_ELECTRONS), std::runtime_error);
    }
}

TEST_CASE("parse VERBOSE")
{
    SECTION("valid input")
    {
        struct TestPair
        {
            bool input;
            elec::Verbose expected;
        };

        const auto pair = GENERATE(TestPair {true, elec::Verbose::TRUE}, TestPair {false, elec::Verbose::FALSE});

        auto input_stream = std::stringstream {};
        input_stream << "verbose = " << std::boolalpha << pair.input << '\n';

        auto parser = elec::InputFileParser {input_stream};
        parser.parse(elec::InputFileKey::VERBOSE);

        const auto& info = parser.parsed_information();
        const auto actual = info.verbose();

        REQUIRE(actual == pair.expected);
    }

    SECTION("invalid throws")
    {
        auto input_stream = std::stringstream {};
        input_stream << R"(verbose = invalid_type\n)";

        auto parser = elec::InputFileParser {input_stream};

        REQUIRE_THROWS_AS(parser.parse(elec::InputFileKey::VERBOSE), std::runtime_error);
    }
}
