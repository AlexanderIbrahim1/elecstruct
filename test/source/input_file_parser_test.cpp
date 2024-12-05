#include <sstream>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "elecstruct/geometry.hpp"
#include "elecstruct/input_file_parser/input_file_parser.hpp"


TEST_CASE("parse ATOM_INFORMATION")
{
    constexpr auto ABS_TOL = double{1.0e-8};

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
