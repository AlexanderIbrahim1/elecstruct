#include <sstream>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

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
        REQUIRE_THAT(atom_information[0].position.x, Catch::Matchers::WithinAbs(0.0, ABS_TOL));
        REQUIRE_THAT(atom_information[0].position.y, Catch::Matchers::WithinAbs(0.0, ABS_TOL));
        REQUIRE_THAT(atom_information[0].position.z, Catch::Matchers::WithinAbs(0.0, ABS_TOL));
    }
}
