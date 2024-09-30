#include <string>

#include <catch2/catch_test_macros.hpp>

#include "elecstruct/elecstruct.hpp"

TEST_CASE("elecstruct test")
{
    const auto exported = exported_class {};
    REQUIRE(exported.name() == std::string {"elecstruct"});
}
