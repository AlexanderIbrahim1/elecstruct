#include <string>

#include "elecstruct/elecstruct.hpp"

auto main() -> int
{
  auto const exported = exported_class {};

  return std::string("elecstruct") == exported.name() ? 0 : 1;
}
