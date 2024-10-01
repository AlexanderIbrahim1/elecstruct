#include <string>

#include "elecstruct/elecstruct.hpp"

exported_class::exported_class()
    : m_name {"elecstruct"}
{}

auto exported_class::name() const -> const char*
{
    return m_name.c_str();
}
