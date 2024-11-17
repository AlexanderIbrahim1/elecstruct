#pragma once

#include <Eigen/Dense>

/*
    NOTE: right now we only have the zero matrix, but I may add other ones later.
*/

namespace elec
{

/*
    A poor and non-physical, but very easy guess to make.
*/
inline auto zero_matrix(std::size_t size) -> Eigen::MatrixXd
{
    return Eigen::MatrixXd::Zero(size, size);
}

}  // namespace elec
