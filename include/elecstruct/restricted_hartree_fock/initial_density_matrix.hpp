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
    const auto size_ = static_cast<Eigen::Index>(size);
    return Eigen::MatrixXd::Zero(size_, size_);
}

inline auto huckel_guess(
    const Eigen::MatrixXd& overlap_mtx,
    const Eigen::MatrixXd& core_hamiltonian_mtx,
    double huckel_constant
)
{
    const auto size = overlap_mtx.cols();

    auto fock_huckel_guess = Eigen::MatrixXd {size, size};
    for (Eigen::Index i0 {0}; i0 < size; ++i0) {
        for (Eigen::Index i1 {0}; i1 < size; ++i1) {
            const auto core_ham_avg = 0.5 * (core_hamiltonian_mtx(i0, i0) + core_hamiltonian_mtx(i1, i1));
            if (i0 == i1) {
                fock_huckel_guess(i0, i1) = core_ham_avg;
            } else {
                fock_huckel_guess(i0, i1) = huckel_constant * overlap_mtx(i0, i1) * core_ham_avg;
            }
        }
    }

    return fock_huckel_guess;
}

}  // namespace elec
