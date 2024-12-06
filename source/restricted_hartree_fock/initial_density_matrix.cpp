#include <Eigen/Dense>

#include "elecstruct/restricted_hartree_fock/initial_density_matrix.hpp"

/*
    NOTE: right now we only have the zero matrix, but I may add other ones later.
*/

namespace
{

/*
    The value of the Extended Huckel constant provided in the reference paper.

    SOURCE:

    title: "Extended Hückel and Slater's rule initial guess for real space grid-based density functional theory"
    authors: M. Lee and K. Leiter and C. Eisner and J. Crone and J. Knap
    journal: Computational and Theoretical Chemistry
    volume: 1062
    year: 2015

    Equation 12, middle of the right column on page 25
*/
constexpr auto DEFAULT_EXTENDED_HUCKEL_CONSTANT = double {1.75};

}  // anonymous namespace

namespace elec
{

auto zero_matrix(std::size_t size) -> Eigen::MatrixXd
{
    const auto size_ = static_cast<Eigen::Index>(size);
    return Eigen::MatrixXd::Zero(size_, size_);
}

auto extended_huckel_guess(const Eigen::MatrixXd& overlap_mtx, const Eigen::MatrixXd& core_hamiltonian_mtx)
    -> Eigen::MatrixXd
{
    const auto size = overlap_mtx.cols();
    const auto constant = DEFAULT_EXTENDED_HUCKEL_CONSTANT;

    auto fock_huckel_guess = Eigen::MatrixXd {size, size};
    for (Eigen::Index i0 {0}; i0 < size; ++i0) {
        for (Eigen::Index i1 {0}; i1 < size; ++i1) {
            const auto core_ham_avg = 0.5 * (core_hamiltonian_mtx(i0, i0) + core_hamiltonian_mtx(i1, i1));
            fock_huckel_guess(i0, i1) = constant * overlap_mtx(i0, i1) * core_ham_avg;
        }
    }

    return fock_huckel_guess;
}

auto core_hamiltonian_guess(const Eigen::MatrixXd& core_hamiltonian_mtx) -> Eigen::MatrixXd
{
    return core_hamiltonian_mtx;
}

}  // namespace elec
