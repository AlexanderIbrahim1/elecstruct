#include <algorithm>
#include <vector>

#include <Eigen/Dense>

#include "elecstruct/atoms.hpp"
#include "elecstruct/cartesian3d.hpp"
#include "elecstruct/geometry.hpp"
#include "elecstruct/matrices.hpp"

#include "elecstruct/restricted_hartree_fock/step.hpp"

namespace elec
{

/*
    Create a new matrix with the sorted columns.
*/
auto matrix_with_sorted_columns(const Eigen::MatrixXd& matrix, const std::vector<Eigen::Index>& indices)
    -> Eigen::MatrixXd
{
    const auto size = matrix.cols();

    if (size != static_cast<Eigen::Index>(indices.size())) {
        throw std::runtime_error {"Mismatch between number of columns in matrix, and number of indices."};
    }

    auto output = Eigen::MatrixXd {size, size};

    for (Eigen::Index i {0}; i < size; ++i) {
        output.col(i) = matrix.col(indices[static_cast<std::size_t>(i)]);
    }

    return output;
}

/*
    Returns a vector of indices such that if the elements in the Eigen::VectorXd `elements` were rearranged
    to these indices, then it would be sorted.
*/
auto indices_to_sort(const Eigen::VectorXd& elements) -> std::vector<Eigen::Index>
{
    auto indices = [&]()
    {
        const auto size = static_cast<std::size_t>(elements.size());

        auto indices_ = std::vector<Eigen::Index> {};
        indices_.reserve(size);
        for (std::size_t i {0}; i < size; ++i) {
            indices_.push_back(static_cast<Eigen::Index>(i));
        }

        return indices_;
    }();

    const auto compare = [&](Eigen::Index i0, Eigen::Index i1) { return elements[i0] < elements[i1]; };

    std::sort(indices.begin(), indices.end(), compare);

    return indices;
}

auto new_density_matrix(
    const Eigen::MatrixXd& fock_mtx,
    const Eigen::MatrixXd& basis_transformation_mtx,
    std::size_t n_electrons
) -> Eigen::MatrixXd
{
    const auto fock_mtx_trans = basis_transformation_mtx.transpose() * fock_mtx * basis_transformation_mtx;

    const auto eigensolver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> {fock_mtx_trans};
    if (eigensolver.info() != Eigen::Success) {
        throw std::runtime_error {"Failed to perform eigenvalue decomposition of Fock matrix."};
    }

    // Get the sorted eigenvalues and eigenvectors; later operations depend on acessing the
    // eigenstates in order of ascending energy, which is why this reordering is required
    const auto eigenvalues = eigensolver.eigenvalues();
    const auto sorted_indices = indices_to_sort(eigenvalues);

    const auto coefficient_mtx_trans = matrix_with_sorted_columns(eigensolver.eigenvectors(), sorted_indices);
    const auto coefficient_mtx = basis_transformation_mtx * coefficient_mtx_trans;

    const auto new_density_mtx = density_matrix_restricted_hartree_fock(coefficient_mtx, n_electrons);

    return new_density_mtx;
}

auto density_matrix_difference(const Eigen::MatrixXd& old_density_mtx, const Eigen::MatrixXd& new_density_mtx) -> double
{
    const auto size = old_density_mtx.cols();

    auto difference = double {0.0};
    for (Eigen::Index i0 {0}; i0 < size; ++i0) {
        for (Eigen::Index i1 {0}; i1 < size; ++i1) {
            const auto element_diff = new_density_mtx(i0, i1) - old_density_mtx(i0, i1);
            difference += element_diff * element_diff;
        }
    }

    return 0.5 * std::sqrt(difference);
}

auto electron_energy(
    const Eigen::MatrixXd& density_mtx,
    const Eigen::MatrixXd& fock_mtx,
    const Eigen::MatrixXd& core_hamiltonain_mtx
) -> double
{
    const auto size = density_mtx.cols();

    auto energy = double {0.0};

    for (Eigen::Index i0 {0}; i0 < size; ++i0) {
        for (Eigen::Index i1 {0}; i1 < size; ++i1) {
            energy += 0.5 * density_mtx(i0, i1) * (fock_mtx(i0, i1) + core_hamiltonain_mtx(i0, i1));
        }
    }

    return energy;
}

auto nuclear_energy(const std::vector<AtomInfo>& atoms) -> double
{
    const auto size = atoms.size();

    auto energy = double {0.0};
    for (std::size_t i0 {0}; i0 < size - 1; ++i0) {
        for (std::size_t i1 {i0 + 1}; i1 < size; ++i1) {
            const auto pos0 = atoms[i0].position;
            const auto pos1 = atoms[i1].position;
            const auto charge0 = nuclear_charge(atoms[i0].label);
            const auto charge1 = nuclear_charge(atoms[i1].label);

            energy += charge0 * charge1 / coord::distance(pos0, pos1);
        }
    }

    return energy;
}

auto total_energy(
    const Eigen::MatrixXd& density_mtx,
    const Eigen::MatrixXd& fock_mtx,
    const Eigen::MatrixXd& core_hamiltonain_mtx,
    const std::vector<AtomInfo>& atoms
) -> double
{
    const auto elec_energy = electron_energy(density_mtx, fock_mtx, core_hamiltonain_mtx);
    const auto nucl_energy = nuclear_energy(atoms);

    return elec_energy + nucl_energy;
}

}  // namespace elec
