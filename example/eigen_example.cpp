#include <iostream>
#include <Eigen/Dense>

auto main() -> int {
    // Define a 2x2 matrix
    Eigen::Matrix2d mat;
    mat(0, 0) = 3;
    mat(1, 0) = 2.5;
    mat(0, 1) = -1;
    mat(1, 1) = mat(1, 0) + mat(0, 1);

    // Print the matrix
    std::cout << "Here is the 2x2 matrix:\n" << mat << std::endl;

    // Perform some simple operations
    Eigen::Vector2d vec(1, 2);
    Eigen::Vector2d result = mat * vec;

    // Print the result of the matrix-vector multiplication
    std::cout << "Matrix-vector product:\n" << result << std::endl;

    return 0;
}

