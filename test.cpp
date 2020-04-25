// Eigen using multithreading to speed up

#include <Eigen/Core>
#include <iostream>
using namespace Eigen;
int main() {
 int n = 4000;
 MatrixXd A = MatrixXd::Ones(n,n);
 MatrixXd B = MatrixXd::Ones(n,n);
 MatrixXd C = MatrixXd::Ones(n,n);
 C += A*B;
 std::cout << C.sum() << "\n";
}
