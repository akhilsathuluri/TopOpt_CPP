/*
This file implements the 88 line matlab code.
Cite: Efficient topology optimization in MATLAB using 88 lines of code, E. Andreassen, A. Clausen, M. Schevenels, B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, Volume 43, Issue 1, p.1 - 16, (2011).

This paper uses the modified SIMP (Solid Isotropic Material with Penalisation) approach for optimisation.
*/

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<Eigen/Core>
#include <Eigen/SparseCore>
#include <unsupported/Eigen/KroneckerProduct>
// #include<math.h>
#include<cmath>
#include<algorithm>

using namespace Eigen;

void top(int nelx, int nely, double volfrac, double penal, double rmin, int ft);

int main(int argc, char const *argv[]) {
  top(4,3,0.5,3,1.5,1);
  return 0;
}

void top(int nelx, int nely, double volfrac, double penal, double rmin, int ft){
  // Define required variables
  double E0, Emin, nu;
  // Define material properties
  E0 = 1;
  Emin = pow(10,-9);
  nu = 0.3;
  // Element stiffness matrix
  MatrixXd A11(4,4), A12(4,4), B11(4,4), B12(4,4);
  A11 << 12,3,-6,-3,3,12,3,0,-6,3,12,-3,-3,0,-3,12;
  A12 << -6,-3,0,3,-3,-6,-3,-6,0,-3,-6,3,3,-6,3,-6;
  B11 << -4,3,-2,9,3,-4,-9,4,-2,-9,-4,-3,9,4,-3,-4;
  B12 << 2,-3,4,-9,-3,2,9,-2,4,9,2,3,-9,-2,3,2;
  // Four node bi-linear finite element
  MatrixXd KE(8,8), tKE1(8,8), tKE2(8,8);
  tKE1 << A11, A12, A12.transpose(), A11;
  tKE2 << B11, B12, B12.transpose(), B11;
  KE << 1/(1-pow(nu,2))/24*tKE1 + nu*tKE2;
  // Preparing indices for assembly
  ArrayXXi nodenrs((nelx+1)*(nely+1), 1);
  nodenrs.col(0) = ArrayXi::LinSpaced((nelx+1)*(nely+1), 0, (nelx+1)*(nely+1)-1);
  nodenrs.resize(nely+1, nelx+1);
  MatrixXi tnodenrs(nely, nelx);
  // Since node numbers start from 0 in this case
  tnodenrs = 2*(nodenrs.block(0,0,nely,nelx)+1);
  VectorXi edofVec(Map<VectorXi>(tnodenrs.data(), tnodenrs.rows()*tnodenrs.cols()));
  MatrixXi tedofMat(1,8), small(1,4), edofMat(nelx*nely,8);
  small << 2,3,0,1;
  tedofMat << 0, 1, 2*nely+small.array(), -2, -1;
  edofMat = tedofMat.replicate(nelx*nely, 1);
  edofMat = edofVec.replicate(1,8)+tedofMat.replicate(nelx*nely, 1);
  // VectorXi iK(64*nelx*nely), jK(64*nelx*nely);
  // kroneckerProduct(edofMat, MatrixXd::Ones(8,1)).eval()
  MatrixXi tempiK(8*nelx*nely,8), tempjK(8*nelx*nely,8);
  tempiK = kroneckerProduct(edofMat,MatrixXi::Constant(8,1,1)).eval();
  tempjK = kroneckerProduct(edofMat,MatrixXi::Constant(1,8,1)).eval();
  VectorXi iK(Map<VectorXi>(tempiK.data(), tempiK.rows()*tempiK.cols()));
  VectorXi jK(Map<VectorXi>(tempjK.data(), tempjK.rows()*tempjK.cols()));
  // Define loads and supports (assumes half MBB-beam )
  ArrayXd F = ArrayXd::Zero(2*(nelx+1)*(nely+1));
  F(1) = -1;
  // Update all these matrices to sparse in the next iteration
  // SparseMatrix<double> sF;
  // sF = (F.matrix()).sparseView();
  // Initialising all U to zero as a MatrixXd
  MatrixXd U = MatrixXd::Zero(2*(nelx+1)*(nely+1), 1);
  // ArrayXd U = ArrayXd::Zero(2*(nelx+1)*(nely+1));
  ArrayXi alldof = ArrayXi::LinSpaced(2*(nelx+1)*(nely+1), 0, 2*(nelx+1)*(nely+1)-1);
  ArrayXi fixeddof(nely+2);
  fixeddof.head(nely+1) = VectorXi::LinSpaced(nely+1, 0, 2*(nely+1));
  fixeddof(nely+1) = 2*(nelx+1)*(nely+1)-1;
  ArrayXi freedof(alldof.size()-fixeddof.size()); std::set_difference(alldof.data(), alldof.data()+alldof.size(), fixeddof.data(), fixeddof.data()+fixeddof.size(), freedof.data());
  // Preparing filter
  // Creating patches of neighbours based on rmin
  // CHANGE THE SIZE TO REQUIRED ONLY
  ArrayXi iH = ArrayXi::Ones(nelx*nely*pow(2*(ceil(rmin)-1)+1, 2), 1);
  ArrayXi jH = ArrayXi::Ones(iH.size(), 1);
  ArrayXd sH = ArrayXd::Zero(iH.size(), 1);
  int k = 0, e1 = 0, e2 = 0;
  for (int i1 = 0; i1 < nelx; i1++) {
    for (int j1 = 0; j1 < nely; j1++) {
      e1 = i1*nely+j1;
      for (int i2 = std::max(i1-(int(ceil(rmin))-1), 0) ; i2 <= std::min(i1+(int(ceil(rmin))-1), nelx-1); i2++) {
        for (int j2 = std::max(j1-(int(ceil(rmin))-1), 0) ; j2 <= std::min(j1+(int(ceil(rmin))-1), nely-1); j2++) {
          e2 = (i2)*nely+j2;
          k++;
          iH(k) = e1;
          jH(k) = e2;
          sH(k) = std::max(0.0, rmin-sqrt(double(pow(i1-i2, 2)+pow(j1-j2, 2))));
        }
      }
    }
  }
  // Similar to the sparse function in MATLAB
  ArrayXXd H = ArrayXXd::Zero(iH.size(), iH.size());
  for (int ii = 0; ii < iH.size(); ii++) {
    H(iH(ii), jH(ii)) +=  sH(ii);
  }
  VectorXd Hs(iH.size());
  Hs = H.rowwise().sum();
  // IS THERE ANY WAY TO VERIFY Hs?
  // Initialise iteration
  ArrayXXd x = ArrayXXd::Constant(nelx, nely, volfrac);
  ArrayXXd xPhys = x;
  ArrayXXd Ex = ArrayXXd::Constant(nelx, nely, 0);
  int loop = 0;
  double change = 1;
  VectorXd vKE(Map<VectorXd>(KE.data(), KE.cols()*KE.rows()));
  ArrayXXd sK = ArrayXXd::Constant(64, nelx*nely, 0);
  // Start iteration
  ArrayXXd K = ArrayXXd::Constant(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1), 0);
  while (change > 0.01) {
    loop++;
    // FE analysis
    Ex = (Emin+(xPhys.array().pow(penal)*(E0-Emin))).matrix();
    VectorXd vEx(Map<VectorXd>(Ex.data(), Ex.cols()*Ex.rows()));
    sK = vKE*vEx.transpose();
    VectorXd vsK(Map<VectorXd>(sK.data(), sK.cols()*sK.rows()));
    for (size_t jj = 0; jj < iK.size(); jj++) {
      K(iK(jj), jK(jj)) = sK(jj);
    }
    K = (K+K.transpose())/2;
  }
 // std::cout<<U(all, freedof)<<std::endl;
 // All keyword is not recognised by Eigen library
}
