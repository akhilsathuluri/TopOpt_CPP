/*
This file implements the 88 line matlab code.
Cite: Efficient topology optimization in MATLAB using 88 lines of code, E. Andreassen, A. Clausen, M. Schevenels, B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, Volume 43, Issue 1, p.1 - 16, (2011).

This paper uses the modified SIMP (solid isotropic material with penalisation) approach for optimisation.

Compilation: time clang++ -Ofast main.cpp
If using openmp parallelisation
Compilation: time g++ -fopenmp -Ofast main.cpp
Execution: OMP_NUM_THREADS=threads ./a.out
Requires Eigen library

Author: Akhil Sathuluri

Next update: Change dense matrices to sparse
*/

#include<stdio.h>
#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<Eigen/Core>
#include<Eigen/SparseCore>
#include<unsupported/Eigen/KroneckerProduct>
#include<cmath>
#include<algorithm>
// For plotting
#include <fstream>
#include<stdlib.h>

using namespace Eigen;

void top(int nelx, int nely, double volfrac, double penal, double rmin, int ft, int plot);

int main(int argc, char const *argv[]) {
  // top(80,30,0.5,3,1.5,1,0);
  top(120,40,0.5,3,1.5,1,1);
  return 0;
}

void top(int nelx, int nely, double volfrac, double penal, double rmin, int ft, int plot){
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
  KE << 1/(1-pow(nu,2))/24*(tKE1 + nu*tKE2);
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
  MatrixXi iK(8*nelx*nely,8), jK(8*nelx*nely,8);
  iK = (kroneckerProduct(edofMat,MatrixXi::Constant(8,1,1)).eval()).transpose();
  iK.resize(64*nelx*nely,1);
  jK = (kroneckerProduct(edofMat,MatrixXi::Constant(1,8,1)).eval()).transpose();
  jK.resize(64*nelx*nely,1);
  // Define loads and supports (assumes half MBB-beam )
  VectorXd F = VectorXd::Zero(2*(nelx+1)*(nely+1));
  F(1) = -1;
  // std::cout << F(freedof) << std::endl;
  // Update all these matrices to sparse in the next iteration
  // SparseMatrix<double> F;
  // F = (sF.matrix()).sparseView();
  VectorXd U = VectorXd::Zero(2*(nelx+1)*(nely+1));
  ArrayXi alldof = ArrayXi::LinSpaced(2*(nelx+1)*(nely+1), 0, 2*(nelx+1)*(nely+1)-1);
  ArrayXi fixeddof(nely+2);
  fixeddof.head(nely+1) = VectorXi::LinSpaced(nely+1, 0, 2*(nely+1));
  fixeddof(nely+1) = 2*(nelx+1)*(nely+1)-1;
  ArrayXi freedof(alldof.size()-fixeddof.size()); std::set_difference(alldof.data(), alldof.data()+alldof.size(), fixeddof.data(), fixeddof.data()+fixeddof.size(), freedof.data());
  // Preparing filter
  // Creating patches of neighbours based on rmin
  VectorXi iH = VectorXi::Ones(nelx*nely*pow(2*(ceil(rmin)-1)+1, 2), 1);
  VectorXi jH = VectorXi::Ones(iH.size(), 1);
  VectorXd sH = VectorXd::Zero(iH.size(), 1);
  int k = -1, e1 = 0, e2 = 0;
  for (int i1 = 0; i1 < nelx; i1++) {
    for (int j1 = 0; j1 < nely; j1++) {
      e1 = i1*nely+j1;
      for (int i2 = std::max(i1-(int(ceil(rmin))-1), 0) ; i2 <= std::min(i1+(int(ceil(rmin))-1), nelx-1); i2++) {
        for (int j2 = std::max(j1-(int(ceil(rmin))-1), 0) ; j2 <= std::min(j1+(int(ceil(rmin))-1), nely-1); j2++) {
          e2 = i2*nely+j2;
          k++;
          iH(k) = e1;
          jH(k) = e2;
          sH(k) = std::max(0.0, rmin-sqrt(double(pow(i1-i2, 2)+pow(j1-j2, 2))));
        }
      }
    }
  }
  // Similar to the sparse function in MATLAB
  MatrixXd H = MatrixXd::Zero(nelx*nely, nelx*nely);
  for (int ii = 0; ii < iH.size(); ii++) {
    H(iH(ii), jH(ii)) +=  sH(ii);
  }
  VectorXd Hs(H.rows());
  Hs = H.rowwise().sum();
  // Initialise iteration
  MatrixXd x = MatrixXd::Constant(nely, nelx, volfrac);
  MatrixXd xPhys = x;
  MatrixXd xnew = MatrixXd::Constant(nely, nelx, 0);
  MatrixXd Ex = MatrixXd::Constant(nely, nelx, 0);
  int loop = 0;
  double change = 1;
  VectorXd vKE(Map<VectorXd>(KE.data(), KE.cols()*KE.rows()));
  MatrixXd sK = MatrixXd::Constant(64, nelx*nely, 0);
  // Start iteration
  MatrixXd ce = MatrixXd::Zero(edofMat.rows(), 1);
  double c = 0.0;
  MatrixXd dc = MatrixXd::Constant(nely, nelx, 0);
  MatrixXd dv = MatrixXd::Ones(nely, nelx);
  // Initialising iterations
  while (change > 0.01) {
    loop++;
    Ex = (Emin+(xPhys.array().pow(penal)*(E0-Emin))).matrix();
    VectorXd vEx(Map<VectorXd>(Ex.data(), Ex.cols()*Ex.rows()));
    sK = vKE*vEx.transpose();
    VectorXd vsK(Map<VectorXd>(sK.data(), sK.cols()*sK.rows()));
    // We need a new K matrix at each iteration
    MatrixXd K = MatrixXd::Constant(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1), 0);
    for (int jj = 0; jj < iK.size(); jj++) {
      K(iK(jj), jK(jj)) += sK(jj);
    }
    K = (K+K.transpose())/2.0;
    // The following step needs Eigen 3.3.9 unstable update
    U(freedof) = (K(freedof, freedof)).colPivHouseholderQr().solve(F(freedof));
    // Resize the ce matrix into vector
    ce.resize(nely*nelx,1);
    // Objective function and sensitivity analysis
    for (int ii = 0; ii < edofMat.rows(); ii++) {
      ce(ii,0) = U(edofMat.row(ii)).transpose()*KE*U(edofMat.row(ii));
    }
    // Resize into a matrix
    ce.resize(nely, nelx);
    // The matrix Ex is already computed above
    c = (Ex.array()*ce.array()).sum();
    dc = (-penal*(E0-Emin)*xPhys.array().pow(penal-1))*ce.array();
    dc.resize(nely*nelx,1);
    dv.resize(nely*nelx,1);
    // Filtering/Modification of sensitivities
    if (ft == 1) {
      VectorXd vx(Map<VectorXd>(x.data(), x.cols()*x.rows()));
      VectorXd mvx = vx;
      // Removing zeros for division
      mvx = (vx.array() < pow(10,-3)).select(pow(10,-3), vx);
      dc = ((H*(vx.array()*dc.array()).matrix()).array())/Hs.array()/mvx.array();
    }
    else if (ft == 2) {
      dc = H*((dc.array()/Hs.array()).matrix());
      dv = H*((dv.array()/Hs.array()).matrix());
    }
    // Resizing dc and dv back to normal
    dc.resize(nely, nelx);
    dv.resize(nely, nelx);
    // Optimality criteria update of design variables
    double l1 = 0.0, l2 = pow(10, 9), move = 0.2, lmid = 0.0;
    while ((l2-l1)/(l1+l2) > pow(10, -3)) {
      lmid = 0.5*(l1+l2);
      // sqrt produces nan values
      xnew = (ArrayXXd::Zero(nely, nelx)).max((x.array()-move).max((ArrayXXd::Ones(nely, nelx)).min((x.array()+move).min(x.array()*(-dc.array()/dv.array()/lmid).sqrt()))));
      if (ft == 1) {
        xPhys = xnew;
      }
      else if (ft == 2) {
        xnew.resize(nely*nelx,1);
        xPhys.resize(nely*nelx,1);
        xPhys = (H*xnew).array()/Hs.array();
        xnew.resize(nely, nelx);
        xPhys.resize(nely, nelx);
      }
      if (xPhys.sum() > volfrac*nely*nelx) {
        l1 = lmid;
      }
      else {
        l2 = lmid;
      }
    }
    change = (xnew - x).array().abs().maxCoeff();
    x = xnew;
    // Print results
    printf("It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n",loop, c, xPhys.mean(), change);
    // save results
    if(plot != 0){
      std::ofstream file("sol.txt", std::ios::trunc);
      file << (1-xPhys.array()).matrix() << std::endl << std::endl;
      // write code to read from this file and plot
      system("./callplot.sh");
      file.close();
    }
  }
  std::cout << "Saving optimal solution ..." << std::endl;
  std::ofstream file("sol.txt", std::ios::trunc);
  file << (1-xPhys.array()).matrix() << std::endl << std::endl;
  // Leave the final result plotted
  std::cout << "Displaying optimal solution ..." << std::endl;
  system("./callplot_final.sh");
  file.close();
}
