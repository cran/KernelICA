/*
 * RcppIncompleteCholesky.cpp
 *
 *
 *  Created on: Apr 07, 2020
 *      Author: Christoph Koesner
 */

#include <RcppEigen.h>
#include "IncompleteCholesky.h"
#include "KernelGauss.h"
#include "KernelHermite.h"

// [[Rcpp::depends(RcppEigen)]]


// computes the incomplete cholesky composition for the Gauss kernel
// [[Rcpp::export]]
Rcpp::List incompleteCholeskyGauss(const Eigen::Map<Eigen::VectorXd> v, double eps, double sigma){
  Eigen::VectorXd vec = v;
  Kernel* ker = new KernelGauss(sigma);
  choleskyReturnVal retval = choleskyLowRank(vec, *ker, eps);
  
  
  Eigen::MatrixXd L = retval.L;
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> pp = retval.perm;
  
  delete ker;
  
  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("L") = retval.L,
    Rcpp::Named("perm") = retval.perm.indices().array() + 1);

  return(output);
}

// computes the incomplete cholesky composition for the Hermite kernel
// [[Rcpp::export]]
Rcpp::List incompleteCholeskyHermite(const Eigen::Map<Eigen::VectorXd> v, double eps, double sigma, int d){
  Eigen::VectorXd vec = v;
  Kernel* ker = new KernelHermite(sigma, d);
  choleskyReturnVal retval = choleskyLowRank(vec, *ker, eps);
  
  
  Eigen::MatrixXd L = retval.L;
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> pp = retval.perm;
  
  delete ker;
  
  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("L") = retval.L,
    Rcpp::Named("perm") = retval.perm.indices().array() + 1);
  
  return(output);
}
