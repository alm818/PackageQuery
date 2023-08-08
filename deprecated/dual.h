#pragma once

#include <Eigen/Dense>
#include <boost/dynamic_bitset.hpp>
#include "utility.h"

using namespace Eigen;

class Dual{
public:
  // Original data
  const RMatrixXd& A;
  const VectorXd& bbl;
  const VectorXd& bbu;
  const VectorXd& c;
  const VectorXd& l;
  const VectorXd& u;
  // Original members
  int n, m, iteration_count, mini_iteration_count, status;
  VectorXd x;
  VectorXi bhead;
  boost::dynamic_bitset<> inv_bhead;
  double exe_solve, score;
  // Internal members
  VectorXd bu, bl, beta, d;
  VectorXd rho_r, alpha_r, alpha_q, tau;
  MatrixXd Binv;
public:
  ~Dual();
  Dual(int core, const RMatrixXd& A, const VectorXd& bl, const VectorXd& bu, const VectorXd& c, const VectorXd& l, const VectorXd& u);
};