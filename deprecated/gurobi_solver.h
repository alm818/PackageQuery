#pragma once

#include <Eigen/Dense>
#include "gurobi_c++.h"

using namespace Eigen;
using namespace std;

class GurobiSolver{
public:
  const VectorXd& c;
  int n, m;
  VectorXd r0, x0;
  double exe_init, exe_relaxed, exe_ilp, relaxed_cscore, ilp_cscore, iteration_count;
  int relaxed_status, ilp_status;
  GRBenv* env;
  GRBmodel* model;
public:
  GurobiSolver(const MatrixXd& A, const VectorXd& bl, const VectorXd& bu, const VectorXd& c, const VectorXd& l, const VectorXd& u, bool is_presolve=true);
  GurobiSolver(const MatrixXd& A, const VectorXd& b, const VectorXd& c, const VectorXd& u, bool is_presolve=true);
  void solveIlp(bool is_concurrent=true);
  void solveRelaxed(bool is_concurrent=true);
};