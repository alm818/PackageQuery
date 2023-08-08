#pragma once

#include <Eigen/Dense>
#include "coo_sparse.h"
#include "gurobi_c++.h"

using namespace Eigen;
using namespace std;

class GurobiLatticeSolver {
public:
  VectorXd b, bq, c, ub, lb, r0, floor_r0, e_score, i_score, best_x, x0;
  MatrixXd A, Aq, tableau;
  CooSparse lattice_dirs;
  vector<int> bhead;
  double c_score, best_c_score, avg_step_count, milp_score;
  double exe_relaxed, exe_tableau, exe_init, exe_milp;
  int n, m, try_count, degenerate_count, status;
  // List of variable violations for each direction
  vector<vector<int>> degenerate_violation;
  GRBenv* env;
  GRBmodel* model;

public:
  ~GurobiLatticeSolver();
  GurobiLatticeSolver(MatrixXd Aq, MatrixXd bq, MatrixXd A, VectorXd b, VectorXd c, VectorXd lb, VectorXd ub);
  GurobiLatticeSolver(MatrixXd A, VectorXd b, VectorXd c);
  GurobiLatticeSolver(MatrixXd A, VectorXd b, VectorXd c, double lb_value, double ub_value);
  GurobiLatticeSolver(string problem);
  double getObjValue(VectorXd sol);
  bool checkFeasibility(VectorXd sol);
  void solve(double time_budget);
  void milpSolve();

private:
  void solveRelaxed();
  bool latticeStep(int d, VectorXd& x, VectorXd& current_e_score, VectorXd& current_i_score, double& current_c_score, unordered_map<int,int>& bound_map);
  bool checkFeasibilityStep(const VectorXd& current_e_score, const VectorXd& current_i_score, const unordered_map<int,int>& bound_map);
};