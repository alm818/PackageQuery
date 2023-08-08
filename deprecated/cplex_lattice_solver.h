#pragma once

#include <Eigen/Dense>
#include "coo_sparse.h"
#include "ilcplex/ilocplex.h"

using namespace Eigen;
using namespace std;

class CplexLatticeSolver {
public:
  VectorXd b, c, ub, lb, r0, floor_r0, feas_scores, best_x;
  MatrixXd A;
  CooSparse lattice_dirs;
  double c_score, best_c_score, avg_step_count;
  double exe_relaxed, exe_presolve, exe_tableau, exe_init;
  int n, m, try_count, degenerate_count;
  // List of variable violations for each direction
  vector<vector<int>> degenerate_violation;

  // Cplex object and attribute for original problem;
  CPXENVptr env;
  CPXLPptr model;
  VectorXd obj, uncrushed_r0, x0;
  int numcols, numrows, sol_status;
  bool is_presolvable;
  double exe_milp, milp_score;

public:
  ~CplexLatticeSolver();
  CplexLatticeSolver(MatrixXd A, VectorXd b, VectorXd c);
  CplexLatticeSolver(MatrixXd A, VectorXd b, VectorXd c, double lb_value, double ub_value);
  CplexLatticeSolver(MatrixXd A, VectorXd b, VectorXd c, VectorXd lb, VectorXd ub);
  CplexLatticeSolver(string problem);
  double getObjValue(VectorXd sol);
  bool checkFeasibility(VectorXd sol);
  void solve(double time_budget);
  void milpSolve();

private:
  VectorXd uncrushX(VectorXd pre_x);
  void initialize(MatrixXd A, VectorXd b, VectorXd c, VectorXd lb, VectorXd ub);
  void extractModel(CPXCLPptr model);
  void solveRelaxed();
  bool latticeStep(int d, VectorXd& x, VectorXd& current_feas_scores, double& current_c_score, unordered_map<int,int>& bound_map);
  bool checkFeasibilityStep(const VectorXd& current_feas_scores, const unordered_map<int,int>& bound_map);
};