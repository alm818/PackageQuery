#pragma once

#include <Eigen/Dense>
#include "coo_sparse.h"

using namespace Eigen;
using namespace std;

class LatticeSolver{
public:
  const MatrixXd& A;
  const VectorXd& b;
  const VectorXd& c;
  const VectorXd& u;
  vector<int> bhead;
  CooSparse lattice_dirs;
  unordered_map<int, int> inv_bhead;
  VectorXd r0, near_r0, fscores, best_x, fracs;
  double relaxed_cscore, cscore, best_cscore, avg_step_count;
  double exe_init, exe_relaxed, exe_tableau, exe_solved, exe_find_dir, exe_init_walk, exe_walk;
  int n, m, core, try_count, status, opt;
public:
  ~LatticeSolver();
  LatticeSolver(int core, const MatrixXd& A, const VectorXd& b, const VectorXd& c, const VectorXd& u);
  double getObjValue(VectorXd sol);
  bool checkFeasibility(VectorXd sol);
  void solve(double time_budget);
  void report();
  void compareReport(VectorXd sol, double solved_time=-1);
private:
  bool latticeStep(int d, VectorXd& x, VectorXd& current_fscores, double& current_cscore, unordered_map<int,int>& bound_map);
  bool checkFeasibilityStep(const VectorXd& current_fscores, const unordered_map<int,int>& bound_map);
};