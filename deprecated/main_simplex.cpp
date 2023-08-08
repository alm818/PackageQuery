#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <bitset>
#include <random>

#include "utility.h"
#include "walker.h"
#include "pseudo_walker.h"
#include "gurobi_lattice_solver.h"
#include "cplex_lattice_solver.h"
#include "fmt/core.h"
#include "simplex.h"

void generateProlem(int n, MatrixXd& A, VectorXd& b, VectorXd& c){
  A.resize(6, n); b.resize(6); c.resize(n);
  //default_random_engine gen;
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  int expected_numvar = 10;
  normal_distribution n_dist(0.5*expected_numvar, 1.0/12*expected_numvar);
  normal_distribution n_dist_c(0.0, 200.0);
  for (int i = 0; i < n; i ++){
    A(0, i) = u_dist(gen);
    A(1, i) = -A(0, i);
    A(2, i) = u_dist(gen);
    A(3, i) = -A(2, i);
    A(4, i) = 1;
    A(5, i) = -1;
    c(i) = n_dist_c(gen);
  }
  double tol_scale = 5;
  double tol1 = u_dist(gen) * tol_scale / 2;
  double tol2 = u_dist(gen) * tol_scale / 2;
  b(0) = n_dist(gen) + tol1;
  b(1) = tol1 - b(0);
  b(2) = n_dist(gen) + tol2;
  b(3) = tol2 - b(2);
  b(4) = 20;
  b(5) = -5;
}

// void test2(){
//   #pragma omp for
//   for (int i = 0; i < 32; i ++){
//     #pragma omp critical
//     cout << "Wtf " << i << endl;
//   }
// }

// void test(){
//   #pragma omp parallel num_threads(32)
//   {
//     while (1){
//       test2();
//       break;
//     }
//     #pragma omp barrier
//   }
// }

void reportGurobi(MatrixXd A, VectorXd b, VectorXd c, double lb_value, double ub_value, double time_budget, bool is_milp){
  GurobiLatticeSolver ls = GurobiLatticeSolver(A, b, c, lb_value, ub_value);
  fmt::print("n:{} m:{} nnz:{}\n", ls.n, ls.m, ls.lattice_dirs.getSize());
  if (!kIgnoreDegenerate){
    cout << ls.degenerate_count << " degenerate direction out of " << ls.n << endl;
  }
  fmt::print("Try count:{} Average step count:{:.2Lf} Relaxed time:{:.5Lf}ms Tableau time:{:.5Lf}ms Init time:{:.5Lf}ms\n", ls.try_count, ls.avg_step_count, ls.exe_relaxed, ls.exe_tableau, ls.exe_init);
  double relaxed_score = ls.getObjValue(ls.r0);
  cout << "Their score: " << relaxed_score << endl;
  return;
  ls.solve(time_budget);
  // print(ls.r0);
  double relax_gap = abs((ls.best_c_score - relaxed_score) / relaxed_score) * 100;
  bool feas = false;
  if (ls.best_c_score == -DBL_MAX){
    ls.best_c_score = nan("");
    relax_gap = nan("");
  } else feas = ls.checkFeasibility(ls.best_x);
  fmt::print("Best_score:{:.2Lf}\n", ls.best_c_score);
  fmt::print("Relaxed_score:{:.2Lf} gap:{:.2Lf}%\n", relaxed_score, relax_gap);
  if (is_milp){
    ls.milpSolve();
    double obj_score = ls.milp_score;
    double l1_dist = (ls.x0-ls.r0).lpNorm<1>();
    double obj_gap = abs((ls.best_c_score - obj_score) / obj_score) * 100;
    if (ls.best_c_score == -DBL_MAX) obj_gap = nan("");
    double lp_gap = abs((obj_score - relaxed_score) / relaxed_score) * 100;
    fmt::print("Milp time:{:.5Lf}ms Obj_score:{:.2Lf} gap:{:.2Lf}% LP_gap:{:.2Lf}% L1_dist:{:.2Lf}\n", ls.exe_milp, obj_score, obj_gap, lp_gap, l1_dist);
  }
  cout << feas << " " << ls.checkFeasibility(ls.r0) << endl;
}

int main(){
  // test();
  int N = 1000000;
  MatrixXd A (6, N);
  VectorXd b (6);
  VectorXd c (N);
  generateProlem(N, A, b, c);
  VectorXd u (N); u.fill(1);
  // MatrixXd A (2, 2); A << -1, 3, -1, -1;
  // VectorXd b (2); b << 6, -1;
  // VectorXd c (2); c << 6, 1;
  // VectorXd u (2); u << 10, 10;
  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
  cout << "START " << endl;
  Simplex simplex = Simplex(16, A, b, c, u);
  cout << simplex.status << endl;
  double t = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count() / 1000000.0;
  // cout << "IN " << t << "ms" << endl;
  // double obj_value = 0;
  // for (int i = 0; i < N; i ++) obj_value += simplex.tableau[simplex.numcols + i]*c(i);
  // cout << "MY OBJ " << obj_value << endl;
  // reportGurobi(A, b, c, 0, 1, 2, false);
}