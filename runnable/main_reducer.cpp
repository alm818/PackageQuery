#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <queue>
#include <random>
#include <chrono>
#include <numeric>

#include "utility.h"
#include "parallel_pq.h"
#include "fmt/core.h"
#include "fitsio.h"
#include "omp.h"
#include "lattice_solver.h"
#include "pseudo_walker.h"
#include "simplex.h"
#include "gurobi_lattice_solver.h"
#include "gurobi_solver.h"
#include "reducer.h"
#include "fmt/core.h"

using namespace std;
using namespace Eigen;

#define tab(row, col) (simplex.tableau[simplex.numcols*(row)+(col)])

void generateProlem(int n, MatrixXd& A, VectorXd& b, VectorXd& c){
  A.resize(6, n); b.resize(6); c.resize(n);
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  int expected_numvar = 500;
  normal_distribution n_dist(0.5*expected_numvar, 1.0/12*expected_numvar);
  normal_distribution n_dist_c(0.0, 500.0);
  #pragma omp parallel for num_threads(CORE_COUNT)
  for (int i = 0; i < n; i ++){
    A(0, i) = u_dist(gen);
    A(1, i) = -A(0, i);
    A(2, i) = u_dist(gen);
    A(3, i) = -A(2, i);
    A(4, i) = 1;
    A(5, i) = -1;
    c(i) = n_dist_c(gen);
    //c(i) = i;
  }
  double tol_scale = 5;
  double tol1 = u_dist(gen) * tol_scale / 2;
  double tol2 = u_dist(gen) * tol_scale / 2;
  b(0) = n_dist(gen) + tol1;
  b(1) = tol1 - b(0);
  b(2) = n_dist(gen) + tol2;
  b(3) = tol2 - b(2);
  b(4) = 1000;
  b(5) = -5;
}

void test_ppw(){
  // unordered_map<int, double> test;
  // test[0] += 1;
  // test[2] += 100;
  // for (auto pi : test){
  //   cout << pi.first << " " << pi.second << endl;
  // }
  // int T = 500000;
  // VectorXd coefs (n); coefs.fill(1.0/n);
  // VectorXd res (n); res.fill(0);
  // LatticeSolver ls = LatticeSolver(80, A, b, c, u);
  // ls.lattice_dirs.inplaceVectorProduct(coefs, res);
  // vector<int> r (T);
  // vector<int> s (T);
  // cout << "START" << endl;
  // chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
  // PseudoWalker pw1 = PseudoWalker(res, true, 1);
  // for (int i = 0; i < T; i ++){
  //   r[i] = pw1.step();
  // }
  // chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
  // double exe1 = chrono::duration_cast<chrono::nanoseconds>(end - start).count() / 1000000.0;
  // cout << "DONE1 " << exe1 << endl;
  // start = chrono::high_resolution_clock::now();
  // PseudoWalker pw2 = PseudoWalker(res, true, 80);
  // for (int i = 0; i < T; i ++){
  //   s[i] = pw2.step();
  // }
  // end = chrono::high_resolution_clock::now();
  // double exe2 = chrono::duration_cast<chrono::nanoseconds>(end - start).count() / 1000000.0;
  // cout << "DONE2 " << exe2  << endl;
  // for (int i = 0; i < T; i ++){
  //   if (r[i] != s[i]) cout << "WTF " << endl;
  // }
  // cout << res.norm() << endl;
}

int main(){
  int n = 1000000;
  int m = 6;
  MatrixXd A (m, n);
  VectorXd b (m);
  VectorXd c (n);
  generateProlem(n, A, b, c);
  VectorXd u (n); u.fill(1);

  // GalaxyDB gdb = GalaxyDB(100000000);
  // MatrixXd A;
  // VectorXd b, c, u;
  // gdb.generateQuery(1, 3, GalaxyDB::Q2, 10, A, b, c, u);
  // int n = c.size();
  // int m = b.size();

  cout << "START REDUCER" << endl;
  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
  Reducer reducer = Reducer(80, &A, &b, &c, &u);
  chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
  double exe = chrono::duration_cast<chrono::nanoseconds>(end - start).count() / 1000000.0;
  cout << exe << endl;

  // start = chrono::high_resolution_clock::now();
  // GurobiSolver gs = GurobiSolver(A, b, c, u);
  // gs.solveRelaxed();
  // cout << gs.iteration_count << endl;
  // end = chrono::high_resolution_clock::now();
  // double exe2 = chrono::duration_cast<chrono::nanoseconds>(end - start).count() / 1000000.0;
  // cout << exe2 << endl;

  // VectorXd r0 (n);
  // for (int i = 0; i < n; i ++) r0(i) = tab(1, i);
  // cout << c.dot(r0) << " " << c.dot(gs.r0) << endl;
}