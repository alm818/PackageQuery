#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <queue>
#include <random>
#include <chrono>
#include <numeric>

#include "fmt/core.h"
#include "fitsio.h"
#include "omp.h"
#include "pb/util/udeclare.h"
#include "pb/util/unumeric.h"
#include "pb/util/udebug.h"
#include "pb/core/map_sort.h"
#include "pb/det/det_prob.h"
#include "pb/core/gurobi_solver.h"
#include "pb/core/pseudo_walker.h"
#include "pb/core/dual.h"

#define CORE_COUNT2 16

using namespace std;
using namespace Eigen;

#define ta(row, col) (primal.tableau[primal.numcols*(row)+(col)])

void generateBoundedProlem2(int expected_n, double outlier_prob, double att_var, int n, RMatrixXd& A, VectorXd& bl, VectorXd& bu, VectorXd& c){
  A.resize(4, n); bl.resize(4); bu.resize(4); c.resize(n); 
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  double multiplicity = att_var * 12;
  //normal_distribution att_dist(10.0, att_var);
  int expected_numvar = expected_n;
  double mean = 0.5*multiplicity*expected_numvar;
  double var = att_var*expected_numvar;
  normal_distribution n_dist(0.0, var);
  normal_distribution n_dist_c(0.0, 200.0);
  #pragma omp parallel for num_threads(CORE_COUNT2)
  for (int i = 0; i < n; i ++){
    A(0, i) = u_dist(gen)*multiplicity;
    A(1, i) = u_dist(gen)*multiplicity;
    A(2, i) = u_dist(gen)*multiplicity;
    A(3, i) = 1;
    c(i) = n_dist_c(gen);
  }
  double tol = var*sqrt(1/outlier_prob);
  bl(0) = mean - tol;
  bu(0) = mean + tol;
  bl(1) = -DBL_MAX;
  bu(1) = mean + tol;
  bl(2) = mean- tol;
  bu(2) = DBL_MAX;
  bl(3) = expected_n;
  bu(3) = expected_n;
}

void generateBoundedProlem(int expected_n, double outlier_prob, double att_var, int n, MatrixXd& A, VectorXd& bl, VectorXd& bu, VectorXd& c){
  A.resize(4, n); bl.resize(4); bu.resize(4); c.resize(n); 
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  // uniform_real_distribution u_dist(0.0, 1.0);
  normal_distribution att_dist(0.0, att_var);
  int expected_numvar = expected_n;
  double mean = 0;
  double var = att_var*expected_numvar;
  normal_distribution n_dist(0.0, var);
  normal_distribution n_dist_c(0.0, 200.0);
  #pragma omp parallel for num_threads(CORE_COUNT2)
  for (int i = 0; i < n; i ++){
    A(0, i) = att_dist(gen);
    A(1, i) = att_dist(gen);
    A(2, i) = att_dist(gen);
    A(3, i) = 1;
    c(i) = n_dist_c(gen);
  }
  double tol = var*sqrt(1/outlier_prob);
  bl(0) = mean - tol;
  bu(0) = mean + tol;
  bl(1) = -DBL_MAX;
  bu(1) = mean + tol;
  bl(2) = mean- tol;
  bu(2) = DBL_MAX;
  bl(3) = expected_n;
  bu(3) = expected_n;
}

bool verify(VectorXd x, const MatrixXd& A, const VectorXd& bl, const VectorXd& bu, const VectorXd& c, const VectorXd& l, const VectorXd& u){
  int n = l.size();
  VectorXd sol = x(seqN(0, n));
  int m = bl.size();
  for (int i = 0; i < n; i ++){
    if (isLess(sol(i), l(i)) || isGreater(sol(i),u(i))){
      cout << "N:" << i << " VIOLATIONS" << endl;
      return false;
    }
  }
  VectorXd b = A * sol;
  for (int i = 0; i < m; i ++){
    if (isLess(b(i), bl(i)) || isGreater(b(i), bu(i))){
      cout << "M:" << i << " VIOLATIONS" << endl;
      return false;
    }
  }
  return true;
}

void convert(const MatrixXd& A, const VectorXd& bl, const VectorXd& bu, MatrixXd& AA, VectorXd& bb){
  int m = A.innerSize();
  int n = A.outerSize();
  int mm = 0;
  for (int i = 0; i < m; i ++){
    if (bl(i) != -DBL_MAX) mm ++;
    if (bu(i) != DBL_MAX) mm ++;
  }
  AA.resize(mm, n);
  bb.resize(mm);
  int index = 0;
  for (int i = 0; i < m; i ++){
    if (bl(i) != -DBL_MAX){
      for (int j = 0; j < n; j ++) AA(index, j) = -A(i, j);
      bb(index) = -bl(i);
      index ++;
    }
    if (bu(i) != DBL_MAX){
      for (int j = 0; j < n; j ++) AA(index, j) = A(i, j);
      bb(index) = bu(i);
      index ++;
    }
  }
}

void testL3Cache(){
  chrono::high_resolution_clock::time_point start, end;
  int N = 100000000; int r = 10;
  VectorXd u (N); u.fill(0.24);
  VectorXd g1 (N); VectorXd g2 (N); VectorXd g3 (N);
  start = chrono::high_resolution_clock::now();
  // for (int i = 0; i < r; i ++){
  //   VectorXd g1 = u;
  //   VectorXd g2 = u+u;
  //   VectorXd g3 = u+u+u;
  // }
  // L3 Effect
  #pragma omp parallel num_threads(16)
  {
    #pragma omp for nowait
    for (int i = 0; i < N; i ++){
      g1(i) = u(i)+u(i);
    }
    // #pragma omp for nowait
    // for (int i = 0; i < N; i ++){
    //   g2(i) = u(i) + u(i);
    // }
    #pragma omp for nowait
    for (int i = 0; i < N; i ++){
      g3(i) = u(i) + u(i) + u(i);
    }
  }
  end = chrono::high_resolution_clock::now();
  double exe = chrono::duration_cast<chrono::nanoseconds>(end - start).count() / 1000000.0;
  cout << exe << endl;

  VectorXd l1 (N); VectorXd l2 (N); VectorXd l3 (N);
  cout << "START" << endl;
  start = chrono::high_resolution_clock::now();
  VectorXd uu = u+u;
  VectorXd uuu = u+u+u;
  // #pragma omp parallel num_threads(16)
  // {
  //   #pragma omp for
  //   for (int i = 0; i < N; i ++){
  //     for (int j = 0; j < r; j ++){
  //       l1(i) = u(i);
  //       l2(i) = u(i)+u(i);
  //       l3(i) = u(i)+u(i)+u(i);
  //     }
  //   }
  // }
  end = chrono::high_resolution_clock::now();
  double exe2 = chrono::duration_cast<chrono::nanoseconds>(end - start).count() / 1000000.0;
  cout << exe2 << endl;
}

void quickRun(){
  int n = 10000000;
  RMatrixXd A;
  VectorXd bl, bu, c;
  generateBoundedProlem2(20, 0.9, 4, n, A, bl, bu, c);
  VectorXd u (n); u.fill(100);
  VectorXd l (n); l.fill(0);
  // testL3Cache();

  // int n = 2;
  // int m = 4;
  // MatrixXd A (m, n); A << 1.0/3, 1, 4.0/3, -1, -1.0/4, 1, 1, -1;
  // VectorXd bl (m); bl << 8.0/3, -DBL_MAX, 1.0/2, 0;
  // VectorXd bu (m); bu << DBL_MAX, 17.0/3, 3, DBL_MAX;
  // VectorXd c (n); c << 1, 1;
  // VectorXd l (n); l << 0, 0;
  // VectorXd u (n); u << 10, 10;

  // MatrixXd AA; VectorXd bb;
  // convert(A, bl, bu, AA, bb);
  // int mm = bb.size();
  // Simplex primal = Simplex(16, AA, bb, c, u);
  // VectorXd primal_sol (n);
  // for (int i = 0; i < n; i ++) primal_sol(i) = ta(1, i);
  // double obj = 0;
  // for (int i = 0; i < n; i ++){
  //   obj += primal_sol(i) * c(i);
  // }
  // fmt::print("{:.10Lf}\n", obj);
  // for (int i = 0; i < mm; i ++){
  //   cout << primal.bhead[i] << " ";
  // }
  // cout << endl;
  
  // Dual dual1 = Dual(1, A, bl, bu, c, l, u);
  // cout << dual1.exe_solve << endl;
  cout << "START" << endl;
  int m = bl.size();
  DetProb db = DetProb(m, n);
  db.A = A;
  db.bl = bl;
  db.bu = bu;
  db.c = c;
  db.l = l;
  db.u = u;
  Dual dual = Dual(8, db);
  cout << solMessage(dual.status) << endl;
  cout << dual.exe_solve << endl;
  fmt::print("{:.10Lf}\n", dual.score);
  shortPrint(dual.sol);
  cout << "----------------------" << endl;
  GurobiSolver gs = GurobiSolver(db);
  gs.solveLp();
  cout << gs.getLpTime() << endl;
  // cout << "AMDA1:" << 8.0/7.0 * (1 - dual.exe_solve / dual1.exe_solve) << endl;
  // cout << "AMDA2:" << 8.0/7.0 * (1 - (gs.exe_relaxed+gs.exe_init) / dual1.exe_solve) << endl;
  fmt::print("{:.10Lf}\n", gs.lp_score);
  shortPrint(gs.lp_sol);
  // // print(dual.x);
  // // print(gs.r0);
  // cout << verify(dual.x, A, bl, bu, c, l, u) << " " << verify(gs., A, bl, bu, c, l, u) << endl;
  if (gs.lp_status == Found && !isEqual(gs.lp_score, dual.score, 1e-4)){
    for (int i = 0; i < n; i ++){
      assert(isEqual(dual.sol(i), gs.lp_sol(i), 1e-6));
    }
  }
}

void mapSort(){
  vector<string> names = {"1", "2"};
  Profiler pro = Profiler(names);
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  int N = 10000000;
  int repeat = 1;
  map_sort::MapSort<pair<double,int>> map_sort;
  map_sort.Init(N);
  for (int j = 0; j < repeat; j ++){
    pair<double,int>* x = new pair<double,int>[N];
    for (int i = 0; i < N ; i ++) x[i] = {u_dist(gen), i};
    pair<double,int>* y = new pair<double,int>[N];
    for (int i = 0; i < N ; i ++) y[i] = {x[i].first, i};
    pro.clock(0, false);
    map_sort.Sort(x, N, 16);
    pro.stop(0, false);
    pro.clock(1, false);
    sort(y, y+N);
    pro.stop(1, false);
    for (int i = 0; i < N; i ++) assert(x[i].first == y[i].first);
    delete[] x;
    delete[] y;
  }
  pro.print();
}

void testCopy(){
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  vector<string> names = {"1", "2", "3", "4", "5"};
  Profiler pro = Profiler(names);
  int N = 10000000;
  MatrixXd A (6, N); 
  for (int i = 0; i < 6; i ++){
    #pragma omp parallel for num_threads(16)
    for (int j = 0; j < N; j ++){
      A(i,j) = u_dist(gen);
    }
  }
  pro.clock(0);
  MatrixXd B = A;
  pro.stop(0);
  pro.clock(1);
  MatrixXd BB (6, N);
  memcpy(&BB(0,0), &A(0,0), 6*N);
  pro.stop(1);
  pro.print();
}

int main(){
  //mapSort();
  // for (int i = 0; i < 10000; i ++){
  //   quickRun();
  //   cout << "--------------------------------" << endl;
  // }
  quickRun();
  //testCopy();
}