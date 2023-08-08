#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <random>

#include "utility.h"
#include "pseudo_walker.h"
#include "lattice_solver.h"
#include "gurobi_lattice_solver.h"
#include "gurobi_solver.h"
#include "fmt/core.h"
#include "CCfits/CCfits.h"

using namespace std;
using namespace Eigen;

#define ta(row, col) (ls.simplex->tableau[ls.simplex->numcols*(row)+(col)])

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

void test(){
  int n = 200000;
  int m = 6;
  MatrixXd A (m, n);
  VectorXd b (m);
  VectorXd c (n);
  generateProlem(n, A, b, c);
  VectorXd u (n); u.fill(1);
  // LatticeSolver ls = LatticeSolver(1, A, b, c, u);
  // int B = 100;
  // VectorXd histo (B); histo.fill(0);
  // for (int i = 0; i < n; i ++){
  //   if (ls.r0(i) < 0.01){
  //     int index = (int)floor(ls.r0(i)*10000);
  //     if (index == B) index --;
  //     histo[index] ++;
  //   }
  // }
  // print(histo);
  LatticeSolver ls2 = LatticeSolver(80, A, b, c, u);
  ls2.solve(0.5);
  ls2.report();
  // ls2.report();
  // cout << "--------------------------------------------" << endl;
  // LatticeSolver ls2 = LatticeSolver(80, A, b, c, u);
  // cout << "___________________" << endl;
  // for (int i = 0; i < m; i ++) cout << ls.bhead[i] << " ";
  // cout << endl;
  // print(ls.c);
  // print(ls.r0);
  // ls.lattice_dirs.print();
  // print(ls.fcv->sample());
  // ls.solve(2);
  // ls.report();
  ////////////////////
  // MatrixXd A(5, 3); A << 0, -4, 4, 2, -4, 4, 2, 4, 0, 0, 4, 0, -4, 0, 4;
  // VectorXd b(5); b << 0, 4, 12, 8, 0;
  // VectorXd c(3); c << 1, 0, 0;
  // VectorXd u (3); u.fill(10);
  // LatticeSolver ls = LatticeSolver(80, A, b, c, u);
  // ls.lattice_dirs.print();
  // cout << "GLS" << endl;
  // GurobiLatticeSolver gls = GurobiLatticeSolver(A, b, c, 0, 1);
  // VectorXd rr = gls.r0 - ls.r0;
  // cout << "WTF " << rr.norm() << endl;
  // cout << "OBJ " << c.dot(gls.r0) << " " << c.dot(ls.r0) << endl;
  // cout << "NNS " << ls.lattice_dirs2.getSize() << " " << gls.lattice_dirs.getSize() << endl;
  // print(gls.r0);
  // cout << endl;
  // gls.lattice_dirs.print();
  // cout << endl;
  // default_random_engine gen {static_cast<long unsigned int>(time(0))};
  // uniform_real_distribution u_dist(0.0, 1.0);
  // VectorXd x(n);
  // for (int i = 0; i < n; i ++) x(i) = u_dist(gen);
  // x /= x.norm();
  // VectorXd res1 (n); res1.fill(0);
  // ls.lattice_dirs.inplaceVectorProduct(x, res1);
  // cout << "NBASIC:";
  // for (int i = 0; i < n; i ++){
  //   cout << ls.nbasic[i] << " ";
  //   int comp = abs(ls.nbasic[i]) - 1;
  //   if (ls.nbasic[i] != 0) res1(comp) = sign(ls.nbasic[i]) * x(i);
  // }
  // cout << endl;
  // VectorXd res2 = gls.lattice_dirs.vectorProduct(x);
  // VectorXd res3 = ls.lattice_dirs2.vectorProduct(x);
  // print(x);
  // print(res1);
  // print(res2);
  // print(res3);
  // double error = 0;
  // int cnt = 0;
  // unordered_map<int, int> gls_inv_bhead;
  // for (int i = 0; i < m; i ++) gls_inv_bhead[gls.bhead[i]] = i;
  // cout << "THEIR NOT EXISTS: ";
  // for (int i = 0; i < m; i ++){
  //   int basic_col = gls.bhead[i];
  //   if (ls.inv_bhead.find(basic_col) == ls.inv_bhead.end()){
  //     cout << basic_col << " ";
  //   }
  // }
  // cout << endl;
  // cout << "MY NOT EXISTS: ";
  // for (int i = 0; i < m; i ++){
  //   int basic_col = ls.bhead[i];
  //   if (gls_inv_bhead.find(basic_col) == gls_inv_bhead.end()){
  //     cout << basic_col << " ";
  //   }
  // }
  // cout << endl;
  // for (int i = 0; i < m; i ++){
  //   int basic_col = gls.bhead[i];
  //   if (ls.inv_bhead.find(basic_col) != ls.inv_bhead.end()){
  //     cnt ++;
  //     int my_row = ls.inv_bhead[basic_col];
  //     for (int j = 0; j < n+m; j ++){
  //       error += abs(ta(my_row+2, j)-gls.tableau(i, j));
  //     }
  //   }
  // }
  // error /= (cnt*(m+n));
  // cout << error << endl;
  // gls.solve(2);
  // gls.milpSolve();
  // if (gls.status == LS_FOUND) ls.compareReport(gls.x0);
  // else ls.report();
  // if (gls.status == LS_FOUND) ls2.compareReport(gls.best_x);
  // else ls2.report();
  // fmt::print("{:.8Lf} {:.8Lf} {:.8Lf} {:.8Lf} {:.8Lf} {}\n", c.dot(ls2.r0), c.dot(gls.r0), gls.getObjValue(gls.r0), ls2.getObjValue(gls.r0), gls.best_c_score, solMessage(gls.status));
  // fmt::print("{} {} {}\n", gls.try_count, gls.avg_step_count, gls.exe_relaxed);
  // cout << "HERE " << c.dot(gls.x0) << endl;
}

void quickRun(){
  int n = 200000;
  int m = 6;
  MatrixXd A (m, n);
  VectorXd b (m);
  VectorXd c (n);
  generateProlem(n, A, b, c);
  VectorXd u (n); u.fill(1);
  LatticeSolver ls = LatticeSolver(80, A, b, c, u);
  // showHistogram(ls.r0, 10, 0, 0);
  // ls.solve(2);
  // ls.report();
  // GurobiSolver gs = GurobiSolver(A, b, c, u);
  // gs.solveIlp();
  // ls.compareReport(gs.x0, gs.exe_ilp);
  // VectorXd len (n); len.fill(0);
  // for (int i = 0; i < n; i ++){
  //   for (const auto& pi : ls.lattice_dirs.rows[i]){
  //     len(i) += pi.second*pi.second;
  //   }
  //   len(i) = sqrt(len(i));
  // }
  // showHistogram(len, 10, 0, 0);
}

int main(){
  // quickRun();
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  VectorXd x (10000000);
  for (int i = 0; i < (int) x.size(); i ++) x(i) = u_dist(gen);
  PseudoWalker pw = PseudoWalker(x, true, 16);
  // showHistogram(x, 10, 0, 0.1);
}