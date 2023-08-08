#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <chrono>
#include <random>

#include "utility.h"
#include "lattice_solver.h"
#include "gurobi_lattice_solver.h"
#include "gurobi_solver.h"
#include "fmt/core.h"
#include "simplex.h"

using namespace std;
using namespace Eigen;

void generateProlem(int n, MatrixXd& A, VectorXd& b, VectorXd& c){
  A.resize(6, n); b.resize(6); c.resize(n);
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  int expected_numvar = 10;
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
  b(4) = 30;
  b(5) = -5;
}

double gurobiILP(const MatrixXd& A, const VectorXd& b, const VectorXd& c, const VectorXd& u, bool is_concurrent){
  GurobiSolver gs = GurobiSolver(A, b, c, u);
  gs.solveIlp(is_concurrent);
  return gs.exe_init + gs.exe_ilp;
}

double gurobiSimplex(const MatrixXd& A, const VectorXd& b, const VectorXd& c, const VectorXd& u, bool is_concurrent){
  GurobiSolver gs = GurobiSolver(A, b, c, u);
  gs.solveRelaxed(is_concurrent);
  return gs.exe_init + gs.exe_relaxed;
}

double parallelSimplex(int core, const MatrixXd& A, const VectorXd& b, const VectorXd& c, const VectorXd& u){
  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
  Simplex simplex = Simplex(core, A, b, c, u);
  double res = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count() / 1000000.0;
  return res;
}

int main(){
  vector<int> cores = {1, 2, 4, 8, 16, 32, 64, 80};
  int low = 2; int high = 7;
  MatrixXd benchmarks (high-low+1, cores.size()+4); benchmarks.fill(-1);
  for (int i = low; i <= high; i ++){
    int n = (int) pow(10, i);  
    int m = 6;
    MatrixXd A (m, n);
    VectorXd b (m);
    VectorXd c (n);
    generateProlem(n, A, b, c);
    VectorXd u (n); u.fill(1);
    if (i <= 6){
      benchmarks(i-2, 0) = gurobiILP(A, b, c, u, false);
      benchmarks(i-2, 1) = gurobiILP(A, b, c, u, true);
    }
    benchmarks(i-2, 2) = gurobiSimplex(A, b, c, u, false);
    benchmarks(i-2, 3) = gurobiSimplex(A, b, c, u, true);
    for (int j = 0; j < (int) cores.size(); j++){
      benchmarks(i-2, j+4) = parallelSimplex(cores[j], A, b, c, u);
    }
    fmt::print("Finished {}x{} problem\n", m, n);
  }
  cout << benchmarks << endl;
}