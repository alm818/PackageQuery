#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <random>
#include <numeric>

#include "utility.h"
#include "lattice_solver.h"
#include "gurobi_lattice_solver.h"
#include "gurobi_solver.h"
#include "fmt/core.h"
#include "fitsio.h"
#include "omp.h"

using namespace std;
using namespace Eigen;

void test_time_utilization(){
  int core = 80;
  GalaxyDB gdb = GalaxyDB(100000000);
  MatrixXd A;
  VectorXd b, c, u;
  int n_partitions = 10;
  int repeat = 3;
  MatrixXd benchmark (n_partitions, 6); benchmark.fill(0);
  for (int i = 0; i < n_partitions; i++){
    double percent = (i+1)/10.0;
    int cur = repeat;
    for (int j = 0; j < repeat; j ++){
      gdb.generateQuery(percent, -1, GalaxyDB::Q1, 10, A, b, c, u);
      LatticeSolver ls = LatticeSolver(core, A, b, c, u);
      ls.solve(2);
      if (ls.status == LS_FOUND){
        benchmark(i, 0) += ls.exe_init;
        benchmark(i, 1) += ls.exe_relaxed;
        benchmark(i, 2) += ls.exe_tableau;
        benchmark(i, 3) += ls.exe_find_dir * ls.try_count / core;
        benchmark(i, 4) += ls.exe_init_walk * ls.try_count / core;
        benchmark(i, 5) += ls.exe_walk * ls.try_count / core;
      } else cur --;
      fmt::print("Finished pair ({},{})\n", i, j);
    }
    for (int j = 0; j < 6; j ++) benchmark(i, j) /= cur;
  }
  cout << benchmark << endl;
  //gdb.generateQuery(1, 3, GalaxyDB::Q2, 10, A, b, c, u);
}

void test_fail_prob(){
  int core = 80;
  GalaxyDB gdb = GalaxyDB(100000000);
  MatrixXd A;
  VectorXd b, c, u;
  vector<double> sz = {1000000, 2000000, 4000000, 6000000, 8000000, 10000000};
  int n_partitions = sz.size();
  vector<int> repeats = {20, 20, 20, 10, 10, 10};
  MatrixXd benchmark (n_partitions, 2); benchmark.fill(0);
  for (int i = 0; i < n_partitions; i++){
    double percent = sz[i] / gdb.num_rows;
    benchmark(i, 1) = percent;
    int cnt = 0;
    for (int j = 0; j < repeats[i]; j ++){
      gdb.generateQuery(percent, -1, GalaxyDB::Q1, 10, A, b, c, u);
      LatticeSolver ls = LatticeSolver(core, A, b, c, u);
      ls.solve(2);
      if (ls.status != LS_FOUND){
        cnt ++;
      }
      fmt::print("Finished pair ({},{})\n", i, j);
    }
    benchmark(i, 0) = cnt/(double)repeats[i];
  }
  cout << benchmark << endl;
}

void compare_ilp(){
  int core = 80;
  GalaxyDB gdb = GalaxyDB(100000000);
  MatrixXd A;
  VectorXd b, c, u;
  vector<double> sz = {200000, 400000, 600000, 800000, 1000000, 2000000, 4000000, 6000000, 8000000, 10000000};
  int n_partitions = sz.size();s  1
  int repeat = 1;
  MatrixXd benchmark (n_partitions, 4); benchmark.fill(0);
  for (int i = 0; i < n_partitions; i ++){
    double percent = sz[i] / gdb.num_rows;
    benchmark(i, 3) = percent;
    if (c.size() <= 1000000) repeat = 1;
    else repeat = 3;
    int cur = repeat;
    for (int j = 0; j < repeat; j ++){
      //gdb.generateQuery(percent, -1, GalaxyDB::Q1, 10, A, b, c, u);
      gdb.generateQuery(percent, 3, GalaxyDB::Q2, 10, A, b, c, u);
      LatticeSolver ls = LatticeSolver(core, A, b, c, u);
      ls.solve(2);
      fmt::print("Finished pair ({},{}) using LS\n", i, j);
      if (ls.status == LS_FOUND){
        benchmark(i, 0) += ls.exe_init + ls.exe_relaxed + ls.exe_tableau + ls.exe_solved;
        if (c.size() <= 1000000){
          GurobiSolver gs = GurobiSolver(A, b, c, u);
          gs.solveIlp();
          benchmark(i, 1) += gs.exe_init + gs.exe_ilp;
          benchmark(i, 2) += gs.ilp_cscore / ls.best_cscore;
        }
      } else cur --;
      fmt::print("Finished pair ({},{})\n", i, j);
    }
    for (int j = 0; j < 3; j ++) benchmark(i, j) /= cur;
  }
  cout << benchmark << endl;
}

int main(){
  //test_time_utilization();
  //compare_ilp();
  test_fail_prob();
}