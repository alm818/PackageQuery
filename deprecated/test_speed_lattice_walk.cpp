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

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

boost::random::mt19937 gen;
boost::random::uniform_int_distribution<> dist(-1000, 1000);

using namespace std;
using namespace Eigen;

void RandomWalk(int N, VectorXd p) {
  Walker w = Walker(p);
  for (int i = 0; i < N; i++) w.step();
}

void RandomPseudoWalk(int N, VectorXd p, bool enable_correction) {
  PseudoWalker w = PseudoWalker(p, enable_correction);
  for (int i = 0; i < N; i++) w.step();
}

void test_lattice_walk(){
  int log_step_start = 2;
  int log_step_end = 5;
  int log_dim_start = 2;
  int log_dim_end = 5;
  for (int log_dim = log_dim_start; log_dim <= log_dim_end; log_dim++){
    int dim = (int)pow(10, log_dim);
    VectorXd p(dim);
    for (int i = 0; i < dim; i++) p(i) = dist(gen);
    for (int log_step = log_step_start; log_step <= log_step_end; log_step++) {
      int step = (int)pow(10, log_step);
      double exe_time = exeTime(RandomWalk, step, p);
      fmt::print("{},{},{:.5Lf}\n",log_step, log_dim, exe_time);
      //fmt::print("log_step={} log_dim={} time={:.5Lf}ms\n",log_step, log_dim, exe_time);
    }
  }
}

void test_pseudo_lattice_walk(){
  int log_step_start = 2;
  int log_step_end = 7;
  int log_dim_start = 2;
  int log_dim_end = 7;
  for (int log_dim = log_dim_start; log_dim <= log_dim_end; log_dim++) {
    int dim = (int)pow(10, log_dim);
    VectorXd p(dim);
    for (int i = 0; i < dim; i++) p(i) = dist(gen);
    for (int log_step = log_step_start; log_step <= log_step_end; log_step++) {
      int step = (int)pow(10, log_step);
      double exe_time = exeTime(RandomPseudoWalk, step, p, false);
      fmt::print("{},{},{:.5Lf}\n",log_step, log_dim, exe_time);
      //Logger::WriteMessage(format("log_step={} log_dim={} time={}ms", log_step, log_dim, time / 1000000.0).c_str());
    }
  }
}

void test_pseudo_lattice_walk_with_correction(){
  int log_step_start = 2;
  int log_step_end = 7;
  int log_dim_start = 2;
  int log_dim_end = 7;
  for (int log_dim = log_dim_start; log_dim <= log_dim_end; log_dim++) {
    int dim = (int)pow(10, log_dim);
    VectorXd p(dim);
    for (int i = 0; i < dim; i++) p(i) = dist(gen);
    for (int log_step = log_step_start; log_step <= log_step_end; log_step++) {
      int step = (int)pow(10, log_step);
      double exe_time = exeTime(RandomPseudoWalk, step, p, true);
      fmt::print("{},{},{:.5Lf}\n",log_step, log_dim, exe_time);
      //Logger::WriteMessage(format("log_step={} log_dim={} time={}ms", log_step, log_dim, time / 1000000.0).c_str());
    }
  }
}

int main(){
  //test_lattice_walk();
  //test_pseudo_lattice_walk();
  //test_pseudo_lattice_walk_with_correction();
}