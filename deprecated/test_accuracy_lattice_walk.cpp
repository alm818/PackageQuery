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

double getProjDistance(VectorXd x, VectorXd p) {
  double v = x.dot(p) / p.norm();
  v *= v;
  return sqrt(abs(x.squaredNorm() - v));
}

void all_three(){
  int log_step_start = 1;
  int log_step_end = 12;
  int log_dim_start = 1;
  int log_dim_end = 12;
  for (int log_dim = log_dim_start; log_dim <= log_dim_end; log_dim++) {
    for (int log_step = log_step_start; log_step <= log_step_end; log_step++) {
      int dim = (int)pow(2, log_dim);
      VectorXd p(dim);
      for (int i = 0; i < dim; i++) p(i) = dist(gen);
      VectorXd pp = p / p.norm();
      double max_abs_val = -1;
      for (int i = 0; i < dim; i++) {
        max_abs_val = max(max_abs_val, abs(pp(i)));
      }
      double conj = 0.5 / max_abs_val;
      int step = (int)pow(2, log_step);
      Walker w = Walker(p);
      PseudoWalker pw = PseudoWalker(p, false);
      PseudoWalker pw2 = PseudoWalker(p, true);
      VectorXd x(p.size()); x.fill(0);
      VectorXd px(p.size()); px.fill(0);
      VectorXd px2(p.size()); px2.fill(0);
      double avg_score = 0;
      double avg_pscore = 0;
      double avg_pscore2 = 0;
      for (int n = 0; n < step; n++) {
        int d = w.step();
        if (d != 0) {
          int i = abs(d) - 1;
          x(i) += sign(d);
        }
        avg_score += getProjDistance(x, p);
        int pd = pw.step();
        if (pd != 0) {
          int i = abs(pd) - 1;
          px(i) += sign(pd);
        }
        avg_pscore += getProjDistance(px, p);
        int pd2 = pw2.step();
        if (pd2 != 0) {
          int i = abs(pd2) - 1;
          px2(i) += sign(pd2);
        }
        avg_pscore2 += getProjDistance(px2, p);
      }
      avg_score /= step;
      avg_pscore /= step;
      avg_pscore2 /= step;
      double advantage = (avg_pscore - avg_score) / avg_score;
      double advantage2 = (avg_pscore2 - avg_score) / avg_score;
      fmt::print("{},{},{:.10Lf},{:.10Lf},{:.2Lf},{:.10Lf},{:.2Lf},{:.5Lf}\n", log_step, log_dim, avg_score, avg_pscore, advantage*100, avg_pscore2, advantage2*100, conj);
      //fmt::print("log_step={} log_dim={} avg_score={:.10Lf} avg_pscore={:.10Lf} advantage={:.2Lf}% avg_pscore2={:.10Lf} advantage2={:.2Lf}% conj={:.5Lf}\n", log_step, log_dim, avg_score, avg_pscore, advantage*100, avg_pscore2, advantage2*100, conj);
    }
  }
}

void only_pseudo(){
  int log_step_start = 2;
  int log_step_end = 5;
  int log_dim_start = 2;
  int log_dim_end = 5;
  for (int log_dim = log_dim_start; log_dim <= log_dim_end; log_dim++) {
    for (int log_step = log_step_start; log_step <= log_step_end; log_step++) {
      int dim = (int)pow(10, log_dim);
      VectorXd p(dim);
      for (int i = 0; i < dim; i++) p(i) = dist(gen);
      VectorXd pp = p / p.norm();
      double max_abs_val = -1;
      for (int i = 0; i < dim; i++) {
        max_abs_val = max(max_abs_val, abs(pp(i)));
      }
      double conj = 0.5 / max_abs_val;
      int step = (int)pow(10, log_step);
      PseudoWalker pw = PseudoWalker(p, true, 2);
      PseudoWalker pw2 = PseudoWalker(p, true, 16);
      VectorXd px(p.size()); px.fill(0);
      VectorXd px2(p.size()); px2.fill(0);
      double avg_pscore = 0;
      double avg_pscore2 = 0;
      for (int n = 0; n < step; n++) {
        int pd = pw.step();
        if (pd != 0) {
          int i = abs(pd) - 1;
          px(i) += sign(pd);
        }
        avg_pscore += getProjDistance(px, p);
        int pd2 = pw2.step();
        if (pd2 != 0) {
          int i = abs(pd2) - 1;
          px2(i) += sign(pd2);
        }
        avg_pscore2 += getProjDistance(px2, p);
      }
      avg_pscore /= step;
      avg_pscore2 /= step;
      double advantage = (avg_pscore - avg_pscore2) / avg_pscore2;
      fmt::print("{},{},{:.10Lf},{:.10Lf},{:.5Lf}\n", log_step, log_dim, avg_pscore, avg_pscore2, advantage*100);
      //fmt::print("log_step={} log_dim={} avg_score={:.10Lf} avg_pscore={:.10Lf} advantage={:.2Lf}% avg_pscore2={:.10Lf} advantage2={:.2Lf}% conj={:.5Lf}\n", log_step, log_dim, avg_score, avg_pscore, advantage*100, avg_pscore2, advantage2*100, conj);
    }
  }
}

void only_pseudo_multicore(){
  vector<string> names = {"1", "2", "3", "4", "5"};
  int log_step_start = 2;
  int log_step_end = 5;
  int log_dim_start = 2;
  int log_dim_end = 7;
  for (int log_dim = log_dim_start; log_dim <= log_dim_end; log_dim++) {
    int dim = (int)pow(10, log_dim);
    Profiler pro = Profiler(names);
    VectorXd p(dim);
    for (int i = 0; i < dim; i++) p(i) = dist(gen);
    VectorXd pp = p / p.norm();
    double max_abs_val = -1;
    for (int i = 0; i < dim; i++) {
      max_abs_val = max(max_abs_val, abs(pp(i)));
    }
    double conj = 0.5 / max_abs_val;
    int step = (int)dim/log2(dim)+1;
    pro.clock(0, false);
    PseudoWalker pw = PseudoWalker(p, true, 1);
    pro.stop(0, false);
    pro.clock(1, false);
    PseudoWalker pw2 = PseudoWalker(p, true, 16);
    pro.stop(1, false);
    VectorXd px(p.size()); px.fill(0);
    VectorXd px2(p.size()); px2.fill(0);
    double avg_pscore = 0;
    double avg_pscore2 = 0;
    cout << "START" << endl;
    for (int n = 0; n < step; n++) {
      pro.clock(2, false);
      int pd = pw.step();
      pro.stop(2, false);
      if (pd != 0) {
        int i = abs(pd) - 1;
        px(i) += sign(pd);
      }
      pro.clock(3, false);
      int pd2 = pw2.step();
      pro.stop(3, false);
      if (pd2 != 0) {
        int i = abs(pd2) - 1;
        px2(i) += sign(pd2);
      }
      assert(pd == pd2);
    }
    cout << step << " " << log_dim << endl;
    pro.print();
    // avg_pscore /= step;
    // avg_pscore2 /= step;
    // double advantage = (avg_pscore - avg_pscore2) / avg_pscore2;
    // fmt::print("{},{},{:.10Lf},{:.10Lf},{:.5Lf}\n", log_step, log_dim, avg_pscore, avg_pscore2, advantage*100);
    //fmt::print("log_step={} log_dim={} avg_score={:.10Lf} avg_pscore={:.10Lf} advantage={:.2Lf}% avg_pscore2={:.10Lf} advantage2={:.2Lf}% conj={:.5Lf}\n", log_step, log_dim, avg_score, avg_pscore, advantage*100, avg_pscore2, advantage2*100, conj);
  }
}

int main(){
  // all_three();
  //only_pseudo();
  only_pseudo_multicore();
}