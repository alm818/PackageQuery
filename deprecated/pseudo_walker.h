#pragma once

#include <queue>
#include <Eigen/Dense>
#include "abstract_walker.h"
#include "parallel_pq.h"

using namespace Eigen;
using namespace std;

class PseudoWalker: public AbstractWalker{

public:
  VectorXd p, x, steps;
  priority_queue<pair<double, int>> pq;
  ParallelPQ* ppq;
  bool enable_correction;
  int n, step_count, correction_count, core;
  double sum_p;

private:
  int executeStep(int i);

public:
  ~PseudoWalker();
  PseudoWalker(VectorXd& p, bool enable_correction=true, int core=1);
  int step();
};
