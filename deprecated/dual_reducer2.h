#pragma once

#include "pb/util/udeclare.h"
#include "pb/core/dual.h"
#include "pb/det/det_prob.h"

using namespace pb;

class DualReducer2{
private:
  int layer_count;
  vector<DetProb*> probs;
  vector<VectorXi*> indices;
  vector<Dual*> duals;
public:
  VectorXd best_sol;
  double best_score, exe_solve;
  int status;
private:
  void filterStay(int core, Dual* dual, DetProb* prob, VectorXi* cur_index, int stay_mode, int stay_count, const VectorXi& stay);
public:
  ~DualReducer2();
  DualReducer2(int core, DetProb *origin_prob, int opt);
  VectorXd getLpSol();
  double getLpScore();
  double getLpTime();
};
