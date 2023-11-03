#pragma once

#include "pb/util/udeclare.h"
#include "pb/core/dual.h"
#include "pb/det/det_prob.h"

using namespace pb;

static const double kMipGap = 1e-4;
static const int kIlpSize = 500;

// All the time is in ms
class DualReducer{
public:
  VectorXd ilp_sol, lp_sol;
  double ilp_score, lp_score, exe_ilp, exe_lp, exe_gb;
  int status, failure_count;
private:
  DetProb* filtering(VectorXi &reduced_index, int core, const DetProb &prob, VectorXd &dual_sol, int stay_count, int stay_mode, vector<int> &stay);
public:
  ~DualReducer();
  DualReducer(int core, const DetProb &prob, bool is_safe=false, double mip_gap=kMipGap, double time_limit=kTimeLimit, int ilp_size=kIlpSize);
};