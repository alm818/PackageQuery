#pragma once

#include "pb/det/det_prob.h"
#include "pb/det/lsr_prob.h"
#include "pb/det/dlv_partition.h"

#include "pb/util/udebug.h"
#include "pb/util/udeclare.h"
#include "pb/util/upostgres.h"
#include "pb/util/uconfig.h"
#include "pb/util/unumeric.h"
#include "libpq-fe.h"

using namespace pb;

// All the time is in ms
class POPLayeredSketchRefine{
public:
  map<long long, long long> ilp_sol;
  map<long long, double> lp_sol;
  double ilp_score, exe_ilp, lp_score, exe_lp;
  int status;
public:
  ~POPLayeredSketchRefine();
  POPLayeredSketchRefine(DetProb &prob);
};