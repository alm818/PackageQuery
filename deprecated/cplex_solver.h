#pragma once

#include "pb/det/det_prob.h"
#include "pb/util/udeclare.h"
#include "ilcplex/cplexx.h"
#include "checker.h"

using namespace pb;

class CplexSolver: public Checker{
private:
  double exe_init, exe_lp, exe_ilp;
  CPXENVptr env;
  CPXLPptr model;
public:
  int lp_status, ilp_status;
  VectorXd lp_sol, ilp_sol;
  double lp_score, ilp_score;
public:
  ~CplexSolver();
  CplexSolver(const DetProb &prob);
  void solveIlp();
  void solveLp();
  double getLpTime();
  double getIlpTime();
  double getScore(const VectorXd &sol);
  int checkLpFeasibility(const VectorXd &sol);
  int checkIlpFeasibility(const VectorXd &sol);
};