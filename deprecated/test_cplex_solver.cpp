#include "pb/util/udeclare.h"
#include "pb/util/udebug.h"
#include "pb/core/cplex_solver.h"

using namespace pb;

void testProb(){
  DetProb prob (2, 2);
  prob.A << 3, 1, -1,  2;
  prob.bl << -DBL_MAX, -DBL_MAX;
  prob.bu << 11, 5;
  prob.u.fill(10);
  prob.c << 6, 5;
  CplexSolver cs = CplexSolver(prob);
  cs.solveLp();
  print(cs.lp_sol);
  fmt::print("Time:{:.2Lf}ms LP_Score:{:.4Lf}\n", cs.getLpTime(), cs.lp_score);
  fmt::print("SolStatus:{} LPStatus:{} ILPStatus:{}\n", solMessage(cs.lp_status), feasMessage(cs.checkLpFeasibility(cs.lp_sol)), feasMessage(cs.checkIlpFeasibility(cs.lp_sol)));
  cs.solveIlp();
  print(cs.ilp_sol);
  fmt::print("Time:{:.2Lf}ms ILP_Score:{:.4Lf}\n", cs.getIlpTime(), cs.ilp_score);
  fmt::print("SolStatus:{} LPStatus:{} ILPStatus:{}\n", solMessage(cs.ilp_status), feasMessage(cs.checkLpFeasibility(cs.ilp_sol)), feasMessage(cs.checkIlpFeasibility(cs.ilp_sol)));
}

int main(){
  testProb();
}