#include "pb/util/udeclare.h"
#include "pb/util/udebug.h"
#include "pb/core/gurobi_solver.h"
#include "pb/core/checker.h"

using namespace pb;

void testProb(){
  DetProb prob (2, 2);
  prob.A << 3, 1, -1,  2;
  prob.bl << -DBL_MAX, -DBL_MAX;
  prob.bu << 11, 5;
  prob.u.fill(10);
  prob.c << 6, 5;
  GurobiSolver gs = GurobiSolver(prob);
  Checker ch = Checker(prob);
  gs.solveLp();
  print(gs.lp_sol);
  fmt::print("Time:{:.2Lf}ms LP_Score:{:.4Lf}\n", gs.getLpTime(), gs.lp_score);
  fmt::print("BC Obj:{:.10Lf} SolStatus:{} LPStatus:{} ILPStatus:{}\n", ch.getScore(gs.lp_sol), solMessage(gs.lp_status), feasMessage(ch.checkLpFeasibility(gs.lp_sol)), feasMessage(ch.checkIlpFeasibility(gs.lp_sol)));
  fmt::print("GS Obj:{:.10Lf} SolStatus:{} LPStatus:{} ILPStatus:{}\n", gs.getScore(gs.lp_sol), solMessage(gs.lp_status), feasMessage(gs.checkLpFeasibility(gs.lp_sol)), feasMessage(gs.checkIlpFeasibility(gs.lp_sol))); 
  gs.solveIlp();
  print(gs.ilp_sol);
  fmt::print("Time:{:.2Lf}ms ILP_Score:{:.4Lf}\n", gs.getIlpTime(), gs.ilp_score);
  fmt::print("BC Obj:{:.10Lf} SolStatus:{} LPStatus:{} ILPStatus:{}\n", ch.getScore(gs.ilp_sol), solMessage(gs.ilp_status), feasMessage(ch.checkLpFeasibility(gs.ilp_sol)), feasMessage(ch.checkIlpFeasibility(gs.ilp_sol)));
  fmt::print("GS Obj:{:.10Lf} SolStatus:{} LPStatus:{} ILPStatus:{}\n", gs.getScore(gs.ilp_sol), solMessage(gs.ilp_status), feasMessage(gs.checkLpFeasibility(gs.ilp_sol)), feasMessage(gs.checkIlpFeasibility(gs.ilp_sol)));
}

int main(){
  testProb();
}