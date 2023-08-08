#include "pb/util/udeclare.h"
#include "pb/util/udebug.h"
#include "pb/util/uconfig.h"
#include "pb/core/gurobi_solver.h"
#include "pb/core/dual_reducer.h"
#include "pb/det/det_prob.h"

using namespace pb;

// void benchmark(){
//   int n = 10000;
//   bool is_positive = true;
//   double outlier_prob = 10000;
//   int expected_n = 50;
//   double variance = 1000;
//   bool is_translate = true;
//   int rep = 100;
//   VectorXd res (14); res.fill(0);
//   for (int r = 0; r < rep; r++){
//     cout << "r " << r << endl;
//     DetProb prob; prob.uniformGenerate(n, expected_n, variance, outlier_prob, false, is_positive, is_translate);
//     GurobiSolver gs = GurobiSolver(prob);
//     DualReducer dr = DualReducer(kPCore, &prob);
//     if (dr.status == Timeout){
//       res(12) ++;
//       continue;
//     } else if (dr.status == Infeasible || dr.status == DualUnbounded){
//       gs.solveIlp(10.0);
//       if (gs.ilp_status == Infeasible) res(13) ++;
//       else if (gs.ilp_status == Timeout) res(13) ++;
//     }
//     int count = 0;
//     VectorXd sol = dr.getLpSol();
//     for (int i = 0; i < n; i ++){
//       if (sol(i) > 0) count ++;
//     }
//     res(0) += count;
//     //cout << "OK1" << endl;
//     double gs_obj = gs.getScore(dr.best_sol);
//     double lp_obj = gs.getScore(sol);
//     double lp_gap = (lp_obj - gs_obj) / lp_obj * 100;
//     res(1) += lp_gap;
//     for (int i = 0; i < 4; i ++){
//       DualReducer2 dr2 = DualReducer2(kPCore, &prob, i);
//       if (gs.checkIlpFeasibility(dr2.best_sol) == Feasibility){
//         double my_gap = (gs_obj - gs.getScore(dr2.best_sol)) / gs_obj * 100;
//         res(i+2) += my_gap;
//       } else{
//         res(i+7) ++;
//       }
//     }
//     //cout << "OK2" << endl;
//     if (gs.checkIlpFeasibility(dr.best_sol) == Feasibility){
//       double gap = (gs_obj - gs.getScore(dr.best_sol)) / gs_obj * 100;
//       res(6) += gap;
//     } else{
//       res(11) ++;
//     }
//     //cout << "OK4" << endl;
//     VectorXd ress = res;
//     double rr = r+1;
//     ress(0) /= rr;
//     ress(1) /= rr;
//     for (int i = 2; i <= 6; i ++){
//       if (rr-ress(12)-ress(i+5) > 0) ress(i) /= (rr -ress(12) -ress(i+5));
//     }
//     for (int i = 7; i < 12; i ++) ress(i) /= ((rr-ress(12))/100.0);
//     ress(12) /= (rr/100.0);
//     ress(13) /= (rr/100.0);
//     print(ress);
//   }
//   res(0) /= rep;
//   res(1) /= rep;
//   for (int i = 2; i <= 6; i ++){
//     if (rep-res(12)-res(i+5) > 0) res(i) /= (rep - res(12) - res(i+5));
//   }
//   for (int i = 7; i < 12; i ++) res(i) /= ((rep-res(12))/100.0);
//   res(12) /= (rep/100.0);
//   res(13) /= (rep/100.0);
//   print(res);
// }

// void compareGs(){
//   DetProb prob; prob.uniformGenerate(10000, 500, 1, 10000, false, true, true);
//   DualReducer dr = DualReducer(kPCore, &prob);
//   GurobiSolver gs = GurobiSolver(prob);
//   CplexSolver cs = CplexSolver(prob);
//   for (int i = 0; i < 4; i ++){
//     DualReducer2 dr2 = DualReducer2(kPCore, &prob, i);
//     cout << dr2.exe_solve << " " << solMessage(dr2.status) << " " << cs.getScore(dr2.best_sol) << " " << feasMessage(cs.checkIlpFeasibility(dr2.best_sol)) << endl;
//   }
//   gs.solveLp();
//   gs.solveIlp();
//   // cs.solveLp();
//   // cs.solveIlp();
//   if (dr.status == Found) showHistogram(dr.best_sol, 10, 0, 0);
//   if (gs.lp_status == Found) showHistogram(gs.lp_sol, 10, 0, 0);
//   if (gs.ilp_status == Found) showHistogram(gs.ilp_sol, 10, 0, 0);
//   // if (cs.lp_status == Found) showHistogram(gs.lp_sol, 10, 0, 0);
//   // if (cs.ilp_status == Found) showHistogram(gs.ilp_sol, 10, 0, 0);
//   // if (dr.status == Found) cout << solCombination(dr.best_sol) << endl;
//   // cout << solCombination(gs.lp_sol) << endl;
//   // cout << solCombination(gs.ilp_sol) << endl;
//   fmt::print("{:.10Lf} {:.10Lf} {:.10Lf} {:.10Lf} {:.10Lf}\n", dr.best_score, gs.lp_score, gs.ilp_score, cs.lp_score, cs.ilp_score);
//   fmt::print("{:.4Lf} {:.4Lf} {:.4Lf} {:.4Lf} {:.4Lf}\n", dr.exe_solve, gs.getLpTime(), gs.getIlpTime(), cs.getLpTime(), cs.getIlpTime());
//   fmt::print("DRG Obj:{:.10Lf} SolStatus:{} LPStatus:{} ILPStatus:{}\n", gs.getScore(dr.best_sol), solMessage(dr.status), feasMessage(gs.checkLpFeasibility(dr.best_sol)), feasMessage(gs.checkIlpFeasibility(dr.best_sol)));
//   fmt::print("DRC Obj:{:.10Lf} SolStatus:{} LPStatus:{} ILPStatus:{}\n", cs.getScore(dr.best_sol), solMessage(dr.status), feasMessage(cs.checkLpFeasibility(dr.best_sol)), feasMessage(cs.checkIlpFeasibility(dr.best_sol)));
//   fmt::print("GLP Obj:{:.10Lf} SolStatus:{} LPStatus:{} ILPStatus:{}\n", gs.getScore(gs.lp_sol), solMessage(gs.lp_status), feasMessage(gs.checkLpFeasibility(gs.lp_sol)), feasMessage(gs.checkIlpFeasibility(gs.lp_sol)));
//   fmt::print("GILP Obj:{:.10Lf} SolStatus:{} LPStatus:{} ILPStatus:{}\n", gs.getScore(gs.ilp_sol), solMessage(gs.ilp_status), feasMessage(gs.checkLpFeasibility(gs.ilp_sol)), feasMessage(gs.checkIlpFeasibility(gs.ilp_sol)));
//   // fmt::print("CLP Obj:{:.10Lf} SolStatus:{} LPStatus:{} ILPStatus:{}\n", cs.getScore(cs.lp_sol), solMessage(cs.lp_status), feasMessage(cs.checkLpFeasibility(cs.lp_sol)), feasMessage(cs.checkIlpFeasibility(cs.lp_sol)));
//   // fmt::print("CILP Obj:{:.10Lf} SolStatus:{} LPStatus:{} ILPStatus:{}\n", cs.getScore(cs.ilp_sol), solMessage(cs.ilp_status), feasMessage(cs.checkLpFeasibility(cs.ilp_sol)), feasMessage(cs.checkIlpFeasibility(cs.ilp_sol)));
// }

// void selfCompare(){
//   DetProb prob; prob.uniformGenerate(10000, 50, 1, false, false);
//   DualReducer dr = DualReducer(kPCore, &prob);
//   VectorXd lp_sol = dr.getLpSol();
//   Checker ch = Checker(prob);
//   showHistogram(lp_sol, 10, 0, 0);
//   fmt::print("{:.10Lf} {:.10Lf}\n", dr.best_score, dr.getLpScore());
//   fmt::print("{:.4Lf} {:.4Lf}\n", dr.exe_solve, dr.getLpTime());
//   fmt::print("SolStatus:{} LPStatus:{} ILPStatus:{}\n", solMessage(dr.status), feasMessage(ch.checkLpFeasibility(dr.best_sol)), feasMessage(ch.checkIlpFeasibility(dr.best_sol)));
//   fmt::print("LPStatus:{} ILPStatus:{}\n", feasMessage(ch.checkLpFeasibility(lp_sol)), feasMessage(ch.checkIlpFeasibility(lp_sol)));
// }

void testRunningTime(){
  fstream out_file; out_file.open(kProjectHome + separator() + "plots" + separator() + "dual_reducer.csv", std::ios::out);
  vector<int> expected_ns = {50, 50, 1000, 1000};
  vector<double> inner_probs;
  vector<double> inner_probs_ssds = {0.8, 0.5, 0.8, 0.5};
  vector<double> inner_probs_tpch = {0.1, 0.01, 0.1, 0.01};
  vector<int> Ns = {2000, 20000, 200000, 2000000, 20000000};
  vector<int> seeds = {1, 2, 3, 4};
  for (int N : Ns){
    for (int q = 0 ; q < 8; q ++){
      int index = q % 4;
      for (int seed : seeds){
        DetProb prob;
        string table_name;
        if (q < 4){
          inner_probs = inner_probs_ssds;
          vector<string> cols = {"tmass_prox", "j", "h", "k"};
          table_name = "ssds";
          prob.tableGenerate(table_name, cols, false, N, expected_ns[index], inner_probs[index], seed);
        } else{
          inner_probs = inner_probs_tpch;
          vector<string> cols = {"price", "quantity", "discount", "tax"};
          table_name = "tpch";
          prob.tableGenerate(table_name, cols, true, N, expected_ns[index], inner_probs[index], seed);
        }
        DualReducer dr1 = DualReducer(1, prob);
        DualReducer dr4 = DualReducer(4, prob);
        DualReducer dr80 = DualReducer(80, prob);
        string s;
        if (N <= 2000000){
          GurobiSolver gs = GurobiSolver(prob);
          gs.solveIlp();
          if (gs.lp_status != Found){
            gs.lp_score = 0.0;
          }
          s = fmt::format("{},{},{},{},{},{:.2Lf},{:.2Lf},{:.2Lf},{:.2Lf},{:.8Lf},{:.8Lf},{:.8Lf},{:.8Lf},{:.8Lf}\n", table_name, N, expected_ns[index], inner_probs[index], seed, dr1.exe_ilp, dr4.exe_ilp, dr80.exe_ilp, gs.getIlpTime(), dr1.lp_score, dr1.ilp_score, dr4.ilp_score, dr80.ilp_score, gs.ilp_score);
        } else{
          s = fmt::format("{},{},{},{},{},{:.2Lf},{:.2Lf},{:.2Lf},{:.2Lf},{:.8Lf},{:.8Lf},{:.8Lf},{:.8Lf},{:.8Lf}\n", table_name, N, expected_ns[index], inner_probs[index], seed, dr1.exe_ilp, dr4.exe_ilp, dr80.exe_ilp, -1.0, dr1.lp_score, dr1.ilp_score, dr4.ilp_score, dr80.ilp_score, -1.0);
        }
        cout << s;
        out_file << s;
        out_file.flush();
      }
    }
  }
  out_file.close();
}

int main(){
  testRunningTime();
}