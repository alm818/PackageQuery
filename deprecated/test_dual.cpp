#include "pb/util/udeclare.h"
#include "pb/util/udebug.h"
#include "pb/util/uconfig.h"
#include "pb/util/umisc.h"
#include "pb/core/dual.h"
#include "pb/core/gurobi_solver.h"

using namespace pb;

void testAcc(){
  int N = 100000;
  int E = 50;
  double var = 10;
  double outlier_prob = 2.0;
  int T = 1000;
  for (int i = 0; i < T; i ++){
    DetProb prob; prob.uniformGenerate(N, E, var, outlier_prob);
    Dual dual = Dual(kPCore, prob);
    GurobiSolver gs = GurobiSolver(prob);
    gs.solveLp();
    fmt::print("{:.10Lf} {:.10Lf}\n", dual.score, gs.lp_score);
    assert(gs.lp_score - dual.score < 1e-6);
  }
}

void testRunningTime(){
  fstream out_file; out_file.open(kProjectHome + separator() + "plots" + separator() + "dual.csv", std::ios::out);
  vector<int> expected_ns = {50, 50, 1000, 1000};
  vector<double> inner_probs_ssds = {0.8, 0.5, 0.8, 0.5};
  vector<double> inner_probs_tpch = {0.1, 0.01, 0.1, 0.01};
  vector<int> Ns = {2000, 20000, 200000, 2000000, 20000000};
  // vector<int> Ns = {100000000};
  vector<int> seeds = {1, 2, 3, 4};
  for (int N : Ns){
    for (int q = 0 ; q < 8; q ++){
      int index = q % 4;
      for (int seed : seeds){
        DetProb prob;
        string table_name;
        if (q < 4){
          vector<string> cols = {"tmass_prox", "j", "h", "k"};
          table_name = "ssds";
          prob.tableGenerate(table_name, cols, false, N, expected_ns[index], inner_probs_ssds[index], seed);
        } else{
          vector<string> cols = {"price", "quantity", "discount", "tax"};
          table_name = "tpch";
          prob.tableGenerate(table_name, cols, true, N, expected_ns[index], inner_probs_tpch[index], seed);
        }
        Dual d1 = Dual(1, prob);
        Dual d4 = Dual(4, prob);
        Dual d80 = Dual(80, prob);
        string s;
        if (N <= 20000000){
          GurobiSolver gs = GurobiSolver(prob);
          gs.solveLp();
          if (gs.lp_status != Found){
            gs.lp_score = 0.0;
          }
          s = fmt::format("{},{},{},{},{},{:.2Lf},{:.2Lf},{:.2Lf},{:.2Lf},{:.8Lf},{:.8Lf},{:.8Lf},{:.8Lf}\n", table_name, N, expected_ns[index], inner_probs[index], seed, d1.exe_solve, d4.exe_solve, d80.exe_solve, gs.getLpTime(), d1.score, d4.score, d80.score, gs.lp_score);
        } else{
          s = fmt::format("{},{},{},{},{},{:.2Lf},{:.2Lf},{:.2Lf},{:.2Lf},{:.8Lf},{:.8Lf},{:.8Lf},{:.8Lf}\n", table_name, N, expected_ns[index], inner_probs[index], seed, d1.exe_solve, d4.exe_solve, d80.exe_solve, -1.0, d1.score, d4.score, d80.score, -1.0);
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
  //testAcc();
  testRunningTime();
}