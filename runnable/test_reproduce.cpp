#include "pb/det/det_sql.h"
#include "pb/det/det_exp.h"
#include "pb/det/det_prob.h"
#include "pb/det/lsr_prob.h"
#include "pb/det/dlv.h"
#include "pb/det/lsr.h"
#include "pb/det/sr.h"
#include "pb/det/kd_tree.h"

#include "pb/tmp/random_dual_reducer.h"
#include "pb/tmp/random_lsr.h"
#include "pb/tmp/test_dual_reducer.h"
#include "pb/tmp/lp_lsr.h"

#include "pb/core/checker.h"
#include "pb/core/dual.h"
#include "pb/core/dual_reducer.h"
#include "pb/core/gurobi_solver.h"

#include "pb/util/udeclare.h"
#include "pb/util/udebug.h"

using namespace pb;

// /******************************************/
void A3(){
  PgManager pg = PgManager();
  string exp_name = "A3";
  DetExp exp = DetExp(exp_name);
  exp.o = 6;
  int R = 20;
  for (int q = 0; q < 2; q ++){
    exp.q = q;
    for (int i = 1; i <= R; i ++){
      exp.seed = i;
      DetSql det_sql = exp.generate();
      exp.kdPartition();
      exp.dlvPartition();
      for (auto H : exp.H8){
        exp.H = H;
        exp.E = exp.Es[exp.q];
        string label = fmt::format("{}", exp.datasets[exp.q]);
        LsrProb lsr_prob = LsrProb(det_sql, "", exp.seed);
        lsr_prob.generateBounds(exp.E, exp.a, exp.H);
        DetProb det_prob = DetProb(det_sql, -1, exp.seed);
        det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
        lsr_prob.partition_name = exp.getDlvPartitionName();
        LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        // cout << "LSR finished\n";
        lsr_prob.partition_name = exp.getKdPartitionName();
        SketchRefine sr = SketchRefine(lsr_prob);
        map<long long, long long> sr_sol;
        bool is_success = sr.sketchAndRefine(sr_sol);
        // cout << "SR finished\n";
        LsrChecker ch = LsrChecker(lsr_prob);
        Checker det_ch = Checker(det_prob);
        double ground = 0;
        // GurobiSolver gs_lp = GurobiSolver(det_prob);
        // gs_lp.solveLp();
        // cout << "GS LP finished\n";
        // ground = gs_lp.lp_score;
        // It is mathematically proved that dr.status or lsr.status == Found means Checker is correct
        
        // if (lsr.status == Found){
        //   is_positive = true;
        // } else{
        //   cout << "Start gs decide has solution or not\n"; 
        //   GurobiSolver gs = GurobiSolver(det_prob);
        //   is_positive = gs.hasIlpSolution();
        //   cout << "gs finished\n";
        // }
        // cout << "Start dual\n";
        Dual dual = Dual(exp.C, det_prob);
        ground = dual.score;
        // gs.solveIlp();
        // if (gs.ilp_status == Found){
        //   is_positive = true;
        //   exp.write("GR_err_" + label, exp.H, pctError(gs.ilp_score, ground));
        // } else {
        //   is_positive = false;
        // }

        if (is_success && ch.checkIlpFeasibility(sr_sol)==Feasibility){
          double sr_score = ch.getScore(sr_sol);
          exp.write(label + "_SR_P", exp.H, 1.0);
          exp.write(label + "_SR_ilp", exp.H, sr_score);
          exp.write(label + "_SR_ground", exp.H, ground);
          exp.write(label + "_SR_igap", exp.H, intGap(sr_score, ground));
        }
        
        if (lsr.status == Found){
          exp.write(label + "_LSR_P", exp.H, 1.0);
          exp.write(label + "_LSR_ilp", exp.H, lsr.ilp_score);
          exp.write(label + "_LSR_ground", exp.H, ground);
          exp.write(label + "_LSR_igap", exp.H, intGap(lsr.ilp_score, ground));
          exp.write(label + "_GDR_P", exp.H, 1.0);
        } else {
          // cout << "Start gs\n"; 
          // GurobiSolver gs = GurobiSolver(det_prob);
          // if (gs.hasIlpSolution()){
          //   exp.write(label + "_GDR_P", exp.H, 1.0);
          // }
        }
      }
    }
  }
}

// /******************************************/
void A4(){
  for (int q = 0; q < 2; q ++){
    string exp_name = "A4_" + DetExp::datasets[q];
    DetExp exp = DetExp(exp_name);
    PgManager pg = PgManager();
    exp.q = q;
    exp.E = exp.Es[exp.q];
    for (auto H : exp.H4){
      exp.H = H;
      for (auto o : exp.o6){
        int R = 5;
        exp.o = o;
        int seed_count = 0;
        exp.seed = 1;
        double N = pow(10, exp.o);
        string label = fmt::format("{}_H{}", exp.datasets[exp.q], exp.H);
        while (seed_count < R){
          DetSql det_sql = exp.generate();

          if (o <= 8){
            exp.kdPartition();
          }

          ////////
          exp.dlvPartition();
          ////////

          LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
          lsr_prob.generateBounds(exp.E, exp.a, exp.H);
          DetProb det_prob;

          if (o <= 8){
            det_prob = DetProb(det_sql, -1, exp.seed);
            det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
          }

          lsr_prob.partition_name = exp.getDlvPartitionName();
          LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
          if (lsr.status == Found){
            seed_count ++;
            LsrChecker ch = LsrChecker(lsr_prob);
            assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
            assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);

            double ground = 0;
            if (o <= 8){
              Dual dual = Dual(exp.C, det_prob);
              ground = dual.score;
            } else ground = lsr.lp_score;

            if (o <= 6){
              GurobiSolver gs = GurobiSolver(det_prob);
              gs.solveIlp(1e-4, kTimeLimit);
              if (gs.ilp_status == Found){
                exp.write(label + "_GDR_time" , N, gs.exe_ilp);
                exp.write(label + "_GDR_ilp", N, gs.ilp_score);
                exp.write(label + "_GDR_ground", N, ground);
                exp.write(label + "_GDR_igap", N, intGap(gs.ilp_score, ground));
              }
            }

            exp.write(label + "_LSR_time", N, lsr.exe_ilp);
            exp.write(label + "_LSR_ilp", N, lsr.ilp_score);
            exp.write(label + "_LSR_ground", N, ground);
            exp.write(label + "_LSR_igap", N, intGap(lsr.ilp_score, ground));

            if (o <= 8){
              lsr_prob.partition_name = exp.getKdPartitionName();
              SketchRefine sr = SketchRefine(lsr_prob);
              map<long long, long long> sr_sol;
              bool is_success = sr.sketchAndRefine(sr_sol);
              LsrChecker ch = LsrChecker(lsr_prob);
              if (is_success && ch.checkIlpFeasibility(sr_sol)==Feasibility){
                double sr_score = ch.getScore(sr_sol);
                exp.write(label + "_SR_time", N, sr.exec_sr);
                exp.write(label + "_SR_ilp", N, sr_score);
                exp.write(label + "_SR_ground", N, ground);
                exp.write(label + "_SR_igap", N, intGap(sr_score, ground));
              }
            }

          }

          exp.seed ++;
        }
      }
    }
  }
}

int main() {
  A4();
  A3();
  return 0;
}