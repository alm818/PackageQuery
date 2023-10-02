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
#include "pb/tmp/ilp_lsr.h"

#include "pb/core/checker.h"
#include "pb/core/dual.h"
#include "pb/core/dual_reducer.h"
#include "pb/core/gurobi_solver.h"

#include "pb/util/udeclare.h"
#include "pb/util/udebug.h"

using namespace pb;

bool is_finished = false;
double max_RAM = 0;
const int kRAMSleep = 10; // In ms

void inspectRAM(){
  is_finished = false;
  max_RAM = 0;
  while (!is_finished){
    max_RAM = max(max_RAM, currentRAM());
    std::this_thread::sleep_for(std::chrono::milliseconds(kRAMSleep));
  }
}

/******************************************/
void Z1(){
  string exp_name = "Z1";
  DetExp exp = DetExp(exp_name);
  exp.o = 7;
  int R = 5;
  int cnt = 0;
  for (int q = 0; q < 2; q ++){
    exp.q = q;
    // if (q == 0) continue;
    for (auto g : exp.g2){
      exp.g = g;
      // if (g < 0.01) continue;
      for (auto S : exp.S3){
        // if (cnt <= 4) continue;
        // if (S <= 100000) continue;
        exp.S = S;
        for (int i = 1; i <= R; i ++){
          cout << "SEED " << i << endl;
          cnt ++;
          // if (S <= 1000000) continue;
          // if (cnt <= 2) continue;
          // if (i <= 1) continue;
          // if (S == 100000) continue;
          // if (S == 1000000 && i <= 3) continue;
          string label = fmt::format("{}_LSR_g{}_s{}", exp.datasets[exp.q], formatFloat(g), exp.S);
          exp.seed = i;
          cout << "det_sql generate\n";
          DetSql det_sql = exp.generate();
          cout << "dlv partition\n";
          exp.dlvPartition();
          // exp.write(label + "_ptime", p_exe);
          cout << "lsr_prob generate\n";
          LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
          cout << "det_prob generate\n";
          DetProb det_prob = DetProb(det_sql, -1, exp.seed);
          for (auto H : exp.H8){
            exp.H = H;
            exp.E = exp.Es[exp.q];
            cout << "lsr_bound_generate\n";
            lsr_prob.generateBounds(exp.E, exp.a, exp.H);
            cout << "det_prob_copy_bound\n";
            det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
            cout << "lsr\n";
            LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
            cout << label << " " << solMessage(lsr.status) << endl;
            if (lsr.status == Found){
              cout << "dual\n";
              Dual dual = Dual(exp.C, det_prob);
              // LsrChecker ch = LsrChecker(lsr_prob);
              // cout << feasMessage(ch.checkIlpFeasibility(lsr.ilp_sol)) << " " << feasMessage(ch.checkLpFeasibility(lsr.lp_sol)) << endl;
              // assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
              // assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
              exp.write(label + "_time", H, lsr.exe_ilp);
              exp.write(label + "_ilp", H, lsr.ilp_score);
              exp.write(label + "_ground", H, dual.score);
              exp.write(label + "_igap", H, intGap(lsr.ilp_score, dual.score));
            }
          }
        }
      }
    }
  }
}

/******************************************/
void Z2(){
  string exp_name = "Z2";
  DetExp exp = DetExp(exp_name);
  exp.o = 6;
  int R = 10;
  for (int q = 0; q < 2; q ++){
    exp.q = q;
    // if (q == 0) continue;
    for (int i = 1; i <= R; i ++){
      exp.seed = i;
      // if (i <= 1) continue;
      DetSql det_sql = exp.generate();
      for (auto H : exp.H8){
        cout << "SEED " << i << endl;
        DetProb det_prob = DetProb(det_sql, -1, exp.seed);
        exp.H = H;
        // if (H <= 9) continue;
        exp.E = exp.Es[exp.q];
        // if (H == 1) continue;
        det_prob.generateBounds(exp.E, exp.a, exp.H);
        cout << "To tdr" << endl;
        DualReducer tdr = DualReducer(exp.C, det_prob, true);
        // GurobiSolver dgs = GurobiSolver(det_prob);
        // dgs.solveLp();
        // cout << solMessage(dgs.lp_status) << " " << solMessage(tdr.status) << endl;
        double ground = tdr.lp_score;
        if (tdr.status == Found){
          exp.write(exp.datasets[exp.q] + "_DR_P", exp.H, 1);
          exp.write(exp.datasets[exp.q] + "_DR_ilp_score", H, tdr.ilp_score);
          exp.write(exp.datasets[exp.q] + "_DR_ground", H, ground);
          exp.write(exp.datasets[exp.q] + "_DR_igap", H, intGap(tdr.ilp_score, ground));
          // exp.write("TDR", tdr.ilp_score, ground);
        } else{
          exp.write(exp.datasets[exp.q] + "_DR_N", exp.H, 1);
        }
        cout << "To rdr" << endl;
        RandomDualReducer rdr = RandomDualReducer(exp.C, det_prob, true, tdr.failure_count);
        if (rdr.status == Found){
          exp.write(exp.datasets[exp.q] + "_RDR_P", exp.H, 1);
          exp.write(exp.datasets[exp.q] + "_RDR_ilp_score", H, rdr.ilp_score);
          exp.write(exp.datasets[exp.q] + "_RDR_ground", H, ground);
          exp.write(exp.datasets[exp.q] + "_RDR_igap", H, intGap(rdr.ilp_score, ground));
          // exp.write("RDR", rdr.ilp_score, ground);
        } else{
          exp.write(exp.datasets[exp.q] + "_RDR_N", exp.H, 1);
        }
        // return;
      }
    }
  }
}

/******************************************/
void Z3(){
  string exp_name = "Z3";
  DetExp exp = DetExp(exp_name);
  exp.o = 7;
  int R = 5;
  for (int q = 0; q < 2; q ++){
    exp.q = q;
    exp.E = exp.Es[exp.q];
    if (q == 0) continue;
    for (int i = 1; i <= R; i ++){
      cout << "SEED " << i << endl;
      exp.seed = i;
      DetSql det_sql = exp.generate();
      exp.dlvPartition();
      cout << "Done partitioning!" << endl;
      LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
      // if (i <= 1) continue;
      for (auto H : exp.H8){ 
        // if (H >= 15) continue;
        exp.H = H;
        lsr_prob.generateBounds(exp.E, exp.a, exp.H);
        cout << "To lsr" << endl;
        LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        double ground = lsr.lp_score;
        cout << "To rlsr" << endl;
        RandomLayeredSketchRefine lplsr = RandomLayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        if (lsr.status == Found){
          exp.write(exp.datasets[exp.q] + "_LSR_P", H, 1);
          exp.write(exp.datasets[exp.q] + "_LSR_ilp", H, lsr.ilp_score);
          exp.write(exp.datasets[exp.q] + "_LSR_ground", H, ground);
          exp.write(exp.datasets[exp.q] + "_LSR_igap", H, intGap(lsr.ilp_score, ground));
        }
        if (lplsr.status == Found){
          exp.write(exp.datasets[exp.q] + "_RLSR_P", H, 1);
          exp.write(exp.datasets[exp.q] + "_RLSR_ilp", H, lplsr.ilp_score);
          exp.write(exp.datasets[exp.q] + "_RLSR_ground", H, ground);
          exp.write(exp.datasets[exp.q] + "_RLSR_igap", H, intGap(lplsr.ilp_score, ground));
        }
      }
    }
  }
}

/******************************************/
void Z4(){
  string exp_name = "Z4";
  DetExp exp = DetExp(exp_name);
  exp.o = 7;
  int R = 5;
  for (int q = 0; q < 2; q ++){
    exp.q = q;
    exp.E = exp.Es[exp.q];
    if (q == 0) continue;
    for (int i = 1; i <= R; i ++){
      // if (i <= 1) continue;
      cout << "SEED " << i << endl;
      exp.seed = i;
      DetSql det_sql = exp.generate();
      exp.dlvPartition();
      cout << "Done partitioning!" << endl;
      LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
      // if (i <= 1) continue;
      for (auto H : exp.H8){ 
        // if (H <= 13) continue;
        exp.H = H;
        lsr_prob.generateBounds(exp.E, exp.a, exp.H);
        cout << "To lsr" << endl;
        LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        double ground = lsr.lp_score;
        cout << "To rlsr" << endl;
        LPLayeredSketchRefine rlsr = LPLayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        if (lsr.status == Found){
          exp.write(exp.datasets[exp.q] + "_LSR_P", H, 1);
          exp.write(exp.datasets[exp.q] + "_LSR_ilp", H, lsr.ilp_score);
          exp.write(exp.datasets[exp.q] + "_LSR_ground", H, ground);
          exp.write(exp.datasets[exp.q] + "_LSR_igap", H, intGap(lsr.ilp_score, ground));
        }
        if (rlsr.status == Found){
          exp.write(exp.datasets[exp.q] + "_LPLSR_P", H, 1);
          exp.write(exp.datasets[exp.q] + "_LPLSR_ilp", H, rlsr.ilp_score);
          exp.write(exp.datasets[exp.q] + "_LPLSR_ground", H, ground);
          exp.write(exp.datasets[exp.q] + "_LPLSR_igap", H, intGap(rlsr.ilp_score, ground));
        }
      }
    }
  }
}

/******************************************/
void Z5(){
  string exp_name = "Z5";
  DetExp exp = DetExp(exp_name);
  exp.o = 7;
  int R = 5;
  for (int q = 1; q >= 0; q --){
    exp.q = q;
    exp.E = exp.Es[exp.q];
    for (int i = 1; i <= R; i ++){
      for (int z = 0; z <= 1; z ++){
        cout << "SEED " << i << endl;
        exp.seed = i;
        if (z == 0) exp.is_max_var = true;
        else exp.is_max_var = false;
        DetSql det_sql = exp.generate();
        DetProb det_prob = DetProb(det_sql, -1, exp.seed);
        exp.dlvPartition();
        cout << "Done partitioning!" << endl;
        LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
        // if (i <= 1) continue;
        for (auto H : exp.H8){
          if (H == 15) continue;
          exp.H = H;
          lsr_prob.generateBounds(exp.E, exp.a, exp.H);
          LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
          det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
          Dual dual = Dual(exp.C, det_prob);
          double ground = dual.score;
          string name = "_LSR" + to_string(z);
          if (lsr.status == Found){
            exp.write(exp.datasets[exp.q] + name + "_P", H, 1);
            exp.write(exp.datasets[exp.q] + name + "_ilp", H, lsr.ilp_score);
            exp.write(exp.datasets[exp.q] + name + "_ground", H, ground);
            exp.write(exp.datasets[exp.q] + name + "_igap", H, intGap(lsr.ilp_score, ground));
          }
        }
      }
    }
  }
}

/******************************************/
void A1(){
  string exp_name = "A1";
  DetExp exp = DetExp(exp_name);
  exp.o = 8;
  int R = 5;
  DetSql det_sql = exp.generate();
  DetProb prob = DetProb(det_sql, -1, exp.seed);
  for (auto C : exp.C7){
    exp.C = C;
    for (int i = 1; i <= R; i ++){
      prob.generateBounds(exp.Es[exp.q], exp.a, exp.H);
      DualReducer dr = DualReducer(exp.C, prob);
      exp.write("D", exp.C, dr.exe_lp);
      exp.write("DR", exp.C, dr.exe_ilp);
    }
  }
  for (auto C : exp.C7){
    exp.C = C;
    LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
    double p_exe = exp.dlvPartition(false);
    exp.write("DLV", exp.C, p_exe);
    for (int i = 1; i <= R; i ++){
      lsr_prob.generateBounds(exp.Es[exp.q], exp.a, exp.H);
      LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
      // cout << solMessage(lsr.status) << endl;
      if (lsr.status == Found){
        LsrChecker ch = LsrChecker(lsr_prob);
        assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
        assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
        exp.write("LSR", exp.C, lsr.exe_ilp);
      }
    }
  }
}
// /******************************************/
// void A2(){
//   string exp_name = "A2";
//   DetExp exp = DetExp(exp_name);
//   int R = 3;
//   PgManager pg = PgManager();
//   for (auto o : exp.o6){
//     // if (o > 7) continue;
//     double N = pow(10, o);
//     exp.o = o;
//     exp.partition_name = fmt::format("{}_{}", exp_name, exp.o);
//     string kd_partition_name = "kd_" + exp.partition_name;
//     DetSql det_sql = exp.generate();
//     DetProb det_prob;
//     vector<LsrProb> lsr_probs;
//     double ground = 0;
//     // cout << "WTF1\n";
//     if (o <= 9){
//       LsrProb lsr_prob = LsrProb(det_sql, exp.partition_name, exp.seed);
//       for (int i = 0; i < R; i ++){
//         exp.H = 1;
//         lsr_prob.generateBounds(exp.E, exp.a, exp.H);
//         lsr_probs.push_back(lsr_prob);
//       }
//       if (o <= 7){
//         //sr partition
//         int tau = (int) (0.001 * N);
//         // cout << o << " " << tau << endl;
//         pg.dropTable("[1P]_" + exp.getTableName() + "_" + kd_partition_name); 
//         pg.dropTable("[1G]_" + exp.getTableName() + "_" + kd_partition_name); 
//         KDTree kt;
//         kt.partitionTable(exp.getTableName(), kd_partition_name, exp.getCols(), tau, DBL_MAX);
//       }
//       if (o >= 6 && o <= 9){
//         //lsr partition
//         exp.partition();
//       }
//     }
//     if (o <= 8){
//       det_prob = DetProb(det_sql, -1, exp.seed);
//       for (int i = 0; i < R; i ++){
//         if (o < 5) det_prob.generateBounds(exp.E, exp.a, exp.H);
//         else {
//           LsrProb lsr_prob = lsr_probs[i];
//           det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
//         }
//         if (o <= 7){
//           GurobiSolver gs = GurobiSolver(det_prob);
//           gs.solveLp();
//           exp.write("GD", N, gs.exe_lp);
//           ground = gs.lp_score;
//           exp.write("GD_aux", N, pctError(gs.lp_score, ground));
//         }
//         if (o <= 6 && 0==1){
//           GurobiSolver gs = GurobiSolver(det_prob);
//           gs.solveIlp(0);
//           exp.write("GDR", N, gs.exe_ilp);
//           exp.write("GDR_aux", N, pctError(gs.ilp_score, ground));
//         }
//         DualReducer dr = DualReducer(exp.C, det_prob);
//         if (dr.status == Found){
//           if (o == 8) ground = dr.lp_score;
//           exp.write("D", N, dr.exe_lp);
//           exp.write("D_aux", N, pctError(dr.lp_score, ground));
//           exp.write("DR", N, dr.exe_ilp);
//           exp.write("DR_aux", N, pctError(dr.ilp_score, ground));
//         }
//         if (o <= 7){
//           LsrProb lsr_prob = lsr_probs[i];
//           lsr_prob.partition_name = kd_partition_name;
//           SketchRefine sr = SketchRefine(lsr_prob);
//           map<long long, long long> sol;
//           bool is_success = sr.sketchAndRefine(sol);
//           LsrChecker ch = LsrChecker(lsr_prob);
//           if (is_success) exp.write("SR", N, sr.exec_sr);
//           if (is_success && ch.checkIlpFeasibility(sol)==Feasibility){
//             double score = ch.getScore(sol);
//             exp.write("SR_aux", N, pctError(score, ground));
//           }
//         }
//         if (o >= 6 && o <= 8){
//           LsrProb lsr_prob = lsr_probs[i];
//           lsr_prob.partition_name = exp.partition_name;
//           LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
//           if (lsr.status == Found){
//             LsrChecker ch = LsrChecker(lsr_prob);
//             assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
//             assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
//             exp.write("LSR", N, lsr.exe_ilp);
//             exp.write("LSR_aux", N, pctError(lsr.ilp_score, ground));
//           }
//         }
//       }
//     }
//     if (o == 9){
//       for (int i = 0; i < R; i ++){
//         LsrProb lsr_prob = lsr_probs[i];
//         LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
//         if (lsr.status == Found){
//           LsrChecker ch = LsrChecker(lsr_prob);
//           assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
//           assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
//           ground = lsr.lp_score;
//           exp.write("LSR", N, lsr.exe_ilp);
//           exp.write("LSR_aux", N, pctError(lsr.ilp_score, ground));
//         }
//       }
//     }
//   }
// }
// /******************************************/
void A3(){
  PgManager pg = PgManager();
  string exp_name = "A3";
  DetExp exp = DetExp(exp_name);
  exp.o = 6;
  int R = 20;
  // double N = pow(10, exp.o);
  for (int q = 0; q < 2; q ++){
    // if (q == 0) continue;
    exp.q = q;
    for (int i = 1; i <= R; i ++){
      // if (q == 0 && i == 21) continue;
      // if (q == 1 && i == 9) continue;
      // if (q == 0 && i <= 13) continue;
      // if (q == 0) continue;
      // if (q == 1 and i <= 9) continue;
      // cout << "Seed: " << i << endl;
      exp.seed = i;
      DetSql det_sql = exp.generate();
      exp.kdPartition();
      exp.dlvPartition();
      for (auto H : exp.H8){
        // if (i == 14 && H <= 13) continue;
        exp.H = H;
        exp.E = exp.Es[exp.q];
        string label = fmt::format("{}", exp.datasets[exp.q]);
        LsrProb lsr_prob = LsrProb(det_sql, "", exp.seed);
        lsr_prob.generateBounds(exp.E, exp.a, exp.H);
        DetProb det_prob = DetProb(det_sql, -1, exp.seed);
        det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
        lsr_prob.partition_name = exp.getDlvPartitionName();
        LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        cout << "LSR finished\n";
        lsr_prob.partition_name = exp.getKdPartitionName();
        SketchRefine sr = SketchRefine(lsr_prob);
        map<long long, long long> sr_sol;
        bool is_success = sr.sketchAndRefine(sr_sol);
        cout << "SR finished\n";
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
        cout << "Start dual\n";
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
// void P1(){
//   string exp_name = "P1";
//   DetExp exp = DetExp(exp_name);
//   exp.o = 8;
//   DetSql det_sql = exp.generate();
//   for (double M : exp.M6){
//     exp.M = M;
//     std::thread ram (inspectRAM);
//     exp.partition_name = fmt::format("{}_M{}", exp_name, exp.M);
//     double p_exe = exp.partition(false);
//     is_finished = true;
//     ram.join();
//     exp.write("DLV", exp.M, max_RAM);
//     exp.write("DLV_aux", exp.M, p_exe);
//   }
// }
// /******************************************/
// void P2(){
//   string exp_name = "P2";
//   DetExp exp = DetExp(exp_name);
//   for (int o : exp.o4){
//     exp.o = o;
//     exp.partition_name = fmt::format("{}_o{}", exp_name, exp.o);
//     DynamicLowVariance dlv = DynamicLowVariance(exp.C, exp.g, exp.M, exp.S);
//     string table_name = exp.getTableName();
//     dlv.dropPartition(table_name, exp.partition_name);
//     dlv.partition(table_name, exp.partition_name);
//     double N = pow(10.0, o);
//     exp.write("DLV_read", N, dlv.pro.time(0));
//     exp.write("DLV_write", N, dlv.pro.time(1));
//     exp.write("DLV_index", N, dlv.pro.time(2));
//     exp.write("DLV_dlv", N, dlv.exe - dlv.pro.time(0) - dlv.pro.time(1) - dlv.pro.time(2));
//     exp.write("DLV_aux", N, dlv.exe);
//   }
// }
// /******************************************/
void L1(){
  string exp_name = "L1";
  DetExp exp = DetExp(exp_name);
  int R = 5;
  exp.o = 8;
  for (int q = 0; q < 2; q ++){
    exp.q = q;
    exp.E = exp.Es[exp.q];
    for (int i = 1; i <= R; i ++){
      exp.seed = i;
      exp.dlvPartition();
      DetSql det_sql = exp.generate();
      LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
      for (auto H : exp.H8){
        lsr_prob.generateBounds(exp.E, exp.a, exp.H);
        LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
        if (lsr.status == Found){
          lsr.pro.print();
          exp.write(exp.datasets[exp.q] + "_LSR_read", H, lsr.pro.time(0));
          exp.write(exp.datasets[exp.q] + "_LSR_dual", H, lsr.exe_dual);
          exp.write(exp.datasets[exp.q] + "_LSR_gb", H, lsr.exe_gb);
          exp.write(exp.datasets[exp.q] + "_LSR_misc", H, lsr.exe_ilp-lsr.exe_gb-lsr.exe_dual-lsr.pro.time(0));
          exp.write(exp.datasets[exp.q] + "_LSR_total", H, lsr.exe_ilp);
        }
      }
    }
  }
}
// /******************************************/
void L2(){
  string exp_name = "L2";
  DetExp exp = DetExp(exp_name);
  exp.o = 7;
  int R = 3;
  for (int q = 0; q < 2; q ++){
    // if (q == 1) continue;
    exp.q = q;
    string label = "LSR_" + exp.datasets[q];
    for (int i = 0; i < R; i ++){
      exp.seed = i+1;
      exp.dlvPartition();
      for (auto F : exp.F5){
        // if (F != 0.7) continue;
        exp.F = F;
        DetSql det_sql = exp.generate();
        det_sql.addFilterWithRatio(exp.filtered_cols[exp.q], exp.F, Bounded);
        // cout << "OK-3\n";
        LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
        // cout << "OK-2\n";
        lsr_prob.generateBounds(exp.E, exp.a, exp.H);
        LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob);
        DetProb det_prob = DetProb(det_sql, -1, exp.seed);
        det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
        DualReducer dr = DualReducer(exp.C, det_prob);
        // cout << label << " " << solMessage(lsr.status) << endl;
        if (lsr.status == Found){
          LsrChecker ch = LsrChecker(lsr_prob);
          assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
          assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
          double ground = dr.lp_score;
          exp.write(label+"_lp", exp.F, pctError(lsr.lp_score, ground));
          exp.write(label+"_ilp", exp.F, pctError(lsr.ilp_score, ground));
          exp.write(label+"_aux", exp.F, lsr.exe_ilp);
        }
      }
    }
  }
}
// /******************************************/
void A4(){
  string exp_name = "A4";
  DetExp exp = DetExp(exp_name);
  PgManager pg = PgManager();
  for (int q = 0; q < 2; q ++){
    if (q == 0) continue;
    // if (q == 1) continue;
    exp.q = q;
    exp.E = exp.Es[exp.q];
    for (auto H : exp.H4){
      exp.H = H;
      for (auto o : exp.o6){
        if (o == 9 && q == 1) continue;
        int R = 5;
        if (o == 9) R = 3;
        // if (q == 0 && H <= 1) continue;
        // if (q == 0 && H <= 3 && o <= 8) continue;
        // if (q == 0 && H == 7 && o <= 6) continue;
        exp.o = o;
        int seed_count = 0;
        exp.seed = 1;
        double N = pow(10, exp.o);
        string label = fmt::format("{}_H{}", exp.datasets[exp.q], exp.H);
        while (seed_count < R){
          DetSql det_sql = exp.generate();

          if (o <= 8){
            cout << "Start kd " + label << endl;
            exp.kdPartition();
          }

          ////////
          cout << "Start dlv "  + label << endl;
          exp.dlvPartition();
          ////////

          cout << "lsr_prob\n";
          LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
          lsr_prob.generateBounds(exp.E, exp.a, exp.H);
          DetProb det_prob;

          if (o <= 8){
            cout << "det_prob\n";
            det_prob = DetProb(det_sql, -1, exp.seed);
            det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
          }

          cout << "Start lsr " + label << endl;
          lsr_prob.partition_name = exp.getDlvPartitionName();
          LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
          if (lsr.status == Found){
            seed_count ++;
            LsrChecker ch = LsrChecker(lsr_prob);
            assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
            assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);

            double ground = 0;
            if (o <= 8){
              cout << "dual\n";
              Dual dual = Dual(exp.C, det_prob);
              ground = dual.score;
            } else ground = lsr.lp_score;

            if (o <= 6){
              cout << "Start gs " + label << endl;
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
              cout << "Start sr " + label << endl;
              lsr_prob.partition_name = exp.getKdPartitionName();
              SketchRefine sr = SketchRefine(lsr_prob);
              map<long long, long long> sr_sol;
              bool is_success = sr.sketchAndRefine(sr_sol);
              LsrChecker ch = LsrChecker(lsr_prob);
              cout << "Sr finished " + label << endl; 
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
/******************************************/
void A5(){
  string exp_name = "A5";
  DetExp exp = DetExp(exp_name);
  int R = 5;
  for (int o = 7; o <= 8; o ++){
    exp.o = o;
    for (int i = 1; i <= R; i ++){
      exp.seed = i;
      DetSql det_sql = exp.generate();
      DetProb prob = DetProb(det_sql, -1, exp.seed);
      prob.generateBounds(exp.Es[exp.q], exp.a, exp.H);
      Dual dual = Dual(exp.C, prob);
      exp.write("D", exp.o, dual.exe_solve);
      GurobiSolver gs = GurobiSolver(prob);
      gs.solveLp();
      exp.write("GD", exp.o, gs.exe_lp);
    }
  }
}
/******************************************/
void test(){
  PgManager pg = PgManager();
  // DynamicLowVariance dlv = DynamicLowVariance();
  // auto ps = dlv.getPartitions();
  // for (auto pi : ps){
  //   string table_name = pi.first;
  //   auto partition_names = pi.second;
  //   if (table_name.find("ssds") != string::npos){
  //     for (auto n : partition_names){
  //       cout << table_name << " " << n << endl;
  //     }
  //   }
  //   for (auto partition_name : partition_names){
  //     dlv.dropPartition(table_name, partition_name);
  //   }
  //   pg.dropTable(table_name); 
  // }

  // PgManager pg = PgManager();
  // auto tables = pg.listTables();
  // for (string table : tables){
  //   if (getLayerIndex(table) > 0){
  //     pg.dropTable(table);
  //   }
  // }

  // PgManager pg = PgManager();
  auto tables = pg.listTables();
  for (string table : tables){
    if (table.find("ssds_") != string::npos) pg.dropTable(table);
    // if (table.find("ssds") != string::npos && table.find("kd") != string::npos) pg.dropTable(table);
  }

  // Inverse hardness test
  // string exp_name = "A4";
  // DetExp exp = DetExp(exp_name);
  // exp.o = 7;
  // DetSql det_sql = exp.generate();
  // DetProb det_prob = DetProb(det_sql, -1, exp.seed);
  // det_prob.generateBounds(exp.E, 0, 11.2);
  // double minE = 0;
  // double h = det_prob.det_bound.minHardness(minE, det_prob.bl, det_prob.bu);
  // cout << h << endl;
}
/******************************************/

int main() {
  // Z1();
  // Z2();
  // Z3();
  // Z4();
  Z5();
  // A1();
  // A2();
  // A3();
  // A4();
  // A5();
  // P1();
  // P2();
  // L1();
  // L2();
  // test();
  return 0;
}