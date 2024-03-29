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
#include "pb/tmp/gdr_lsr.h"
#include "pb/tmp/tmp_sr.h"
#include "pb/tmp/pop_lsr.h"

#include "pb/core/checker.h"
#include "pb/core/dual.h"
#include "pb/core/dual_reducer.h"
#include "pb/core/gurobi_solver.h"

#include "pb/util/udeclare.h"
#include "pb/util/udebug.h"

using namespace pb;

// /******************************************/
// Grid-search downscale factor d_f and augmenting size \alpha
void G1(){ 
  string exp_name = "G1";
  DetExp exp = DetExp(exp_name);
  exp.o = 7;
  int R = 5;
  int cnt = 0;
  for (int q = 0; q <= 3; q ++){
    exp.q = q;
    // if (q == 0) continue;
    for (auto g : exp.g3){
      // if (g == 0.001) continue;
      exp.g = g;
      for (auto S : exp.S3){
        // if (S == 10000) continue;
        exp.S = S;
        for (int i = 1; i <= R; i ++){
          cnt ++;
          string label = fmt::format("{}{}_df{}_alpha{}", exp.datasets[exp.q], q/2, formatFloat(g), exp.S);
          exp.seed = i;
          DetSql det_sql = exp.generate();
          double p_exe = exp.dlvPartition();
          exp.write(label + "_ptime", p_exe);
          LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
          DetProb det_prob = DetProb(det_sql, -1, exp.seed);
          for (auto H : exp.H4){
            exp.H = H;
            exp.E = exp.Es[exp.q];
            lsr_prob.generateBounds(exp.E, exp.a, exp.H);
            det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
            LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
            // lsr.pro.print();
            if (lsr.status == Found){
              Dual dual = Dual(exp.C, det_prob);
              LsrChecker ch = LsrChecker(lsr_prob);
              assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
              assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
              exp.write(label + "_found", H, 1);
              exp.write(label + "_time", H, lsr.exe_ilp);
              exp.write(label + "_ilp", H, lsr.ilp_score);
              exp.write(label + "_lp", H, dual.score);
              exp.write(label + "_igap", H, intGap(lsr.ilp_score, dual.score));
            } else{
              exp.write(label + "_nofound", H, 1);
            }
          }
        }
      }
    }
  }
}

// /******************************************/
// The effect of q on the performance of Dual Reducer
void G2(){
  string exp_name = "G2";
  DetExp exp = DetExp(exp_name);
  exp.o = 5;
  int R = 10;
  for (int q = 0; q <= 3; q ++){
    exp.q = q;
    for (int i = 1; i <= R; i ++){
      exp.seed = i;
      DetSql det_sql = exp.generate();
      DetProb det_prob = DetProb(det_sql, -1, exp.seed);
      for (auto H : exp.H7){
        exp.H = H;
        exp.E = exp.Es[exp.q];
        det_prob.generateBounds(exp.E, exp.a, exp.H);
        for (auto Q : exp.Q4){
          string label = fmt::format("{}{}_Q{}", exp.datasets[exp.q], q/2, Q);
          DualReducer dr = DualReducer(exp.C, det_prob, true, 1e-4, kTimeLimit, Q);
          if (dr.status == Found){
            Checker ch = Checker(det_prob);
            assert(ch.checkLpFeasibility(dr.lp_sol) == Feasibility);
            assert(ch.checkIlpFeasibility(dr.ilp_sol) == Feasibility);
            exp.write(label + "_found", H, 1);
            exp.write(label + "_time", H, dr.exe_ilp);
            exp.write(label + "_ilp", H, dr.ilp_score);
            exp.write(label + "_lp", H, dr.lp_score);
            exp.write(label + "_igap", H, intGap(dr.ilp_score, dr.lp_score));
          } else{
            exp.write(label + "_nofound", H, 1);
          }
        }
      }
    }
  }
}

// /******************************************/
// Mini-Experiment 1. Comparing Progressive Shading with LP solution and with ILP solution
void M1(){
  string exp_name = "M1";
  DetExp exp = DetExp(exp_name);
  exp.o = 7;
  int R = 5;
  for (int q = 0; q <= 3; q ++){
    exp.q = q;
    exp.E = exp.Es[exp.q];
    string label = fmt::format("{}{}", exp.datasets[exp.q], q/2);
    for (int i = 1; i <= R; i ++){
      exp.seed = i;
      DetSql det_sql = exp.generate();
      exp.dlvPartition();
      LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
      for (auto H : exp.H7){ 
        exp.H = H;
        lsr_prob.generateBounds(exp.E, exp.a, exp.H);
        LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        double ground = lsr.lp_score;
        ILPLayeredSketchRefine rlsr = ILPLayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        if (lsr.status == Found){
          exp.write(label + "_LPLSR_found", H, 1);
          exp.write(label + "_LPLSR_time", H, lsr.exe_ilp);
          exp.write(label + "_LPLSR_ilp", H, lsr.ilp_score);
          exp.write(label + "_LPLSR_lp", H, ground);
          exp.write(label + "_LPLSR_igap", H, intGap(lsr.ilp_score, ground));
        } else{
          exp.write(label + "_LPLSR_nofound", H, 1);
        }
        if (rlsr.status == Found){
          exp.write(label + "_ILPLSR_found", H, 1);
          exp.write(label + "_ILPLSR_time", H, rlsr.exe_ilp);
          exp.write(label + "_ILPLSR_ilp", H, rlsr.ilp_score);
          exp.write(label + "_ILPLSR_lp", H, ground);
          exp.write(label + "_ILPLSR_igap", H, intGap(rlsr.ilp_score, ground));
        } else{
          exp.write(label + "_ILPLSR_nofound", H, 1);
        }
      }
    }
  }
}

// /******************************************/
// Mini-Experiment 2. Comparing Progressive Shading with Neighbor Sampling and with Random Sampling
void M2(){
  string exp_name = "M2";
  DetExp exp = DetExp(exp_name);
  exp.o = 7;
  int R = 5;
  for (int q = 0; q <= 3; q ++){
    exp.q = q;
    exp.E = exp.Es[exp.q];
    string label = fmt::format("{}{}", exp.datasets[exp.q], q/2);
    for (int i = 1; i <= R; i ++){
      exp.seed = i;
      DetSql det_sql = exp.generate();
      exp.dlvPartition();
      LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
      for (auto H : exp.H7){ 
        exp.H = H;
        lsr_prob.generateBounds(exp.E, exp.a, exp.H);
        // cout << "OK 1\n"; 
        LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        double ground = lsr.lp_score;
        // cout << "OK 2\n"; 
        RandomLayeredSketchRefine lplsr = RandomLayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        if (lsr.status == Found){
          exp.write(label + "_LSR_found", H, 1);
          exp.write(label + "_LSR_time", H, lsr.exe_ilp);
          exp.write(label + "_LSR_ilp", H, lsr.ilp_score);
          exp.write(label + "_LSR_lp", H, ground);
          exp.write(label + "_LSR_igap", H, intGap(lsr.ilp_score, ground));
        } else{
          exp.write(label + "_LSR_nofound", H, 1);
        }
        if (lplsr.status == Found){
          exp.write(label + "_RLSR_found", H, 1);
          exp.write(label + "_RLSR_time", H, lplsr.exe_ilp);
          exp.write(label + "_RLSR_ilp", H, lplsr.ilp_score);
          exp.write(label + "_RLSR_lp", H, ground);
          exp.write(label + "_RLSR_igap", H, intGap(lplsr.ilp_score, ground));
        } else{
          exp.write(label + "_RLSR_nofound", H, 1);
        }
      }
    }
  }
}

// /******************************************/
// Mini-Experiment 3. Parallel Dual Simplex performance as the number of cores increases
void M3(){
  string exp_name = "M3";
  DetExp exp = DetExp(exp_name);
  int R = 1;
  for (int o = 7; o <= 8; o ++){
    exp.o = o;
    string label = fmt::format("D_o{}", exp.o);
    for (int i = 1; i <= R; i ++){
      exp.seed = i;
      DetSql det_sql = exp.generate();
      DetProb prob = DetProb(det_sql, -1, exp.seed);
      prob.generateBounds(exp.E, exp.a, exp.H);
      for (auto C : exp.C7){
        exp.C = C;
        Dual dual = Dual(exp.C, prob);
        exp.write(label, exp.C, dual.exe_solve);
      }
    }
  }
}

// /******************************************/
// Mini-Experiment 4. Comparing Dual Reducer with Auxiliary LP and with random sampling
void M4(){
  string exp_name = "M4";
  DetExp exp = DetExp(exp_name);
  exp.o = 6;
  int R = 5;
  for (int q = 0; q <= 3; q ++){
    exp.q = q;
    string label = fmt::format("{}{}", exp.datasets[exp.q], q/2);
    for (int i = 1; i <= R; i ++){
      exp.seed = i;
      DetSql det_sql = exp.generate();
      for (auto H : exp.H7){
        DetProb det_prob = DetProb(det_sql, -1, exp.seed);
        exp.H = H;
        exp.E = exp.Es[exp.q];
        det_prob.generateBounds(exp.E, exp.a, exp.H);
        DualReducer tdr = DualReducer(exp.C, det_prob, true);
        double ground = tdr.lp_score;
        if (tdr.status == Found){
          exp.write(label + "_DR_found", exp.H, 1);
          exp.write(label + "_DR_time", H, tdr.exe_ilp);
          exp.write(label + "_DR_ilp", H, tdr.ilp_score);
          exp.write(label + "_DR_lp", H, ground);
          exp.write(label + "_DR_igap", H, intGap(tdr.ilp_score, ground));
        } else{
          exp.write(label + "_DR_nofound", exp.H, 1);
        }
        RandomDualReducer rdr = RandomDualReducer(exp.C, det_prob, true, tdr.failure_count);
        if (rdr.status == Found){
          exp.write(label + "_RDR_found", exp.H, 1);
          exp.write(label + "_RDR_time", H, rdr.exe_ilp);
          exp.write(label + "_RDR_ilp", H, rdr.ilp_score);
          exp.write(label + "_RDR_lp", H, ground);
          exp.write(label + "_RDR_igap", H, intGap(rdr.ilp_score, ground));
        } else{
          exp.write(label + "_RDR_nofound", exp.H, 1);
        }
      }
    }
  }
}

// /******************************************/
// Mini-Experiment 5. Comparing DLV and KD-Tree
void M5(){
  string exp_name = "M5";
  DetExp exp = DetExp(exp_name);
  exp.q = 0;
  for (int o = 8; o <= 9; o ++){
    exp.o = o;
    DetSql det_sql = exp.generate();
    // Kd cannot run at 10^9 due to memory issue
    if (o == 8){
      double kd_exe = exp.kdPartition(false);
      exp.write("KD", o, kd_exe);
    }
    double dlv_exe = exp.dlvPartition(false);
    exp.write("DLV", o, dlv_exe);
  }
}

// /******************************************/
// Mini-Experiment 6. Comparing Progressive Shading with Dual Reducer and with Gurobi
void M6(){
  string exp_name = "M6";
  DetExp exp = DetExp(exp_name);
  exp.o = 8;
  int R = 5;
  for (int q = 0; q <= 3; q ++){
    exp.q = q;
    exp.E = exp.Es[exp.q];
    string label = fmt::format("{}{}", exp.datasets[exp.q], q/2);
    for (int i = 1; i <= R; i ++){
      exp.seed = i;
      DetSql det_sql = exp.generate();
      exp.dlvPartition();
      LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
      for (auto H : exp.H7){ 
        exp.H = H;
        lsr_prob.generateBounds(exp.E, exp.a, exp.H);
        LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        double ground = lsr.lp_score;
        // RandomLayeredSketchRefine lplsr = RandomLayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        GDRLayeredSketchRefine gdrlsr = GDRLayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        if (lsr.status == Found){
          exp.write(label + "_LSR_found", H, 1);
          exp.write(label + "_LSR_time", H, lsr.exe_ilp);
          exp.write(label + "_LSR_ilp", H, lsr.ilp_score);
          exp.write(label + "_LSR_lp", H, ground);
          exp.write(label + "_LSR_igap", H, intGap(lsr.ilp_score, ground));
        } else{
          exp.write(label + "_LSR_nofound", H, 1);
        }
        if (gdrlsr.status == Found){
          exp.write(label + "_GDRLSR_found", H, 1);
          exp.write(label + "_GDRLSR_time", H, gdrlsr.exe_ilp);
          exp.write(label + "_GDRLSR_ilp", H, gdrlsr.ilp_score);
          exp.write(label + "_GDRLSR_lp", H, ground);
          exp.write(label + "_GDRLSR_igap", H, intGap(gdrlsr.ilp_score, ground));
        } else{
          exp.write(label + "_GDRLSR_nofound", H, 1);
        }
      }
    }
  }
}

// /******************************************/
// Mini-Experiment 7. Progressive Shading performance as the number of cores increases
void M7(){
  string exp_name = "M7";
  DetExp exp = DetExp(exp_name);
  int R = 1;
  exp.q = 0;
  for (int o = 8; o <= 9; o ++){
    exp.o = o;
    string label = fmt::format("o{}", exp.o);
    for (int i = 1; i <= R; i ++){
      exp.seed = i;
      DetSql det_sql = exp.generate();
      cout << "dlv:" << exp.dlvPartition() << "ms\n";
      LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
      lsr_prob.generateBounds(exp.Es[exp.q], exp.a, exp.H);
      for (auto C : exp.C7){
        exp.C = C;
        LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        exp.write(label + "_time", exp.C, lsr.exe_ilp);
      }
    }
  }
}

// /******************************************/
// Mini-Experiment 8. Progressive Shading performance using random partition
void M8(){
  string exp_name = "M8";
  DetExp exp = DetExp(exp_name);
  exp.o = 7;
  int R = 5;
  for (int q = 0; q <= 3; q ++){
    exp.q = q;
    exp.E = exp.Es[exp.q];
    for (int i = 1; i <= R; i ++){
      exp.seed = i;
      DetSql det_sql = exp.generate();
      exp.dlvPartition();
      DetProb det_prob = DetProb(det_sql, -1, exp.seed);
      for (auto H : exp.H7){
        exp.H = H;
        string label = fmt::format("{}{}", exp.datasets[exp.q], q/2);
        LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
        lsr_prob.generateBounds(exp.E, exp.a, exp.H);
        det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);

        lsr_prob.partition_name = exp.getDlvPartitionName();
        LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
      
        POPLayeredSketchRefine poplsr = POPLayeredSketchRefine(det_prob);
        double ground = 0;
        if (lsr.status == Found){
          ground = lsr.lp_score;
          exp.write(label + "_LSR_found", H, 1);
          exp.write(label + "_LSR_time", H, lsr.exe_ilp);
          exp.write(label + "_LSR_ilp", H, lsr.ilp_score);
          exp.write(label + "_LSR_lp", H, ground);
          exp.write(label + "_LSR_igap", H, intGap(lsr.ilp_score, ground));
        } else{
          exp.write(label + "_LSR_nofound", H, 1);
        }

        if (poplsr.status == Found){
          exp.write(label + "_POPLSR_found", H, 1);
          exp.write(label + "_POPLSR_time", H, poplsr.exe_ilp);
          exp.write(label + "_POPLSR_ilp", H, poplsr.ilp_score);
          exp.write(label + "_POPLSR_lp", H, ground);
          exp.write(label + "_POPLSR_igap", H, intGap(poplsr.ilp_score, ground));
        } else{
          exp.write(label + "_POPLSR_nofound", H, 1);
        }
      }
    }
  }
}

// /******************************************/
// Query performance as relation size increases for Q1,Q2,Q3,Q4 as 4 cores
void E1(){
  string exp_name = "E1";
  DetExp exp = DetExp(exp_name);
  int R = 10;
  for (int q = 3; q >= 0; q --){
    // if (q == 3) continue;
    exp.q = q;
    exp.E = exp.Es[exp.q];
    for (auto o : exp.o6){
      // if (q == 2 && o <= 7) continue;
      if (exp.datasets[exp.q] == "ssds" && o == 9) continue;
      exp.o = o;
      double N = pow(10, exp.o);
      for (int i = 1; i <= R; i ++){
        exp.seed = i;
        DetSql det_sql = exp.generate();
        if (o <= 8){
          cout << "kd:" << exp.kdPartition() << "ms\n";
        }
        cout << "dlv:" << exp.dlvPartition() << "ms\n";
        LsrProb lsr_prob = LsrProb(det_sql, exp.getDlvPartitionName(), exp.seed);
        DetProb det_prob;
        if (o <= 6) det_prob = DetProb(det_sql, -1, exp.seed);

        for (auto H : exp.H4){
          exp.H = H;
          string label = fmt::format("{}{}_H{}", exp.datasets[exp.q], q/2, exp.H);

          lsr_prob.generateBounds(exp.E, exp.a, exp.H);
          if (o <= 6){
            det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
          }

          lsr_prob.partition_name = exp.getDlvPartitionName();
          LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
          cout << "lsr!!\n";
          lsr_prob.partition_name = exp.getDlvPartitionName();
          LayeredSketchRefine lsr4 = LayeredSketchRefine(4, lsr_prob, exp.S, true);
          cout << "lsr4!!\n";
          double ground = -DBL_MAX;
          if (o <= 6){
            Dual dual = Dual(exp.C, det_prob);
            ground = dual.score;
          } else{
            if (lsr.status == Found) ground = lsr.lp_score;
          }
          if (ground == -DBL_MAX) continue;

          LsrChecker ch = LsrChecker(lsr_prob);

          if (lsr.status == Found){
            assert(ch.checkLpFeasibility(lsr.lp_sol) == Feasibility);
            assert(ch.checkIlpFeasibility(lsr.ilp_sol) == Feasibility);
            exp.write(label + "_LSR_found", N, 1);
            exp.write(label + "_LSR_time", N, lsr.exe_ilp);
            exp.write(label + "_LSR_ilp", N, lsr.ilp_score);
            exp.write(label + "_LSR_ground", N, ground);
            exp.write(label + "_LSR_igap", N, intGap(lsr.ilp_score, ground));
          } else {
            exp.write(label + "_LSR_nofound", N, 1);
          }

          if (lsr4.status == Found){
            assert(ch.checkLpFeasibility(lsr4.lp_sol) == Feasibility);
            assert(ch.checkIlpFeasibility(lsr4.ilp_sol) == Feasibility);
            exp.write(label + "_LSR4_found", N, 1);
            exp.write(label + "_LSR4_time", N, lsr4.exe_ilp);
            exp.write(label + "_LSR4_ilp", N, lsr4.ilp_score);
            exp.write(label + "_LSR4_ground", N, ground);
            exp.write(label + "_LSR4_igap", N, intGap(lsr4.ilp_score, ground));
          } else {
            exp.write(label + "_LSR4_nofound", N, 1);
          }

          if (o <= 6){
            GurobiSolver gs = GurobiSolver(det_prob);
            gs.solveIlp(1e-4, kTimeLimit);
            if (gs.ilp_status == Found){
              exp.write(label + "_GDR_found", N, 1);
              exp.write(label + "_GDR_time" , N, gs.exe_ilp);
              exp.write(label + "_GDR_ilp", N, gs.ilp_score);
              exp.write(label + "_GDR_ground", N, ground);
              exp.write(label + "_GDR_igap", N, intGap(gs.ilp_score, ground));
            } else{
              exp.write(label + "_GDR_nofound", N, 1);
            }
          }
          cout << "gs!!\n";

          if (o <= 8){
            lsr_prob.partition_name = exp.getKdPartitionName();
            TmpSketchRefine sr = TmpSketchRefine(lsr_prob, 80);
            map<long long, long long> sr_sol;
            bool is_success = sr.sketchAndRefine(sr_sol);
            if (is_success && ch.checkIlpFeasibility(sr_sol)==Feasibility){
              double sr_score = ch.getScore(sr_sol);
              exp.write(label + "_SR_found", N, 1);
              exp.write(label + "_SR_time", N, sr.exec_sr);
              exp.write(label + "_SR_ilp", N, sr_score);
              exp.write(label + "_SR_ground", N, ground);
              exp.write(label + "_SR_igap", N, intGap(sr_score, ground));
            } else{
              exp.write(label + "_SR_nofound", N, 1);
            }
          }
          cout << "sr!!\n";
        }
        if (o == 9){
          exp.dropDlvPartition();
          exp.dropGeneratedTable();
        }
      }
    }
  }
}

// /******************************************/
// False infeasibility as hardness increases for Q1,Q2,Q3,Q4 as 4 cores
void E2(){
  string exp_name = "E2";
  DetExp exp = DetExp(exp_name);
  exp.o = 6;
  int R = 20;
  for (int q = 3; q >= 0; q --){
    exp.q = q;
    exp.E = exp.Es[exp.q];
    for (int i = 1; i <= R; i ++){
      exp.seed = i;
      DetSql det_sql = exp.generate();
      exp.kdPartition();
      exp.dlvPartition();
      for (auto H : exp.H8){
        exp.H = H;
        string label = fmt::format("{}{}_H{}", exp.datasets[exp.q], q/2, exp.H);
        LsrProb lsr_prob = LsrProb(det_sql, "", exp.seed);
        lsr_prob.generateBounds(exp.E, exp.a, exp.H);
        DetProb det_prob = DetProb(det_sql, -1, exp.seed);
        det_prob.copyBounds(lsr_prob.bl, lsr_prob.bu, lsr_prob.cl, lsr_prob.cu);
        lsr_prob.partition_name = exp.getDlvPartitionName();
        LayeredSketchRefine lsr = LayeredSketchRefine(exp.C, lsr_prob, exp.S, true);
        lsr_prob.partition_name = exp.getKdPartitionName();
        SketchRefine sr = SketchRefine(lsr_prob);
        map<long long, long long> sr_sol;
        bool is_success = sr.sketchAndRefine(sr_sol);
        bool is_positive = false;
        LsrChecker ch = LsrChecker(lsr_prob);
        Checker det_ch = Checker(det_prob);

        if (lsr.status == Found){
          is_positive = true;
        } else{
          GurobiSolver gs = GurobiSolver(det_prob);
          is_positive = gs.hasIlpSolution(1800);
        }

        if (is_success && ch.checkIlpFeasibility(sr_sol)==Feasibility){
          exp.write(label + "_SR_found", exp.H, 1.0);
        } else{
          exp.write(label + "_SR_nofound", exp.H, 1.0);
        }
        
        if (lsr.status == Found){
          exp.write(label + "_LSR_found", exp.H, 1.0);
          exp.write(label + "_GDR_found", exp.H, 1.0);
        } else {
          exp.write(label + "_LSR_nofound", exp.H, 1.0);
          if (is_positive){
            exp.write(label + "_GDR_found", exp.H, 1.0);
          } else{
            exp.write(label + "_GDR_nofound", exp.H, 1.0);
          }
        }
      }
    }
  }
}

/******************************************/
// Reset database
void test(){
  PgManager pg = PgManager();
  auto tables = pg.listTables();
  for (string table : tables){
    if (table != "ssds" && table != "tpch"){
      cout << table << endl;
      pg.dropTable(table);
    }
  }
}
/******************************************/

void T1(){
  string exp_name = "T1";
  DetExp exp = DetExp(exp_name);
  exp.o = 7;
  int R = 5;
  for (int q = 0; q < 2; q ++){
    exp.q = q;
    exp.E = exp.Es[exp.q];
    for (int i = 1; i <= R; i ++){
      exp.seed = i;
      DetSql det_sql = exp.generate();
      exp.kdPartition();
      LsrProb lsr_prob = LsrProb(det_sql, exp.getKdPartitionName(), exp.seed);
      string label = exp.datasets[exp.q];
      for (auto H : exp.H7){ 
        exp.H = H;
        lsr_prob.generateBounds(exp.E, exp.a, exp.H);
        SketchRefine sr = SketchRefine(lsr_prob, 1);
        map<long long, long long> sr_sol;
        bool is_success = sr.sketchAndRefine(sr_sol);    
        TmpSketchRefine tmp_sr = TmpSketchRefine(lsr_prob, 80);
        map<long long, long long> tmp_sr_sol;
        bool tmp_is_success = tmp_sr.sketchAndRefine(tmp_sr_sol); 
        LsrChecker ch = LsrChecker(lsr_prob);

        double ground = 0;

        if (is_success && ch.checkIlpFeasibility(sr_sol)==Feasibility){
          double sr_score = ch.getScore(sr_sol);
          ground = sr_score;
          exp.write(label + "_SR_found", H, 1);
          exp.write(label + "_SR_time", H, sr.exec_sr);
          exp.write(label + "_SR_ilp", H, sr_score);
          exp.write(label + "_SR_ground", H, ground);
          exp.write(label + "_SR_igap", H, intGap(sr_score, ground));
        } else{
          exp.write(label + "_SR_nofound", H, 1);
        }

        if (tmp_is_success && ch.checkIlpFeasibility(tmp_sr_sol)==Feasibility){
          double tmp_sr_score = ch.getScore(tmp_sr_sol);
          exp.write(label + "_TMPSR_found", H, 1);
          exp.write(label + "_TMPSR_time", H, tmp_sr.exec_sr);
          exp.write(label + "_TMPSR_ilp", H, tmp_sr_score);
          exp.write(label + "_TMPSR_ground", H, ground);
          exp.write(label + "_TMPSR_igap", H, intGap(tmp_sr_score, ground));
        } else{
          exp.write(label + "_TMPSR_nofound", H, 1);
        }
      }
    }
  }
}

/******************************************/

int main() {
  // test();
  G1();
  G2();
  M1();
  M2();
  M3();
  M4();
  M5();
  M6();
  // M7();
  M8();
  E1();
  E2();
  // T1();
  return 0;
}