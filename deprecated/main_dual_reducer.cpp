#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <queue>
#include <random>
#include <chrono>
#include <numeric>

#include "map_sort.h"
#include "utility.h"
#include "parallel_pq.h"
#include "fmt/core.h"
#include "fitsio.h"
#include "omp.h"
#include "lattice_solver.h"
#include "pseudo_walker.h"
#include "simplex.h"
#include "gurobi_lattice_solver.h"
#include "gurobi_solver.h"
#include "dual.h"
#include "dual_reducer.h"
#include "reducer.h"
#include "fmt/core.h"

using namespace std;
using namespace Eigen;

void generateNormalProlem(int expected_n, double outlier_prob, double att_var, int n, MatrixXd& A, VectorXd& bl, VectorXd& bu, VectorXd& c){
  A.resize(4, n); bl.resize(4); bu.resize(4); c.resize(n); 
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  // double multiplicity = att_var * 12;
  double one_mean = 100;
  normal_distribution att_dist(one_mean, att_var);
  int expected_numvar = expected_n;
  double mean = one_mean*expected_numvar;
  double var = att_var*expected_numvar;
  normal_distribution n_dist(0.0, var);
  normal_distribution n_dist_c(0.0, 1.0);
  #pragma omp parallel for num_threads(CORE_COUNT)
  for (int i = 0; i < n; i ++){
    A(0, i) = att_dist(gen);
    A(1, i) = att_dist(gen);
    A(2, i) = att_dist(gen);
    A(3, i) = 1;
    c(i) = n_dist_c(gen);
  }
  double tol = var*sqrt(1/outlier_prob);
  bl(0) = mean - tol;
  bu(0) = mean + tol;
  bl(1) = -DBL_MAX;
  bu(1) = mean + tol;
  bl(2) = mean- tol;
  bu(2) = DBL_MAX;
  bl(3) = expected_n/2;
  bu(3) = expected_n*3/2;
}

void generateBoundedProlem(int expected_n, double outlier_prob, double att_var, int n, MatrixXd& A, VectorXd& bl, VectorXd& bu, VectorXd& c){
  A.resize(4, n); bl.resize(4); bu.resize(4); c.resize(n); 
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(-1.0, 1.0);
  double multiplicity = att_var * 6;
  //normal_distribution att_dist(10.0, att_var);
  int expected_numvar = expected_n;
  double mean = 0.5*multiplicity*expected_numvar;
  double var = att_var*expected_numvar;
  normal_distribution n_dist(0.0, var);
  //normal_distribution n_dist_c(0.0, 1.0/12);
  uniform_real_distribution n_dist_c(-1.0, 1.0);
  #pragma omp parallel for num_threads(CORE_COUNT)
  for (int i = 0; i < n; i ++){
    A(0, i) = u_dist(gen)*multiplicity;
    A(1, i) = u_dist(gen)*multiplicity;
    A(2, i) = u_dist(gen)*multiplicity;
    A(3, i) = 1;
    c(i) = n_dist_c(gen)*multiplicity;
    //c(i) = A(0, i) + A(1, i) + A(2, i);
  }
  double tol = var*sqrt(1/outlier_prob);
  bl(0) = mean - tol;
  bu(0) = mean + tol;
  bl(1) = -DBL_MAX;
  bu(1) = mean + tol;
  bl(2) = mean- tol;
  bu(2) = DBL_MAX;
  bl(3) = expected_n/2;
  bu(3) = expected_n*3/2;
}

void centroid_test(){
  int core = 16;
  int n = 2;
  int m = 4;
  MatrixXd A (m, n); A << 1.0/3, 1, 4.0/3, -1, -1.0/4, 1, 1, -1;
  VectorXd bl (m); bl << 8.0/3, -DBL_MAX, 1.0/2, 0;
  VectorXd bu (m); bu << DBL_MAX, 17.0/3, 6, DBL_MAX;
  VectorXd c (n); c << 1, 1;
  VectorXd l (n); l << 0, 0;
  VectorXd u (n); u << 10, 10;

  Dual* dual = new Dual(core, A, bl, bu, c, l, u);
  print(dual->x);

  MatrixXd BA = -dual->Binv * dual->A;
  cout << BA << endl;
  cout << endl;
  cout << dual->Binv << endl;
  cout << "BASIS" << endl;
  print(dual->bhead);
  cout << endl;

  if (dual->status == LS_FOUND){
    VectorXd norms (m+n); norms.fill(0);
    VectorXd centroid_dir (n); centroid_dir.fill(0);
    VectorXd r0 (n);
    VectorXd stay_r0 (n);
    VectorXd sum_row(m); sum_row.fill(0);
    int stay_count = 0;
    #pragma omp parallel num_threads(core)
    {
      double whole, frac;
      int local_stay_count = 0;
      #pragma omp for nowait
      for (int i = 0; i < n; i ++){
        r0(i) = dual->x(i);
        frac = modf(dual->x(i), &whole);
        if (isEqual(frac, 0)) stay_r0(i) = whole;
        else{
          stay_r0(i) = -1;
          local_stay_count ++;
        }
      }
      #pragma omp atomic
      stay_count += local_stay_count;
      #pragma omp for nowait
      for (int i = 0; i < m+n; i ++){
        if (!dual->inv_bhead[i]){
          double norm = 0;
          if (i < n){
            norm ++;
            for (int j = 0; j < m; j ++){
              if (dual->bhead(j) < n){
                double val = 0;
                for (int k = 0; k < m; k ++){
                  val += dual->Binv(j, k) * dual->A(k, i);
                }
                norm += val*val;
              }
            }
          } else{
            int index = i - n;
            for (int j = 0; j < m; j ++){
              if (dual->bhead(j) < n){
                norm += dual->Binv(j, index)*dual->Binv(j, index);
              }
            }
          }
          norms(i) = sqrt(norm);
        }
      }
      #pragma omp barrier
      for (int j = 0; j < m; j ++){
        double local_sum_row = 0;
        #pragma omp for nowait
        for (int i = 0; i < m+n; i ++){
          if (!dual->inv_bhead[i]){
            // Non-basic
            if (i < n){
              if (isEqual(dual->x(i), dual->l(i))) local_sum_row += dual->A(j, i) / norms(i);
              else local_sum_row -= dual->A(j, i) / norms(i);
            } else{
              int index = i - n;
              if (index == j){
                if (isEqual(dual->x(i), dual->bl(index))) local_sum_row -= 1.0 / norms(i);
                else local_sum_row += 1.0 / norms(i);
              }
            }
          }
        }
        #pragma omp atomic
        sum_row(j) += local_sum_row;
      }
      #pragma omp for nowait
      for (int i = 0; i < n; i ++){
        if (!dual->inv_bhead[i]){
          // Non-basic
          if (isEqual(dual->x(i), dual->l(i))) centroid_dir(i) = 1.0/norms(i)/n;
          else centroid_dir(i) = -1.0/norms(i)/n;
        }
      }
      #pragma omp barrier
      #pragma omp master
      {
        // cout << "SUM ROW" << endl;
        // print(sum_row);
        for (int i = 0; i < m; i ++){
          if (dual->bhead(i) < n){
            double val = 0;
            for (int j = 0; j < m; j ++){
              val += dual->Binv(i, j) * sum_row(j);
            }
            centroid_dir(dual->bhead(i)) = val / n;
          }
        }
      }
      #pragma omp barrier
      #pragma omp single
      {
        print(centroid_dir);
      }
    }
  }
}

bool verify(VectorXd x, const MatrixXd& A, const VectorXd& bl, const VectorXd& bu, const VectorXd& c, const VectorXd& l, const VectorXd& u){
  int n = l.size();
  VectorXd sol = x(seqN(0, n));
  int m = bl.size();
  for (int i = 0; i < n; i ++){
    if (isLess(sol(i), l(i)) || isGreater(sol(i),u(i))){
      cout << "N:" << i << " VIOLATIONS" << endl;
      return false;
    }
  }
  VectorXd b = A * sol;
  for (int i = 0; i < m; i ++){
    if (isLess(b(i), bl(i)) || isGreater(b(i), bu(i))){
      cout << "M:" << i << " VIOLATIONS" << endl;
      return false;
    }
  }
  return true;
}

void quickRun(){
  int n = 10000;
  MatrixXd A;
  VectorXd bl, bu, c;
  generateBoundedProlem(500, 0.7, 1, n, A, bl, bu, c);
  //generateNormalProlem(15, 0.8, 4, n, A, bl, bu, c);
  VectorXd u (n); u.fill(1);
  VectorXd l (n); l.fill(0);

  // VectorXd* a = new VectorXd(10);
  // a->fill(10);
  // cout << (*a)(1) << endl;
  // VectorXd* b = new VectorXd(10);
  // memcpy(&(*b)(0), &(*a)(0), 10*8);
  // (*a)(1) = 12;
  // print(*a);
  // print(*b);

  GurobiSolver gs = GurobiSolver(A, bl, bu, c, l, u);
  gs.solveRelaxed();
  // print(gs.r0);
  cout << gs.exe_relaxed << " " << gs.iteration_count << endl;
  cout << gs.relaxed_cscore << endl;

  // cout << "-------------------" << endl;
  DualReducer dr2 = DualReducer(8, &A, &bl, &bu, &c, &l, &u, 1);
  for (auto oi : dr2.original_indices){
    cout << oi->size() << " ";
  }
  cout << endl;
  cout << solMessage(dr2.status) << " " << dr2.exe_solve << " " << dr2.duals[0]->exe_solve << endl;
  
  DualReducer dr = DualReducer(8, &A, &bl, &bu, &c, &l, &u, 0);
  for (auto oi : dr.original_indices){
    cout << oi->size() << " ";
  }
  cout << endl;
  cout << solMessage(dr.status) << " " << dr.exe_solve << " " << dr.duals[0]->exe_solve << endl;
  double r0_obj = dr.duals[0]->score;
  double gap = (r0_obj - dr.best_score) / r0_obj * 100;
  cout << dr.best_score << " " << r0_obj << " " << gap << "%" << endl;
  //cout << solCombination(dr.best_x) << endl;
  verify(dr.best_x, A, bl, bu, c, l, u);

  double gap2 = (r0_obj - dr2.best_score) / r0_obj * 100;
  cout << dr2.best_score << " " << r0_obj << " " << gap2 << "%" << endl;
  //cout << solCombination(dr2.best_x) << endl;
  verify(dr2.best_x, A, bl, bu, c, l, u);
  GurobiSolver igs = GurobiSolver(A, bl, bu, c, l, u);
  igs.solveIlp();
  double gap3 = (igs.ilp_cscore - dr.best_score) / igs.ilp_cscore * 100;
  double gap4 = (igs.ilp_cscore - dr2.best_score) / igs.ilp_cscore * 100;
  cout << igs.ilp_cscore << " " << gap3 << "%" << " " << gap4 << "%" << endl;
  cout << igs.exe_ilp + igs.exe_init << endl;

  // int n = 2;
  // int m = 4;
  // MatrixXd A (m, n); A << 1.0/3, 1, 4.0/3, -1, -1.0/4, 1, 1, -1;
  // VectorXd bl (m); bl << 8.0/3, -DBL_MAX, 1.0/2, 0;
  // VectorXd bu (m); bu << DBL_MAX, 17.0/3, 3, DBL_MAX;
  // VectorXd c (n); c << 1, 1;
  // VectorXd l (n); l << 0, 0;
  // VectorXd u (n); u << 10, 10;

  // Dual dual = Dual(16, A, bl, bu, c, l, u);
  // cout << solMessage(dual.status) << endl;
  // cout << dual.exe_solve << " " << dual.iteration_count << " " << dual.mini_iteration_count << endl;
  // cout << dual.score << endl;

  // GurobiSolver gs = GurobiSolver(A, bl, bu, c, l, u);
  // gs.solveRelaxed();
  // cout << gs.exe_relaxed << " " << gs.iteration_count << endl;
  // cout << gs.relaxed_cscore << endl;
  // // print(dual.x);
  // // print(gs.r0);
  // if (gs.relaxed_status == LS_FOUND && !isEqual(gs.relaxed_cscore, dual.score, 1e-8)){
  //   for (int i = 0; i < n; i ++){
  //     assert(isEqual(dual.x(i), gs.r0(i), 1e-6));
  //   }
  // }
}

void L3cache(){
  int core = 8;
  int L3 = 8 * 1024 * 1024;
  int max_double = (int)floor(L3 / 8);
  vector<string> names = {"0", "1", "2", "3", "4", "5", "6", "7"};
  Profiler pro = Profiler(names);
  int m = 4;
  int n = max_double * 10;
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  vector<VectorXd> arr (m);
  for (int i = 0; i < m; i ++) arr[i].resize(n);
  for (int i = 0; i < m; i ++){
    #pragma omp parallel for num_threads(core)
    for (int j = 0; j < n; j ++) arr[i](j) = u_dist(gen);
  }
  pro.clock(0);
  #pragma omp parallel for num_threads(core)
  for (int i = 0; i < n; i ++){
    for (int j = 0; j < m; j ++){
      arr[j](i) *= 2;
      arr[j](i) /= 3;
    }
  }
  pro.stop(0);
  pro.clock(1);
  #pragma omp parallel num_threads(core)
  {
    for (int j = 0; j < m; j ++){
      #pragma omp for nowait
      for (int i = 0; i < n; i ++){
        arr[j](i) *= 2;
        arr[j](i) /= 3;
      }
    }
  }
  pro.stop(1);
  pro.clock(2);
  #pragma omp parallel for num_threads(core)
  for (int i = 0; i < n; i ++){
    for (int j = 1; j < m; j ++){
      arr[j](i) *= 2;
      arr[j](i) -= arr[0](i);
    }
  }
  pro.stop(2);
  pro.clock(3);
  #pragma omp parallel num_threads(core)
  {
    for (int j = 1; j < m; j ++){
      #pragma omp for nowait
      for (int i = 0; i < n; i ++){
        arr[j](i) *= 2;
        arr[j](i) -= arr[0](i);
      }
    }
  }
  pro.stop(3);
  pro.print();
}

void testL3Cache(){
  chrono::high_resolution_clock::time_point start, end;
  int N = 100000000; int r = 10;
  VectorXd u (N); u.fill(0.24);
  VectorXd g1 (N); VectorXd g2 (N); VectorXd g3 (N);
  start = chrono::high_resolution_clock::now();
  // for (int i = 0; i < r; i ++){
  //   VectorXd g1 = u;
  //   VectorXd g2 = u+u;
  //   VectorXd g3 = u+u+u;
  // }
  // L3 Effect
  #pragma omp parallel num_threads(16)
  {
    #pragma omp for nowait
    for (int i = 0; i < N; i ++){
      g1(i) = u(i)+u(i);
    }
    // #pragma omp for nowait
    // for (int i = 0; i < N; i ++){
    //   g2(i) = u(i) + u(i);
    // }
    #pragma omp for nowait
    for (int i = 0; i < N; i ++){
      g3(i) = u(i) + u(i) + u(i);
    }
  }
  end = chrono::high_resolution_clock::now();
  double exe = chrono::duration_cast<chrono::nanoseconds>(end - start).count() / 1000000.0;
  cout << exe << endl;

  VectorXd l1 (N); VectorXd l2 (N); VectorXd l3 (N);
  cout << "START" << endl;
  start = chrono::high_resolution_clock::now();
  VectorXd uu = u+u;
  VectorXd uuu = u+u+u;
  // #pragma omp parallel num_threads(16)
  // {
  //   #pragma omp for
  //   for (int i = 0; i < N; i ++){
  //     for (int j = 0; j < r; j ++){
  //       l1(i) = u(i);
  //       l2(i) = u(i)+u(i);
  //       l3(i) = u(i)+u(i)+u(i);
  //     }
  //   }
  // }
  end = chrono::high_resolution_clock::now();
  double exe2 = chrono::duration_cast<chrono::nanoseconds>(end - start).count() / 1000000.0;
  cout << exe2 << endl;
}

int main(){
  quickRun();
  //centroid_test();
  //L3cache();
  //testL3Cache();
}