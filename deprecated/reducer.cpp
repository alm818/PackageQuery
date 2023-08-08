#define FMT_HEADER_ONLY

#include <iostream>
#include <float.h>
#include <chrono>
#include <random>
#include <unordered_map>
#include <numeric>

#include "fmt/core.h"
#include "pseudo_walker.h"
#include "gurobi_solver.h"
#include "simplex.h"
#include "utility.h"
#include "reducer.h"
#include "omp.h"

using namespace Eigen;
using namespace std;

#define tab(row, col) (simplex.tableau[simplex.numcols*(row)+(col)])

Reducer::~Reducer(){
}

Reducer::Reducer(int core, MatrixXd* AA, VectorXd* bb, VectorXd* cc, VectorXd* uu): A(AA), b(bb), c(cc), u(uu){
  MatrixXd* nextA;
  VectorXd*nextb, *nextc, *nextu;
  int nn = cc->size();
  vector<int> original_index (nn);
  vector<int> next_original_index;
  //iota(original_index.begin(), original_index.end(), 0);
  best_x.resize(nn);
  best_score = 0;
  int m = b->size();
  reduce_count = 0;
  status = LS_NOT_FOUND;
  while (true){
    int n = c->size();
    if (n < kReduce) break;
    Simplex simplex = Simplex(core, *A, *b, *c, *u);
    if (simplex.status == LS_FOUND){
      VectorXd centroid_dir (n); centroid_dir.fill(0);
      VectorXd r0 (n);
      unordered_map<int, int> inv_bhead;
      for (int i = 0; i < m; i ++) inv_bhead[simplex.bhead[i]] = i;
      int stay_count = 0;
      #pragma omp parallel num_threads(core)
      {
        if (reduce_count == 0){
          #pragma omp for nowait
          for (int i = 0; i < nn; i ++) original_index[i] = i;
        }
        double whole;
        int local_stay_count = 0;
        #pragma omp for nowait
        for (int i = 0; i < n; i ++){
          whole = round(tab(1, i));
          if (isEqual(whole-tab(1, i), 0)) r0(i) = whole;
          else{
            r0(i) = -1;
            local_stay_count ++;
          }
        }
        #pragma omp atomic
        stay_count += local_stay_count;
        #pragma omp barrier
        unordered_map<int, double> local_centroid_dir;
        #pragma omp for nowait
        for (int i = 0; i < m+n; i ++){
          if (inv_bhead.find(i) == inv_bhead.end()){
            // Non-basic
            int index = i;
            for (int j = 0; j < m; j ++){
              if (i > simplex.bhead[j]) index --;
            }
            double dir = -1.0;
            if (isEqual(tab(1, i), 0)) dir = 1.0;
            if (i < n) local_centroid_dir[i] += (dir/n);
            for (int j = 2; j <= m+1; j ++){
              int basic_col = simplex.bhead[j-2];
              if (basic_col < n) local_centroid_dir[basic_col] -= (tab(j, i) * dir / n);
            }
          }
        }
        for (const auto& pi : local_centroid_dir){
          #pragma omp atomic
          centroid_dir(pi.first) += pi.second;
        }
      }
      PseudoWalker walker = PseudoWalker(centroid_dir, true, core);
      int step_horizon = (int)ceil(n / log2(n) * kRedundant);
      for (int i = 0; i < step_horizon; i ++){
        int index = abs(walker.step())-1;
        if (r0(index) >= 0){
          r0(index) = -1;
          stay_count ++;
        }
      }
      vector<int> stay_index (stay_count);
      next_original_index.resize(stay_count);
      nextA = new MatrixXd(m, stay_count);
      nextc = new VectorXd(stay_count);
      nextu = new VectorXd(stay_count);
      nextb = new VectorXd(m);
      for (int i = 0; i < m; i ++){
        (*nextb)(i) = (*b)(i);
      }
      #pragma omp parallel num_threads(core)
      {
        int local_start_index = -1;
        int local_stay_count = 0;
        vector<int> local_stay_index (stay_count);
        VectorXd local_bl (m); local_bl.fill(0);
        #pragma omp master
        {
          stay_count = 0;
        }
        #pragma omp barrier
        #pragma omp for nowait
        for (int i = 0; i < n; i ++){
          if (r0(i) < 0){
            local_stay_index[local_stay_count] = i;
            local_stay_count ++;
          } else{
            best_x(original_index[i]) = r0(i);
            #pragma omp atomic
            best_score += r0(i)*((*c)(i));
            for (int j = 0; j < m; j ++){
              local_bl(j) += (*A)(j, i) * r0(i);
            }
          }
        }
        #pragma omp critical
        {
          local_start_index = stay_count;
          stay_count += local_stay_count;
        }
        for (int i = 0; i < m; i ++){
          #pragma omp atomic
          (*nextb)(i) -= local_bl(i);
        }
        #pragma omp barrier
        for (int i = local_start_index; i < local_start_index + local_stay_count; i ++){
          stay_index[i] = local_stay_index[i-local_start_index];
          next_original_index[i] = original_index[stay_index[i]];
        }
        #pragma omp barrier
        #pragma omp for
        for (int i = 0; i < stay_count; i ++){
          original_index[i] = next_original_index[i];
          for (int j = 0; j < m; j ++){
            (*nextA)(j, i) = (*A)(j, stay_index[i]);
          }
          (*nextc)(i) = (*c)(stay_index[i]);
          (*nextu)(i) = (*u)(stay_index[i]);
        }
        #pragma omp master
        {
          if (reduce_count > 0){
            delete A;
            delete b;
            delete c;
            delete u;
          }
          A = nextA;
          b = nextb;
          c = nextc;
          u = nextu;
        }
        #pragma omp barrier
      }
      reduce_count ++;
    } else{
      status = simplex.status;
      break;
    }
  }
  if (status == LS_NOT_FOUND){
    GurobiSolver gs = GurobiSolver(*A, *b, *c, *u);
    gs.solveIlp();
    if (gs.ilp_status == LS_FOUND){
      for (int i = 0; i < gs.n; i ++){
        best_x(original_index[i]) = gs.x0(i);
        best_score += gs.x0(i) * ((*c)(i));
      }
    } else status = gs.ilp_status;
  }
}
