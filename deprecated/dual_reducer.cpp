#include "pb/core/gurobi_solver.h"
#include "pb/core/pseudo_walker.h"
#include "pb/util/unumeric.h"
#include "dual_reducer.h"

int DualReducer::kMinIlp = 500;
int DualReducer::kMaxIlp = 5000;
double DualReducer::kEpsilon = 1e-8;
enum StayStatus {NotStay, LikelyStay, Stay};

DualReducer::~DualReducer(){
  for (int i = 1; i < (int) probs.size(); i ++) delete probs[i];
  for (int i = 0; i < (int) duals.size(); i ++) delete duals[i];
  for (int i = 0; i < (int) indices.size(); i ++) delete indices[i];
}

void DualReducer::filterStay(int core, Dual* dual, DetProb* prob, VectorXi* cur_index, int stay_mode, int stay_count, const VectorXi& stay){
  // Stage 2
  int m = (int) prob->bl.size();
  int n = (int) prob->l.size();
  VectorXi stay_index (stay_count);
  VectorXi *next_index = new VectorXi(stay_count);
  indices.push_back(next_index);
  DetProb *nextProb = new DetProb(m, stay_count);
  nextProb->bl = prob->bl;
  nextProb->bu = prob->bu;
  probs.push_back(nextProb);
  int start_stay_count = 0;
  // pro.stop(5, false);
  // cout << "B " << layer_count << endl;
  #pragma omp parallel num_threads(core)
  {
    // pro.clock(6);
    int local_start_index = -1;
    int local_stay_count = 0;
    VectorXi local_stay_index (stay_count);
    VectorXd local_change_b (m); local_change_b.fill(0);
    #pragma omp for nowait
    for (int i = 0; i < n; i ++){
      // #pragma omp critical
      // cout << stay.size() << " " << dual->sol.size() << " " << best_sol.size() << " " << local_change_b.size() << " " << m << " " << n;
      if (stay(i) == stay_mode || stay(i) == Stay){
        local_stay_index(local_stay_count) = i;
        local_stay_count ++;
      } else{
        best_sol((*cur_index)(i)) = dual->sol(i);
        for (int j = 0; j < m; j ++){
          local_change_b(j) += prob->A(j, i) * dual->sol(i);
        }
      }
    }
    // #pragma omp single
    // cout << "B1 " << layer_count << endl;
    #pragma omp critical
    {
      // cout << start_stay_count << " " << local_stay_count << " " << stay_count << endl;
      local_start_index = start_stay_count;
      start_stay_count += local_stay_count;
    }
    // #pragma omp single
    // cout << "B11 " << layer_count << endl;
    for (int i = local_start_index; i < local_start_index + local_stay_count; i ++){
      // assert(local_stay_index(i-local_start_index) >= 0 && local_stay_index(i-local_start_index) < n);
      stay_index(i) = local_stay_index(i-local_start_index);
      (*next_index)(i) = (*cur_index)(stay_index(i));
    }
    // #pragma omp single
    // cout << "B12 " << layer_count << endl;
    #pragma omp barrier
    for (int i = 0; i < m; i ++){
      if (prob->bl(i) != -DBL_MAX){
        #pragma omp atomic
        nextProb->bl(i) -= local_change_b(i);
      }
      if (prob->bu(i) != DBL_MAX){
        #pragma omp atomic
        nextProb->bu(i) -= local_change_b(i);
      }
    }
    // #pragma omp single
    // cout << "B13 " << layer_count << endl;
    // #pragma omp single
    // {
    //   cout << "FIRST " << stay_count << " " << m << endl;
    //   cout << nextProb->A.innerSize() << " " << nextProb->A.outerSize() << endl;
    //   cout << nextProb->c.size() << " " << nextProb->l.size() << " " <<  nextProb->u.size()<< endl;
    //   cout << stay_index
    // }
    // pro.stop(6);
    // pro.clock(7);
    #pragma omp for nowait
    for (int i = 0; i < start_stay_count; i ++){
      // if (!(stay_index(i) < n && stay_index(i) >= 0)){
      //   cout << stay_index(i) << " " << i << " " << n << endl;
      //   assert(0==1);
      // }
      for (int j = 0; j < m; j ++){
        nextProb->A(j, i) = prob->A(j, stay_index(i));
      }
      nextProb->c(i) = prob->c(stay_index(i));
      nextProb->l(i) = prob->l(stay_index(i));
      nextProb->u(i) = prob->u(stay_index(i));
    }
    // pro.stop(7);
    // #pragma omp barrier
    // #pragma omp single
    // cout << "B2 " << layer_count << endl;
  }
  layer_count ++;
}

DualReducer::DualReducer(int core, DetProb *origin_prob){
  // vector<string> names = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
  // Profiler pro = Profiler(names);
  std::chrono::high_resolution_clock::time_point start;
  start = std::chrono::high_resolution_clock::now();
  // pro.clock(0, false);
  probs.push_back(origin_prob);

  int nn = (int) origin_prob->c.size();
  int m = (int) origin_prob->bl.size();

  indices.push_back(new VectorXi(nn));
  best_sol.resize(nn);
  best_score = 0;
  layer_count = 0;
  status = NotFound;

  DetProb* prob;
  VectorXi* cur_index;
  // pro.stop(0, false);
  while (true){
    //cout << "Z " << layer_count << endl;
    prob = probs[layer_count];
    cur_index = indices[layer_count];
    //cout << "Z1 " << layer_count << endl;
    int n = (int) prob->c.size();
    if (n < kMinIlp){
      if (layer_count == 0){
        // cout << "Z11 " << layer_count << endl;
        for (int i = 0; i < nn; i ++) (*cur_index)(i) = i;
        // pro.clock(1, false);
        Dual* dual = new Dual(core, *prob);
        // pro.stop(1, false);
        duals.push_back(dual);
        // cout << "Z12 " << layer_count << endl;
      }
      break;
    }
    //cout << "Z2 " << layer_count << endl;
    // pro.clock(1, false);
    Dual* dual = new Dual(core, *prob);
    // pro.stop(1, false);
    duals.push_back(dual);
    //cout << "A " << layer_count << endl;

    if (dual->status == Found){
      VectorXd norms (n+m); norms.fill(0);
      VectorXd centroid_dir (n); centroid_dir.fill(0);
      VectorXd sum_row(m); sum_row.fill(0);
      VectorXi stay (n);
      int stay_count = 0;
      int frac_count = 0;
      #pragma omp parallel num_threads(core)
      {
        // pro.clock(2);
        if (layer_count == 0){
          #pragma omp for nowait
          for (int i = 0; i < nn; i ++) (*cur_index)(i) = i;
        }
        double whole, frac;
        int local_stay_count = 0;
        int local_frac_count = 0;
        #pragma omp for nowait
        for (int i = 0; i < n; i ++){
          frac = modf(dual->sol(i), &whole);
          if (isEqual(frac, 0, kEpsilon)){
            if (whole == 0) stay(i) = NotStay;
            else{
              stay(i) = LikelyStay;
              local_stay_count ++;
            }
          } else{
            stay(i) = Stay;
            local_stay_count ++;
            local_frac_count ++;
          }
        }
        #pragma omp atomic
        stay_count += local_stay_count;
        #pragma omp atomic
        frac_count += local_frac_count;
        #pragma omp barrier
        if (stay_count < kMinIlp || stay_count > kMaxIlp){
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
                      val += dual->Binv(j, k) * prob->A(k, i);
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
                  if (isEqual(dual->sol(i), prob->l(i), kEpsilon)) local_sum_row += prob->A(j, i) / norms(i);
                  else local_sum_row -= prob->A(j, i) / norms(i);
                } else{
                  int index = i - n;
                  if (index == j){
                    if (isEqual(dual->sol(i), prob->bl(index), kEpsilon)) local_sum_row -= 1.0 / norms(i);
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
              if (isEqual(dual->sol(i), prob->l(i), kEpsilon)) centroid_dir(i) = 1.0/norms(i)/n;
              else centroid_dir(i) = -1.0/norms(i)/n;
            }
          }
          #pragma omp barrier
          #pragma omp master
          {
            //cout << "A1 " << layer_count << endl;
            for (int i = 0; i < m; i ++){
              if (dual->bhead(i) < n){
                double val = 0;
                for (int j = 0; j < m; j ++){
                  val += dual->Binv(i, j) * sum_row(j);
                }
                //cout << "HERE " << i << " " << dual->bhead(i) << " VAL:" << val << endl;
                centroid_dir(dual->bhead(i)) = val / n;
              }
            }
          }
          // pro.stop(2);
        }
      }
      int stay_mode = LikelyStay;
      VectorXi max_stay;
      int max_stay_count = stay_count;
      if (stay_count < kMinIlp){
        max_stay = stay;
        PseudoWalker walker = PseudoWalker(centroid_dir, true, core);
        while (true){
          int step = walker.step();
          int index = abs(step)-1;
          if (stay(index) == NotStay){
            if (stay_count < kMinIlp){
              stay(index) = LikelyStay;
              stay_count ++;
            }
            max_stay(index) = LikelyStay;
            max_stay_count ++;
            if (max_stay_count == kMaxIlp) break;
          }
        }
      } else if (stay_count > kMaxIlp){
        stay_mode = Stay;
        PseudoWalker walker = PseudoWalker(centroid_dir, true, core);
        stay_count = frac_count;
        while (true){
          int step = walker.step();
          int index = abs(step)-1;
          if (stay(index) == LikelyStay){
            stay(index) = Stay;
            stay_count ++;
            if (stay_count == kMaxIlp) break;
          }
        }
      }
      // pro.stop(4, false);
      // pro.clock(5, false);
      if (max_stay.size() == n) filterStay(core, dual, prob, cur_index, stay_mode, max_stay_count, max_stay);
      filterStay(core, dual, prob, cur_index, stay_mode, stay_count, stay);
      if (stay_count == n || (kMinIlp <= stay_count && stay_count <= kMaxIlp)){
        break;
      }
      //cout << "B3 " << layer_count << endl;
    } else{
      status = dual->status;
      break;
    }
  }
  //cout << "C " << layer_count << endl;
  //cout << bl->size() << " " << A->innerSize() << " " << A->outerSize() << endl;
  // pro.clock(8, false);

  // if (prob->c.size() == 0){
  //   status = Found;
  //   #pragma omp parallel num_threads(core)
  //   {
  //     double local_best_score = 0;
  //     #pragma omp for nowait
  //     for (int i = 0; i < nn; i ++) local_best_score += best_sol(i) * origin_prob->c(i);
  //     #pragma omp atomic
  //     best_score += local_best_score;
  //   }
  // }
  if (status == NotFound){
    // cout << "HERE " << prob->c.size() << endl;
    //cout << "C02 " << layer_count << endl;
    for (int j = layer_count; j >= 0; j --){
      prob = probs[j];
      cur_index = indices[j];
      //cout << "HERE " << prob->c.size() << endl;
      GurobiSolver gs = GurobiSolver(*prob);
      gs.solveIlp(10.0);
      //cout << "C03" << layer_count << endl;
      status = gs.ilp_status;
      if (gs.ilp_status == Found){
        for (int i = 0; i < gs.ilp_sol.size(); i ++) best_sol((*cur_index)(i)) = gs.ilp_sol(i);
        //cout << "C1 " << layer_count << endl;
        #pragma omp parallel num_threads(core)
        {
          double local_best_score = 0;
          #pragma omp for nowait
          for (int i = 0; i < nn; i ++) local_best_score += best_sol(i) * origin_prob->c(i);
          #pragma omp atomic
          best_score += local_best_score;
        }
        break;
        //cout << "C2 " << layer_count << endl;
      } else if (gs.ilp_status == Timeout) break;
    }
  }
  // pro.stop(8, false);
  // pro.print();
  exe_solve = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000000.0;
}

VectorXd DualReducer::getLpSol(){
  return duals[0]->sol(Eigen::seqN(0, duals[0]->n));
}

double DualReducer::getLpScore(){
  return duals[0]->score;
}

double DualReducer::getLpTime(){
  return duals[0]->exe_solve;
}