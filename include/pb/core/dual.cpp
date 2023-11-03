#include "pb/core/dual.h"
#include "pb/core/gurobi_solver.h"
#include "pb/core/checker.h"
#include "pb/lib/map_sort.h"
#include "pb/util/unumeric.h"
#include "pb/util/udebug.h"

using std::upper_bound;

// Primal feasibility tolerance
const double kE_P = 1e-6;
// Relative primal feasibility tolerance
const double kE_r = 1e-12;
// Pivot tolerance
const double kE_ap = 1e-6;
// Equality tolerance
const double kE_Eq = 1e-6;
// Inequality tolerance
const double kE_Ieq = 1e-6;

double primalInfeasibilities(double x, double l, double u){
  double left = l - l * kE_r - kE_P;
  if (x < left) return x-l;
  double right = u + u * kE_r + kE_P;
  if (x > right) return x-u;
  return 0;
}

double getSlope(const VectorXd &u, const VectorXd &l, const VectorXd &bu, const VectorXd &bl, const VectorXd &alpha_r, int i){
  int n = (int) u.size();
  if (i < n) return (u(i)-l(i)) * fabs(alpha_r(i));
  else{
    int index = i-n;
    return (bu(index)-bl(index)) * fabs(alpha_r(i));
  }
}

// All slopes are non-negative
// We want to keep k smallest scores such that the sum is at least max_slope
void pushHeap(const VectorXd &u, const VectorXd &l, const VectorXd &bu, const VectorXd &bl, const VectorXd &alpha_r, const VectorXd &d, pair<double, int> *heap, int &sz, double &sum_slope, double max_slope, int i){
  double score = d(i) / alpha_r(i);
  double current_slope = getSlope(u, l, bu, bl, alpha_r, i);
  if (sum_slope <= max_slope){
    //cout << "PURE ADD1 " << score << " " << i << endl;
    heap[sz] = {score, i};
    sz ++;
    push_heap(heap, heap+sz);
    sum_slope += current_slope;
  } 
  else{
    // For sure n > 0 since sum_slope > max_slope >= 0
    if (score <= heap[0].first){
      double replaced_slope = getSlope(u, l, bu, bl, alpha_r, heap[0].second);
      double next_sum_slope = sum_slope - replaced_slope + current_slope;
      if (next_sum_slope < max_slope){
        //cout << "PURE ADD2 " << score << " " << i << endl;
        heap[sz] = {score, i};
        sz ++;
        push_heap(heap, heap+sz);
        sum_slope += current_slope;
      } else{
        //cout << "REPLACEMENT " << score << " " << i << " " << heap[0].first << " " << heap[0].second << endl;
        pop_heap(heap, heap+sz);
        heap[sz-1] = {score, i};
        push_heap(heap, heap+sz);
        sum_slope = next_sum_slope;
      }
    }
  }
}

Dual::~Dual(){
}

// Maximizing prob.c 
// Dual minimizing c
Dual::Dual(int core, const DetProb &prob){
  // vector<string> names = {"Copy", "RestInit", "BoundStricten", "InitXBound", "Compute XB", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"};
  // Profiler pro = Profiler(names);
  std::chrono::high_resolution_clock::time_point start;
  start = std::chrono::high_resolution_clock::now();

  status = NotFound;
  iteration_count = 0;
  mini_iteration_count = 0;
  score = 0;

  n = (int) prob.c.size();
  m = (int) prob.bl.size();

  // pro.clock(0);
  bl = prob.bl;
  bu = prob.bu;
  // pro.stop(0);
  // pro.clock(1);
  // Norm square steepest edges
  beta.resize(m); beta.fill(1);

  // Reduced cost
  d.resize(n+m); d.fill(0); 

  // Basis indices
  bhead.resize(m);
  for (int i = 0; i < m; i ++) bhead(i) = i+n;
  // Inv Basis indices
  inv_bhead = boost::dynamic_bitset<>(n+m, 0);
  for (int i = 0; i < m; i ++) inv_bhead[bhead(i)].flip();

  // Primal solution
  sol.resize(n+m); sol.fill(0);

  // Leaving variable index in the basis
  int r;
  // Leaving variable index
  int p;
  // Entering variable index
  int q, q_index;

  // Binv
  Binv = MatrixXd::Identity(m, m);
  // rho_r
  rho_r.resize(m);
  // alpha_r
  alpha_r.resize(n+m);
  // delta of alpha_r
  double max_delta = 0;
  int sign_max_delta = 0;
  
  // init for BRFT PPQ
  pair<double, int> *init;
  // init size
  int init_size = 0;
  // Dual step
  double dual_step = 0;
  // Primal step
  double primal_step = 0;

  // alpha_q
  alpha_q.resize(m);
  // tau for DSE
  tau.resize(m);
  
  // slope_partition
  VectorXi slope_partition (core+1);
  // slope_sums_partition
  VectorXd slope_sums_parition (core);
  // slope_sums
  VectorXd slope_sums;

  // Tilde a
  VectorXd tilde_a (m);

  // cout << "A" << endl;
  // cout << A << endl;
  // cout << "d" << endl;
  // print(d);

  // pro.stop(1);
  #pragma omp parallel num_threads(core)
  {
    // Phase-1 Dual
    // pro.clock(2);
    // Stricten bound
    #pragma omp master
    {
      for (int i = 0; i < m; i ++){
        if (prob.bl(i) == -DBL_MAX) bl(i) = 0;
        if (prob.bu(i) == DBL_MAX) bu(i) = 0;
      }
    }
    #pragma omp barrier
    for (int i = 0; i < m; i ++){
      if (prob.bl(i) == -DBL_MAX){
        double local_bound = 0;
        #pragma omp for nowait
        for (int j = 0; j < n; j ++){
          if (prob.A(i, j) > 0) local_bound += prob.l(j) * prob.A(i, j);
          else local_bound += prob.u(j) * prob.A(i, j);
        }
        #pragma omp atomic
        bl(i) += local_bound;
      }
      if (prob.bu(i) == DBL_MAX){
        double local_bound = 0;
        #pragma omp for nowait
        for (int j = 0; j < n; j ++){
          if (prob.A(i, j) > 0) local_bound += prob.u(j) * prob.A(i, j);
          else local_bound += prob.l(j) * prob.A(i, j);
        }
        #pragma omp atomic
        bu(i) += local_bound;
      }
    }
    // pro.stop(2);
    // Initialize x
    // pro.clock(3);
    #pragma omp for nowait
    for (int i = 0; i < n; i ++){
      d(i) = -prob.c(i);
      if (prob.c(i) > 0) sol(i) = prob.u(i);
      else sol(i) = prob.l(i);
    }
    // pro.stop(3);
    // pro.clock(4);
    #pragma omp barrier
    for (int i = 0; i < m; i ++){
      double local_anxn = 0;
      #pragma omp for nowait
      for (int j = 0; j < n; j ++){
        local_anxn += prob.A(i, j) * sol(j);
      }
      #pragma omp atomic
      sol(i+n) += local_anxn;
    }
    // pro.stop(4);
  }

  // Phase-2 Dual
  while (true){
    #pragma omp parallel num_threads(core)
    {
      // Pricing
      // pro.clock(5);
      #pragma omp master
      {
        r = -1;
        double max_dse = -DBL_MAX;
        double delta = 0;

        // cout << "BASIS" << endl;
        // print(bhead);
        // cout << "X" << endl;
        // shortPrint(x);
        // cout << "VERIFY: " << verify2(x, A, bbl, bbu, c, l, u) << endl;
        // // cout << "BOUNDS" << endl;
        // // print(l);
        // // print(u);
        // cout << "SLACK BOUNDS" << endl;
        // print(bl);
        // print(bu);
        // cout << "BETA" << endl;
        // print(beta);

        for (int i = 0; i < m; i ++){
          if (bhead(i) < n){
            // n structural variables
            int index = bhead(i);
            delta = primalInfeasibilities(sol(bhead(i)), prob.l(index), prob.u(index));
          } else{
            // m slack variables
            int index = bhead(i) - n;
            delta = primalInfeasibilities(sol(bhead(i)), bl(index), bu(index));
          }
          if (!isEqual(delta, 0)){
            double dse = delta*delta/beta(i);
            if (max_dse < dse){
              max_dse = dse;
              max_delta = delta;
              r = i;
            }
          }
        }
        if (r == -1){
          // Optimal basis found
          status = Found;
        } else{
          sign_max_delta = (int) sign(max_delta);
          p = bhead(r);
          // BTran
          rho_r = Binv.row(r);
        }
      }
      // pro.stop(5);
      #pragma omp barrier
      if (status != Found){
        // pro.clock(6);
        // Pivot row
        #pragma omp for nowait
        for (int i = 0; i < n+m; i ++){
          double val = 0;
          if (i < n){
            for (int j = 0; j < m; j ++) 
              val -= prob.A(j, i) * rho_r(j);
          } else{
            val += rho_r(i-n);
          }
          alpha_r(i) = val * sign_max_delta;
        }
        // pro.stop(6);
        // pro.clock(7);
        // Initialize Ratio test
        #pragma omp master
        {
          init_size = 0;
        }
        #pragma omp barrier
        int local_init_size = 0;
        int local_init_start = 0;
        if (iteration_count == 0){
          vector<int> local_init ((n+m)/core+1);
          #pragma omp for nowait
          for (int i = 0; i < n+m; i ++){
            if (!inv_bhead[i]){
              // Non-basic
              if (i < n){
                if ((isGreater(alpha_r(i), 0, kE_Ieq) && isEqual(sol(i), prob.l(i), kE_Eq)) || (isLess(alpha_r(i), 0, kE_Ieq) && isEqual(sol(i), prob.u(i), kE_Eq))){
                  local_init[local_init_size] = i;
                  local_init_size ++;
                }
              } else{
                if ((isGreater(alpha_r(i), 0, kE_Ieq) && isEqual(sol(i), bl(i-n), kE_Eq)) || (isLess(alpha_r(i), 0, kE_Ieq) && isEqual(sol(i), bu(i-n), kE_Eq))){
                  local_init[local_init_size] = i;
                  local_init_size ++;
                }
              }
            }
          }
          #pragma omp critical (c1)
          {
            local_init_start = init_size;
            init_size += local_init_size;
          }
          // pro.stop(7);
          // pro.clock(8);
          #pragma omp barrier
          #pragma omp master
          {
            init = new pair<double, int>[init_size];
          }
          #pragma omp barrier
          for (int i = local_init_start; i < local_init_start + local_init_size; i ++){
            int j = local_init[i-local_init_start];
            init[i] = {d(j)/alpha_r(j), j};
          }
          // pro.stop(8);
        } else{
          double sum_slope = 0;
          pair<double, int>* local_init = new pair<double, int>[(n+m)/core+1];

          // cout << "MAX SLOPE: " << fabs(max_delta) << endl;

          #pragma omp for nowait
          for (int i = 0; i < n+m; i ++){
            if (!inv_bhead[i]){
              // Non-basic
              if (i < n){
                if ((isGreater(alpha_r(i), 0, kE_Ieq) && isEqual(sol(i), prob.l(i), kE_Eq)) || (isLess(alpha_r(i), 0, kE_Ieq) && isEqual(sol(i), prob.u(i), kE_Eq))){
                  pushHeap(prob.u, prob.l, bu, bl, alpha_r, d, local_init, local_init_size, sum_slope, fabs(max_delta), i);
                  // cout << "Adding " << i << " " << sum_slope << " " << local_init_size << endl;
                  // for (int i = 0; i < local_init_size; i ++){
                  //   cout << local_init[i].first << "," << local_init[i].second << "," << getSlope(u, l, bu, bl, alpha_r, i) << " ";
                  // }
                  // cout << endl;
                }
              } else{
                if ((isGreater(alpha_r(i), 0, kE_Ieq) && isEqual(sol(i), bl(i-n), kE_Eq)) || (isLess(alpha_r(i), 0, kE_Ieq) && isEqual(sol(i), bu(i-n), kE_Eq))){
                  pushHeap(prob.u, prob.l, bu, bl, alpha_r, d, local_init, local_init_size, sum_slope, fabs(max_delta), i);
                  // cout << "Adding " << i << " " << sum_slope << " " << local_init_size << endl;
                  // for (int i = 0; i < local_init_size; i ++){
                  //   cout << local_init[i].first << "," << local_init[i].second << "," << getSlope(u, l, bu, bl, alpha_r, i) << " ";
                  // }
                  // cout << endl;
                }
              }
            }
          }
          #pragma omp critical (c1)
          {
            local_init_start = init_size;
            init_size += local_init_size;
          }
          // pro.stop(7);
          // pro.clock(8);
          #pragma omp barrier
          #pragma omp master
          {
            init = new pair<double, int>[init_size];
          }
          #pragma omp barrier
          for (int i = local_init_start; i < local_init_start + local_init_size; i ++){
            init[i] = local_init[i - local_init_start];
          }
          delete[] local_init;
          // pro.stop(8);
        }
      } else{
        double local_score = 0;
        #pragma omp for nowait
        for (int i = 0; i < n; i ++) local_score += sol(i) * prob.c(i);
        #pragma omp atomic
        score += local_score;
      }
    }
    if (!init_size){
      status = Infeasible;
      break;
    }
    if (status != Found){
      // Ratio test using Algorithm 4 combining Algorithm 5
      // cout << "A " << iteration_count << " " << init_size << endl;
      // pro.clock(9, false);
      map_sort::Sort(init, init_size, core);
      // pro.stop(9, false);
      // cout << "A1" << endl;
      slope_partition.fill(init_size);
      slope_sums.resize(init_size);
      // double ss = 0;
      #pragma omp parallel num_threads(core)
      {
        // pro.clock(10);
        double sum_slope = 0;
        bool is_first = false;
        #pragma omp for
        for (int i = 0; i < init_size; i ++){
          if (!is_first){
            is_first = true;
            slope_partition(omp_get_thread_num()) = i;
          }
          sum_slope += getSlope(prob.u, prob.l, bu, bl, alpha_r, init[i].second);
          slope_sums(i) = sum_slope;
          // #pragma omp atomic
          // ss += getSlope(prob.u, prob.l, bu, bl, alpha_r, init[i].second);
        }
        #pragma omp master
        {
          // cout << "B " << iteration_count << " " << init_size << endl;
          sort(slope_partition.begin(), slope_partition.end());
          slope_sums_parition(0) = 0;
          for (int i = 1; i < core; i ++){
            slope_sums_parition(i) = slope_sums_parition(i-1) + slope_sums(slope_partition(i)-1);
          }
        }
        #pragma omp barrier
        int thread_num = omp_get_thread_num();
        for (int i = slope_partition(thread_num); i < slope_partition(thread_num+1); i ++){
          slope_sums(i) += slope_sums_parition(thread_num);
        }
        #pragma omp barrier
        // pro.stop(10);
        // pro.clock(11);
        #pragma omp master
        {
          // cout << "ALPHA_R" << endl;
          // print(alpha_r);
          // for (int i = 0; i < init_size; i ++) cout << init[i].first << " " << init[i].second << endl;
          // cout << "MAX_SLOPE: " << fabs(max_delta) << endl;
          // print(slope_sums);
          //cout << "C " << iteration_count << " " << init_size << endl;
          if (isLess(slope_sums(init_size-1), fabs(max_delta), kE_ap)){
            // cout << ss << endl;
            // cout << fabs(max_delta) << " " << slope_sums(init_size-1) << endl;
            status = DualUnbounded;
          } else{
            q_index = (int) (upper_bound(slope_sums.begin(), slope_sums.end(), fabs(max_delta)) - slope_sums.begin());
            // cout << q_index << " " << init_size << " " << slope_sums.size() << endl;
            q_index = min(q_index, init_size-1);
            q = init[q_index].second;

            // cout << "PIVOT " << fabs(alpha_r(q)) << endl;
            // cout << "SCORE" << endl;
            // for (int i = 0; i < init_size; i ++) cout << init[i].first << " " << init[i].second << endl;
            // cout << "PRE_MAX_DELTA " << max_delta << endl;

            if (q_index > 0){
              mini_iteration_count += q_index - 1;
              max_delta = sign_max_delta * (fabs(max_delta) - slope_sums(q_index-1));
            }
            // cout << alpha_r(q) << endl;
            dual_step = d(q) / alpha_r(q) * sign_max_delta;
            d(p) = -dual_step;
            tilde_a.fill(0);

            // cout << "AFT_MAX_DELTA " << max_delta << endl;
            // cout << "Q" << endl;
            // cout << q_index << " " << q << endl;
            // cout << "SLOPE SUMS" << endl;
            // print(slope_sums);

            // FTran
            if (q < n) alpha_q = -Binv * prob.A.col(q);
            else alpha_q = Binv.col(q-n);

            // DSE FTran
            tau = Binv * rho_r;
          }
          //cout << "C1 " << iteration_count << " " << init_size << endl;
        }
        // pro.stop(11);
        #pragma omp barrier
        if (status != DualUnbounded){
          // pro.clock(12);
          // Update d
          #pragma omp for nowait
          for (int i = 0; i < n+m; i ++){
            if (!inv_bhead[i]){
              // Non-basic
              d(i) -= dual_step * alpha_r(i) * sign_max_delta;
            }
          }
          // pro.stop(12);
          // pro.clock(13);
          // Update tilde_a
          VectorXd local_tilde_a (m); local_tilde_a.fill(0);
          #pragma omp for nowait
          for (int i = 0; i < q_index; i ++){
            int j = init[i].second;
            if (j < n){
              if (isEqual(sol(j), prob.l(j), kE_Eq)){
                local_tilde_a += (prob.l(j) - prob.u(j)) * prob.A.col(j);
                sol(j) = prob.u(j);
              } else{
                local_tilde_a += (prob.u(j) - prob.l(j)) * prob.A.col(j);
                sol(j) = prob.l(j);
              }
            } else{
              int index = j - n;
              if (isEqual(sol(j), bl(index), kE_Eq)){
                local_tilde_a(index) += (bu(index) - bl(index));
                sol(j) = bu(index);
              } else{
                local_tilde_a(index) += (bl(index) - bu(index));
                sol(j) = bl(index);
              }
            }
          }
          for (int i = 0; i < m; i ++){
            #pragma omp atomic
            tilde_a(i) += local_tilde_a(i);
          }
          #pragma omp barrier
          // pro.stop(13);
          // pro.clock(14);
          #pragma omp master
          {
            //cout << "D " << iteration_count << " " << init_size << endl;
            // Update xb based on tilde_a
            VectorXd delta_xB = Binv * tilde_a;
            for (int i = 0; i < m; i ++){
              sol(bhead(i)) -= delta_xB(i);
            }

            // cout << "BINV" << endl;
            // cout << Binv << endl;
            // // cout << "TILDE_A" << endl;
            // // print(tilde_a);
            // cout << "AFT TILDE_A" << endl;
            // shortPrint(x);

            // Update primal solution
            primal_step = max_delta / alpha_q(r);
            for (int i = 0; i < m; i ++) sol(bhead(i)) -= primal_step * alpha_q(i);
            sol(q) += primal_step;

            // cout << "AFT PRIMAL" << endl;
            // shortPrint(x);

            // Update edge norm beta
            for (int i = 0; i < m; i ++){
              if (i != r){
                double ratio = alpha_q(i) / alpha_q(r);
                beta(i) += ratio * (ratio * beta(r) - 2 * tau(i));
              }
            }
            beta(r) /= (alpha_q(r) * alpha_q(r));

            // Update bhead and inv_bhead
            inv_bhead[p].flip();
            bhead(r) = q;
            inv_bhead[q].flip();

            // for (int i = 0; i < n+m; i ++){
            //   if (!inv_bhead[i]){
            //     if (i < n){
            //       bool t = (isEqual(x(i),l(i)) && isGreaterEqual(d(i), 0)) || (isEqual(x(i),u(i)) && isLessEqual(d(i), 0));
            //       if (!t){
            //         fmt::print("{} {:.10Lf} {:.10Lf} {:.10Lf} {:.10Lf}\n", i, x(i), l(i), u(i), d(i));
            //         assert(0==1);
            //       }
            //     } else{
            //       int index = i - n;
            //       bool t = (isEqual(x(i),bl(index)) && isGreaterEqual(d(i), 0)) || (isEqual(x(i),bu(index)) && isLessEqual(d(i), 0));
            //       if (!t){
            //         fmt::print("{} {:.10Lf} {:.10Lf} {:.10Lf} {:.10Lf}\n", i, x(i), bl(index), bu(index), d(i));
            //         assert(0==1);
            //       }
            //     }
            //   }
            // }

            // Update Binv
            for (int i = 0; i < m; i ++){
              if (i != r){
                double ratio = alpha_q(i) / alpha_q(r);
                Binv.row(i) -= ratio * Binv.row(r);
              }
            }
            Binv.row(r) /= alpha_q(r);

            // Iteration count for pivot
            iteration_count ++;
          }
          // pro.stop(14);
        }
      }
      if (status == DualUnbounded) break;
    } else break;
    delete[] init; init = nullptr;
  }
  if (init) delete[] init;

  bool require_gurobi = false;
  if (status == Found){
    Checker ch = Checker(prob);
    VectorXd ch_sol = sol;
    ch_sol.conservativeResize(n);
    if (ch.checkLpFeasibility(ch_sol) != Feasibility){
      require_gurobi = true;
    }
  }
  if (require_gurobi){
    // Numerical issue is huge
    GurobiSolver gs = GurobiSolver(prob);
    gs.solveLp();
    status = gs.lp_status;
    if (status == Found){
      sol = gs.lp_sol;
      score = gs.lp_score;
    }
  }
  // pro.print();
  exe_solve = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000000.0;
}