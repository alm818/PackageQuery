#include "pb/util/udeclare.h"
#include "pb/util/uconfig.h"
#include "pb/util/udebug.h"
#include "pb/det/det_prob.h"
#include "pb/core/dual_reducer.h"
#include "pb/core/gurobi_solver.h"
#include "pb/det/det_bound.h"
#include "pb/det/dlv.h"
#include "pb/det/synthetic.h"
#include "pb/det/lsr_prob.h"
#include "pb/util/upostgres.h"
#include "pb/det/lsr.h"
#include "pb/det/det_sql.h"
#include "pb/core/dist.h"
#include "pb/det/det_exp.h"

#include "pb/lib/praxis.hpp"
#include "pb/lib/normal.hpp"
#include "pb/lib/log_normal_truncated_ab.hpp"
#include "pb/lib/log_normal.hpp"
#include "pb/lib/truncated_normal.hpp"
#include "pb/lib/toms178.hpp"
#include "pb/lib/brent.hpp"
#include "pb/core/vg.h"

#include "pb/lib/random_quantile.h"
#include "pb/lib/common.h"

using namespace pb;

class MedianFinder {
public:
    // Adds a new number to the data structure
    void addNum(double num) {
        // If max heap is empty or the new number is less than or equal to the root of max heap
        if (maxHeap.empty() || num <= maxHeap.top()) {
            maxHeap.push(num);
        } else {
            minHeap.push(num);
        }

        // Balance the heaps
        if (maxHeap.size() > minHeap.size() + 1) {
            minHeap.push(maxHeap.top());
            maxHeap.pop();
        } else if (minHeap.size() > maxHeap.size()) {
            maxHeap.push(minHeap.top());
            minHeap.pop();
        }
    }

    // Returns the median of the current data stream
    double findMedian() {
        if (maxHeap.size() == minHeap.size()) {
            return (maxHeap.top() + minHeap.top()) / 2.0;
        } else {
            return maxHeap.top();
        }
    }

private:
    // Max heap to store the lower half of the numbers
    priority_queue<double> maxHeap;

    // Min heap to store the upper half of the numbers
    priority_queue<double, vector<double>, std::greater<double>> minHeap;
};

void dlv_binary_test(){
  // double left = 0;
  // double dif = samples[N-1] - samples[0];
  // double right = 0.25*dif*dif; // Popoviciu inequality
  // while (abs(left-right) > 1e-8){
  //   double mid = (left+right)/2;
  //   int cnt = 1;
  //   ScalarMeanVar smv = ScalarMeanVar();
  //   for (int i = 0; i < N; i ++){
  //     smv.add(samples[i]);
  //     if (smv.getBiasedVar() > mid){
  //       smv.reset();
  //       smv.add(samples[i]);
  //       cnt ++;
  //     }
  //   }
  //   // cout << cnt << " " << g << " " << left << " " << mid << " " << right << endl;
  //   if (cnt < g) right = mid;
  //   else left = mid;
  // }
  // ScalarMeanVar smv_opt_var = ScalarMeanVar();
  // double opt_var = (left+right) / 2;
  // vector<int> opt_sizes (g, 0);
  // int index = 0;
  // ScalarMeanVar smv = ScalarMeanVar();
  // for (int i = 0; i < N; i ++){
  //   ScalarMeanVar tmp_smv = smv;
  //   smv.add(samples[i]);
  //   if (smv.getBiasedVar() > opt_var){
  //     smv_opt_var.add(tmp_smv.getBiasedVar());
  //     smv.reset();
  //     smv.add(samples[i]);
  //     index ++;
  //   }
  //   opt_sizes[index] ++;
  // }
  // smv_opt_var.add(smv.getBiasedVar());
}

void dlv_graph_test(){
  DetExp exp = DetExp("DLV2");
  int N = 100000;
  VectorXd samples;
  double mean = 0;
  double var = 1;
  double a = mean - sqrt(3*var);
  double b = mean + sqrt(3*var);
  UNUSED(a);
  UNUSED(b);
  Normal dist = Normal(mean, var);
  // Uniform dist = Uniform(a, b);
  dist.sample(samples, N);
  sort(samples.begin(), samples.end());
  ScalarMeanVar whole;
  for (int i = 0; i < (int) samples.size(); i ++){
    whole.add(samples[i]);
  }
  double sigmas = whole.getBiasedVar();
  double upper = 0.25*(whole.getMax() - whole.getMin())*(whole.getMax() - whole.getMin());
  double lower = 0;
  int cur_count = N;

  while (cur_count > 1){
    double left = lower;
    double right = upper;
    double my_cnt, my_vRAC, my_vRWC;
    while (abs(right-left)>1e-6){
      double mid = (left+right)/2;
      long long cnt = 1;
      ScalarMeanVar smv = ScalarMeanVar();
      ScalarMeanVar smv_var = ScalarMeanVar();
      for (int i = 0; i < N; i ++){
        ScalarMeanVar pre_smv = smv;
        smv.add(samples[i]);
        if (smv.getBiasedVar() > mid){
          smv_var.add(pre_smv.getBiasedVar());
          smv.reset();
          smv.add(samples[i]);
          cnt ++;
        }
      }

      if (cnt < cur_count){
        right = mid;
        smv_var.add(smv.getBiasedVar());
        my_vRAC = smv_var.getMax()/sigmas;
        my_vRWC = cnt*cnt*smv_var.getMax()/sigmas;
        my_cnt = cnt;
      } else left = mid;
    }
    double mid = (left+right)/2;
    exp.write("p", mid, my_cnt);
    exp.write("rac", mid, my_vRAC);
    exp.write("rwc", mid, my_vRWC);
    cur_count = my_cnt;
  }
}

VectorXd samples;
double sigmas;
double pmax, pmin;

double fv (double x[], int nvars){
  UNUSED(nvars);
  double beta = x[0];
  int cnt = 1;
  ScalarMeanVar smv = ScalarMeanVar();
  ScalarMeanVar smv_var = ScalarMeanVar();
  for (int i = 0; i < (int) samples.size(); i ++){
    ScalarMeanVar pre_smv = smv;
    smv.add(samples[i]);
    if (smv.getBiasedVar() > beta){
      smv_var.add(pre_smv.getBiasedVar());
      smv.reset();
      smv.add(samples[i]);
      cnt ++;
    }
  }
  if (cnt > pmax) return 24;
  if (cnt < pmin) return 24;
  smv_var.add(smv.getBiasedVar());
  double vRWC = cnt*cnt*smv_var.getMax()/sigmas;
  return vRWC;
}

double f (double beta){
  int cnt = 1;
  ScalarMeanVar smv = ScalarMeanVar();
  ScalarMeanVar smv_var = ScalarMeanVar();
  for (int i = 0; i < (int) samples.size(); i ++){
    ScalarMeanVar pre_smv = smv;
    smv.add(samples[i]);
    if (smv.getBiasedVar() > beta){
      smv_var.add(pre_smv.getBiasedVar());
      smv.reset();
      smv.add(samples[i]);
      cnt ++;
    }
  }
  if (cnt > pmax) return 24;
  if (cnt < pmin) return 24;
  smv_var.add(smv.getBiasedVar());
  cout << cnt << " " << smv_var.getMax() << endl;
  double vRWC = cnt*cnt*smv_var.getMax()/sigmas;
  return vRWC;
}

void dlv_iterative_test(){
  DetExp exp = DetExp("DLV3");
  int N = 100000;
  int L = 62914;
  int extra_cnt = 200;
  VectorXd samples;
  double mean = 0;
  double var = 1;
  double a = mean - sqrt(3*var);
  double b = mean + sqrt(3*var);
  UNUSED(a);
  UNUSED(b);
  Normal dist = Normal(mean, var);
  // Uniform dist = Uniform(a, b);
  dist.sample(samples, N);
  sort(samples.begin(), samples.end());
  ScalarMeanVar whole;
  for (int i = 0; i < (int) samples.size(); i ++){
    whole.add(samples[i]);
  }
  double sigmas = whole.getBiasedVar();
  double upper = 0.25*(whole.getMax() - whole.getMin())*(whole.getMax() - whole.getMin());
  double lower = 0;
  int cur_count = N;
  
  int fi_ind = 0;
  double max_count = 0;
  bool is_done = false;

  while (true){
    double left = lower;
    double right = upper;
    double base_cnt = 1;
    double base_RAC = 0;
    if (!is_done){
      while (abs(right-left)>1e-6){
        double mid = (left+right)/2;
        long long cnt = 1;
        ScalarMeanVar smv = ScalarMeanVar();
        ScalarMeanVar smv_var = ScalarMeanVar();
        for (int i = 0; i < N; i ++){
          ScalarMeanVar pre_smv = smv;
          smv.add(samples[i]);
          if (smv.getBiasedVar() > mid){
            smv_var.add(pre_smv.getBiasedVar());
            smv.reset();
            smv.add(samples[i]);
            cnt ++;
          }
        }

        if (cnt < cur_count){
          right = mid;
          smv_var.add(smv.getBiasedVar());
          base_RAC = cnt*cnt*smv_var.getMean()/sigmas;
          base_cnt = cnt;
        } else left = mid;
      }
      max_count = max(max_count, base_cnt);
    }
    if (base_cnt == 1){
      is_done = true;
      fi_ind ++;
      if (fi_ind > extra_cnt) break;
      base_cnt = max_count + fi_ind * (L - max_count) / extra_cnt;
    }

    priority_queue<pair<double, pair<int, int>>> pq;
    pq.push({0, {0, N-1}});

    double g_ratio = 0.1;

    // bool is_total_kd = false;
    bool is_total_iterative = true;

    while (pq.size() < base_cnt){
      auto p = pq.top().second;
      pq.pop();
      int l = p.first;
      int r = p.second;
      int n = r-l+1;
      ScalarMeanVar smv = ScalarMeanVar();
      for (int i = l; i <= r; i ++){
        smv.add(samples(i));
      }
      double beta = smv.getBiasedVar() * g_ratio * g_ratio * 13.5;
      smv.reset();
      int start = l;
      for (int i = l; i <= r; i ++){
        double total_var;
        if (is_total_iterative){
          total_var = smv.getM2();
        } else{
          total_var = smv.getBiasedVar();
        }
        smv.add(samples(i));
        if (smv.getBiasedVar() > beta){
          pq.push({total_var, {start, i-1}});
          start = i;
          smv.reset();
          smv.add(samples(i));
        }
      }
      if (is_total_iterative) pq.push({smv.getM2(), {start, r}});
      else pq.push({smv.getBiasedVar(), {start, r}});
    }
    int this_cnt = pq.size();
    ScalarMeanVar smv_var = ScalarMeanVar();
    while (pq.size() > 0){
      auto p = pq.top().second;
      int n = p.second - p.first + 1;
      if (is_total_iterative){
        smv_var.add(pq.top().first/n);
      } else{
        smv_var.add(pq.top().first);
      }
      pq.pop();
    }
    double iterative_RAC = this_cnt*smv_var.getMean()/sigmas*this_cnt;

    // kd unsort
    pq.push({0, {0, N-1}});
    while (pq.size() < base_cnt){
      double ind = pq.top().first;
      auto p = pq.top().second;
      pq.pop();
      int l = p.first;
      int r = p.second;
      int n = r-l+1;
      ScalarMeanVar smv = ScalarMeanVar();
      for (int i = l; i <= r; i ++){
        smv.add(samples(i));
      }
      double mean = smv.getMean();
      for (int i = l; i <= r; i ++){
        if (samples(i) > mean){
          // cout << mean << " " << l << " " << i-0.5 << " " << r << endl;
          pq.push({ind-1, {l, i-1}});
          pq.push({ind-1, {i, r}});
          break;
        }
      }
    }
    this_cnt = pq.size();
    smv_var.reset();
    while (pq.size() > 0){
      auto p = pq.top().second;
      int n = p.second - p.first + 1;
      ScalarMeanVar smv = ScalarMeanVar();
      for (int i = p.first; i <= p.second; i ++){
        smv.add(samples(i));
      }
      smv_var.add(smv.getBiasedVar());
      pq.pop();
    }
    double kd_unsort_RAC = this_cnt*smv_var.getMean()/sigmas*this_cnt;

    // kd sort
    pq.push({0, {0, N-1}});
    while (pq.size() < base_cnt){
      double ind = pq.top().first;
      auto p = pq.top().second;
      pq.pop();
      int l = p.first;
      int r = p.second;
      int n = r-l+1;

      ScalarMeanVar smv = ScalarMeanVar();
      for (int i = l; i <= r; i ++){
        smv.add(samples(i));
      }
      double mean = smv.getMean();
      smv.reset();
      for (int i = l; i <= r; i ++){
        if (samples(i) > mean){
          // cout << mean << " " << l << " " << i-0.5 << " " << r << endl;
          pq.push({smv.getM2(), {l, i-1}});
          smv.reset();
          for (int j = i; j <= r; j ++){
            smv.add(samples(j));
          }
          pq.push({smv.getM2(), {i, r}});
          break;
        }
        smv.add(samples(i));
      }
    }
    this_cnt = pq.size();
    smv_var.reset();
    while (pq.size() > 0){
      auto p = pq.top().second;
      int n = p.second - p.first + 1;
      ScalarMeanVar smv = ScalarMeanVar();
      for (int i = p.first; i <= p.second; i ++){
        smv.add(samples(i));
      }
      smv_var.add(smv.getBiasedVar());
      pq.pop();
    }
    double kd_sort_RAC = this_cnt*smv_var.getMean()/sigmas*this_cnt;

    double mid = (left+right)/2;
    exp.write("p", mid, base_cnt);
    exp.write("binary_rac", mid, base_RAC);
    exp.write("iterative_rac", mid, iterative_RAC);
    exp.write("kd_unsort_rac", mid, kd_unsort_RAC);
    exp.write("kd_sort_rac", mid, kd_sort_RAC);
    cur_count = base_cnt;
  }
}

void dlv_scale_test(){
  DetExp exp = DetExp("DLV4");
  int N = 100000;
  int L = (int)(N*0.1);
  VectorXd samples;
  double mean = 0;
  double var = 1;
  double a = mean - sqrt(3*var);
  double b = mean + sqrt(3*var);
  UNUSED(a);
  UNUSED(b);
  Normal dist = Normal(mean, var);
  // Uniform dist = Uniform(a, b);
  dist.sample(samples, N);
  sort(samples.begin(), samples.end());
  ScalarMeanVar whole;
  for (int i = 0; i < (int) samples.size(); i ++){
    whole.add(samples[i]);
  }
  double sigmas = whole.getBiasedVar();
  double upper = 0.25*(whole.getMax() - whole.getMin())*(whole.getMax() - whole.getMin());
  double lower = 0;
  int cur_count = N;
  double left = lower;
  double right = upper;
  double base_cnt = 1;
  double base_RAC = 0;
  while (abs(right-left)>1e-6){
    double mid = (left+right)/2;
    long long cnt = 1;
    ScalarMeanVar smv = ScalarMeanVar();
    ScalarMeanVar smv_var = ScalarMeanVar();
    for (int i = 0; i < N; i ++){
      ScalarMeanVar pre_smv = smv;
      smv.add(samples[i]);
      if (smv.getBiasedVar() > mid){
        smv_var.add(pre_smv.getBiasedVar());
        smv.reset();
        smv.add(samples[i]);
        cnt ++;
      }
    }

    if (cnt < cur_count){
      right = mid;
      smv_var.add(smv.getBiasedVar());
      base_RAC = cnt*cnt*smv_var.getMean()/sigmas;
      base_cnt = cnt;
    } else left = mid;
  }
  
  int extra_cnt = 20;
  int fi_ind = 0;
  double max_count = base_cnt;
  double start_x = 1.0;
  double end_x = 90;
  int x_count = 200;
  for (int i = 0; i < x_count; i ++){
    double x = start_x + (end_x-start_x)/x_count*i;
    for (int j = 1; j <= extra_cnt; j ++){
      double base_cnt = max_count + (L - max_count) / extra_cnt * j;
      base_cnt = N*0.1;
      priority_queue<pair<double, pair<int, int>>> pq;
      pq.push({0, {0, N-1}});
      double g_ratio = 0.1;
      // bool is_total_kd = false;
      bool is_total_iterative = true;
      double compute = 0;
      while (pq.size() < base_cnt){
        auto p = pq.top().second;
        pq.pop();
        int l = p.first;
        int r = p.second;
        int n = r-l+1;
        ScalarMeanVar smv = ScalarMeanVar();
        for (int i = l; i <= r; i ++){
          smv.add(samples(i));
        }
        compute += n;
        double beta = smv.getBiasedVar() * g_ratio * g_ratio * x;
        smv.reset();
        int start = l;
        for (int i = l; i <= r; i ++){
          double total_var;
          if (is_total_iterative){
            total_var = smv.getM2();
          } else{
            total_var = smv.getBiasedVar();
          }
          smv.add(samples(i));
          if (smv.getBiasedVar() > beta){
            pq.push({total_var, {start, i-1}});
            start = i;
            smv.reset();
            smv.add(samples(i));
          }
        }
        if (is_total_iterative) pq.push({smv.getM2(), {start, r}});
        else pq.push({smv.getBiasedVar(), {start, r}});
      }
      int this_cnt = pq.size();
      ScalarMeanVar smv_var = ScalarMeanVar();
      while (pq.size() > 0){
        auto p = pq.top().second;
        int n = p.second - p.first + 1;
        if (is_total_iterative){
          smv_var.add(pq.top().first/n);
        } else{
          smv_var.add(pq.top().first);
        }
        pq.pop();
      }
      double iterative_RAC = this_cnt*this_cnt*smv_var.getMean()/sigmas;
      exp.write("rac", x, iterative_RAC);
      exp.write("compute", x, compute/N);
      break;
    }
  }
  base_cnt = N*0.1;
  priority_queue<pair<double, pair<int, int>>> pq;
  pq.push({0, {0, N-1}});
  double kd_compute = 0;
  while (pq.size() < base_cnt){
    double ind = pq.top().first;
    auto p = pq.top().second;
    pq.pop();
    int l = p.first;
    int r = p.second;
    int n = r-l+1;
    kd_compute += n;
    ScalarMeanVar smv = ScalarMeanVar();
    for (int i = l; i <= r; i ++){
      smv.add(samples(i));
    }
    double mean = smv.getMean();
    smv.reset();
    for (int i = l; i <= r; i ++){
      if (samples(i) > mean){
        // cout << mean << " " << l << " " << i-0.5 << " " << r << endl;
        pq.push({smv.getM2(), {l, i-1}});
        smv.reset();
        for (int j = i; j <= r; j ++){
          smv.add(samples(j));
        }
        pq.push({smv.getM2(), {i, r}});
        break;
      }
      smv.add(samples(i));
    }
  }
  double this_cnt = pq.size();
  ScalarMeanVar smv_var = ScalarMeanVar();
  smv_var.reset();
  while (pq.size() > 0){
    auto p = pq.top().second;
    int n = p.second - p.first + 1;
    ScalarMeanVar smv = ScalarMeanVar();
    for (int i = p.first; i <= p.second; i ++){
      smv.add(samples(i));
    }
    smv_var.add(smv.getBiasedVar());
    pq.pop();
  }
  double kd_sort_RAC = this_cnt*this_cnt*smv_var.getMean()/sigmas;
  cout << kd_sort_RAC << " " << kd_compute/N << endl;
}

void dlv_opt(){
  DetExp exp = DetExp("DLV_OPT");
  int N = 1000000;
  double group_ratio = 0.001;
  pmax = (int) ceil(N*group_ratio);
  pmin = (int) ceil(pmax / 13.5);
  double mean = 0;
  double var = 10000;
  double a = mean - sqrt(3*var);
  double b = mean + sqrt(3*var);
  UNUSED(a);
  UNUSED(b);
  Normal dist = Normal(mean, var);
  // Uniform dist = Uniform(a, b);
  dist.sample(samples, N);
  sort(samples.begin(), samples.end());
  ScalarMeanVar whole;
  for (int i = 0; i < (int) samples.size(); i ++){
    whole.add(samples[i]);
  }
  sigmas = whole.getBiasedVar();

  // double ibeta[1] = {1};
  // double ebeta[1];
  // double rho = 0.5;
  // double eps = 1e-8;
  // double itermax = 1000;
  // int iter = hooke(1, ibeta, ebeta, rho, eps, itermax, f);
  // cout << ebeta[0] << " " << iter << endl;

  double x1, x2;
  double dif = samples[N-1] - samples[0];
  double left = 0;
  double eps = 1e-5;
  double right = 0.25*dif*dif; // Popoviciu inequality
  // while (abs(left-right) > 1e-8){
  //   double mid1 = (2*left+right)/3;
  //   double mid2 = (left+2*right)/3;
  //   double fx1 = brent::local_min(0, mid1, eps, f, x1);
  //   double fx2 = brent::local_min(0, mid2, eps, f, x2);
  //   cout << left << " " << mid1 << " " << mid2 << " " << right << endl;
  //   cout << fx1 << " " << fx2 << endl;
  //   if (fx1 <= fx2) right = mid2;
  //   else left = mid1;
  // }
  // double mid = (left+right)/2;
  double x;
  double fx = brent::local_min(0, sigmas/(pmin*pmin), 1e-12, f, x);
  cout << x << " " << fx << endl;
}

int main(){
  // dlv_graph_test();
  dlv_iterative_test();
  // dlv_scale_test();
  // dlv_opt();
}