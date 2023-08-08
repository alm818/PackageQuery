#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <queue>
#include <random>
#include <chrono>
#include <numeric>

#include "utility.h"
#include "parallel_pq.h"
#include "fmt/core.h"
#include "fitsio.h"
#include "omp.h"

using namespace std;
using namespace Eigen;

int main(){
  int n = 10000000;
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  vector<pair<double, int>>* arr = new vector<pair<double, int>>(n);
  vector<pair<double, int>> arr2 (n);
  for (int i = 0; i < n; i ++){
    double  val = u_dist(gen);
    (*arr)[i] = {val, i};
    arr2[i] = {val, i};
  }
  int T = 50000;
  vector<int> a;
  vector<int> b;
  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
  ParallelPQ ppq = ParallelPQ(16, arr);
  for (int i = 0; i < T; i ++){
    auto pi = ppq.peak();
    a.push_back(pi.second);
    ppq.pop();
    //ppq.subtractPeak(0.01);
  }
  chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();
  double exe1 = chrono::duration_cast<chrono::nanoseconds>(end - start).count() / 1000000.0;
  
  // Beat by 10 times
  start = chrono::high_resolution_clock::now();
  priority_queue<pair<double, int>> pq = priority_queue(arr2.begin(), arr2.end());
  for (int i = 0; i < T; i ++){
    auto pi = pq.top();
    b.push_back(pi.second);
    pq.pop();
    //pq.push({pi.first-0.1, pi.second});
  }
  end = chrono::high_resolution_clock::now();
  double exe2 = chrono::duration_cast<chrono::nanoseconds>(end - start).count() / 1000000.0;
  fmt::print("ppq:{}ms pq:{}ms\n", exe1, exe2);
  for (int i = 0; i < T; i ++){
    if (a[i] != b[i]){
      cout << "NOT OKAY" << endl;
    }
  }
}