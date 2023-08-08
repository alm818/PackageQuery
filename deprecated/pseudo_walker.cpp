#include <iostream>


#include "pseudo_walker.h"
#include "utility.h"
#include "omp.h"

using namespace std;

PseudoWalker::~PseudoWalker(){
  if (ppq) delete ppq;
}

PseudoWalker::PseudoWalker(VectorXd& p, bool enable_correction, int core){
  ppq = NULL; // Very important to avoid deleting uninitialized non-null ptr
  n = (int) p.size();
  assert(n > 1);
  steps.resize(n);
  x.resize(n); x.fill(0);
  this->enable_correction = enable_correction;
  this->core = core;
  step_count = 0;
  correction_count = 0;
  sum_p = 0;
  if (core > 1){
    this->p.resize(n);
    vector<pair<double, int>>* init = new vector<pair<double, int>>(n);
    double p_norm = 0;
    double S = 0;
    VectorXd pp = p / p.norm();
    // cout << "THEIR NORM:" << p.norm() << endl;
    // cout << "THEIR S:" << pp.cwiseAbs().sum() << endl; 
    #pragma omp parallel num_threads(core)
    {
      double local_p_norm = 0;
      double local_S = 0;
      #pragma omp for nowait
      for (int i = 0; i < n; i ++){
        local_p_norm += p(i)*p(i);
        local_S += abs(p(i));
      }
      #pragma omp atomic
      p_norm += local_p_norm;
      #pragma omp atomic
      S += local_S;
      #pragma omp barrier
      #pragma omp master
      {
        p_norm = sqrt(p_norm);
        S /= p_norm;
      }
      #pragma omp barrier
      #pragma omp for nowait
      for (int i = 0; i < n; i ++){
        this->p(i) = p(i) / p_norm;
        double abs_pi = abs(this->p(i));
        steps(i) = 2 - 2 * abs_pi * (abs_pi - (S - abs_pi) / (n - 1.0));
        (*init)[i] = { this->p(i) * this->p(i), i };
      }
    }
    // cout << "MY S: " << S << endl;
    // cout << "MY P_NORM: " << p_norm << endl;
    ppq = new ParallelPQ(core, init);
  } else{
    this->p = p / p.norm();
    double S = this->p.cwiseAbs().sum();
    vector<pair<double, int>> init(n);
    for (int i = 0; i < n; i++) {
      double abs_pi = abs(this->p(i));
      steps(i) = 2 - 2 * abs_pi * (abs_pi - (S - abs_pi) / (n - 1.0));
      init[i] = { this->p(i) * this->p(i), i };
    }
    pq = priority_queue(init.begin(), init.end());
  }
}

int PseudoWalker::executeStep(int i){
  step_count ++;
  double s = sign(p(i));
  x(i) += s;
  sum_p += abs(p(i));
  if (enable_correction && step_count >= p.size() / log2(p.size()) * (correction_count+1)){
    correction_count ++;
    if (core > 1){
      vector<pair<double, int>>* init = new vector<pair<double, int>>(n);
      #pragma omp parallel for num_threads(core)
      for (int i = 0; i < n; i++){
        double abs_pi = abs(this->p(i));
        (*init)[i] = { abs_pi * (abs_pi + 2*sum_p) - 2*abs(x(i)), i };
      }
      delete ppq;
      ppq = new ParallelPQ(core, init);
    } else{
      vector<pair<double, int>> init (n);
      for (int i = 0; i < n; i++) {
        double abs_pi = abs(this->p(i));
        init[i] = { abs_pi * (abs_pi + 2*sum_p) - 2*abs(x(i)), i };
      }
      pq = priority_queue(init.begin(), init.end());
    }
  }
  return s*(i+1);
}

int PseudoWalker::step() {
  if (core > 1){
    const auto& pi = ppq->peak();
    // cout << "SUBTRACING STEP AT " << pi.second << " " << steps(pi.second) << endl;
    // double abs_pi = abs(p(pi.second));
    // cout << p(pi.second) << " " << abs_pi << " " << n << endl;
    ppq->subtractPeak(steps(pi.second));
    return executeStep(pi.second);
  } else{
    const auto& pi = pq.top();
    int i = pi.second;
    double new_pv = pi.first - steps(i);
    pq.pop();
    pq.emplace(new_pv, i);
    return executeStep(i);
  }
}