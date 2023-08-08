#include <iostream>
#include <random>

#include "fast_cv.h"
#include "utility.h"
#include  "omp.h"

FastCV::~FastCV(){
}

FastCV::FastCV(int core, CooSparse& mat, const vector<int>& nbasic): mat(mat), nbasic(nbasic){
  m = mat.m;
  n = mat.n;
  if (n >= kFastCV){
    c = (int)ceil(n/(log2(n)+m));
    col_norms.resize(n);
    #pragma omp parallel for num_threads(core)
    for (int i = 0; i < n; i ++){
      double norm = 0;
      for (const auto& pi : mat.rows[i]) norm += pi.second*pi.second;
      col_norms(i) = sqrt(norm);
    }
  }
}

// O(7n+clogn+cm) if non-blind probability
// By the choice of c, this is O(8n)
VectorXd FastCV::sample(){
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  VectorXd coefs (n);
  coefs(0) = 1;
  for (int i = 1; i < n; i ++) coefs(i) = u_dist(gen);
  sort(coefs.begin(), coefs.end());
  for (int i = n-1; i >= 1; i--) coefs(i) -= coefs(i-1);
  VectorXd res (n); res.fill(0);
  if (n >= kFastCV){
    VectorXd p (n);
    for (int i = 0; i < n; i ++){
      if (i == 0) p(i) = col_norms(i) * coefs(i);
      else p(i) = p(i-1) + col_norms(i) * coefs(i);
    }
    for (int i = 0; i < c; i ++){
      double sampled_value = u_dist(gen) * p(n-1);
      int left = 0; int right = n-1;
      while (left < right){
        int mid = (left+right)/2;
        if (p(mid) < sampled_value) left = mid + 1;
        else right = mid;
      }
      for (const auto& pi : mat.rows[left]){
        res(pi.first) += coefs(left) * pi.second;
      }
    }
  } else {
    mat.inplaceVectorProduct(coefs, res);
  }
  for (int i = 0; i < n; i ++){
    int comp = abs(nbasic[i]) - 1;
    if (nbasic[i] != 0) res(comp) = sign(nbasic[i]) * coefs(i);
  }
  // VectorXd base (m); base.fill(0);
  // for (int i = 0; i < m; i ++){
  //   for (int j = 0; j < n; j ++){
  //     base(i) += mat(i, j) * coefs(j);
  //   }
  // }
  // print(base);
  // print(tmp);
  return res;
}
