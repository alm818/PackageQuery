#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <bitset>
#include <random>

#include "utility.h"
#include "fmt/core.h"
#include "fast_cv.h"

using namespace std;
using namespace Eigen;

int main(){
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  normal_distribution n_dist(0.0, 1.0);
  int m = 6;
  int N = 100;
  MatrixXd A (m, N);
  for (int i = 0; i < m; i ++){
    for (int j = 0; j < N; j ++){
      A(i, j) = n_dist(gen);
    }
  }
  FastCV fcv = FastCV(40, A);
  VectorXd res = fcv.sample();
  print(res);
}