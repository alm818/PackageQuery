#pragma once

#include <Eigen/Dense>
#include "coo_sparse.h"

using namespace Eigen;
using namespace std;

class FastCV{
public:
  int n, m, c;
  CooSparse& mat;
  const VectorXd& nbasic;
  VectorXd col_norms;
public:
  ~FastCV();
  FastCV(int core, CooSparse& mat, const VectorXd& nbasic);
  VectorXd sample();
};