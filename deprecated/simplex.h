#pragma once

#include <Eigen/Dense>

using namespace Eigen;

class Simplex{
public:
  int n, m, ma, numcols, npm, status, lea_row, iteration_count;
  VectorXd v_cache;
  double min_ratio, d_cache, sum;
  double* tableau;
  double* bhead;
public:
  ~Simplex();
  Simplex(int core, const MatrixXd& A, const VectorXd& b, const VectorXd& c, const VectorXd& u);
private:
  void pivot(int horizon, const VectorXd& u, int ent_col);
  void selectEnteringColumn(int horizon, const VectorXd& c, const VectorXd& u, double& ent_value, int& ent_col);
};