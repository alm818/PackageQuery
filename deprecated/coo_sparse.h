#pragma once

#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class CooSparse{

public:
  int n, m;
  vector<vector<pair<int,double>>> rows;

public:
  CooSparse();
  CooSparse(int n, int m);
  int getSize();
  void addEntry(int r, int c, double v);
  VectorXd vectorProduct(const VectorXd& x);
  void inplaceVectorProduct(const VectorXd& x, VectorXd& res);
  void print();
};