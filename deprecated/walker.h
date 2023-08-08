#pragma once

#include <Eigen/Dense>
#include "abstract_walker.h"

using namespace Eigen;
using namespace std;

class Walker: public AbstractWalker{

public:
  VectorXd p, x;
  double xp;
  vector<int> nz_dim;

public:
  ~Walker();
  Walker(VectorXd p);
  int step();

private:
  int executeStep(int i);
  double tryStep(int i);

};