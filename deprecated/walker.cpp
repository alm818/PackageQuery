#include <iostream>
#include <float.h>

#include "walker.h"
#include "utility.h"

using namespace std;

Walker::~Walker(){
}

Walker::Walker(VectorXd p) {
  int dim = (int) p.size();
  this->p = p / p.norm();
  x.resize(dim);
  x.fill(0);
  xp = 0;
  for (int i = 0; i < dim; i++) {
    if (abs(this->p(i)) > kFloatEps)
      nz_dim.push_back(i);
  }
}

double Walker::tryStep(int i) {
  double s = sign(p(i));
  return 2 * s * (x(i) - p(i) * xp) - p(i) * p(i);
}

int Walker::executeStep(int i) {
  double s = sign(p(i));
  x(i) += s;
  xp += s * p(i);
  return s*(i+1);
}

int Walker::step() {
  double min_score = DBL_MAX;
  int min_i = -1;
  for (const auto &i : nz_dim) {
    double score = tryStep(i);
    if (min_score > score) {
      min_score = score;
      min_i = i;
    }
  }
  return executeStep(min_i);
}