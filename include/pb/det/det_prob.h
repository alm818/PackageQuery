#pragma once

#include "det_bound.h"

#include "pb/util/udeclare.h"

using namespace pb;

// bl <= Ax <= bu
// l <= x <= u
// Maximizing cx
class DetProb{
public:
  RMatrixXd A;
  VectorXd bl, bu, c, l, u;
  vector<long long> ids;
  DetBound detBound;
public:
  ~DetProb();
  DetProb();
  DetProb(int m, int n);
  void resize(int m, int n);
  void uniformGenerate(int n, int expected_n, double att_var, double outlier_prob, bool restrict_count=true, bool is_positive=false, bool is_translate=false, int seed=-1);
  void normalGenerate(int n, int expected_n, double att_var, double outlier_prob, bool restrict_count=true, int seed=-1);
  void tableGenerate(string table_name, vector<string>& cols, bool is_maximize, int n, int seed=-1);
  double boundGenerate(double E, double alpha, double hardness);
  void normalizeObjective();
};