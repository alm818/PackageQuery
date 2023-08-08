#pragma once

#include <vector>

using namespace std;

// Max priority queue
class ParallelPQ {
public:
  vector<pair<double, int>>* q;
  double lazy;
  int n;
public:
  ~ParallelPQ();
  ParallelPQ(int core, vector<pair<double, int>>* arr);
  pair<double, int> peak();
  void pop();
  void subtractPeak(double pos_val);
  int size();
private:
  void heapify(int i);
};