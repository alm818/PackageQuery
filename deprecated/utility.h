#pragma once

#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <thread>

using namespace Eigen;
using namespace std;

#define CORE_COUNT 80

constexpr char CPX_EQUAL = 'E';
constexpr char CPX_LESS_EQUAL = 'L';
constexpr char CPX_GREATER_EQUAL = 'G';
constexpr char CPX_RANGE = 'R';

constexpr double kFloatEps = 1e-4;
constexpr double kFeasTolerance = 1e-6;
constexpr double kStubbornMultiplier = 2;
constexpr double kSizeTol = 1e5;
constexpr double kRedundant = 1.05;

constexpr int kReduce = 1000;
constexpr int kFastCV = 100000000;
constexpr int kRetargetUB = 1 << 10;
constexpr int kLatticeUB = 1 << 6;
constexpr int kBoundTol = 1;
constexpr int LS_NOT_FOUND = 0;
constexpr int LS_FOUND = 1;
constexpr int LS_INFEASIBLE = 2;
constexpr int LS_UNBOUNDED = 3;
constexpr int LS_FEASIBLE = 4;
constexpr int LS_DUAL_UNBOUNDED = 5;

constexpr bool kIgnoreDegenerate = true;
constexpr bool kDoPresolve = false;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> RMatrixXd;

//const string kProjectHome = "/home/alm818/package_query";
const string kProjectHome = "C:/Users/xuana/Desktop/VisualStudioCode/PackageQuery";

bool isEqual(double x, double y, double eps=kFloatEps);
bool isLess(double x, double y, double eps=kFloatEps);
bool isGreater(double x, double y, double eps=kFloatEps);
bool isLessEqual(double x, double y, double eps=kFloatEps);
bool isGreaterEqual(double x, double y, double eps=kFloatEps);
long long ceilDiv(long long x, long long q);

void showHistogram(vector<int> x, int bucket_count, double start, double end);
void showHistogram(VectorXd x, int bucket_count, double start, double end);
VectorXd readSolution(string problem);
vector<string> getAllFiles(string root, string ext);
unordered_map<string, int> getProblemSizes();
VectorXd bucketSort(vector<double> array);
string solMessage(int sol_status);
string solCombination(VectorXd sol);
string join(vector<string> names, string delim=",");

//ccfits
MatrixXd readTable(string file_name, const vector<int>& cols, vector<string>& column_names);
int getNumRows(string file_name);

template <typename T>
double sign(T value) {
  if (value < 0) return -1.0;
  if (value > 0) return 1.0;
  return 0.0;
}

template <typename T>
double nonNegativeSign(T value) {
  if (value < 0) return -1.0;
  return 1.0;
}

template <typename T>
string print(Matrix<T, Dynamic, 1> v) {
  string s = "[";
  for (int i = 0; i < v.size() - 1; i++) {
    s += to_string(v(i)) + " ";
  }
  if (v.size() > 0) {
    s += to_string(v(v.size() - 1));
  }
  s += "]\n";
  cout << s;
  return s;
}

template <typename T>
string shortPrint(Matrix<T, Dynamic, 1> v) {
  string s = "[";
  for (int i = 0; i < v.size() - 1; i++) {
    if (isEqual(v(i), 0)) continue;
    s += to_string(i) + ":" + to_string(v(i)) + " ";
  }
  if (v.size() > 0 && !isEqual(v(v.size() - 1), 0)) {
    s += to_string(v.size() - 1) + ":" + to_string(v(v.size() - 1));
  }
  s += "]\n";
  cout << s;
  return s;
}

template<typename F, typename... Args>
double exeTime(F func, Args&&... args) {
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  func(forward<Args>(args)...);
  return chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - t1).count() / 1000000.0;
}

// C++ implementation of Welford's running mean-var algorithm
class MeanVar{
public:
  VectorXd mean, M2;
  int sample_count, attr_count;
public:
  MeanVar();
  MeanVar(int attr_count);
  void add(const VectorXd& x);
  void add(double* start, int cycle=1);
  VectorXd getMean();
  VectorXd getVar();
};

class ScalarMeanVar{
public:
  double mean, M2;
  int sample_count;
public:
  ScalarMeanVar();
  void add(double x);
  void reset();
  double getMean();
  double getVar();
};

constexpr double kGalaxyOutlierProb = 0.6;

class GalaxyDB{
public:
  static vector<int> galaxy_cols, Q1, Q2;
public:
  MatrixXd data;
  MeanVar mv;
  int num_rows, num_cols;
  vector<string> column_names;
public:
  GalaxyDB(int max_rows);
  void generateQuery(double percent, int c_att, vector<int> atts, int expected_sol_size, MatrixXd& A, VectorXd& b, VectorXd& c, VectorXd& u);
};

class Profiler{
public:
  int n;
  VectorXd time;
  VectorXi count;
  vector<chrono::high_resolution_clock::time_point> eps;
  vector<string> names;
  Profiler();
  Profiler(vector<string> names);
  void clock(int i, bool is_parallel=true);
  void stop(int i, bool is_parallel=true);
  void print();
};

// Return the number of threads that would be executed in parallel regions
int GetMaxThreads();
// Set the number of threads that would be executed in parallel regions
void SetNumThreads(int num_threads);
// Return the thread number, which lies in [0, the number of threads)
int GetThreadId();