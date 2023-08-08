#define FMT_HEADER_ONLY

#include <sstream>
#include <string>
#include <numeric>
#include <fstream>
#include <random>
#include <float.h>
#include <limits.h>
#include <omp.h>
#include <filesystem>

#include "ilcplex/ilocplex.h"
#include "utility.h"
#include "fitsio.h"
#include <fmt/core.h>

namespace fs = filesystem;

void showHistogram(vector<int> x, int bucket_count, double start, double end){
  VectorXd dx (x.size());
  for (int i = 0; i < (int)x.size(); i ++) dx(i) = x[i];
  showHistogram(dx, bucket_count, start, end);
}

void showHistogram(VectorXd x, int bucket_count, double start, double end){
  assert(bucket_count >= 1);
  double min_x = x.minCoeff();
  double max_x = x.maxCoeff();
  if (start <= min_x) start = min_x;
  if (end >= max_x) end = max_x;
  if (end <= start){
    start = min_x;
    end = max_x;
  }
  int outlier = 0;
  vector<int> buckets (bucket_count, 0);
  double dist = (end-start) / bucket_count;
  for (double x_val : x){
    if (start <= x_val && x_val <= end){
      int index = (int) min(floor((x_val - start) / dist), bucket_count-1.0);
      buckets[index] ++;
    } else outlier ++;
  }
  for (int i = 0; i <= bucket_count; i ++){
    double val = start + i * dist;
    fmt::print("{: 8.2Lf}", val);
  }
  cout << endl;
  for (int i = 0; i < bucket_count; i ++){
    fmt::print("{: 8d}", buckets[i]);
  }
  fmt::print("{: 8d}", outlier);
  cout << endl;
}

long long ceilDiv(long long x, long long q){
  lldiv_t d = lldiv(x, q);
  return d.quot + (d.rem != 0LL);
}

bool isEqual(double x, double y, double eps){
  return fabs(x-y) < eps;
}

bool isLessEqual(double x, double y, double eps){
  return x-y < eps;
}

bool isGreaterEqual(double x, double y, double eps){
  return x-y > -eps;
}

bool isLess(double x, double y, double eps){
  return x-y < -eps;
}

bool isGreater(double x, double y, double eps){
  return x-y > eps;
}

vector<string> sol_messages = {
  "NOT FOUND",
  "FOUND",
  "INFEASIBLE",
  "UNBOUNDED",
  "FEASIBLE",
  "DUAL UNBOUNDED"
};

string solMessage(int sol_status){
  return sol_messages[sol_status];
}

string solCombination(VectorXd sol){
  string combs = "";
  double frac, whole;
  for (int i = 0; i < sol.size(); i ++){
    frac = modf(sol(i), &whole);
    assert(isEqual(frac, 0));
    int count = sol(i);
    if (count == 1) combs += to_string(i) + " ";
    else if (count > 1) combs += fmt::format("{}x{} ", to_string(i), count);
  }
  return combs;
}

string join(vector<string> names, string delim){
  string res = "";
  for (int i = 0; i < names.size()-1; i++) res += names[i] + delim;
  res += names[names.size()-1];
  return res;
}

const char* space = " \t\n\r\f\v";

// trim from end of string (right)
inline string& rtrim(string& s, const char* t = space){
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

// trim from beginning of string (left)
inline string& ltrim(string& s, const char* t = space){
  s.erase(0, s.find_first_not_of(t));
  return s;
}

// trim from both ends of string (right then left)
inline string& trim(string& s, const char* t = space){
  return ltrim(rtrim(s, t), t);
}

VectorXd readSolution(string problem) {
  int status = 0;
  CPXENVptr env = CPXXopenCPLEX(&status); assert(!status);
  CPXLPptr model = CPXXcreateprob(env, &status, problem.c_str()); assert(!status);
  assert(!CPXXreadcopyprob(env, model, fmt::format("{}/resource/test_cases/{}.mps", kProjectHome, problem).c_str(), NULL));
  int n = CPXXgetnumcols(env, model);
  long long surplus;
  assert(CPXXgetcolname(env, model, NULL, NULL, 0, &surplus, 0, n-1) == CPXERR_NEGATIVE_SURPLUS);
  long long storespace = -surplus;
  char** name = new char*[n];
  char* namestore = new char[storespace];
  assert(!CPXXgetcolname(env, model, name, namestore, storespace, &surplus, 0, n-1));
  assert(!CPXXfreeprob(env, &model));
  assert(!CPXXcloseCPLEX(&env));

  unordered_map<string, int> name_map;
  for (int i = 0; i < n; i ++){
    string var_name (name[i]);
    name_map[trim(var_name)] = i;
  }
  delete[] namestore;

  string line;
  VectorXd s (n); s.fill(0);
  ifstream infile(fmt::format("{}/resource/solutions/{}.sol", kProjectHome, problem));
  while (getline(infile, line)) {
    if (line.length() > 0 && line.at(0) != '='){
      auto pos = line.find_first_of(space);
      if (pos != string::npos){
        string val_str = line.substr(pos, string::npos);
        string name_str = line.substr(0, pos);
        double val = stod(trim(val_str));
        s(name_map[trim(name_str)]) = val;
      }
    }
  }
  infile.close();
  return s;
}

/**
 * \brief   Return the filenames of all files that have the specified extension
 *          in the specified directory and all subdirectories.
 */
vector<string> getAllFiles(string root, string ext){
  vector<string> paths;
  for (auto &p : fs::recursive_directory_iterator(root)){
    if (p.path().extension() == ext)
      paths.push_back(p.path().stem().string());
  }
  return paths;
}

unordered_map<string, int> getProblemSizes(){
  unordered_map<string, int> problems;
  auto paths = getAllFiles(fmt::format("{}/resource/solutions", kProjectHome), ".sol");
  for (auto path : paths){
    VectorXd sol = readSolution(path);
    problems[path] = sol.size();
  }
  return problems;
}

// Assuming uniform distribution [0,1]
// Return sorted array with 1 appended at the end
// O(n) on average
VectorXd bucketSort(vector<double> array){
  int n = array.size() + 1;
  VectorXd res (n); res(n-1) = 1.0;
  vector<vector<double>> buckets (n);
  for (int i = 0; i < (int)array.size(); i ++){
    int bucket_index = min((int)floor(array[i]*n), n-1);
    buckets[bucket_index].push_back(array[i]);
  }
  int index = 0;
  for (int i = 0; i < n; i ++){
    for (auto v : buckets[i]){
      res(index) = v;
      index ++;
    }
  }
  // Insertion sort
  for (int i = 1; i < n - 1; i ++){
    double key = res(i);
    int j = i - 1;
    while (j >= 0 && res(j) > key){
      res(j+1) = res(j);
      j --;
    }
    res(j+1) = key;
  }
  return res;
}

MeanVar::MeanVar(){
}


MeanVar::MeanVar(int attr_count){
  mean.resize(attr_count); mean.fill(0);
  M2.resize(attr_count); M2.fill(0);
  this->attr_count = attr_count;
  sample_count = 0;
}

void MeanVar::add(const VectorXd& x){
  sample_count ++;
  VectorXd delta = x - mean;
  mean += delta / sample_count;
  M2 += delta.cwiseProduct(x - mean);
}

void MeanVar::add(double* start, int cycle){
  sample_count ++;
  for (int i = 0; i < attr_count; i ++){
    double x = start[i*cycle];
    double delta = x - mean(i);
    mean(i) += delta / sample_count;
    M2(i) += delta * (x - mean(i));
  }
}

VectorXd MeanVar::getMean(){
  return mean;
}

VectorXd MeanVar::getVar(){
  if (sample_count == 0) return M2;
  return M2 / sample_count;
}

ScalarMeanVar::ScalarMeanVar(){
  mean = 0;
  M2 = 0;
  sample_count = 0;
}

void ScalarMeanVar::reset(){
  mean = 0;
  M2 = 0;
  sample_count = 0;
}

void ScalarMeanVar::add(double x){
  sample_count ++;
  double delta = x - mean;
  mean += delta / sample_count;
  M2 += delta * (x - mean);
}

double ScalarMeanVar::getMean(){
  return mean;
}

double ScalarMeanVar::getVar(){
  if (sample_count == 0) return 0;
  return M2 / sample_count;
}

int getNumRows(string file_name){
  fitsfile *fptr;
  int status = 0;
  int hdunum, hdutype;
  long nrows;
  assert(!fits_open_file(&fptr, file_name.c_str(), READONLY, &status));
  if (fits_get_hdu_num(fptr, &hdunum) == 1) fits_movabs_hdu(fptr, 2, &hdutype, &status);
  else fits_get_hdu_type(fptr, &hdutype, &status);
  assert(hdutype != IMAGE_HDU);
  fits_get_num_rows(fptr, &nrows, &status);
  assert(!fits_close_file(fptr, &status));
  return nrows;
}

MatrixXd readTable(string file_name, const vector<int>& cols, vector<string>& column_names){
  fitsfile *fptr;
  int status = 0;
  int hdunum, hdutype, ncols, anynul;
  long nrows;
  assert(!fits_open_file(&fptr, file_name.c_str(), READONLY, &status));
  if (fits_get_hdu_num(fptr, &hdunum) == 1) fits_movabs_hdu(fptr, 2, &hdutype, &status);
  else fits_get_hdu_type(fptr, &hdutype, &status);
  assert(hdutype != IMAGE_HDU);
  fits_get_num_rows(fptr, &nrows, &status);
  fits_get_num_cols(fptr, &ncols, &status);
  MatrixXd table (nrows, cols.size()); // Column major
  char colname[100];
  char tpl[] = "*";
  int colnum = -1;
  int cur = 0;
  column_names.clear();
  while (status != COL_NOT_FOUND && cur < (int)cols.size()){
    fits_get_colname(fptr, CASEINSEN, tpl, colname, &colnum, &status);
    if (colnum == cols[cur]){
      string name (colname);
      column_names.push_back(name);
      cur ++;
    }
  }
  status = 0;
  for (int i = 0; i < (int)cols.size(); i ++){
    assert(!ffgcvd(fptr, cols[i]+1, 1, 1, nrows, -DBL_MAX, &table(0, i), &anynul, &status));
  }
  assert(!fits_close_file(fptr, &status));
  return table;
}

vector<int> GalaxyDB::galaxy_cols = {2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 55, 56, 57, 58, 59, 60, 61, 62, 63, 70, 74, 75, 76, 79};
vector<int> GalaxyDB::Q1 = {9, -9, 16, -20};
vector<int> GalaxyDB::Q2 = {14, -14, 4, -26};

GalaxyDB::GalaxyDB(int max_rows){
  vector<string> paths = getAllFiles(fmt::format("{}/resource/galaxy", kProjectHome), ".fits");
  vector<MatrixXd> tables;
  num_rows = 0;
  for (int i = 0; i < (int)paths.size(); i ++){
    string name = paths[i];
    string file_name = fmt::format("{}/resource/galaxy/{}.fits", kProjectHome, name);
    int cur_rows = getNumRows(file_name);
    if (num_rows + cur_rows < max_rows){
      MatrixXd table = readTable(file_name, GalaxyDB::galaxy_cols, column_names);
      num_rows += table.rows();
      tables.push_back(table);
    }
    fmt::print("[{}/{}]Finished reading {}\n", i+1, paths.size(), file_name);
  }
  num_cols = column_names.size();
  data.resize(num_rows, num_cols);
  int row_offset = 0;
  for (MatrixXd table : tables){
    int cur_rows = table.rows();
    data.middleRows(row_offset, cur_rows) = table;
    row_offset += cur_rows;
  }
  mv = MeanVar(num_cols);
  for (int i = 0; i < num_rows; i ++) mv.add(data.row(i).transpose());
  fmt::print("Finished loading galaxy database with {} rows and {} cols\n", num_rows, num_cols);
}

void GalaxyDB::generateQuery(double percent, int c_att, vector<int> atts, int expected_sol_size, MatrixXd& A, VectorXd& b, VectorXd& c, VectorXd& u){
  int n = (int) num_rows*percent;
  int m = atts.size() + 1;
  A.resize(m, n); b.resize(m); c.resize(n); u.resize(n); u.fill(1);
  vector<int> index (num_rows);
  iota(index.begin(), index.end(), 0);
  shuffle(index.begin(), index.end(), std::mt19937{std::random_device{}()});
  VectorXd chebyshev (m);
  VectorXd mean = mv.getMean();
  VectorXd var = mv.getVar();
  double k = sqrt(1/kGalaxyOutlierProb);
  for (int i = 0; i < m-1; i ++){
    int att_index = abs(atts[i]) - 1;
    double mu = expected_sol_size * mean(att_index);
    double dist = k * expected_sol_size * var(att_index);
    if (atts[i] > 0) b(i) = mu + dist;
    else b(i) = dist - mu;
  }
  b(m-1) = -1;
  #pragma omp parallel for num_threads(CORE_COUNT)
  for (int j = 0; j < n; j ++){
    if (c_att >= 0) c(j) = data(index[j], c_att);
    else c(j) = 1.0;
    for (int i = 0; i < m; i ++){
      if (i < m-1){
        int att_index = abs(atts[i]) - 1;
        A(i, j) = sign(atts[i]) * data(index[j], att_index);
      } else {
        // i = m-1
        A(i, j) = -1.0;
      }
    }
  }
}

Profiler::Profiler(){
}

Profiler::Profiler(vector<string> names): names(names){
  n = names.size();
  time.resize(n); time.fill(0);
  count.resize(n); count.fill(0);
  eps.resize(n);
}

void Profiler::clock(int i, bool is_parallel){
  if (is_parallel){
    #pragma omp barrier
    #pragma omp master
    {
      if (0 <= i && i < n) eps[i] = chrono::high_resolution_clock::now();
    }
    #pragma omp barrier
  } else{
    if (0 <= i && i < n) eps[i] = chrono::high_resolution_clock::now();
  }
}

void Profiler::stop(int i, bool is_parallel){
  if (is_parallel){
    #pragma omp barrier
    #pragma omp master
    {
      if (0 <= i && i < n){
        time(i) += chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - eps[i]).count() / 1000000.0;
        count(i) ++;
      }
    }
    #pragma omp barrier
  } else{
    if (0 <= i && i < n){
      time(i) += chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - eps[i]).count() / 1000000.0;
      count(i) ++;
    }
  }
}

void Profiler::print(){
  for (int i = 0; i < n; i ++){
    if (count(i) > 0){
      cout << names[i] << " Total:" << time(i) << "ms Count:" << count(i) << " Average:" << time(i)/count(i) << "ms" << endl;
    }
  }
}

// Return the number of threads that would be executed in parallel regions
int GetMaxThreads() {
  #ifdef _OPENMP
    return omp_get_max_threads();
  #else
    return 1;
  #endif
}

// Set the number of threads that would be executed in parallel regions
void SetNumThreads(int num_threads) {
  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
  #else
    if (num_threads != 1) {
      assert(!"compile with -fopenmp");
    }
  #endif
}

// Return the thread number, which lies in [0, the number of threads)
int GetThreadId() {
  #ifdef _OPENMP
    return omp_get_thread_num();
  #else
    return 0;
  #endif
}