#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <queue>
#include <random>
#include <chrono>
#include <typeinfo>
#include <numeric>

#include "map_sort.h"
#include "utility.h"
#include "parallel_pq.h"
#include "fmt/core.h"
#include "fitsio.h"
#include "omp.h"
#include "pseudo_walker.h"
#include "gurobi_solver.h"
#include "dual.h"
#include "dual_reducer.h"
#include "config.h"
#include "libpq-fe.h"
#include "problem.h"

using namespace std;
using namespace Eigen;

MeanVar readDataFromPG(string dbname, string table_name, string obj_att, vector<string> constr_atts, RMatrixXd& A, VectorXd& c, int n, VectorXi& indices){
  int nn;
  {
    string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", user, hostaddr, port, dbname, password);
    PGconn* conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;
    string count_sql = fmt::format("SELECT COUNT(*) FROM {}", table_name);
    res = PQexec(conn, count_sql.c_str());
    nn = atoi(PQgetvalue(res, 0, 0));
    PQclear(res);
    PQfinish(conn);
  }
  vector<int> sub_indices (nn); iota(sub_indices.begin(), sub_indices.end(), 0);
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  shuffle(sub_indices.begin(), sub_indices.end(), default_random_engine(seed));
  vector<int> selected (nn, -1);
  indices.resize(n);
  memcpy(&indices(0), &sub_indices[0], sizeof(int)*n);
  for (int i = 0; i < n; i ++) selected[indices(i)] = i;
  int m = constr_atts.size();
  A.resize(m, n);
  c.resize(n);
  vector<bool> is_count (m, false);
  vector<string> cols;
  for (int i = 0; i < m; i ++){
    if (constr_atts[i].length() == 0) is_count[i] = true;
    else cols.push_back(constr_atts[i]);
  }
  cols.push_back(obj_att);
  string col_names = join(cols);
  int nS = 0;
  VectorXd mean (m); mean.fill(0);
  VectorXd M2 (m); M2.fill(0);
  #pragma omp parallel num_threads(db_physical_core)
  {
    int db_chunk = nn / omp_get_num_threads();
    string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", user, hostaddr, port, dbname, password);
    PGconn* conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    MeanVar mv = MeanVar(m);
    int tn = omp_get_thread_num();
    int left = tn * db_chunk;
    int right = (tn+1)*db_chunk;
    PGresult *res = NULL;
    string select_sql = fmt::format("SELECT {} FROM {} WHERE {}.id BETWEEN {} AND {}", col_names, table_name, table_name, left+1, right);
    res = PQexec(conn, select_sql.c_str());
    for (int i = 0; i < PQntuples(res); i ++){
      if (selected[i+left] >= 0){
        int index = selected[i+left];
        int cnt = 0;
        for (int j = 0; j < m; j ++){
          if (!is_count[j]){
            A(j, index) = atof(PQgetvalue(res, i, cnt));
            cnt ++;
          } else A(j, index) = 1;
        }
        c(index) = atof(PQgetvalue(res, i, cols.size()-1));
        mv.add(A.col(index));
      }
    }
    PQclear(res);
    PQfinish(conn);
    #pragma omp critical
    {
      int tmp_nS = nS + mv.sample_count;
      VectorXd delta = mv.mean - mean;
      mean = (nS*mean + mv.sample_count*mv.mean) / tmp_nS;
      M2 += mv.M2 + delta.cwiseProduct(delta) * (nS * mv.sample_count / (double) (tmp_nS));
      nS = tmp_nS;
    }
  }
  MeanVar mv = MeanVar();
  mv.sample_count = nS;
  mv.M2 = M2;
  mv.mean = mean;
  mv.var = mv.M2 / (nS - 1);
  return mv;
}

void generateQuery(string dbname, string table_name, vector<string> constr_atts, string obj_att, vector<int> bounds, int k, RMatrixXd& A, VectorXd& bl, VectorXd& bu, VectorXd& c, int n, VectorXi& indices){
  MeanVar mv = readDataFromPG(dbname, table_name, obj_att, constr_atts, A, c, n, indices);
  int m = A.outerSize();
  bl.resize(m);
  bu.resize(m);
  VectorXd mean = mv.mean * k;
  VectorXd var = mv.var * k;
  double outlier_prob = 0.9;
  VectorXd tol = var*sqrt(1/outlier_prob);
  for (int i = 0; i < m; i ++){
    if (constr_atts[i].length() > 0){
      if (bounds[i] <= 0) bl(i) = mean(i) - tol(i);
      else bl(i) = -DBL_MAX;
      if (bounds[i] >= 0) bu(i) = mean(i) + tol(i);
      else bu(i) = DBL_MAX;
    } else{
      if (bounds[i] <= 0) bl(i) = k*1.0/2;
      else bl(i) = -DBL_MAX;
      if (bounds[i] >= 0) bu(i) = k*3.0/2;
      else bu(i) = DBL_MAX;
    }
  }
}

double partition_reducer(string dbname, string table_name, int pindex, RMatrixXd& A, VectorXd& bl, VectorXd& bu, VectorXd& c, VectorXd& l, VectorXd& u, VectorXi& indices){
  string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", user, hostaddr, port, dbname, password);
  PGconn* conn = PQconnectdb(conninfo.c_str());
  assert(PQstatus(conn) == CONNECTION_OK);
  double group_ratio = 0.01;
  PGresult *res = NULL;
  string count_sql = fmt::format("SELECT COUNT(*) FROM {}", table_name);
  res = PQexec(conn, count_sql.c_str());
  int nn = atoi(PQgetvalue(res, 0, 0));
  PQclear(res);
  string gc_sql = fmt::format("SELECT COUNT(*) FROM {}_G{}", table_name, pindex);
  res = PQexec(conn, gc_sql.c_str());
  int group_count = atoi(PQgetvalue(res, 0, 0));
  PQclear(res);
  vector<vector<int>> groups (group_count, vector<int>());
  int n = indices.size();
  vector<int> selected (nn, -1);
  for (int i = 0; i < n; i ++) selected[indices(i)] = i;
  string select_sql = fmt::format("SELECT * FROM {}_P{}", table_name, pindex);
  res = PQexec(conn, select_sql.c_str());
  for (int i = 0; i < PQntuples(res); i ++){
    int tid = atoi(PQgetvalue(res, i, 0)) - 1;
    if (selected[tid] >= 0){
      int gid = atoi(PQgetvalue(res, i, 1)) - 1;
      groups[gid].push_back(tid);
    }
  }
  PQclear(res);
  int m = bl.size();
  RMatrixXd Ar (m, group_count);
  VectorXd cr (group_count);
  VectorXd lr (group_count); lr.fill(0);
  VectorXd ur (group_count);

  // string mean_sql = fmt::format("SELECT id,mean_a1,mean_a2,mean_a3,mean_a4 FROM {}_G{}", table_name, pindex);
  // res = PQexec(conn, mean_sql.c_str());
  // for (int i = 0; i < group_count; i ++){
  //   int gid = atoi(PQgetvalue(res, i, 0)) - 1;
  //   for (int j = 0; j < m; j ++){
  //     if (j == m-1) Ar(j, gid) = 1;
  //     else{
  //       Ar(j, gid) = atof(PQgetvalue(res, i, j+1));
  //     }
  //   }
  //   cr(gid) = atof(PQgetvalue(res, i, m));
  //   ur(gid) = groups[gid].size();
  // }
  // PQclear(res);

  #pragma omp parallel num_threads(db_physical_core)
  {
    #pragma omp for nowait
    for (int i = 0; i < group_count; i ++){
      MeanVar mv = MeanVar(m);
      ScalarMeanVar smv = ScalarMeanVar();
      ur(i) = groups[i].size();
      for (int j = 0; j < groups[i].size(); j++){
        int index = selected[groups[i][j]];
        mv.add(A.col(index));
        smv.add(c(index));
      }
      Ar.col(i) = mv.mean;
      cr(i) = smv.mean;
    }
  }

  Dual dual = Dual(db_physical_core, Ar, bl, bu, cr, lr, ur);
  //DualReducer dual = DualReducer(db_physical_core, &Ar, &bl, &bu, &cr, &lr, &ur, 1);
  int cnt = 0;
  for (int i = 0; i < group_count; i ++){
    if (dual.x(i) > 0){
      cout << groups[i].size() << " ";
      cnt += groups[i].size();
    }
    // if (dual.best_x(i) > 0){
    //   cout << groups[i].size() << " ";
    //   cnt += groups[i].size();
    // }
  }
  cout << endl;
  RMatrixXd Arr (m, cnt);
  VectorXd crr (cnt);
  VectorXd lrr (cnt); lrr.fill(0);
  VectorXd urr (cnt); urr.fill(1);
  int ind = 0;
  for (int i = 0; i < group_count; i ++){
    if (dual.x(i) > 0){
      auto g = groups[i];
      for (auto j : g){
        int index = selected[j];
        Arr.col(ind) = A.col(index);
        crr(ind) = c(index);
        ind ++;
      }
    }
    // if (dual.best_x(i) > 0){
    //   auto g = groups[i];
    //   for (auto j : g){
    //     int index = selected[j];
    //     Arr.col(ind) = A.col(index);
    //     crr(ind) = c(index);
    //     ind ++;
    //   }
    // }
  }
  DualReducer dr = DualReducer(db_physical_core, &Arr, &bl, &bu, &crr, &lrr, &urr, 0);
  return dr.best_score;
}

void testReader(){
  string dbname = "synthetic";
  string table_name = "n100000000_t1_v100";
  vector<string> constr_atts = {"a2", "a3", "a4", ""};
  vector<int> bounds = {0, -1, 1, 0};
  string obj_att = "a1";
  RMatrixXd A;
  VectorXd c, bl, bu;
  VectorXi indices;
  int n = 10000000;
  generateQuery(dbname, table_name, constr_atts, obj_att, bounds, 20, A, bl, bu, c, n, indices);
  int m = A.outerSize();
  VectorXd l (n); l.fill(0);
  VectorXd u (n); u.fill(1);
  double my_obj = partition_reducer(dbname, table_name, 0, A, bl, bu, c, l, u, indices);
  DualReducer dr = DualReducer(db_physical_core, &A, &bl, &bu, &c, &l, &u, 1);
  double r0_obj = dr.duals[0]->score;
  double gap = (r0_obj - dr.best_score) / r0_obj * 100;
  cout << solMessage(dr.status) << endl;
  // cout << solCombination(dr.best_x) << endl;
  cout << dr.best_score << " " << r0_obj << " " << gap << "%" << endl;
  double gap2 = (r0_obj - my_obj) / r0_obj * 100;
  cout << my_obj << " " << gap2 << "%" << " " << r0_obj / my_obj << endl;
  // GurobiSolver gs = GurobiSolver(A, bl, bu, c, l, u);
  // gs.solveRelaxed();
  // cout << gs.exe_relaxed << " " << gs.iteration_count << endl;
  // cout << gs.relaxed_cscore << endl;
  // GurobiSolver igs = GurobiSolver(A, bl, bu, c, l, u);
  // igs.solveIlp();
  // double gap3 = (igs.ilp_cscore - dr.best_score) / igs.ilp_cscore * 100;
  // cout << igs.ilp_cscore << " " << gap3 << "%" << endl;
}

int main(){
  testReader();
}