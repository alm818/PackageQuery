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
#include "lattice_solver.h"
#include "pseudo_walker.h"
#include "simplex.h"
#include "gurobi_lattice_solver.h"
#include "gurobi_solver.h"
#include "dual.h"
#include "dual_reducer.h"
#include "reducer.h"
#include "fmt/core.h"
#include "config.h"
#include "libpq-fe.h"

using namespace std;
using namespace Eigen;

void generateQuery(int qi, int n, RMatrixXd& A, VectorXd& bl, VectorXd& bu, VectorXd& c){
  double outlier_prob = 0.6;
  vector<double> lbs = {-1, -1, -1, -1, 0, 0, 0, 0};
  vector<int> expected_ns = {20, 20, 500, 500, 20, 20, 500, 500};
  vector<double> att_vars = {1, 100, 1, 100, 1, 100, 1, 100};
  int expected_n = expected_ns[qi];
  double att_var = att_vars[qi];
  A.resize(4, n); bl.resize(4); bu.resize(4); c.resize(n); 
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  double left = lbs[qi];
  double right = 1.0;
  uniform_real_distribution u_dist(left, right);
  double multiplicity = att_var / (right-left) * 12;
  //normal_distribution att_dist(10.0, att_var);
  int expected_numvar = expected_n;
  double mean = 0.5*multiplicity*expected_numvar;
  double var = att_var*expected_numvar;
  normal_distribution n_dist(0.0, var);
  //normal_distribution n_dist_c(0.0, 1.0/12);
  uniform_real_distribution n_dist_c(left, right);
  #pragma omp parallel for num_threads(CORE_COUNT)
  for (int i = 0; i < n; i ++){
    A(0, i) = u_dist(gen)*multiplicity;
    A(1, i) = u_dist(gen)*multiplicity;
    A(2, i) = u_dist(gen)*multiplicity;
    A(3, i) = 1;
    c(i) = n_dist_c(gen)*multiplicity;
    //c(i) = A(0, i) + A(1, i) + A(2, i);
  }
  double tol = var*sqrt(1/outlier_prob);
  bl(0) = mean - tol;
  bu(0) = mean + tol;
  bl(1) = -DBL_MAX;
  bu(1) = mean + tol;
  bl(2) = mean- tol;
  bu(2) = DBL_MAX;
  bl(3) = expected_n*1.0/2;
  bu(3) = expected_n*3.0/2;
}

void mapSort(){
  vector<string> names = {"1", "2"};
  Profiler pro = Profiler(names);
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  uniform_real_distribution u_dist(0.0, 1.0);
  int N = 100000000;
  int repeat = 1;
  map_sort::MapSort<pair<double,int>> map_sort;
  map_sort.Init(N);
  for (int j = 0; j < repeat; j ++){
    pair<double,int>* x = new pair<double,int>[N];
    for (int i = 0; i < N ; i ++) x[i] = {u_dist(gen), i};
    pair<double,int>* y = new pair<double,int>[N];
    for (int i = 0; i < N ; i ++) y[i] = {x[i].first, i};
    pro.clock(0, false);
    map_sort.Sort(x, N, 8);
    pro.stop(0, false);
    pro.clock(1, false);
    sort(y, y+N);
    pro.stop(1, false);
    for (int i = 0; i < N; i ++) assert(x[i].first == y[i].first);
    delete[] x;
    delete[] y;
  }
  pro.print();
}

void partitionTest(){
  default_random_engine gen {static_cast<long unsigned int>(time(0))};
  double var = 100;
  normal_distribution n_dist(0.0, var);
  int N = 100000;
  VectorXd a (N);
  for (int i = 0; i < N; i ++){
    a(i) = n_dist(gen);
  }
  double reduce_factor = 2;
  double lim_var = var/reduce_factor;
  vector<int> delim;
  sort(a.begin(), a.end());
  ScalarMeanVar mv = ScalarMeanVar();
  for (int i = 0; i < N; i ++){
    mv.add(a(i));
    if (mv.var > lim_var){
      mv.reset();
      delim.push_back(i);
    }
  }
  for (int i = 0; i < delim.size(); i ++) cout << delim[i] << " ";
  cout << endl;
  int cur = 0;
  mv.reset();
  for (int i = 0; i < delim.size(); i ++){
    while (cur < delim[i]){
      mv.add(a(cur));
      cur ++;
    }
    cout << i << " "<< mv.mean << " " << mv.var << endl;
    mv.reset();
  }
  showHistogram(a, 20, 0, 0);
}

// void createSyntheticData(connection& C, int n, int type, double var){
//   try{
//     work W(C);
//     string table_name = fmt::format("N{}_T{}_V{}", n, type, (int)var);

//     string clear_sql = fmt::format("DROP TABLE IF EXISTS {};", table_name);
//     W.exec(clear_sql);

//     string create_sql = fmt::format("CREATE TABLE IF NOT EXISTS {} ("\
//     "id SERIAL PRIMARY KEY,"\
//     "a1 DOUBLE PRECISION,"\
//     "a2 DOUBLE PRECISION,"\
//     "a3 DOUBLE PRECISION,"\
//     "a4 DOUBLE PRECISION"\
//     ");", table_name);
//     W.exec(create_sql);

//     double left, right;
//     if (type == 1){
//       left = -sqrt(3*var);
//       right = sqrt(3*var);
//     } else if (type == 2){
//       left = 0;
//       right = 2*sqrt(3*var);
//     }
//     W.commit();

//     default_random_engine gen {static_cast<long unsigned int>(time(0))};
//     uniform_real_distribution u_dist(left, right);
//     string insert_sql = fmt::format("INSERT INTO {}(a1,a2,a3,a4) VALUES ", table_name);
//     int bulk = 1000000;
//     for (int i = 0; i < n/bulk; i ++){
//       work Wi (C);
//       string bulk_sql = insert_sql;
//       for (int j = 0; j < bulk-1; j ++){
//         bulk_sql += fmt::format("({:.16Lf},{:.16Lf},{:.16Lf},{:.16Lf}),", u_dist(gen), u_dist(gen), u_dist(gen), u_dist(gen));
//       }
//       bulk_sql += fmt::format("({:.16Lf},{:.16Lf},{:.16Lf},{:.16Lf});", u_dist(gen), u_dist(gen), u_dist(gen), u_dist(gen));
//       Wi.exec(bulk_sql);
//       Wi.commit();
//       fmt::print("{}%\n", 100.0*(i+1)*bulk/n);
//     }
//   } catch (const exception& e){
//     cerr << e.what() << endl;
//   }
// }

// void postgresql(){
//   try{
//     string dbname = "synthetic";
//     connection C (fmt::format("dbname={} user={} password={} hostaddr={} port={}", dbname, user, password, hostaddr, port));
//     assert(C.is_open());
//     createSyntheticData(C, 10000000, 1, 100.0);
//   } catch (const exception& e){
//     cerr << e.what() << endl;
//   }
// }

struct IndexComp{
  const RMatrixXd& mat;
  int j;
  IndexComp(const RMatrixXd& mat, int j): mat(mat), j(j){
  }
  inline bool operator()(int i1, int i2){
    return mat(j, i1) < mat(j, i2);
  }
};

double kGroupTolerance = 10;

// work* getTransaction(connection* C, string dbname){
//   if (C) delete C;
//   C = new connection(fmt::format("dbname={} user={} password={} hostaddr={} port={}", dbname, user, password, hostaddr, port));
//   work* W = new work(*C);
//   return W;
// }

// ID in POSTGRES is BASED-1
// ID-1 corresponding to the index in C++
void dynamic_kd_partition(int pindex, string dbname, string table_name, int n, vector<string> cols, double group_ratio, double var_ratio){
//   vector<string> names = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
  Profiler pro = Profiler(names);
  pro.clock(0);
  int db_chunk = n / db_logical_core;
  int m = cols.size();
  RMatrixXd A (m, n);
  int nS = 0;
  VectorXd mean (m); mean.fill(0);
  VectorXd M2 (m); M2.fill(0);
  string col_names = "";
  for (int i = 0; i < cols.size()-1; i++) col_names += cols[i] + ",";
  col_names += cols[cols.size()-1];
  string g_name = fmt::format("{}_G{}", table_name, pindex);
  string p_name = fmt::format("{}_P{}", table_name, pindex);
  bool is_error = false;
  #pragma omp parallel num_threads(db_logical_core)
  {
    try{
      connection C (fmt::format("dbname={} user={} password={} hostaddr={} port={}", dbname, user, password, hostaddr, port));
      assert(C.is_open());
      nontransaction N (C);
      MeanVar mv = MeanVar(m);
      #pragma omp for nowait
      for (int i = 0; i < db_logical_core; i ++){
        int left = i * db_chunk;
        int right = (i+1)*db_chunk;
        string select_sql = fmt::format("SELECT {} FROM {} WHERE id BETWEEN {} and {};", col_names, table_name, left+1, right);
        result R (N.exec(select_sql));
        C.close();
        for (int j = 0; j < R.size(); j ++){
          const auto& tup = R[j];
          for (int k = 0; k < m; k ++){
            string str_val (tup[k].view());
            A(k, j+left) = stod(str_val);
          }
          mv.add(A.col(j+left));
        }
      }
      #pragma omp critical
      {
        int tmp_nS = nS + mv.sample_count;
        VectorXd delta = mv.mean - mean;
        mean = (nS*mean + mv.sample_count*mv.mean) / tmp_nS;
        M2 += mv.M2 + delta.cwiseProduct(delta) * (nS * mv.sample_count / (double) (tmp_nS));
        nS = tmp_nS;
      }
    } catch (const exception& e){
      cerr << e.what() << endl;
      is_error = true;
    }
  }
  if (is_error) return;
  pro.stop(0);
//   pro.clock(1);
//   // Phase-1
//   int partition_index = -1;
//   double max_var = -1;
//   for (int i = 0; i < m; i ++){
//     double var = M2(i) / (nS - 1);
//     if (max_var < var){
//       max_var = var;
//       partition_index = i;
//     }
//   }
//   pair<double, int>* pis = new pair<double, int>[n];
//   #pragma omp parallel for num_threads(db_physical_core)
//   for (int i = 0; i < n; i ++) pis[i] = {A(partition_index, i), i};
//   map_sort::Sort(pis, n, db_physical_core);
//   pro.stop(1);
//   int soft_group_lim = (int) ceil(n * group_ratio);
//   int soft_partition_lim = (int) ceil(kGroupTolerance / var_ratio);
//   int hard_group_lim = soft_group_lim + db_physical_core * soft_partition_lim;
//   int group_count = db_physical_core;
//   int heap_length = db_physical_core;
//   vector<vector<double>> lows (m, vector<double>(hard_group_lim, -DBL_MAX));
//   vector<vector<double>> highs (m, vector<double>(hard_group_lim, DBL_MAX));
//   vector<tuple<double, int, int>> max_heap (hard_group_lim);
//   vector<VectorXi*> groups (hard_group_lim, NULL);
//   pro.clock(2);
//   int chunk = n / db_physical_core + (n % db_physical_core != 0);
//   string insert_sql, g_col_names;

//   #pragma omp parallel num_threads(db_physical_core)
//   {
//     #pragma omp for nowait
//     for (int i = 0; i < db_physical_core; i ++){
//       // pro.clock(6);
//       int left = i * chunk;
//       int right = min((i+1) * chunk, n);
//       groups[i] = new VectorXi(right-left);
//       if (i > 0) lows[partition_index][i] = A(partition_index, pis[left].second);
//       if (i < db_physical_core-1) highs[partition_index][i] = A(partition_index, pis[right].second);
//       MeanVar mv = MeanVar(m);
//       for (int j = left; j < right; j ++){
//         int col_index = pis[j].second;
//         (*groups[i])(j-left) = col_index;
//         mv.add(A.col(col_index));
//       }
//       // pro.stop(6);
//       // pro.clock(7);
//       int local_pindex = -1;
//       double max_var = -1;
//       for (int j = 0; j < m; j ++){
//         double total_var = mv.var(j) * mv.sample_count;
//         if (max_var < total_var){
//           max_var = total_var;
//           local_pindex = j;
//         }
//       }
//       max_heap[i] = {max_var, local_pindex, i};
//       // pro.stop(7);
//     }
//     #pragma omp barrier
//     pro.stop(2);
//     pro.clock(3);
//     #pragma omp master
//     {
//       delete[] pis;
//       make_heap(max_heap.begin(), max_heap.begin() + heap_length);
//     }
//     pro.stop(3);
//     pro.clock(4);
//     #pragma omp barrier
//     // Phase-2
//     while (group_count < soft_group_lim){
//       int mi = -1;
//       int gi = -1; 
//       double max_var = 0;
//       #pragma omp critical
//       {
//         while (heap_length > 0 && max_var == 0){
//           tie(max_var, mi, gi) = max_heap[0];
//           pop_heap(max_heap.begin(), max_heap.begin() + heap_length);
//           heap_length --;
//         }
//       }
//       if (gi == -1) break;
//       auto& g = *groups[gi];
//       sort(g.begin(), g.end(), IndexComp(A, mi));
//       int delim_sz = soft_partition_lim;
//       vector<int> delims (delim_sz);
//       int delim_count = 0;
//       ScalarMeanVar smv = ScalarMeanVar();
//       double reduced_var = max_var / g.size() * var_ratio;
//       for (int i = 0; i < g.size(); i ++){
//         smv.add(A(mi, g(i)));
//         if (smv.var > reduced_var){
//           if (delim_count < delim_sz) delims[delim_count] = i;
//           else{
//             delim_sz += soft_partition_lim;
//             delims.resize(delim_sz);
//             delims[delim_count] = i;
//           }
//           delim_count ++;
//           smv.reset();
//           smv.add(A(mi, g(i)));
//         }
//       }
//       if (delim_count < delim_sz) delims[delim_count] = g.size();
//       else delims.push_back(g.size());
//       delim_count ++;

//       int g_start_index;
//       #pragma omp critical
//       {
//         g_start_index = group_count;
//         group_count += delim_count-1;
//         if (group_count > hard_group_lim){
//           hard_group_lim += db_physical_core * soft_partition_lim;
//           for (int i = 0; i < m; i ++){
//             lows[i].resize(hard_group_lim, -DBL_MAX);
//             highs[i].resize(hard_group_lim, DBL_MAX);
//           }
//           max_heap.resize(hard_group_lim);
//           groups.resize(hard_group_lim, NULL);
//         }
//       }
//       for (int i = 0; i < delim_count-1; i ++){
//         int g_index = g_start_index + i;
//         int group_sz = delims[i+1] - delims[i];
//         VectorXi* gptr = new VectorXi(group_sz);
//         memcpy(&(*gptr)(0), &g(delims[i]), group_sz*sizeof(int));

//         groups[g_index] = gptr;
//         for (int j = 0; j < m; j ++){
//           lows[j][g_index] = lows[j][gi];
//           highs[j][g_index] = highs[j][gi];
//         }
//         lows[mi][g_index] = A(mi, g(delims[i]));
//         if (i < delim_count-2) highs[mi][g_index] = A(mi, g(delims[i+1]));
//         else highs[mi][g_index] = highs[mi][gi];

//         // Begin Repeatable code 
//         MeanVar mv = MeanVar(m);
//         for (int j = 0; j < group_sz; j ++){
//           mv.add(A.col((*gptr)(j)));
//         }
//         int m_index = -1;
//         double m_var = 0;
//         for (int j = 0; j < m; j ++){
//           double total_var = mv.var(j) * mv.sample_count;
//           if (m_var < total_var){
//             m_var = total_var;
//             m_index = j;
//           }
//         }
//         if (m_index != -1){
//           #pragma omp critical
//           {
//             max_heap[heap_length] = {m_var, m_index, g_index};
//             heap_length ++;
//             push_heap(max_heap.begin(), max_heap.begin() + heap_length);
//           }
//         }
//         // End Repeatable code 
//       }

//       {
//         // lows is good no need to change
//         if (delim_count != 1) highs[mi][gi] = A(mi, g(delims[0]));
//         g.conservativeResize(delims[0]);

//         int group_sz = delims[0];
//         // Begin Repeatable code 
//         MeanVar mv = MeanVar(m);
//         for (int j = 0; j < group_sz; j ++){
//           mv.add(A.col(g(j)));
//         }
//         int m_index = -1;
//         double m_var = 0;
//         for (int j = 0; j < m; j ++){
//           double total_var = mv.var(j) * mv.sample_count;
//           if (m_var < total_var){
//             m_var = total_var;
//             m_index = j;
//           }
//         }
//         if (m_index != -1){
//           #pragma omp critical
//           {
//             max_heap[heap_length] = {m_var, m_index, gi};
//             heap_length ++;
//             push_heap(max_heap.begin(), max_heap.begin() + heap_length);
//           }
//         }
//         // End Repeatable code 
//       }
//     }
//   }
//   #pragma omp parallel num_threads(db_physical_core)
//   {
//     pro.stop(4);
//     pro.clock(5);
//     // Phase-3
//     #pragma omp master
//     {
//       try{
//         // int cnt = 0;
//         // for (int i = 0; i < m; i ++){
//         //   for (int j = 0; j < group_count; j ++){
//         //     if (lows[i][j] == -DBL_MAX) cnt ++;
//         //     if (highs[i][j] == DBL_MAX) cnt ++;
//         //   }
//         // }
//         // cout << m << " " << group_count << " " << cnt << endl;
//         // vector<bool> a (n, false);
//         // for (int j = 0; j < group_count; j ++){
//         //   auto& g = *groups[j];
//         //   for (auto k : g){
//         //     if (!a[k]){
//         //       a[k] = true;
//         //     } else{
//         //       cout << "WTF " << j << " " << k << endl;
//         //     }
//         //   }
//         // }
//         connection* wC = NULL;
//         connection C (fmt::format("dbname={} user={} password={} hostaddr={} port={}", dbname, user, password, hostaddr, port));
//         assert(C.is_open());
//         // Create kd_partition table
//         work* create_kdW = getTransaction(wC, dbname);
//         string create_kd = fmt::format("CREATE TABLE IF NOT EXISTS {} ("\
//         "pindex INT NOT NULL,"\
//         "table_name VARCHAR(64),"\
//         "col_names VARCHAR(512),"\
//         "group_ratio DOUBLE PRECISION,"\
//         "var_ratio DOUBLE PRECISION"\
//         ");", kd_summary_table);
//         create_kdW->exec(create_kd);
//         create_kdW->commit();
//         delete create_kdW;

//         // Check if partition exists
//         nontransaction N(C);
//         string exist_sql = fmt::format("SELECT COUNT(*) FROM {} as kd WHERE kd.pindex={} AND kd.table_name='{}';", kd_summary_table, pindex, table_name);
//         result r (N.exec(exist_sql));
//         int exists = stoi(string(r[0][0].view()));
//         if (exists){
//           work* drop_kdW = getTransaction(wC, dbname);
//           string drop_sql = fmt::format("DROP TABLE IF EXISTS {};", p_name);
//           drop_kdW->exec(drop_sql);
//           drop_sql = fmt::format("DROP TABLE IF EXISTS {};", g_name);
//           drop_kdW->exec(drop_sql);
//           string delete_sql = fmt::format("DELETE FROM {} as kd WHERE kd.pindex={} AND kd.table_name='{}';", kd_summary_table, pindex, table_name);
//           drop_kdW->exec(delete_sql);
//           drop_kdW->commit();
//           delete drop_kdW;
//         }
//         string insert_kd_sql = fmt::format("INSERT INTO {} (pindex, table_name, col_names, group_ratio, var_ratio) VALUES ({}, '{}', '{}', {:.16Lf}, {:.16Lf})", kd_summary_table, pindex, table_name, col_names, group_ratio, var_ratio);
//         work* insert_kdW = getTransaction(wC, dbname);
//         insert_kdW->exec(insert_kd_sql);
//         insert_kdW->commit();
//         delete insert_kdW;

//         // Create G table
//         g_col_names = "id,";
//         work* create_gW = getTransaction(wC, dbname);
//         string atts_names = "";
//         for (auto col : cols){
//           atts_names += fmt::format("low_{} DOUBLE PRECISION,"\
//           "high_{} DOUBLE PRECISION,"\
//           "mean_{} DOUBLE PRECISION,"\
//           "var_{} DOUBLE PRECISION,", col, col, col, col);
//           g_col_names += fmt::format("low_{},high_{},mean_{},var_{},", col, col, col, col);
//         }
//         insert_sql = fmt::format("INSERT INTO {}({}) VALUES ", g_name, g_col_names.substr(0, g_col_names.length()-1));

//         string create_g = fmt::format("CREATE TABLE IF NOT EXISTS {} ("\
//         "id INT PRIMARY KEY,"\
//         "{}"\
//         ");", g_name, atts_names.substr(0, atts_names.length()-1));
//         create_gW->exec(create_g);
//         create_gW->commit();
//         delete create_gW;
//         // Create P table
//         work* create_pW = getTransaction(wC, dbname);
//         string create_p = fmt::format("CREATE TABLE IF NOT EXISTS {} ("\
//         "tid INT,"\
//         "gid INT"\
//         ");", p_name, table_name, g_name);
//         create_pW->exec(create_p);
//         create_pW->commit();
//         delete create_pW;
//       } catch (const exception& e){
//         cerr << e.what() << endl;
//         is_error = true;
//       }
//     }
//     pro.stop(5);
//     pro.clock(6);
//     #pragma omp barrier
//     if (!is_error){
//       // Populate G table
//       try{
//         connection CC (fmt::format("dbname={} user={} password={} hostaddr={} port={}", dbname, user, password, hostaddr, port));
//         assert(CC.is_open());
//         work WW (CC);
//         string local_sql = insert_sql;
//         #pragma omp for nowait
//         for (int i = 0; i < group_count; i ++){
//           const auto& g = *groups[i];
//           MeanVar mv = MeanVar(m);
          
//           string value_insert = to_string(i+1) + ",";
//           for (int j = 0; j < g.size(); j++){
//             mv.add(A.col(g(j)));
//           }
//           for (int j = 0; j < m; j ++){
//             value_insert += fmt::format("{:.20Lf},{:.20Lf},{:.20Lf},{:.20Lf},", lows[j][i], highs[j][i], mv.mean(j), mv.var(j));
//           }
//           local_sql += fmt::format("({}),", value_insert.substr(0, value_insert.length()-1));
//         }
//         WW.exec(local_sql.substr(0, local_sql.length()-1) + ";");
//         WW.commit();
//       } catch (const exception& e){
//         cerr << e.what() << endl;
//         is_error = true;
//       }
//     }
//     pro.stop(6);
//     pro.clock(7);
//     #pragma omp barrier
//     // Populate P table
//     if (!is_error){
//       try{
//         connection CC (fmt::format("dbname={} user={} password={} hostaddr={} port={}", dbname, user, password, hostaddr, port));
//         assert(CC.is_open());
//         work WW (CC);
//         string local_sql = fmt::format("INSERT INTO {}(tid, gid) VALUES ", p_name);
//         #pragma omp for nowait
//         for (int i = 0; i < group_count; i ++){
//           const auto& g = *groups[i];
//           for (const auto& j : g){
//             local_sql += fmt::format("({},{}),", j+1, i+1);
//           }
//         }
//         WW.exec(local_sql.substr(0, local_sql.length()-1) + ";");
//         WW.commit();
//       } catch (const exception& e){
//         cerr << e.what() << endl;
//         is_error = true;
//       }
//     }
//     #pragma omp barrier
//     #pragma omp for nowait
//     for (int i = 0; i < group_count; i ++) delete groups[i];
//   }
//   pro.stop(7);
  pro.print();
}

void customSort(){
  RMatrixXd a (2,5); a << 3, 1, 0, 5, -1, -3, -1, -0, -5, 1;
  VectorXi b (5); b << 0,1,2,3,4;
  cout << a << endl;
  cout << endl;
  // print(b);
  // b.conservativeResize(3);
  print(b);
  sort(b.begin(), b.end(), IndexComp(a,1));
  print(b);
  // tuple<double, double, double> a = {3,4,3};
  // tuple<double, double, double> b = {3,2,3};
  // double c1,c2,c3;
  // tie(c1,c2,c3) = a;
  // cout << c1 << " " << c2 << " " << c3 << endl;
}

void smvTest(){
  // ScalarMeanVar smv = ScalarMeanVar();
  // cout << smv.mean << " " << smv.var << endl;
  // smv.add(1);
  // cout << smv.mean << " " << smv.var << endl;
  // smv.add(4);
  // cout << smv.mean << " " << smv.var << endl;
  int n = 20;
  VectorXi a (n);
  for (int i = 0; i < n; i ++) a(i) = i;
  VectorXi b (5); b.fill(0);
  int start = 7;
  print(b);
  memcpy(&b(0), &a(start), 5*sizeof(double));
  print(b);
}

void GPtest(){
}

int main(){
  //GPtest();
  //partitionTest();
  //postgresql();
  //mapSort();
  //smvTest();
  vector<string> cols = {"a1", "a2", "a3", "a4"};
  dynamic_kd_partition(0, "synthetic", "n10000000_t1_v100", 10000000, cols, 0.01, 0.1);
  //dynamic_kd_partition(0, "synthetic", "n10000000_t1_v100", 10000000, cols, 0.01, 0.1);
  //customSort();
  //fmt::print("{:.16Lf}\n", -DBL_MAX);
}