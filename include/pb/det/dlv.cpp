#define FMT_HEADER_ONLY

#include "fmt/core.h"
#include "dlv.h"
#include "pb/util/uconfig.h"
#include "pb/util/umisc.h"
#include "pb/util/unumeric.h"
#include "pb/util/udeclare.h"
#include "pb/core/map_sort.h"

using namespace pb;

const string kIntervalType = "floatrange";
const string kTempTable = "tmp";
const double kSizeBias = 0.5;
const int kTempReserveSize = 1000;
// const string kTraverseFunction = "traverse";
// const string kClearLockFunction = "clear_lock";
// const string kProcessReserveFunction = "process_reserve";
// const int kReserveSize = (int) ceilDiv(kInMemorySize, kPCore);

#define at(row, col) (m*(col)+(row))

struct IndexComp{
  const double* mat;
  int j, m;
  IndexComp(const double* mat, int j, int m): mat(mat), j(j), m(m){
  }
  inline bool operator()(int i1, int i2){
    return mat[at(j, i1)] < mat[at(j, i2)];
  }
};

double DynamicLowVariance::kGroupRatio = 0.01;
double DynamicLowVariance::kVarScale = 2.5;
long long DynamicLowVariance::kLpSize = 100000;

DynamicLowVariance::~DynamicLowVariance(){
  PQfinish(_conn);
  delete pg;
}

DynamicLowVariance::DynamicLowVariance(string dbname): dbname(dbname){
  vector<string> names = {"Init", "ComputeStat", "All", "FetchData", "ProcessData", "WriteData", "CreateIndex", "CreateTable", "gist", "id", "tid", "gid"};
  pro = Profiler(names);
  init();
}

void DynamicLowVariance::init(){
  pro.clock(0);
  pg = new PgManager(dbname);

  _conn = PQconnectdb(pg->conninfo.c_str());
  assert(PQstatus(_conn) == CONNECTION_OK);
  _res = NULL;

  _sql = fmt::format("SET client_min_messages = warning;");
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);

  _sql = fmt::format(""\
    "DO $$ BEGIN "\
      "CREATE TYPE {} AS RANGE ("\
      "subtype = float8,"\
      "subtype_diff = float8mi"\
      ");"
    "EXCEPTION "\
      " WHEN duplicate_object THEN null;"\
    "END $$;", kIntervalType);
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);

  _sql = fmt::format(""
    "CREATE TABLE IF NOT EXISTS {}("
    "	table_name VARCHAR(63) UNIQUE NOT NULL,"
    "	size BIGINT,"
    " cols TEXT[],"
    "	mean DOUBLE PRECISION[],"
    "	M2 DOUBLE PRECISION[],"
    "	amin DOUBLE PRECISION[],"
    "	amax DOUBLE PRECISION[]"
    ")", kStatTable);
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);

  _sql = fmt::format(""
    "CREATE TABLE IF NOT EXISTS {}("
    "	table_name VARCHAR(63) NOT NULL,"
    "	partition_name VARCHAR(31) NOT NULL,"
    "	cols TEXT[],"
    "	group_ratio DOUBLE PRECISION,"
    "	max_size INTEGER,"
    " layer_count INTEGER,"
    " UNIQUE (table_name, partition_name)"
    ");", kPartitionTable);
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);
  pro.stop(0);
}

bool DynamicLowVariance::checkStats(string table_name){
  _sql = fmt::format("SELECT COUNT(*) FROM {} WHERE table_name='{}';", kStatTable, table_name);
  _res = PQexec(_conn, _sql.c_str());
  int count = atoi(PQgetvalue(_res, 0, 0));
  PQclear(_res);
  return count;
}

void DynamicLowVariance::writeStats(string table_name, Stat *stat){
  _sql = fmt::format(""
    "INSERT INTO {}(table_name, size, cols, mean, M2, amin, amax) "
    "VALUES ('{}', {}, {}, {}, {}, {}, {}) "
    "ON CONFLICT (table_name) "
    "DO UPDATE SET size=EXCLUDED.size,cols=EXCLUDED.cols,mean=EXCLUDED.mean,M2=EXCLUDED.M2,amin=EXCLUDED.amin,amax=EXCLUDED.amax;", 
    kStatTable, table_name, stat->size, pgJoin(stat->cols), pgJoin(stat->mean, kPrecision), pgJoin(stat->M2, kPrecision), pgJoin(stat->amin, kPrecision), pgJoin(stat->amax, kPrecision));
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);
}

Stat* DynamicLowVariance::readStats(string table_name){
  _sql = fmt::format("SELECT size, cols, mean, M2, amin, amax FROM \"{}\" WHERE table_name='{}';", kStatTable, table_name);
  _res = PQexec(_conn, _sql.c_str());
  Stat* stat = new Stat(pgStringSplit(PQgetvalue(_res, 0, 1)));
  stat->add(atoll(PQgetvalue(_res, 0, 0)), pgValueSplit(PQgetvalue(_res, 0, 2)), pgValueSplit(PQgetvalue(_res, 0, 3)), pgValueSplit(PQgetvalue(_res, 0, 4)), pgValueSplit(PQgetvalue(_res, 0, 5)));
  PQclear(_res);
  return stat;
}

void DynamicLowVariance::partition(string table_name, string partition_name){
  vector<string> cols = pg->getNumericCols(table_name);
  partition(table_name, partition_name, cols);
}

void DynamicLowVariance::partition(string table_name, string partition_name, const vector<string> &cols){
  pro.clock(2);
  long long size = doPartition(table_name, partition_name, cols);
  string g_name = nextGName(table_name + "_" + partition_name);
  while (size){
    size = doPartition(g_name, "", cols);
    g_name = nextGName(g_name);
  }
  pro.stop(2);
}

long long DynamicLowVariance::doPartition(string table_name, string suffix, const vector<string> &cols){
  if (!checkStats(table_name)){
    vector<string> cols = pg->getNumericCols(table_name);
    Stat *stat = pg->computeStats(table_name, cols);
    writeStats(table_name, stat);
    delete stat;
  }
  Stat *stat = readStats(table_name);
  long long size = stat->size;
  if (size <= kLpSize) return 0;

  int stat_max_var_index = -1;
  int max_var_index = -1;
  double max_var = -1;
  int m = (int) cols.size();
  for (int i = 0; i < m; i ++){
    int index = stat->getIndex(cols[i]);
    double var = stat->getVar(index);
    if (max_var < var){
      max_var = var;
      max_var_index = i;
      stat_max_var_index = index;
    }
  }
  double min_att = stat->amin(stat_max_var_index);
  double max_att = stat->amax(stat_max_var_index);

  long long chunk = ceilDiv(size, kPCore);
  long long bucket = ceilDiv(size, kInMemorySize);
  long long partition_count = bucket;
  vector<MeanVar> bucket_stat (bucket, MeanVar(m));
  vector<pair<double, double>> intervals (bucket);
  map<double, long long> key_indices;
  vector<double> keys (bucket);
  vector<PGconn*> wait_conns;

  string create_tmp_table;
  string drop_tmp_table = "DROP TABLE IF EXISTS {}{}";
  {
    vector<string> att_names (m);
    for (int i = 0; i < m; i ++) att_names[i] = fmt::format("{} DOUBLE PRECISION", cols[i]);
    
    create_tmp_table = fmt::format(""
      "CREATE TABLE IF NOT EXISTS {}("
      "	tid BIGINT,"
      " {}"
      ");", "{}{}", join(att_names, ","));
  }

  string symbolic_name = table_name;
  if (suffix.length()) symbolic_name += "_" + suffix;
  string g_name = nextGName(symbolic_name);
  string p_name = nextPName(symbolic_name);
  string col_names = join(cols, ",");

  cout << "Begin 1a" << endl;
  // Phase-1a: Initial quick-partition for #bucket partitions
  double bucket_width = (max_att - min_att) / bucket;
  for (long long i = 0; i < bucket; i ++){
    _sql = fmt::format(drop_tmp_table, kTempTable, i);
    _res = PQexec(_conn, _sql.c_str());
    PQclear(_res);
    _sql = fmt::format(create_tmp_table, kTempTable, i);
    _res = PQexec(_conn, _sql.c_str());
    PQclear(_res);

    intervals[i] = {min_att + i*bucket_width, min_att + (i+1)*bucket_width};
    key_indices[(double) i] = i;
    keys[i] = i;
  }

  #pragma omp parallel num_threads(kPCore)
  {
    string sql;
    PGconn *conn = PQconnectdb(pg->conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    PGconn *_conn = PQconnectdb(pg->conninfo.c_str());
    assert(PQstatus(_conn) == CONNECTION_OK);
    PGresult *_res = NULL;

    int seg = omp_get_thread_num();
    long long start_id = seg * chunk + 1;
    long long end_id = min((seg + 1) * chunk, size);
    vector<int> sz (bucket, 0);
    vector<vector<pair<long long, VectorXd>>> cache (bucket, vector<pair<long long, VectorXd>>(kTempReserveSize));
    vector<MeanVar> cache_stat (bucket, MeanVar(m));

    sql = fmt::format("SELECT {},{} FROM \"{}\" WHERE {} BETWEEN {} AND {};", kId, col_names, table_name, kId, start_id, end_id);
    res = PQexec(conn, sql.c_str());

    for (int i = 0; i < PQntuples(res); i++){
      long long tid = atol(PQgetvalue(res, i, 0));
      VectorXd vs (m);
      for (int j = 0; j < m; j ++) vs(j) = atof(PQgetvalue(res, i, j+1));
      long long bucket_ind = min((long long) floor((vs(max_var_index) - min_att) / (max_att - min_att) * bucket), bucket-1);
      cache[bucket_ind][sz[bucket_ind]] = {tid, vs};
      cache_stat[bucket_ind].add(vs);
      sz[bucket_ind] ++;
      if (sz[bucket_ind] == kTempReserveSize){
        sql = fmt::format("COPY {}{} FROM STDIN with(delimiter ',');", kTempTable, bucket_ind);
        _res = PQexec(_conn, sql.c_str());
        assert(PQresultStatus(_res) == PGRES_COPY_IN);
        PQclear(_res);
        string data = "";
        for (int j = 0; j < sz[bucket_ind]; j ++) data += fmt::format("{},{}\n", cache[bucket_ind][j].first, join(cache[bucket_ind][j].second, kPrecision));
        assert(PQputCopyData(_conn, data.c_str(), (int) data.length()) == 1);
        assert(PQputCopyEnd(_conn, NULL) == 1);
        _res = PQgetResult(_conn);
        assert(PQresultStatus(_res) == PGRES_COMMAND_OK);
        sz[bucket_ind] = 0;
      }
    }
    PQclear(res);

    for (long long i = 0; i < bucket; i ++){
      #pragma omp critical
      {
        bucket_stat[i].add(cache_stat[i]);
      }
      if (sz[i]){
        sql = fmt::format("COPY {}{} FROM STDIN with(delimiter ',');", kTempTable, i);
        _res = PQexec(_conn, sql.c_str());
        assert(PQresultStatus(_res) == PGRES_COPY_IN);
        PQclear(_res);
        string data = "";
        for (int j = 0; j < sz[i]; j ++) data += fmt::format("{},{}\n", cache[i][j].first, join(cache[i][j].second, kPrecision));
        assert(PQputCopyData(_conn, data.c_str(), (int) data.length()) == 1);
        assert(PQputCopyEnd(_conn, NULL) == 1);
        _res = PQgetResult(_conn);
        assert(PQresultStatus(_res) == PGRES_COMMAND_OK);
      }
    }
    
    PQfinish(conn);
    PQfinish(_conn);
  }

  cout << "End 1a" << endl;

  pro.clock(3);
  pg->dropTable(g_name);
  string atts_names = fmt::format("{} BIGINT,", kId);
  for (auto col : cols){
    atts_names += fmt::format("interval_{} {},{} DOUBLE PRECISION, m2_{} DOUBLE PRECISION,", col, kIntervalType, col, col);
  }
  atts_names += "size BIGINT";
  _sql = fmt::format("CREATE TABLE IF NOT EXISTS \"{}\"({});", g_name, atts_names);
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);

  pg->dropTable(p_name);
  _sql = fmt::format("CREATE TABLE IF NOT EXISTS \"{}\"("\
    "tid BIGINT,"\
    "gid BIGINT"\
    ");", p_name);
  _res = PQexec(_conn, _sql.c_str());
  PQclear(_res);
  pro.stop(3);

  cout << "Begin 1b" << endl;

  // Phase-1b: Recursive quick-partition for partition with size > #kInMemorySize
  long long global_group_count = 0;
  for (long long p = 0; p < partition_count; p ++){
    long long recurse_size = bucket_stat[p].sample_count;
    if (recurse_size > kInMemorySize){
      double start_key = keys[p];
      auto next_it = key_indices.upper_bound(start_key);
      double end_key = start_key + 1;
      if (next_it != key_indices.end()) end_key = next_it->first;
      double key_width = (end_key - start_key) / bucket;
      double bucket_width = (intervals[p].second - intervals[p].first) / bucket;
      for (long long i = 0; i < bucket; i ++){
        _sql = fmt::format(drop_tmp_table, kTempTable, partition_count + i);
        _res = PQexec(_conn, _sql.c_str());
        PQclear(_res);
        _sql = fmt::format(create_tmp_table, kTempTable, partition_count + i);
        _res = PQexec(_conn, _sql.c_str());
        PQclear(_res);
        intervals.emplace_back(intervals[p].first + i*bucket_width, intervals[p].first + (i+1)*bucket_width);
        bucket_stat.emplace_back(m);
        double key = start_key + key_width * i;
        keys.emplace_back(key);
        key_indices[key] = partition_count + i;
      }

      long long recurse_chunk = ceilDiv(recurse_size, kPCore);
      _sql = fmt::format("SELECT * FROM {}{};", kTempTable, p);
      _res = PQexec(_conn, _sql.c_str());

      #pragma omp parallel num_threads(kPCore)
      {
        string sql;
        PGconn* conn = PQconnectdb(pg->conninfo.c_str());
        assert(PQstatus(conn) == CONNECTION_OK);
        PGresult *res = NULL;

        int seg = omp_get_thread_num();
        long long start_tuple_id = seg * recurse_chunk;
        long long end_tuple_id = min((seg + 1) * recurse_chunk - 1, recurse_size - 1);
        vector<int> sz (bucket, 0);
        vector<vector<pair<long long, VectorXd>>> cache (bucket, vector<pair<long long, VectorXd>>(kTempReserveSize));
        vector<MeanVar> cache_stat (bucket, MeanVar(m));

        for (int i = start_tuple_id; i <= end_tuple_id; i ++){
          long long tid = atol(PQgetvalue(_res, i, 0));
          VectorXd vs (m);
          for (int j = 0; j < m; j ++) vs(j) = atof(PQgetvalue(_res, i, j+1));
          long long bucket_ind = min((long long) floor((vs[max_var_index] - intervals[p].first) / bucket_width), bucket-1);
          cache[bucket_ind][sz[bucket_ind]] = {tid, vs};
          cache_stat[bucket_ind].add(vs);
          sz[bucket_ind] ++;
          if (sz[bucket_ind] == kTempReserveSize){
            sql = fmt::format("COPY {}{} FROM STDIN with(delimiter ',');", kTempTable, partition_count + bucket_ind);
            res = PQexec(conn, sql.c_str());
            assert(PQresultStatus(res) == PGRES_COPY_IN);
            PQclear(res);
            string data = "";
            for (int j = 0; j < sz[bucket_ind]; j ++) data += fmt::format("{},{}\n", cache[bucket_ind][j].first, join(cache[bucket_ind][j].second, kPrecision));
            assert(PQputCopyData(conn, data.c_str(), (int) data.length()) == 1);
            assert(PQputCopyEnd(conn, NULL) == 1);
            res = PQgetResult(conn);
            assert(PQresultStatus(res) == PGRES_COMMAND_OK);
            sz[bucket_ind] = 0;
          }
        }

        for (long long i = 0; i < bucket; i ++){
          #pragma omp critical
          {
            bucket_stat[partition_count + i].add(cache_stat[i]);
          }
          if (sz[i]){
            sql = fmt::format("COPY {}{} FROM STDIN with(delimiter ',');", kTempTable, partition_count + i);
            res = PQexec(conn, sql.c_str());
            assert(PQresultStatus(res) == PGRES_COPY_IN);
            PQclear(res);
            string data = "";
            for (int j = 0; j < sz[i]; j ++) data += fmt::format("{},{}\n", cache[i][j].first, join(cache[i][j].second, kPrecision));
            assert(PQputCopyData(conn, data.c_str(), (int) data.length()) == 1);
            assert(PQputCopyEnd(conn, NULL) == 1);
            res = PQgetResult(conn);
            assert(PQresultStatus(res) == PGRES_COMMAND_OK);
          }
        }

        PQfinish(conn);
      }

      PQclear(_res);
      PGconn* conn = PQconnectdb(pg->conninfo.c_str());
      assert(PQstatus(conn) == CONNECTION_OK);
      wait_conns.push_back(conn);
      _sql = fmt::format("DROP TABLE {}{};", kTempTable, p);
      assert(PQsendQuery(conn, _sql.c_str()));
      partition_count += bucket;
      bucket_stat[p].reset();
    }
  }

  cout << "End 1b" << endl;


  // Debug
  // for (auto it : key_indices){
  //   cout << it.second << " " << bucket_stat[it.second].sample_count << endl;
  // }

  // Phase-2a: Aggregate and compute group count for each local partition

  cout << "Begin 2a" << endl;
  vector<MeanVar> agg_bucket_stat;
  vector<pair<double, double>> agg_intervals;
  vector<vector<long long>> agg_groups;
  MeanVar mv = MeanVar(m);
  vector<long long> agg_group;
  long long current_size = 0;
  for (auto it : key_indices){
    long long sz = bucket_stat[it.second].sample_count;
    if (current_size + sz > kInMemorySize){
      agg_bucket_stat.push_back(mv);
      agg_groups.push_back(agg_group);
      double left_bound = -DBL_MAX;
      if (agg_intervals.size() > 0) left_bound = agg_intervals[agg_intervals.size() - 1].second;
      agg_intervals.emplace_back(left_bound, intervals[it.second].first);
      mv = MeanVar(m);
      agg_group = vector<long long>();
      current_size = 0;
    }
    mv.add(bucket_stat[it.second]);
    agg_group.emplace_back(it.second);
    current_size += sz;
  }
  if (agg_group.size() > 0){
    agg_bucket_stat.push_back(mv);
    agg_groups.push_back(agg_group);
    double left_bound = -DBL_MAX;
    if (agg_intervals.size() > 0) left_bound = agg_intervals[agg_intervals.size() - 1].second;
    agg_intervals.emplace_back(left_bound, DBL_MAX);
  }
  // Debug
  // long long s = 0;
  // for (int i = 0; i < agg_intervals.size(); i ++){
  //   s += agg_bucket_stat[i].sample_count;
  //   cout << agg_intervals[i].first << " " << agg_intervals[i].second << " " << agg_bucket_stat[i].sample_count << endl;
  //   print(agg_bucket_stat[i].getVar());
  //   for (long long ind : agg_groups[i]) cout << ind << " ";
  //   cout << endl;
  // }
  // cout << s << endl;

  // Phase-2b: Size assignment problem
  cout << "Begin 2b" << endl;
  long long max_size = -1;
  long long agg_count = (long long) agg_bucket_stat.size();
  max_var = -1;
  for (long long p = 0; p < agg_count; p ++){
    if (agg_bucket_stat[p].sample_count){
      max_size = max(max_size, agg_bucket_stat[p].sample_count);
      max_var = max(max_var, agg_bucket_stat[p].getVar().sum());
    }
  }
  VectorXd scores (agg_count); scores.fill(0);
  vector<long long> group_count (agg_count, 0);
  double current_bias = kSizeBias;
  while (true){
    double score_sum = 0;
    for (long long p = 0; p < agg_count; p ++){
      if (agg_bucket_stat[p].sample_count){
        scores(p) = agg_bucket_stat[p].sample_count / (double) max_size * current_bias + agg_bucket_stat[p].getVar().sum() / max_var * (1 - current_bias);
        score_sum += scores(p);
      }
    }
    bool is_feasible = true;
    for (long long p = 0; p < agg_count; p ++){
      if (agg_bucket_stat[p].sample_count){
        // cout << stat->size << " " << kGroupRatio << endl;
        group_count[p] = (long long) ceil(stat->size * kGroupRatio * scores(p) / score_sum);
        if (group_count[p] > agg_bucket_stat[p].sample_count) is_feasible = false;
      }
    }
    if (is_feasible) break;
    else current_bias = (1 + current_bias) / 2;
  }
  delete stat;

  // Debug
  // cout << "BIAS " << current_bias << endl;
  // for (long long p = 0; p < agg_count; p ++){
  //   if (agg_bucket_stat[p].sample_count){
  //     cout << agg_bucket_stat[p].sample_count << " "  << group_count[p] << endl;
  //   }
  // }

  cout << "Begin 2c" << endl;
  // Phase-2c: Local DLV
  string condition_names, rec_names, t_names;
  {
    vector<string> conditions;
    vector<string> recs;
    vector<string> ts;
    for (int i = 0; i < m; i ++){
      conditions.push_back(fmt::format("t.interval_{} @> %L::float8", cols[i]));
      recs.push_back(fmt::format("rec.{}", cols[i]));
      ts.push_back(fmt::format("t.{}", cols[i]));
    }
    condition_names = join(conditions, " AND ");
    rec_names = join(recs, ", ");
    t_names = join(ts, ",");
  }

  for (long long p = 0; p < agg_count; p ++){
    long long n = agg_bucket_stat[p].sample_count;
    if (n){
      vector<PGconn*> conns (agg_groups[p].size());
      vector<PGresult*> ress (agg_groups[p].size());
      for (long long i = 0; i < (long long) agg_groups[p].size(); i ++){
        conns[i] = PQconnectdb(pg->conninfo.c_str());
        assert(PQstatus(conns[i]) == CONNECTION_OK);
        string select_sql = fmt::format("SELECT * FROM {}{};", kTempTable, agg_groups[p][i]);
        assert(PQsendQuery(conns[i], select_sql.c_str()));
      }
      // A is column-wise
      double *A = new double [m*n];
      vector<long long> tids (n);
      pair<double, long long>* pis = new pair<double, long long>[n];
      int partition_index = -1;
      double max_var = -1;
      VectorXd var = agg_bucket_stat[p].getVar();
      for (int i = 0; i < m; i++){
        if (max_var < var(i)){
          max_var = var(i);
          partition_index = i;
        }
      }
      long long current_group_index = 0;
      long long current_tuple_index = 0;
      long long current_col_index = 0;
      #pragma omp parallel num_threads(kPCore)
      {
        int local_group_index = -1;
        int n_tuple = 0;
        long long local_tuple_index = -1;
        long long local_col_index = -1;
        bool is_done = false;
        while (true){
          #pragma omp critical
          {
            if (current_group_index < (long long) agg_groups[p].size()){
              local_group_index = current_group_index;
              local_tuple_index = current_tuple_index;
              local_col_index = current_col_index;
              long long group_size = bucket_stat[agg_groups[p][current_group_index]].sample_count;
              if (current_tuple_index + kTempReserveSize < group_size){
                current_col_index += kTempReserveSize;
                current_tuple_index += kTempReserveSize;
              } else {
                current_group_index ++;
                current_col_index += group_size - current_tuple_index;
                current_tuple_index = 0;
              }
              n_tuple = (int) (current_col_index - local_col_index);
            } else is_done = true;
            // cout << is_done << " " << p << " " << local_group_index << " " << local_tuple_index << " " << local_col_index << " " << omp_get_thread_num() << endl;
          }
          if (is_done) break;
          #pragma omp critical
          {
            if (ress[local_group_index] == nullptr){
              ress[local_group_index] = PQgetResult(conns[local_group_index]);
              assert(PQgetResult(conns[local_group_index]) == nullptr);
            }
            // cout << PQntuples(ress[local_group_index]) << " " << local_tuple_index << " " << n_tuple << endl;
          }
          for (int i = 0; i < n_tuple; i ++){
            long long col_index = local_col_index + i;
            long long tuple_index = local_tuple_index + i;
            long long tid = atoll(PQgetvalue(ress[local_group_index], tuple_index, 0));
            tids[col_index] = tid;
            for (int j = 0; j < m; j ++){
              A[at(j, col_index)] = atof(PQgetvalue(ress[local_group_index], tuple_index, j+1));
            }
            pis[col_index] = {A[at(partition_index, col_index)], col_index};
          }
        }
      }
      for (long long i = 0; i < (long long) agg_groups[p].size(); i ++){
        PQclear(ress[i]);
        PQfinish(conns[i]);
      }
      for (auto group_ind : agg_groups[p]){
        PGconn* conn = PQconnectdb(pg->conninfo.c_str());
        assert(PQstatus(conn) == CONNECTION_OK);
        wait_conns.push_back(conn);
        _sql = fmt::format("DROP TABLE {}{};", kTempTable, group_ind);
        assert(PQsendQuery(conn, _sql.c_str()));
      }

      map_sort::Sort(pis, n, kPCore);
      long long soft_group_lim = (long long) ceil(n * kGroupRatio);
      double var_ratio = kGroupRatio * kGroupRatio * kVarScale;
      int soft_partition_lim = (int) ceil(kVarScale / kGroupRatio);
      long long hard_group_lim = soft_group_lim + kPCore + soft_partition_lim;
      long long group_count = kPCore;
      long long heap_length = kPCore;
      vector<vector<double>> lows (m, vector<double>(hard_group_lim, -DBL_MAX));
      vector<vector<double>> highs (m, vector<double>(hard_group_lim, DBL_MAX));
      vector<tuple<double, int, long long>> max_heap (hard_group_lim);
      vector<vector<long long>*> groups (hard_group_lim, nullptr);
      long long local_chunk = ceilDiv(n, kPCore);

      #pragma omp parallel num_threads(kPCore)
      {
        {
          int i = omp_get_thread_num();
          lows[max_var_index][i] = agg_intervals[p].first;
          highs[max_var_index][i] = agg_intervals[p].second;
          long long left = i * local_chunk;
          long long right = min((i+1)*local_chunk, n);
          groups[i] = new vector<long long>(right - left);
          if (i > 0) lows[partition_index][i] = A[at(partition_index, pis[left].second)];
          if (i < kPCore - 1) highs[partition_index][i] = A[at(partition_index, pis[right].second)];

          MeanVar mv = MeanVar(m);
          for (int j = left; j < right; j ++){
            long long col_index = pis[j].second;
            (*groups[i])[j-left] = col_index;
            mv.add(A+at(0, col_index));
          }
          int local_pindex = -1;
          VectorXd total_vars = mv.getM2();
          double max_total_var = -1;
          for (int j = 0; j < m; j ++){
            double total_var = total_vars(j);
            if (max_total_var < total_var){
              max_total_var = total_var;
              local_pindex = j;
            }
          }
          max_heap[i] = {max_total_var, local_pindex, i};
        }
        #pragma omp barrier
        #pragma omp master
        {
          delete[] pis;
          make_heap(max_heap.begin(), max_heap.begin() + heap_length);
        }
        #pragma omp barrier
        while (group_count < soft_group_lim){
          int mi = -1;
          long long gi = -1;
          double max_total_var = 0;
          #pragma omp critical
          {
            while (heap_length > 0 && max_total_var == 0){
              tie(max_total_var, mi, gi) = max_heap[0];
              pop_heap(max_heap.begin(), max_heap.begin() + heap_length);
              heap_length --;
            }
          }
          if (max_total_var == 0) break;
          auto& g = *groups[gi];
          sort(g.begin(), g.end(), IndexComp(A, mi, m));
          int delim_sz = soft_partition_lim;
          vector<long long> delims (delim_sz);
          int delim_count = 0;
          ScalarMeanVar smv = ScalarMeanVar();
          double reduced_var = max_total_var / g.size() * var_ratio;
          for (int i = 0; i < (int) g.size(); i ++){
            smv.add(A[at(mi, g[i])]);
            if (smv.getVar() > reduced_var){
              if (delim_count < delim_sz) delims[delim_count] = i;
              else{
                delim_sz += soft_partition_lim;
                delims.resize(delim_sz);
                delims[delim_count] = i;
              }
              delim_count ++;
              smv.reset();
              smv.add(A[at(mi, g[i])]);
            }
          }
          if (delim_count < delim_sz) delims[delim_count] = (long long) g.size();
          else delims.emplace_back(g.size());
          delim_count ++;

          long long g_start_index;
          #pragma omp critical
          {
            g_start_index = group_count;
            group_count += delim_count - 1;
            if (group_count > hard_group_lim){
              hard_group_lim += kPCore * soft_partition_lim;
              for (int i = 0; i < m; i ++){
                if (i == max_var_index){
                  lows[i].resize(hard_group_lim, agg_intervals[p].first);
                  highs[i].resize(hard_group_lim, agg_intervals[p].second);
                } else{
                  lows[i].resize(hard_group_lim, -DBL_MAX);
                  highs[i].resize(hard_group_lim, DBL_MAX);
                }
              }
              max_heap.resize(hard_group_lim);
              groups.resize(hard_group_lim, nullptr);
            }
          }

          for (int i = 0; i < delim_count-1; i ++){
            long long g_index = g_start_index + i;
            long long group_sz = delims[i+1] - delims[i];
            vector<long long>* gptr = new vector<long long>(group_sz);
            memcpy(&(*gptr)[0], &g[delims[i]], group_sz*sizeof(long long));

            groups[g_index] = gptr;
            for (int j = 0; j < m; j ++){
              lows[j][g_index] = lows[j][gi];
              highs[j][g_index] = highs[j][gi];
            }
            lows[mi][g_index] = A[at(mi, g[delims[i]])];
            if (i < delim_count-2) highs[mi][g_index] = A[at(mi, g[delims[i+1]])];
            else highs[mi][g_index] = highs[mi][gi];

            MeanVar mv = MeanVar(m);
            for (long long j = 0; j < group_sz; j ++){
              mv.add(A+at(0, (*gptr)[j]));
            }
            int m_index = -1;
            VectorXd total_vars = mv.getM2();
            double m_var = 0;
            for (int j = 0; j < m; j ++){
              if (m_var < total_vars(j)){
                m_var = total_vars(j);
                m_index = j;
              }
            }
            if (m_index != -1){
              #pragma omp critical
              {
                max_heap[heap_length] = {m_var, m_index, g_index};
                heap_length ++;
                push_heap(max_heap.begin(), max_heap.begin() + heap_length);
              }
            }
          }

          {
            if (delim_count != 1) highs[mi][gi] = A[at(mi, g[delims[0]])];
            g.resize(delims[0]);

            long long group_sz = delims[0];
            MeanVar mv = MeanVar(m);
            for (long long j = 0; j < group_sz; j ++){
              mv.add(A+at(0, g[j]));
            }
            int m_index = -1;
            VectorXd total_vars = mv.getM2();
            double m_var = 0;
            for (int j = 0; j < m; j ++){
              if (m_var < total_vars(j)){
                m_var = total_vars(j);
                m_index = j;
              }
            }
            if (m_index != -1){
              #pragma omp critical
              {
                max_heap[heap_length] = {m_var, m_index, gi};
                heap_length ++;
                push_heap(max_heap.begin(), max_heap.begin() + heap_length);
              }
            }
          }
        }
        #pragma omp barrier
        PGconn *conn = PQconnectdb(pg->conninfo.c_str());
        assert(PQstatus(conn) == CONNECTION_OK);
        PGresult *res = NULL;

        #pragma omp master
        cout << "Populate G/P" << endl;
        // Populate G table
        string sql = fmt::format("COPY \"{}\" FROM STDIN with (delimiter '|', null '{}');", g_name, kNullLiteral);
        res = PQexec(conn, sql.c_str());
        assert(PQresultStatus(res) == PGRES_COPY_IN);
        PQclear(res);
        #pragma omp for nowait
        for (long long i = 0; i < group_count; i ++){
          const auto& g = *groups[i];
          MeanVar mv = MeanVar(m);
          string data = fmt::format("{}|", global_group_count+i+1);
          for (long long j : g) mv.add(A + at(0, j));
          VectorXd mean = mv.getMean();
          VectorXd M2 = mv.getM2();
          for (int j = 0; j < m; j ++){
            data += fmt::format("[{},{}]|{:.{}Lf}|{:.{}Lf}|", infAlias(lows[j][i], kPrecision), infAlias(highs[j][i], kPrecision), mean(j), kPrecision, M2(j), kPrecision);
          }
          data += fmt::format("{}\n", mv.sample_count);
          assert(PQputCopyData(conn, data.c_str(), (int) data.length()) == 1);
        }
        assert(PQputCopyEnd(conn, NULL) == 1);
        res = PQgetResult(conn);
        assert(PQresultStatus(res) == PGRES_COMMAND_OK);
        PQclear(res);

        // Populate P table
        sql = fmt::format("COPY \"{}\" FROM STDIN with(delimiter ',');", p_name);
        res = PQexec(conn, sql.c_str());
        assert(PQresultStatus(res) == PGRES_COPY_IN);
        PQclear(res);
        #pragma omp for nowait
        for (long long i = 0; i < group_count; i ++){
          const auto& g = *groups[i];
          string data = "";
          for (long long j : g){
            data += fmt::format("{},{}\n", tids[j], global_group_count+i+1);
          }
          assert(PQputCopyData(conn, data.c_str(), (int) data.length()) == 1);
        }
        assert(PQputCopyEnd(conn, NULL) == 1);
        res = PQgetResult(conn);
        assert(PQresultStatus(res) == PGRES_COMMAND_OK);
        PQclear(res);
        PQfinish(conn);
        #pragma omp barrier
        #pragma omp for nowait
        for (long long i = 0; i < group_count; i ++) delete groups[i];
      }
      global_group_count += group_count;
      delete[] A;
    }
  }  
  PGconn *conn;

  conn = PQconnectdb(pg->conninfo.c_str());
  vector<string> interval_names;
  for (int i = 0; i < min(m, kMaxMultiColumnIndexes); i ++) interval_names.push_back("interval_" + cols[i]);
  _sql = fmt::format("CREATE INDEX \"{}_group_interval\" ON \"{}\" USING gist ({});", g_name, g_name, join(interval_names, ","));
  assert(PQsendQuery(conn, _sql.c_str()));
  wait_conns.push_back(conn);

  conn = PQconnectdb(pg->conninfo.c_str());
  _sql = fmt::format("ALTER TABLE \"{}\" ADD PRIMARY KEY ({});", g_name, kId);
  assert(PQsendQuery(conn, _sql.c_str()));
  wait_conns.push_back(conn);

  conn = PQconnectdb(pg->conninfo.c_str());
  _sql = fmt::format("CREATE INDEX \"{}_gid_index\" ON \"{}\" USING btree (gid);", p_name, p_name);
  assert(PQsendQuery(conn, _sql.c_str()));
  wait_conns.push_back(conn);

  for (PGconn* conn : wait_conns){
    while (PQgetResult(conn));
    PQfinish(conn);
  }

  return global_group_count;
}