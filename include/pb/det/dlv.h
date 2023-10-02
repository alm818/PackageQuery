#pragma once

#include "pb/util/udebug.h"
#include "pb/util/upostgres.h"
#include "libpq-fe.h"

static const double kGroupRatio = 0.1;
static const long long kTps = (long long) pow(10, 5.5);

class DynamicLowVariance{
private:
  int core;
  double original_group_ratio, main_memory;
  long long tps;
  bool is_max_var;
  string _sql;
  PgManager *pg;
  PGconn *_conn;
  PGresult *_res;
public:
  double exe;
  #if DEBUG
    Profiler pro;
  #endif
private:
  void init();
  long long doPartition(string table_name, string suffix, const vector<string> &cols, double group_ratio);
public:
  ~DynamicLowVariance();
  DynamicLowVariance(int core=kPCore, double group_ratio=kGroupRatio, double main_memory=kMainMemorySize, long long tps=kTps, bool is_max_var=true);
  void dropAllPartitions();
  void dropTempTables();
  bool existPartition(string table_name, string partition_name);
  void dropPartition(string table_name, string partition_name);
  unordered_map<string, vector<string>> getPartitions();
  void partition(string table_name, string partition_name);
  void partition(string table_name, string partition_name, const vector<string> &cols);
};