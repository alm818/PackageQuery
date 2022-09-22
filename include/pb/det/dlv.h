#pragma once

#include "pb/util/udebug.h"
#include "pb/util/upostgres.h"
#include "libpq-fe.h"

class DynamicLowVariance{
private:
  string _sql;
  PgManager *pg;
  PGconn *_conn;
  PGresult *_res;
public:
  static double kGroupRatio, kVarScale;
  static long long kLpSize;
  Profiler pro;
private:
  void init();
  long long doPartition(string table_name, string suffix, const vector<string> &cols);
public:
  ~DynamicLowVariance();
  DynamicLowVariance();
  void partition(string table_name, string partition_name);
  void partition(string table_name, string partition_name, const vector<string> &cols);
};