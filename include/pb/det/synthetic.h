#pragma once

#include "pb/util/upostgres.h"
#include "pb/util/udebug.h"

string getSubtableName(string table_name, double order, int seed);

class Synthetic{
public:
  static string table_name;
  Profiler pro;
  PgManager *pg;
private:
  void init();
  void create(long long N, int ucount, int ncount, vector<double> means, vector<double> vars);
public:
  ~Synthetic();
  Synthetic();
  void createSubtable(string table_name, double order, vector<string> cols, int seed=-1);
  void createSubtable(string table_name, double order, int seed=-1);
  void createUniform(long long N, int count, double mean, double var);
  void createNormal(long long N, int count, double mean, double var);
  void createMixed(long long N, int ucount, int ncount, double mean, double var);
};