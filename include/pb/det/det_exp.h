#pragma once

#include <fstream>

#include "pb/util/upostgres.h"
#include "pb/det/det_sql.h"

using std::fstream;

static string kOutFolder = "plots";
static string kBackupFolder = "backup";

class DetExp{
private:
  fstream out, backup;
  PgManager *pg;
  bool verbose;
  vector<string> lines;
  string backup_path, out_file;
public:
  // Expected number of tuples in solution package
  double E;
  // Variness of means of constraint bounds
  double a;
  // Hardness of the query
  double H;
  // Filter ratio of the base predicate
  double F;
  // Group ratio of the partition
  double g;
  // Main memory used in partition in GB
  double M;
  // Tau ratio for KD Tree
  double tau_ratio;
  // TPS for DLV
  long long tps;
  // Targeted size of LP in LSR
  long long S;
  // Number of cores used
  int C;
  // Order of table
  double o;
  // Query index
  int q;
  // Seed
  int seed;
  // is_max_var for dlv;
  bool is_max_var;

  // Ranges
  static vector<double> H8, H7, H4, E2, M6, F5, g3;
  static vector<int> Q4, C7, o6, o4;
  static vector<long long> S3;

  // Queries
  static vector<string> datasets;
  static vector<string> obj_cols;
  static vector<bool> is_maximizes;
  static vector<vector<string>> arr_att_cols;
  static vector<vector<int>> arr_att_senses;
  static vector<bool> has_count_constraints;
  static vector<long long> us;
  static vector<string> filtered_cols;
  static vector<double> Es;
  
public:
  ~DetExp();
  DetExp(string out_file, bool verbose=true);
  void reset();
  string getTableName();
  vector<string> getCols();
  DetSql generate(bool is_lazy=true);
  double dlvPartition(bool is_lazy=true);
  double kdPartition(bool is_lazy=true);
  string getKdPartitionName();
  string getDlvPartitionName();
  void write(string id, double x, double y);
  void write(string label, double v);
};
