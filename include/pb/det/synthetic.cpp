#define FMT_HEADER_ONLY

#include "fmt/core.h"
#include "synthetic.h"
#include "libpq-fe.h"
#include "pb/util/uconfig.h"
#include "pb/util/unumeric.h"
#include "pb/util/udeclare.h"
#include "pb/util/umisc.h"

using namespace pb;

string Synthetic::table_name = "test_data";

string getSubtableName(string table_name, double order, int seed){
  return fmt::format("{}_{}_{}", table_name, formatFloat(order), seed);
}

Synthetic::~Synthetic(){
  delete pg;
}

Synthetic::Synthetic(){
  vector<string> names = {"Init", "Create table", "Generate data", "Create index"};
  pro = Profiler(names);
  pro.clock(0);
  init();
  pro.stop(0);
}

void Synthetic::init(){
  pg = new PgManager();
}

void Synthetic::createSubtable(string table_name, double order, int seed){
  auto cols = pg->getNumericCols(table_name);
  createSubtable(table_name, order, cols, seed);
}

void Synthetic::createSubtable(string table_name, double order, vector<string> cols, int seed){
  pro.clock(1);
  string subtable_name = getSubtableName(table_name, order, seed);
  {
    string sql;
    PGconn* conn = PQconnectdb(pg->conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    sql = fmt::format("DROP TABLE IF EXISTS \"{}\" CASCADE;", subtable_name);
    res = PQexec(conn, sql.c_str());
    PQclear(res);

    string inject = ", {} DOUBLE PRECISION";
    string declare_inject = "";
    for (int i = 0; i < (int) cols.size(); i ++) declare_inject += fmt::format(inject, cols[i]);

    sql = fmt::format(""
      "CREATE TABLE IF NOT EXISTS \"{}\" ("
      "	{} BIGINT"
      " {}"
      ");", subtable_name, kId, declare_inject);
    res = PQexec(conn, sql.c_str());
    PQclear(res);
    PQfinish(conn);
  }
  pro.stop(1);
  pro.clock(2);
  long long size = pg->getSize(table_name);
  long long chunk = ceilDiv(size, (long long) kPCore);
  double n = pow(10.0, order);
  double probability = n / size;
  unsigned int local_seed;
  if (seed < 0){
    random_device rd;
    local_seed = rd();
  } else{
    local_seed = seed;
  }
  seed_seq seq{local_seed};
  vector<unsigned int> local_seeds (kPCore);
  seq.generate(local_seeds.begin(), local_seeds.end());
  string col_name = join(cols, ",");
  long long global_size = 0;

  #pragma omp parallel num_threads(kPCore)
  {
    PGconn *conn = PQconnectdb(pg->conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGconn *_conn = PQconnectdb(pg->conninfo.c_str());
    assert(PQstatus(_conn) == CONNECTION_OK);

    int seg = omp_get_thread_num();
    default_random_engine gen (local_seeds[seg]);
    uniform_real_distribution dist (0.0, 1.0);
    long long start_id = seg * chunk + 1;
    long long end_id = (seg + 1) * chunk;
    long long slide_size = getTupleCount(2*cols.size()+4);
    long long slide_count = ceilDiv(end_id - start_id + 1, slide_size);
    for (long long sl = 0; sl < slide_count; sl ++){
      // #pragma omp critical
      // {
      //   cout << "Thread: " << omp_get_thread_num() << " starts slide " << sl << " with size " << slide_size << endl;
      // }
      long long slide_start_id = start_id + sl * slide_size;
      long long slide_end_id = min(start_id + (sl + 1) * slide_size - 1, end_id);
      vector<long long> selected_ids; selected_ids.reserve((int)(slide_size * probability));
      for (long long id = slide_start_id; id <= slide_end_id; id ++){
        if (dist(gen) <= probability) selected_ids.push_back(id);
      }
      long long start_index;
      #pragma omp critical
      {
        start_index = global_size;
        global_size += selected_ids.size();
      }
      #pragma omp barrier
      vector<string> sids;
      for (auto id : selected_ids) sids.push_back(to_string(id));
      
      PGresult *res = NULL;
      PGresult *_res = NULL;

      string sql = fmt::format("SELECT {} FROM \"{}\" WHERE {} IN ({})", col_name, table_name, kId, join(sids, ","));
      res = PQexec(conn, sql.c_str());

      string _sql = fmt::format("COPY \"{}\" FROM STDIN with(delimiter ',');", subtable_name);
      _res = PQexec(_conn, _sql.c_str());
      assert(PQresultStatus(_res) == PGRES_COPY_IN);
      PQclear(_res);

      VectorXd vals (cols.size());
      for (long long i = 0; i < PQntuples(res); i++){
        long long index = i + start_index;
        for (int j = 0; j < (int) cols.size(); j ++) vals(j) = atof(PQgetvalue(res, i, j));
        string data = fmt::format("{},{}\n", index, join(vals, kPrecision));
        assert(PQputCopyData(_conn, data.c_str(), (int) data.length()) == 1);
      }
      assert(PQputCopyEnd(_conn, NULL) == 1);
      _res = PQgetResult(_conn);
      assert(PQresultStatus(_res) == PGRES_COMMAND_OK);
      PQclear(_res);
      PQclear(res);
      // #pragma omp critical
      // {
      //   cout << "Thread: " << omp_get_thread_num() << " finishs slide " << sl << endl;
      // }
    }
    PQfinish(_conn);
    PQfinish(conn);
  }
  pro.stop(2);
  pro.clock(3);
  {
    string sql;
    PGconn* conn = PQconnectdb(pg->conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    sql = fmt::format("ALTER TABLE \"{}\" ADD PRIMARY KEY ({});", subtable_name, kId);
    res = PQexec(conn, sql.c_str());
    PQclear(res);

    PQfinish(conn);
  }
  pro.stop(3);
}

void Synthetic::createUniform(long long N, int count, double mean, double var){
  vector<double> means (count, mean);  
  vector<double> vars (count, var);
  create(N, count, 0, means, vars);
}

void Synthetic::createNormal(long long N, int count, double mean, double var){
  vector<double> means (count, mean);  
  vector<double> vars (count, var);
  create(N, 0, count, means, vars);
}

void Synthetic::createMixed(long long N, int ucount, int ncount, double mean, double var){
  int count = ucount + ncount;
  vector<double> means (count, mean);  
  vector<double> vars (count, var);
  create(N, ucount, ncount, means, vars);
}

void Synthetic::create(long long N, int ucount, int ncount, vector<double> means, vector<double> vars){
  pro.clock(1);
  string conninfo = fmt::format("postgresql://{}@{}?port={}&dbname={}&password={}", kPgUser, kPgHostaddr, kPgPort, kPgDatabase, kPgPassword);
  {
    string sql;
    PGconn* conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    sql = fmt::format("DROP TABLE IF EXISTS \"{}\" CASCADE;", Synthetic::table_name);
    res = PQexec(conn, sql.c_str());
    PQclear(res);

    string inject1 = ", u{} DOUBLE PRECISION";
    string inject2 = ", n{} DOUBLE PRECISION";
    string declare_inject = "";
    for (int i = 1; i <= ucount; i ++) declare_inject += fmt::format(inject1, i);
    for (int i = 1; i <= ncount; i ++) declare_inject += fmt::format(inject2, i);

    sql = fmt::format(""
      "CREATE TABLE IF NOT EXISTS \"{}\" ("
      "	{} BIGINT"
      " {}"
      ");", Synthetic::table_name, kId, declare_inject);
    res = PQexec(conn, sql.c_str());
    PQclear(res);
    PQfinish(conn);
  }
  pro.stop(1);

  pro.clock(2);
  long long chunk = ceilDiv(N, (long long) kPCore);
  #pragma omp parallel num_threads(kPCore)
  {
    int seg = omp_get_thread_num();
    std::default_random_engine gen;
    gen.seed(seg);
    std::uniform_real_distribution<double> u_dist(-sqrt(3), sqrt(3));
    std::normal_distribution<double> n_dist(0.0, 1.0);
    long long start_id = seg * chunk + 1;
    long long end_id = min((seg + 1) * chunk, N);

    PGconn *conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    string sql = fmt::format("COPY \"{}\" FROM STDIN with(delimiter ',');", Synthetic::table_name);
    res = PQexec(conn, sql.c_str());
    assert(PQresultStatus(res) == PGRES_COPY_IN);
    PQclear(res);
    VectorXd x (ucount + ncount);
    for (long long id=start_id; id <= end_id; id ++){
      for (int i = 0; i < ucount; i ++) x(i) = u_dist(gen) * sqrt(vars[i]) + means[i];
      for (int i = ucount; i < ucount + ncount; i ++) x(i) = n_dist(gen) * sqrt(vars[i]) + means[i];
      string data = fmt::format("{},{}\n", id, join(x, kPrecision));
      assert(PQputCopyData(conn, data.c_str(), (int) data.length()) == 1);
    }
    assert(PQputCopyEnd(conn, NULL) == 1);
    res = PQgetResult(conn);
    assert(PQresultStatus(res) == PGRES_COMMAND_OK);
    PQclear(res);
    PQfinish(conn);
  }
  pro.stop(2);

  pro.clock(3);
  {
    string sql;
    PGconn* conn = PQconnectdb(conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;

    sql = fmt::format("ALTER TABLE \"{}\" ADD PRIMARY KEY ({});", Synthetic::table_name, kId);
    res = PQexec(conn, sql.c_str());
    PQclear(res);

    PQfinish(conn);
  }
  pro.stop(3);
}