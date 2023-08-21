#include "lp_lsr.h"

#include "pb/core/dual_reducer.h"
#include "pb/core/checker.h"

const int kSleepPeriod = 25; // In ms

const double kMinGapOpt = 1e-4;
const double kMinGap = 1e-1;

void LPLayeredSketchRefine::init(){
  pg = new PgManager();
  _conn = PQconnectdb(pg->conninfo.c_str());
  assert(PQstatus(_conn) == CONNECTION_OK);
  _res = NULL;
  status = NotFound;
  INIT_CLOCK(pro);
  exe_gb = exe_dual = 0;
}

DLVPartition* LPLayeredSketchRefine::getDLVPartition(LsrProb* prob){
  _sql = fmt::format("SELECT cols, group_ratio, tps, layer_count FROM {} WHERE table_name='{}' AND partition_name='{}';", kPartitionTable, prob->det_sql.table_name, prob->partition_name);
  _res = PQexec(_conn, _sql.c_str());
  if (PQntuples(_res)){
    // cout << "OKPRE\n"; 
    DLVPartition* partition = new DLVPartition(prob, pgStringSplit(PQgetvalue(_res, 0, 0)), 
      atof(PQgetvalue(_res, 0, 1)), atoll(PQgetvalue(_res, 0, 2)), atoi(PQgetvalue(_res, 0, 3)), true);
    // cout << "OKAFT\n"; 
    PQclear(_res);
    return partition;
  }
  return NULL; 
}

LPLayeredSketchRefine::~LPLayeredSketchRefine(){
  PQfinish(_conn);
  if (partition) delete partition;
  delete pg;
}

void LPLayeredSketchRefine::formulateDetProb(int core, LsrProb &prob, DetProb &det_prob, string current_gtable, const vector<long long> &ids){
  int m = prob.det_sql.att_cols.size() + 1;
  int n = (int) ids.size();
  det_prob.resize(m, n);
  det_prob.u.fill((double) prob.det_sql.u);
  for (int i = 0; i < m; i ++){
    if (i < m - 1){
      det_prob.bl(i) = prob.bl(i);
      det_prob.bu(i) = prob.bu(i);
    } else{
      det_prob.bl(i) = prob.cl;
      det_prob.bu(i) = prob.cu;
    }
  }
  int chunk = ceilDiv(n, core);
  int layer = getLayerIndex(current_gtable);
  string col_names = join(prob.det_sql.att_cols, ",");
  #pragma omp parallel num_threads(core)
  {
    CREATE_CLOCK(local_pro);
    string sql;
    PGconn* conn = PQconnectdb(pg->conninfo.c_str());
    assert(PQstatus(conn) == CONNECTION_OK);
    PGresult *res = NULL;
    int seg = omp_get_thread_num();
    int start_id = seg * chunk;
    int end_id = min((seg + 1) * chunk - 1, n - 1);
    string string_ids = "";
    for (int i = start_id; i < end_id; i ++) string_ids += to_string(ids[i]) + ",";
    if (end_id >= start_id) {
      string_ids += to_string(ids[end_id]);
      if (layer > 0){
        sql = fmt::format("SELECT {},{},{},size FROM \"{}\" WHERE {} IN ({});",
          kId, prob.det_sql.obj_col, col_names, current_gtable, kId, string_ids);
      } else{
        sql = fmt::format("SELECT {},{},{} FROM \"{}\" WHERE {} IN ({});",
          kId, prob.det_sql.obj_col, col_names, current_gtable, kId, string_ids);
      }
      START_CLOCK(local_pro, 0);
      res = PQexec(conn, sql.c_str());
      END_CLOCK(local_pro, 0);
      for (int i = 0; i < PQntuples(res); i++){
        int index = i + start_id;
        det_prob.ids[index] = atoll(PQgetvalue(res, i, 0));
        det_prob.c(index) = atof(PQgetvalue(res, i, 1));
        for (int j = 0; j < m-1; j ++){
          det_prob.A(j, index) = atof(PQgetvalue(res, i, 2+j)); 
        }
        if (layer > 0) det_prob.u(index) = atof(PQgetvalue(res, i, m+1));
        det_prob.A(m-1, index) = 1.0;
      }
      ADD_CLOCK(local_pro, core);
    }
    PQclear(res);
    PQfinish(conn);
  }
  if (!prob.det_sql.is_maximize) det_prob.c = -det_prob.c;
  // Remove -inf <= Ax <= inf if exists
  det_prob.truncate();
}

LPLayeredSketchRefine::LPLayeredSketchRefine(int core, LsrProb &prob, long long lp_size, bool is_safe): lp_size(lp_size){
  init();
  std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
  partition = getDLVPartition(&prob);
  // cout << "OKWTF1\n";
  if (!partition){
    if (prob.det_sql.table_size <= lp_size){
      DetProb det_prob = DetProb(prob.det_sql, -1, kGlobalSeed);
      det_prob.copyBounds(prob.bl, prob.bu, prob.cl, prob.cu);

      DualReducer dr = DualReducer(core, det_prob, is_safe, kMinGapOpt, kTimeLimit);
      exe_gb += dr.exe_gb;
      exe_dual += dr.exe_lp;
      exe_lp = dr.exe_lp;
      exe_ilp = dr.exe_ilp;
      status = dr.status;
      if (status == Found){
        ilp_score = dr.ilp_score;
        lp_score = dr.lp_score;
        for (int i = 0; i < (int) det_prob.c.size(); i ++){
          if (isGreater(dr.ilp_sol(i), 0)){
            ilp_sol[det_prob.ids[i]] = (long long) dr.ilp_sol(i);
          }
          if (isGreater(dr.lp_sol(i), 0)){
            lp_sol[det_prob.ids[i]] = dr.lp_sol(i);
          }
        }
      }
      exe_ilp = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000000.0;
      exe_lp = exe_ilp - (dr.exe_ilp - dr.exe_lp);
      return;
    }
    status = NoPartitionFound;
    return;
  }
  if (!partition->isCompatible()){
    status = IncompatiblePartition;
    return;
  }

  // cout << "OKWTF\n";

  double limit_size_per_group = DBL_MAX;
  if (isLess(kOutlierPercentage, 1.0)){
    limit_size_per_group = 
      (1 - kOutlierPercentage*partition->group_ratio)
      /(partition->group_ratio*(1-kOutlierPercentage));
  }

  // Phase-0: Filtering
  string current_gtable = partition->getInitialGName(partition->layer_count);
  int n = pg->getSize(current_gtable);
  vector<long long> last_layer_ids;
  // cout << "OK-1\n";
  if (partition->is_filtering){
    // Filtering
    #if DEBUG
      MeanVar global_mv = MeanVar(partition->query_cols.size());
    #endif
    string gname;
    long long group_count, chunk;
    int start_index = 0;
    // cout << "OK0\n";
    #pragma omp parallel num_threads(core)
    {
      string sql;
      PGconn *conn = PQconnectdb(pg->conninfo.c_str());
      assert(PQstatus(conn) == CONNECTION_OK);
      PGresult *res = NULL;
      DLVPartition *loc_partition = new DLVPartition(&prob, partition->cols, partition->group_ratio, partition->tps, partition->layer_count);
      
      int local_start_index = 0;
      vector<long long> local_ids;
      local_ids.reserve(n);

      for (int layer = 1; layer <= partition->layer_count; layer ++){
        #pragma omp master
        {
          gname = partition->getGName(layer);
          group_count = pg->getSize(partition->getInitialGName(layer));
          chunk = ceilDiv(group_count, (long long) core);
        }
        #pragma omp barrier
        int seg = omp_get_thread_num();
        long long start_id = seg * chunk + 1;
        long long end_id = min((seg + 1) * chunk, group_count);
        if (end_id >= start_id){
          sql = fmt::format("COPY \"{}\" FROM STDIN with(delimiter ',');", gname);
          res = PQexec(conn, sql.c_str());
          assert(PQresultStatus(res) == PGRES_COPY_IN);
          PQclear(res);
          if (layer == 1){
            // Fast filtering if the intervals lie outside the filtering conditions
            unordered_map<string, int> cols_map;
            for (int i = 0; i < (int) partition->cols.size(); i ++){
              cols_map[partition->cols[i]] = i;
            }
            for (long long group_id = start_id; group_id <= end_id; group_id ++){
              vector<pair<double, double>> intervals = loc_partition->getGroupIntervals(layer, group_id);
              bool is_filtered = false;
              for (int j = 0; j < (int) prob.det_sql.filter_cols.size(); j++){
                int col_index = cols_map[prob.det_sql.filter_cols[j]];
                if (!doesIntersect(intervals[col_index], prob.det_sql.filter_intervals[j])){
                  is_filtered = true;
                  break;
                }
              }
              if (!is_filtered){
                // Go through group and check individual original tuples
                MeanVar mv = MeanVar(loc_partition->query_cols.size());
                loc_partition->getGroupFilteredStat(mv, layer, group_id);
                if (mv.sample_count){
                  string data = fmt::format("{},{}\n", group_id, join(mv.getVar(), kPrecision));
                  // cout << data << " " << sql << endl;
                  assert(PQputCopyData(conn, data.c_str(), (int) data.length()) == 1);
                }
                #if DEBUG
                #pragma omp critical (c6)
                {
                  global_mv.add(mv);
                }
                #endif
                if (layer == partition->layer_count){
                  local_ids.push_back(group_id);
                }
              }
            }
          } else{
            // Normal filtering since the intervals are no longer correct for original tuples
            for (long long group_id = start_id; group_id <= end_id; group_id ++){
              MeanVar mv = MeanVar(loc_partition->query_cols.size());
              loc_partition->getGroupFilteredStat(mv, layer, group_id);
              if (mv.sample_count){
                string data = fmt::format("{},{}\n", group_id, join(mv.getVar(), kPrecision));
                assert(PQputCopyData(conn, data.c_str(), (int) data.length()) == 1);
              }
              if (layer == partition->layer_count){
                local_ids.push_back(group_id);
              }
            }
          }
          assert(PQputCopyEnd(conn, NULL) == 1);
          res = PQgetResult(conn);
          assert(PQresultStatus(res) == PGRES_COMMAND_OK);
          PQclear(res);
        }
        #pragma omp barrier
        #pragma omp master
        {
          _sql = fmt::format("CREATE INDEX \"{}_id\" ON \"{}\" USING btree ({});", gname, gname, kId);
          _res = PQexec(_conn, _sql.c_str());
          PQclear(_res);
        }
        #pragma omp barrier
        if (layer < partition->layer_count) loc_partition->prepareGroupFilteredStatSql(layer+1);
        #pragma omp barrier
      }
      PQfinish(conn);
      delete loc_partition;
      #pragma omp critical (c7)
      {
        local_start_index = start_index;
        start_index += local_ids.size();
      }
      #pragma omp barrier
      #pragma omp master
      {
        last_layer_ids.resize(start_index);
      }
      #pragma omp barrier
      if (local_ids.size()){
        memcpy(&(last_layer_ids[local_start_index]), &(local_ids[0]), local_ids.size()*sizeof(long long));
      }
    }
    #if DEBUG
      VectorXd means = global_mv.getMean();
      means.conservativeResize(prob.det_sql.att_cols.size());
      VectorXd vars = global_mv.getVar();
      vars.conservativeResize(prob.det_sql.att_cols.size());
      prob.rebound(means, vars);
    #endif
  } else{
    last_layer_ids.resize(n);
    iota(last_layer_ids.begin(), last_layer_ids.end(), 1);
  }
  // cout << "OK1\n"; 
  DetProb det_prob;
  // Phase-1: Retrieve last layer
  formulateDetProb(core, prob, det_prob, current_gtable, last_layer_ids);
  // Phase-2: Iteratively go to next layer
  for (int layer = partition->layer_count; layer >= 1; layer --){
    // // Make sure that we augment less tuples for the original layer
    // if (layer == 1) lp_size = kOptLpSize;

    int n = (int) det_prob.c.size();
    // Phase-2a: Sketch
    Dual dual = Dual(core, det_prob);
    exe_dual += dual.exe_solve;
    if (dual.status != Found){
      status = dual.status;
      return;
    }
    
    vector<long long> ids;
    double E = 0;
    for (int i = 0; i < n; i ++){
      E += dual.sol(i);
    }
    double target_size = lp_size * partition->group_ratio * 2;
    DetProb new_prob = det_prob;
    new_prob.u.fill(E / target_size);
    Dual scaled_dual = Dual(core, new_prob);

    unordered_set<long long> group_ids_set;
    long long total_size = 0;
    long long global_start_index = 0;
    #pragma omp parallel num_threads(core)
    {
      DLVPartition *loc_partition = new DLVPartition(&prob, partition->cols, partition->group_ratio, partition->tps, partition->layer_count);
      vector<long long> loc_ids; loc_ids.reserve(lp_size);
      #pragma omp for
      for (int i = 0; i < n; i ++){
        // Condition for sketch
        if (isGreater(scaled_dual.sol(i), 0)){
          long long group_id = det_prob.ids[i];
          if (total_size <= lp_size){
            long long size = loc_partition->getGroupComp(loc_ids, layer, group_id, limit_size_per_group);
            #pragma omp atomic
            total_size += size;
          }
        }
      }

      // Debug
      // #pragma omp critical
      // {
      //   cout << "SIZE: " << loc_ids.size() << endl;
      //   for (auto id : loc_ids){
      //     cout << id << " ";
      //   }
      //   cout << endl;
      // }

      // Phase-2c: Refine
      long long start_index = 0;
      #pragma omp critical (c5)
      {
        start_index = global_start_index;
        global_start_index += loc_ids.size();
      }
      #pragma omp barrier
      #pragma omp master
      {
        ids.resize(total_size);
      }
      #pragma omp barrier
      if (loc_ids.size()) memcpy(&(ids[start_index]), &(loc_ids[0]), loc_ids.size()*sizeof(long long));
      ADD_CLOCK(loc_partition->pro, core);
      delete loc_partition;
    }

    string next_gtable;
    if (layer > 1) next_gtable = partition->getGName(layer-1);
    else next_gtable = prob.det_sql.table_name;
    formulateDetProb(core, prob, det_prob, next_gtable, ids);
    // Debug
    // for (auto v : ids) cout << v << " ";
    // cout << endl;
    // cout << ids.size() << endl;
  }
  if (partition->is_filtering){
    for (int layer = 1; layer <= partition->layer_count; layer ++){
      string gname = partition->getGName(layer);
      pg->dropTable(gname);
    }
  }
  // Phase-3: Final solution
  DualReducer dr = DualReducer(core, det_prob, is_safe, kMinGapOpt, kTimeLimit);
  exe_gb += dr.exe_gb;
  exe_dual += dr.exe_lp;
  if (dr.status != Found){
    status = dr.status;
    return;
  }
  // Checker ch = Checker(det_prob);
  // cout << solMessage(dr.status) << " " << feasMessage(ch.checkIlpFeasibility(dr.ilp_sol)) << " " << feasMessage(ch.checkLpFeasibility(dr.lp_sol)) << endl;
  status = Found;
  ilp_score = lp_score = 0;
  for (int i = 0; i < (int) det_prob.c.size(); i ++){
    if (isGreater(dr.ilp_sol(i), 0)){
      ilp_sol[det_prob.ids[i]] = (long long) dr.ilp_sol(i);
      ilp_score += ilp_sol[det_prob.ids[i]]*det_prob.c(i);
    }
    if (isGreater(dr.lp_sol(i), 0)){
      lp_sol[det_prob.ids[i]] = dr.lp_sol(i);
      lp_score += lp_sol[det_prob.ids[i]]*det_prob.c(i);
    }
  }
  LsrChecker ch = LsrChecker(prob);
  if (ch.checkLpFeasibility(lp_sol) != Feasibility || ch.checkIlpFeasibility(ilp_sol) != Feasibility){
    // cout << feasMessage(ch.checkLpFeasibility(lp_sol)) << endl;
    // cout << feasMessage(ch.checkIlpFeasibility(ilp_sol)) << endl;
    status = NotFound;
  }
  exe_ilp = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000000.0;
  exe_lp = exe_ilp - (dr.exe_ilp - dr.exe_lp);
}