#include "pop_lsr.h"

#include "pb/core/dual_reducer.h"
#include "pb/core/checker.h"

POPLayeredSketchRefine::~POPLayeredSketchRefine(){
}

POPLayeredSketchRefine::POPLayeredSketchRefine(DetProb &prob){
  auto t0 = std::chrono::high_resolution_clock::now();
  int n = (int) prob.c.size();
  int m = (int) prob.bl.size();
  int df = 100;
  vector<long long> indices (n);
  random_device rd;
  default_random_engine gen (rd());
  iota(indices.begin(), indices.end(), 0);
  shuffle(indices.begin(), indices.end(), gen);

  int nn = ceilDiv(n, df);
  DetProb sketch = DetProb(m, nn);
  sketch.bl = prob.bl;
  sketch.bu = prob.bu;
  sketch.l.fill(0);
  sketch.u.fill(df);
  for (int i = 0; i < nn; i ++){
    int st_index = i*df;
    int fn_index = min((i+1)*df, n)-1;
    for (int j = 0; j < m; j ++){
      ScalarMeanVar smv = ScalarMeanVar();
      for (int k = st_index; k <= fn_index; k ++) smv.add(prob.A(j, indices[k]));
      sketch.A(j, i) = smv.getMean();
    }
    ScalarMeanVar smv = ScalarMeanVar();
    for (int k = st_index; k <= fn_index; k ++) smv.add(prob.c(indices[k]));
    sketch.c(i) = smv.getMean();
  }

  DualReducer dr = DualReducer(4, sketch, true);
  vector<long long> sketch_indices (nn);
  iota(sketch_indices.begin(), sketch_indices.end(), 0);
  shuffle(sketch_indices.begin(), sketch_indices.end(), gen);
  vector<long long> sols;
  for (int i = 0; i < nn; i ++){
    if (dr.lp_sol(i) > 0 || dr.ilp_sol(i) > 0) sols.push_back(i);
  }
  int need_size = ceilDiv(nn, df);
  for (int i = 0; i < nn; i ++){
    if (dr.lp_sol(sketch_indices[i]) > 0 || dr.ilp_sol(sketch_indices[i]) > 0) continue;
    sols.push_back(sketch_indices[i]);
    if ((int) sols.size() == need_size) break;
  }

  int nnn = 0;
  for (int i = 0; i < (int) sols.size(); i ++){
    int st_index = sols[i]*df;
    int fn_index = min((int)((sols[i]+1)*df), n)-1;
    nnn += fn_index - st_index + 1;
  }
  DetProb refine = DetProb(m, nnn);
  refine.bl = prob.bl;
  refine.bu = prob.bu;
  int cur_index = 0;
  for (int i = 0; i < (int) sols.size(); i ++){
    int st_index = sols[i]*df;
    int fn_index = min((int)((sols[i]+1)*df), n)-1;
    for (int k = st_index; k <= fn_index; k ++){
      for (int j = 0; j < m; j ++){
        refine.A(j, cur_index) = prob.A(j, indices[k]);
      }
      refine.c(cur_index) = prob.c(indices[k]);
      refine.ids[cur_index] = indices[k];
      cur_index ++;
    }
  }
  DualReducer dr2 = DualReducer(4, refine, true);
  status = dr2.status;
  if (dr2.status == Found){
    lp_score = dr2.lp_score;
    ilp_score = dr2.ilp_score;
    for (int i = 0; i < nnn; i ++){
      if (dr2.lp_sol(i) > 0){
        lp_sol[refine.ids[i]] = dr2.lp_sol(i);
      }
      if (dr2.ilp_sol(i) > 0){
        ilp_sol[refine.ids[i]] = dr2.ilp_sol(i);
      }
    }
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  exe_ilp = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() / 1000000.0;
}