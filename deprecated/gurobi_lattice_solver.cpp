#define FMT_HEADER_ONLY

#include <iostream>
#include <random>
#include <float.h>

#include "fmt/core.h"
#include "gurobi_lattice_solver.h"
#include "pseudo_walker.h"
#include "walker.h"
#include "utility.h"
#include "omp.h"

using namespace Eigen;
using namespace std;

GurobiLatticeSolver::~GurobiLatticeSolver(){
  GRBfreemodel(model);
  GRBfreeenv(env);
}

GurobiLatticeSolver::GurobiLatticeSolver(MatrixXd Aq, MatrixXd bq, MatrixXd A, VectorXd b, VectorXd c, VectorXd lb, VectorXd ub) {
  this->Aq = Aq;
  this->bq = bq;
  this->A = A;
  this->b = b;
  this->c = c;
  this->lb = lb;
  this->ub = ub;
  n = (int)c.size();
  m = (int)(b.size() + bq.size());
  best_c_score = -DBL_MAX;
  solveRelaxed();
}

GurobiLatticeSolver::GurobiLatticeSolver(MatrixXd A, VectorXd b, VectorXd c, double lb_value, double ub_value) {
  this->A = A;
  this->b = b;
  this->c = c;
  n = (int)c.size();
  m = (int)b.size();
  lb.resize(n); lb.fill(lb_value);
  ub.resize(n); ub.fill(ub_value);
  best_c_score = -DBL_MAX;
  // C++ API For Gurobi
  //GRBEnv env = GRBEnv();
  //GRBModel model = GRBModel(env);
  //for (int i = 0; i < n; i++) {
  //  model.addVar(lb(i), ub(i), c(i), GRB_INTEGER);
  //}
  //model.update();
  //GRBVar* vars = model.getVars();
  //for (int i = 0; i < m; i++) {
  //  GRBLinExpr lhs = GRBLinExpr();
  //  double* coeffs = new double[n];
  //  VectorXd::Map(coeffs, n) = A.row(i);
  //  lhs.addTerms(coeffs, vars, n);
  //  model.addConstr(lhs, GRB_LESS_EQUAL, b(i));
  //}
  //model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
  //model.update();
  solveRelaxed();
}

GurobiLatticeSolver::GurobiLatticeSolver(MatrixXd A, VectorXd b, VectorXd c): GurobiLatticeSolver(A, b, c, 0, DBL_MAX) {
}

// Order of business: Aq first then A
void GurobiLatticeSolver::solveRelaxed() {
  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
  env = NULL;
  model = NULL;
  assert(!GRBemptyenv(&env));
  assert(!GRBsetparam(env, "OutputFlag", "0"));
  assert(!GRBstartenv(env));
  assert(!GRBnewmodel(env, &model, "origin", 0, NULL, NULL, NULL, NULL, NULL));
  for (int i = 0; i < n; i++) {
    assert(!GRBaddvar(model, 0, NULL, NULL, c(i), lb(i), ub(i), GRB_INTEGER, NULL));
  }
  assert(!GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE));
  for (int i = 0; i < Aq.innerSize(); i++) {
    int* ind = new int[n];
    double* val = new double[n];
    for (int j = 0; j < n; j++) {
      ind[j] = j;
      val[j] = Aq(i, j);
    }
    assert(!GRBaddconstr(model, n, ind, val, GRB_EQUAL, bq(i), NULL));
    delete ind;
    delete val;
  }
  for (int i = 0; i < A.innerSize(); i++) {
    int* ind = new int[n];
    double* val = new double[n];
    #pragma omp parallel for schedule(dynamic) num_threads(CORE_COUNT)
    for (int j = 0; j < n; j++) {
      ind[j] = j;
      val[j] = A(i, j);
    }
    assert(!GRBaddconstr(model, n, ind, val, GRB_LESS_EQUAL, b(i), NULL));
    delete ind;
    delete val;
  }
  assert(!GRBsetintparam(env, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF));
  assert(!GRBupdatemodel(model));
  exe_init = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count() / 1000000.0;
  GRBmodel* relaxed;
  //cout << "Start relaxed" << endl;
  assert(!GRBrelaxmodel(model, &relaxed));
  exe_relaxed = exeTime([](GRBmodel* relaxed) {
    assert(!GRBoptimize(relaxed));
    }, relaxed);
  r0.resize(n);
  assert(!GRBgetdblattrarray(relaxed, GRB_DBL_ATTR_X, 0, n, &r0(0)));
  //cout << "Finished relaxed" << endl;

  // Process relaxed solution
  floor_r0.resize(n);
  for (int i = 0; i < n; i ++) floor_r0(i) = floor(r0(i));
  e_score.resize(Aq.innerSize());
  for (int i = 0; i < Aq.innerSize(); i ++)e_score(i) = bq(i) - Aq.row(i).dot(floor_r0);
  i_score.resize(A.innerSize());
  for (int i = 0; i < A.innerSize(); i ++) i_score(i) = b(i) - A.row(i).dot(floor_r0);
  c_score = c.dot(floor_r0);

  start = chrono::high_resolution_clock::now();
  // Process simplex tableau
  int npm = n + m;
  bhead.resize(m);
  assert(!GRBgetBasisHead(relaxed, &bhead[0]));
  vector<int> column_index (npm, -1);
  for (int i = 0; i < m; i ++){
    if (bhead[i] >= 0 && bhead[i] < npm){
      column_index[bhead[i]] = i;
    }
  }
  int dir_count = 0;
  for (int i = 0; i < npm; i ++){
    if (column_index[i] == -1){
      dir_count ++;
      column_index[i] = -dir_count;
    }
  }
  lattice_dirs = CooSparse(n, m);
  degenerate_count = 0;
  //////////////////////////
  tableau.resize(m, m+n); tableau.fill(0);
  for (int i = 0; i < m; i ++){
    GRBsvec v = {npm, new int[npm], new double[npm]};
    assert(!GRBBinvRowi(relaxed, i, &v));
    for (int j = 0; j < v.len; j++){
      tableau(i, v.ind[j]) = v.val[j];
    }
    delete v.ind;
    delete v.val;
  }
  // cout << "---------------------------------" << endl;
  // cout << tableau << endl;
  // cout << "---------------------------------" << endl;
  // for (int i = 0; i < m; i ++) cout << bhead[i] << " ";
  // cout << endl;
  // for (int i = 0; i < npm; i ++) cout << column_index[i] << " ";
  // cout << endl;
  //////////////////////////
  if (kIgnoreDegenerate){
    #pragma omp parallel for schedule(dynamic) num_threads(CORE_COUNT)
    for (int i = 0; i < n; i ++){
      if (column_index[i] >= 0){ // i is basic column and column_index[i] is its row
        GRBsvec v = {npm, new int[npm], new double[npm]};
        #pragma omp critical
        {
          assert(!GRBBinvRowi(relaxed, column_index[i], &v));
        }
        for (int j = 0; j < v.len; j++){
          int ci = column_index[v.ind[j]];
          if (ci < 0) lattice_dirs.addEntry(i, abs(ci)-1, -v.val[j]);
        }
        delete v.ind;
        delete v.val;
      } else{
        lattice_dirs.addEntry(i, abs(column_index[i])-1, 1);
      }
    }
  } else{
    // For degenerate case
    VectorXd r0_vertex(m + n);
    for (int i = 0; i < n; i++) r0_vertex(i) = r0(i);
    for (int i = 0; i < Aq.innerSize(); i++) {
      int index = i + n;
      r0_vertex(index) = bq(i) - Aq.row(i).dot(r0);
    }
    for (int i = 0; i < A.innerSize(); i++) {
      int index = i + n + Aq.innerSize();
      r0_vertex(index) = b(i) - A.row(i).dot(r0);
    }
    for (int i = 0; i < n; i++){
      if (r0_vertex(i) - lb(i) < kFloatEps) r0_vertex(i) = -1;
      else if (ub(i) - r0_vertex(i) < kFloatEps) r0_vertex(i) = 1;
      else r0_vertex(i) = 0;
    }
    for (int i = 0; i < A.innerSize(); i ++){
      int index = i + n + Aq.innerSize();
      if (r0_vertex(index) < kFloatEps) r0_vertex(index) = -1;
      else r0_vertex(index) = 0;
    }
    degenerate_violation.resize(n);
    #pragma omp parallel for schedule(dynamic) num_threads(CORE_COUNT)
    for (int i = 0; i < npm; i ++){
      if (column_index[i] >= 0){
        GRBsvec v = {npm, new int[npm], new double[npm]};
        #pragma omp critical
        {
          assert(!GRBBinvRowi(relaxed, column_index[i], &v));
        }
        for (int j = 0; j < v.len; j++){
          int ci = column_index[v.ind[j]];
          if (ci < 0){
            if (i < n) lattice_dirs.addEntry(i, abs(ci)-1, -v.val[j]);
            if (-v.val[j] * r0_vertex[i] > 0){
              #pragma omp critical
              {
                degenerate_violation[abs(ci) - 1].push_back(i);
              }
            }
          }
        }
        delete v.ind;
        delete v.val;
      } else{
        if (i < n) lattice_dirs.addEntry(i, abs(column_index[i])-1, 1);
      }
    }
    for (int i = 0; i < n; i ++){
      if (degenerate_violation[i].size() > 0) degenerate_count ++;
    }
  }
  exe_tableau = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count() / 1000000.0;
}

GurobiLatticeSolver::GurobiLatticeSolver(string problem) {
  GRBEnv env = GRBEnv();
  env.set("OutputFlag", "0");
  GRBModel model = GRBModel(env, fmt::format("{}/resource/test_cases/{}.mps", kProjectHome, problem));
  GRBConstr* constraints = model.getConstrs();
  GRBVar* vars = model.getVars();
  n = model.get(GRB_IntAttr_NumVars);
  m = model.get(GRB_IntAttr_NumConstrs);
  lb.resize(n); lb.fill(0);
  ub.resize(n); ub.fill(DBL_MAX);
  for (int i = 0; i < n; i++) {
    GRBVar vi = vars[i];
    lb(i) = vi.get(GRB_DoubleAttr_LB);
    ub(i) = vi.get(GRB_DoubleAttr_UB);
  }
  int q_count = 0;
  for (int i = 0; i < m; i++) {
    GRBConstr ci = constraints[i];
    if (ci.get(GRB_CharAttr_Sense) == GRB_EQUAL) q_count++;
  }
  Aq.resize(q_count, n); bq.resize(q_count);
  int count = m - q_count;
  A.resize(count, n); b.resize(count);
  int q_index = 0;
  int index = 0;
  for (int i = 0; i < m; i++) {
    GRBConstr ci = constraints[i];
    if (ci.get(GRB_CharAttr_Sense) == GRB_EQUAL) {
      for (int j = 0; j < n; j++) {
        GRBVar vj = vars[j];
        Aq(q_index, j) = model.getCoeff(ci, vj);
      }
      bq(q_index) = ci.get(GRB_DoubleAttr_RHS);
      q_index++;
    }
    else {
      for (int j = 0; j < n; j++) {
        GRBVar vj = vars[j];
        A(index, j) = model.getCoeff(ci, vj);
      }
      b(index) = ci.get(GRB_DoubleAttr_RHS);
      if (ci.get(GRB_CharAttr_Sense) == GRB_GREATER_EQUAL) {
        A.row(index) *= -1;
        b(index) *= -1;
      }
      index++;
    }
  }
  c.resize(n); c.fill(0);
  auto obj = model.getObjective().getLinExpr();
  for (int i = 0; i < (int)obj.size(); i++) {
    c(obj.getVar(i).index()) = obj.getCoeff(i);
  }
  int sense = model.get(GRB_IntAttr_ModelSense);
  if (sense == GRB_MINIMIZE) c *= -1;
  best_c_score = -DBL_MAX;
  solveRelaxed();
}

double GurobiLatticeSolver::getObjValue(VectorXd sol) {
  return c.dot(sol);
}

bool GurobiLatticeSolver::checkFeasibility(VectorXd sol) {
  for (int i = 0; i < n; i++) {
    if (sol(i) > ub(i) + kFeasTolerance || sol(i) < lb(i) - kFeasTolerance) {
      return false;
    }
  }
  for (int i = 0; i < Aq.innerSize(); i++) {
    if (abs(Aq.row(i).dot(sol) - bq(i)) > kFeasTolerance) {
      return false;
    }
  }
  for (int i = 0; i < A.innerSize(); i++) {
    if (A.row(i).dot(sol) > b(i) + kFeasTolerance) {
      return false;
    }
  }
  return true;
}

bool GurobiLatticeSolver::checkFeasibilityStep(const VectorXd& current_e_score, const VectorXd& current_i_score, const unordered_map<int,int>& bound_map){
  if (bound_map.size() > 0) return false;
  for (int j = 0; j < Aq.innerSize(); j ++){
    if (abs(current_e_score(j)) > kFeasTolerance) return false;
  }
  for (int j = 0; j < A.innerSize(); j ++){
    if (current_i_score(j) < 0) return false;
  }
  return true;
}

bool GurobiLatticeSolver::latticeStep(int d, VectorXd& x, VectorXd& current_e_score, VectorXd& current_i_score, double& current_c_score, unordered_map<int,int>& bound_map){
  int i = abs(d) - 1;
  int s = sign(d);
  x(i) += s;
  if (bound_map.find(i) == bound_map.end()){
    int lb_violation = floor(x(i) - lb(i));
    int ub_violation = ceil(x(i) - ub(i));
    if (lb_violation < 0) bound_map[i] = lb_violation;
    else if (ub_violation > 0) bound_map[i] = ub_violation;
  } else {
    bound_map[i] += s;
    if (bound_map[i] == 0) bound_map.erase(i);
    if (abs(bound_map[i]) > kBoundTol) return false;
  }
  for (int j = 0; j < Aq.innerSize(); j ++) current_e_score(j) -= s*Aq(j, i);
  for (int j = 0; j < A.innerSize(); j ++) current_i_score(j) -= s*A(j, i);
  current_c_score += s*c(i);
  return true;
}

void GurobiLatticeSolver::solve(double time_budget) {
  chrono::high_resolution_clock::time_point ep = chrono::high_resolution_clock::now();
  try_count = 0;
  avg_step_count = 0;
  VectorXd ts (20); ts.fill(0);
  #pragma omp parallel num_threads(80)
  {
    // Local declaration
    int local_try_count = 0;
    int local_step_count = 0;
    double local_best_c_score = -DBL_MAX;
    VectorXd local_ts (20); local_ts.fill(0);
    VectorXd local_best_x;
    //default_random_engine gen {0};
    default_random_engine gen {static_cast<long unsigned int>(time(0))};
    uniform_real_distribution u_dist(0.0, 1.0);

    while (true) {
      if (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - ep).count() / 1000000000.0 > time_budget)
        break;
      chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
      local_try_count ++;
      chrono::high_resolution_clock::time_point t11 = chrono::high_resolution_clock::now();
      vector<double> splits(n-1);
      chrono::high_resolution_clock::time_point t12 = chrono::high_resolution_clock::now();
      for (int i = 0; i < n-1; i++) splits[i] = u_dist(gen);
      chrono::high_resolution_clock::time_point t13 = chrono::high_resolution_clock::now();
      VectorXd coefs = bucketSort(splits);
      chrono::high_resolution_clock::time_point t14 = chrono::high_resolution_clock::now();
      for (int i = n-1; i >= 1; i--) coefs(i) -= coefs(i-1);
      chrono::high_resolution_clock::time_point t15 = chrono::high_resolution_clock::now();
      // if (!kIgnoreDegenerate && degenerate_count < n){
      //   for (int i = 0; i < n; i ++){
      //     // Weird bug when coefs(i) = 0
      //     if (degenerate_violation[i].size() > 0) coefs(i) = kFloatEps;
      //   }
      //   double sum_coefs = coefs.sum();
      //   for (int i = 0; i < n; i ++) coefs(i) /= sum_coefs;
      // }
      VectorXd dir = lattice_dirs.vectorProduct(coefs);
      chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
      AbstractWalker* walker;
      if (n < kLatticeUB) walker = new Walker(dir);
      else if (n < kRetargetUB) walker = new PseudoWalker(dir, true);
      else walker = new PseudoWalker(dir, false);
      chrono::high_resolution_clock::time_point t3 = chrono::high_resolution_clock::now();
      VectorXd x = floor_r0;
      VectorXd current_e_score = e_score;
      VectorXd current_i_score = i_score;
      double current_c_score = c_score;
      unordered_map<int, int> bound_map;
      for (int i = 0; i < n; i ++){
        if (u_dist(gen) < r0(i) - floor_r0(i)){
          latticeStep(i+1, x, current_e_score, current_i_score, current_c_score, bound_map);
        }
      }
      chrono::high_resolution_clock::time_point t4 = chrono::high_resolution_clock::now();
      if (checkFeasibilityStep(current_e_score, current_i_score, bound_map)){
        if (local_best_c_score < current_c_score){
          local_best_c_score = current_c_score;
          local_best_x = x;
        }
      }
      int step_count = 0;
      double stubborn_scale = 1.0;
      chrono::high_resolution_clock::time_point t5 = chrono::high_resolution_clock::now();
      while (true){
        int d = walker->step();
        if (!latticeStep(d, x, current_e_score, current_i_score, current_c_score, bound_map)) break;
        if (checkFeasibilityStep(current_e_score, current_i_score, bound_map)){
          if (local_best_c_score < current_c_score){
            local_best_c_score = current_c_score;
            local_best_x = x;
            stubborn_scale = 1.0;
          }
          if (u_dist(gen) < stubborn_scale * (best_c_score - current_c_score) / best_c_score) break;
        }
        step_count ++;
        stubborn_scale *= kStubbornMultiplier;
        if (step_count > kSizeTol) break;
      }
      delete walker;
      chrono::high_resolution_clock::time_point t6 = chrono::high_resolution_clock::now();
      local_step_count += step_count;
      local_ts(6) += chrono::duration_cast<chrono::nanoseconds>(t2-t1).count()/1000000.0;
      local_ts(7) += chrono::duration_cast<chrono::nanoseconds>(t3-t2).count()/1000000.0;
      local_ts(8) += chrono::duration_cast<chrono::nanoseconds>(t4-t3).count()/1000000.0;
      local_ts(9) += chrono::duration_cast<chrono::nanoseconds>(t5-t4).count()/1000000.0;
      local_ts(10) += chrono::duration_cast<chrono::nanoseconds>(t6-t5).count()/1000000.0;
      local_ts(0) += chrono::duration_cast<chrono::nanoseconds>(t11-t1).count()/1000000.0;
      local_ts(1) += chrono::duration_cast<chrono::nanoseconds>(t12-t11).count()/1000000.0;
      local_ts(2) += chrono::duration_cast<chrono::nanoseconds>(t13-t12).count()/1000000.0;
      local_ts(3) += chrono::duration_cast<chrono::nanoseconds>(t14-t13).count()/1000000.0;
      local_ts(4) += chrono::duration_cast<chrono::nanoseconds>(t15-t14).count()/1000000.0;
      local_ts(5) += chrono::duration_cast<chrono::nanoseconds>(t2-t15).count()/1000000.0;
    }
    #pragma omp critical 
    {
      ts += local_ts;
      if (best_c_score < local_best_c_score){
        best_c_score = local_best_c_score;
        best_x = local_best_x;
      }
      try_count += local_try_count;
      avg_step_count += local_step_count;
    }
  }
  ts /= try_count;
  // cout << "GLSTIME ";
  // print(ts);
  avg_step_count /= try_count;
  status = LS_NOT_FOUND;
  if (best_c_score != -DBL_MAX) status = LS_FOUND;
}

void GurobiLatticeSolver::milpSolve(){
  exe_milp = exeTime([](GRBmodel* model){
    assert(!GRBoptimize(model));
  }, model);
  x0.resize(n);
  assert(!GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, n, &x0(0)));
  milp_score = getObjValue(x0);
}