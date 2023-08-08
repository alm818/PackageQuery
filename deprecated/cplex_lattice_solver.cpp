#define FMT_HEADER_ONLY

#include <iostream>
#include <random>
#include <float.h>

#include "fmt/core.h"
#include "cplex_lattice_solver.h"
#include "pseudo_walker.h"
#include "walker.h"
#include "utility.h"
#include "omp.h"

using namespace Eigen;
using namespace std;

CplexLatticeSolver::~CplexLatticeSolver(){
  assert(!CPXfreeprob(env, &model));
  assert(!CPXcloseCPLEX(&env));
}

VectorXd CplexLatticeSolver::uncrushX(VectorXd pre_x){
  if (!is_presolvable) return pre_x;
  VectorXd x (numcols);
  assert(!CPXXuncrushx(env, model, &x(0), &pre_x(0)));
  return x;
}

double CplexLatticeSolver::getObjValue(VectorXd sol){
  return obj.dot(sol);
}

bool CplexLatticeSolver::checkFeasibility(VectorXd sol){
  VectorXd feas (numrows);
  assert(!CPXXgetrowinfeas(env, model, &sol(0), &feas(0), 0, numrows-1));
  for (int i = 0; i < numrows; i ++){
    if (abs(feas(i)) > kFeasTolerance) return false;
  }
  return true;
}

void CplexLatticeSolver::initialize(MatrixXd A, VectorXd b, VectorXd c, VectorXd lb, VectorXd ub){
  int status = 0;
  env = CPXXopenCPLEX(&status); assert(!status);
  string problem = "ls";
  model = CPXXcreateprob (env, &status, problem.c_str()); assert(!status);
  numcols = c.size();
  numrows = b.size();
  int objsen = CPX_MAX;
  obj = c;
  char* sense = new char[numrows];
  fill(sense, sense + numrows, CPX_LESS_EQUAL);

  long long* matbeg = new long long[numcols];
  int* matcnt = new int[numcols];
  for (int i = 0; i < numcols; i ++){
    matbeg[i] = i * (long long)numrows;
    matcnt[i] = numrows;
  }

  int* matind = new int[numcols * numrows];
  double* matval = new double[numcols * numrows];
  for (int i = 0; i < numrows; i ++){
    for (int j = 0; j < numcols; j ++){
      int index = j*numrows + i;
      matind[index] = i;
      matval[index] = A(i, j);
    }
  }

  assert(!CPXXcopylp(env, model, numcols, numrows, objsen, &obj(0), &b(0), sense, matbeg, matcnt, matind, matval, &lb(0), &ub(0), NULL));

  char* ctype = new char[numcols];
  fill(ctype, ctype + numcols, CPX_INTEGER);
  assert(!CPXXcopyctype(env, model, ctype));

  delete[] sense; delete[] ctype;
  delete[] matbeg; delete[] matcnt; delete[] matind; delete[] matval;
  
  solveRelaxed();
}

void CplexLatticeSolver::extractModel(CPXCLPptr model){
  // We need to get info for our n, m, A, b, c, lb, ub
  n = CPXXgetnumcols(env, model);

  c.resize(n);
  assert(!CPXXgetobj(env, model, &c(0), 0, n-1));
  if (CPXXgetobjsen(env, model) == CPX_MIN) c *= -1;

  lb.resize(n);
  assert(!CPXXgetlb(env, model, &lb(0), 0, n-1));

  ub.resize(n);
  assert(!CPXXgetub(env, model, &ub(0), 0, n-1));

  int model_m = CPXXgetnumrows(env, model);
  m = model_m;
  char* sense = new char[model_m];
  assert(!CPXXgetsense(env, model, sense, 0, model_m-1));
  for (int i = 0; i < model_m; i ++){
    if (sense[i] == CPX_EQUAL || sense[i] == CPX_RANGE) m ++;
  }

  b.resize(m);
  A.resize(m, n);
  double* b_ptr = new double[model_m];
  assert(!CPXXgetrhs(env, model, b_ptr, 0, model_m-1));
  double* rngval = new double[model_m];
  assert(!CPXXgetrngval(env, model, rngval, 0, model_m-1));
  int index = 0;
  for (int i = 0; i < model_m; i ++){
    if (sense[i] == CPX_LESS_EQUAL || sense[i] == CPX_GREATER_EQUAL){
      b(index) = b_ptr[i];
      for (int j = 0; j < n; j ++) assert(!CPXXgetcoef(env, model, i, j, &A(index, j)));
      if (sense[i] == CPX_GREATER_EQUAL){
        b(index) *= -1;
        A.row(index) *= -1;
      }
      index ++;
    } 
  }
  for (int i = 0; i < model_m; i ++){
    if (sense[i] == CPX_EQUAL || sense[i] == CPX_RANGE){
      if (sense[i] == CPX_EQUAL){
        b(index) = b_ptr[i];
        b(index+1) = -b_ptr[i];
        for (int j = 0; j < n; j ++){
          double val = 0;
          assert(!CPXXgetcoef(env, model, i, j, &val));
          A(index, j) = val;
          A(index+1, j) = -val;
        }
        index += 2;
      } else{
        b(index) = -min(b_ptr[i], b_ptr[i] + rngval[i]);
        b(index+1) = max(b_ptr[i], b_ptr[i] + rngval[i]);
        for (int j = 0; j < n; j ++){
          double val = 0;
          assert(!CPXXgetcoef(env, model, i, j, &val));
          A(index, j) = -val;
          A(index+1, j) = val;
        }
        index += 2;
      }
    }
  }
  delete[] b_ptr; delete[] rngval; delete[] sense;
}

void CplexLatticeSolver::solveRelaxed(){
  exe_tableau = 0;
  exe_presolve = 0;
  exe_relaxed = 0;
  // Debug
  // assert(!CPXXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON));

  // Presolve
  int status;
  CPXCLPptr redmodel = NULL;
  if (kDoPresolve){
    sol_status = LS_NOT_FOUND;
    assert(!CPXXsetintparam(env, CPXPARAM_Preprocessing_Dependency, 3));
    assert(!CPXXsetintparam(env, CPXPARAM_Preprocessing_CoeffReduce, 0));
    assert(!CPXXsetintparam(env, CPXPARAM_Preprocessing_BoundStrength, 1));
    exe_presolve = exeTime([](CPXENVptr& env, CPXLPptr& model, int& sol_status){
      int CPX_status = CPXXpresolve(env, model, CPX_ALG_PRIMAL);
      if (CPX_status == CPXERR_PRESLV_INF) sol_status = LS_INFEASIBLE;
      else if (CPX_status == CPXERR_PRESLV_UNBD) sol_status = LS_UNBOUNDED;
    }, env, model, sol_status);
    CPXXgetredlp(env, model, &redmodel);
  }

  // Relaxed Presolve
  is_presolvable = redmodel != NULL;
  CPXLPptr presolved_model;
  if (is_presolvable) presolved_model = CPXXcloneprob(env, redmodel, &status);
  else presolved_model = CPXXcloneprob(env, model, &status);
  assert(!status);
  extractModel(presolved_model);

  assert(!CPXXchgprobtype(env, presolved_model, CPXPROB_LP));
  exe_relaxed = exeTime([](CPXENVptr& env, CPXLPptr& presolved_model){
    assert(!CPXXlpopt(env, presolved_model));
  }, env, presolved_model);
  // Process relaxed solution
  best_c_score = -DBL_MAX;
  r0.resize(n);
  assert(!CPXXgetx(env, presolved_model, &r0(0), 0, n-1));
  uncrushed_r0 = uncrushX(r0);
  floor_r0.resize(n);
  for (int i = 0; i < n; i ++) floor_r0(i) = floor(r0(i));
  feas_scores.resize(m);
  for (int i = 0; i < m; i ++) feas_scores(i) = b(i) - A.row(i).dot(floor_r0);
  c_score = c.dot(floor_r0);

  // Process simplex tableau
  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
  int pre_m = CPXXgetnumrows(env, presolved_model);
  int* bhead = new int[pre_m];
  assert(!CPXXgetbhead(env, presolved_model, bhead, NULL));
  int npm = n + pre_m;
  vector<int> column_index (npm, -1);
  for (int i = 0; i < pre_m; i ++){
    if  (bhead[i] >= 0) column_index[bhead[i]] = i;
    else column_index[abs(bhead[i]) - 1 + n] = i;
  }
  delete[] bhead;
  int dir_count = 0;
  for (int i = 0; i < npm; i ++){
    if (column_index[i] == -1){
      dir_count ++;
      column_index[i] = -dir_count;
    }
  }
  lattice_dirs = CooSparse(n, pre_m);
  degenerate_count = 0;
  if (kIgnoreDegenerate){
    for (int i = 0; i < npm; i ++){
      if (column_index[i] < 0){
        if (i < n){
          double* col = new double[pre_m];
          assert(!CPXXbinvacol(env, presolved_model, i, col));
          for (int j = 0; j < n; j ++){
            if (column_index[j] >= 0) lattice_dirs.addEntry(j, abs(column_index[i])-1, -col[column_index[j]]);
            else if (j == i) lattice_dirs.addEntry(j, abs(column_index[i])-1, 1);
          }
          delete[] col;
        } else{
          double* col = new double[pre_m];
          assert(!CPXXbinvcol(env, presolved_model, i-n, col));
          for (int j = 0; j < n; j ++){
            if (column_index[j] >= 0) lattice_dirs.addEntry(j, abs(column_index[i])-1, -col[column_index[j]]);
          }
          delete[] col;
        }
      }
    }
  }
  exe_tableau = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count() / 1000000000.0;
}

CplexLatticeSolver::CplexLatticeSolver(MatrixXd A, VectorXd b, VectorXd c, VectorXd lb, VectorXd ub){
  exe_init = exeTime([](CplexLatticeSolver* ls, MatrixXd& A, VectorXd& b, VectorXd& c, VectorXd& lb, VectorXd& ub){
    ls->initialize(A, b, c, lb, ub);
  }, this, A, b, c, lb, ub);
}

CplexLatticeSolver::CplexLatticeSolver(MatrixXd A, VectorXd b, VectorXd c, double lb_value, double ub_value){
  VectorXd lb (c.size()); lb.fill(lb_value);
  VectorXd ub (c.size()); ub.fill(ub_value);
  exe_init = exeTime([](CplexLatticeSolver* ls, MatrixXd& A, VectorXd& b, VectorXd& c, VectorXd& lb, VectorXd& ub){
    ls->initialize(A, b, c, lb, ub);
  }, this, A, b, c, lb, ub);
}

CplexLatticeSolver::CplexLatticeSolver(MatrixXd A, VectorXd b, VectorXd c): CplexLatticeSolver(A, b, c, 0, DBL_MAX){
}

CplexLatticeSolver::CplexLatticeSolver(string problem){
  int status = 0;
  env = CPXXopenCPLEX(&status); assert(!status);
  model = CPXXcreateprob(env, &status, problem.c_str()); assert(!status);
  assert(!CPXXreadcopyprob(env, model, fmt::format("{}/resource/test_cases/{}.mps", kProjectHome, problem).c_str(), NULL));
  numcols = CPXXgetnumcols(env, model);
  numrows = CPXXgetnumrows(env, model);
  obj.resize(numcols);
  assert(!CPXXgetobj(env, model, &obj(0), 0, numcols-1));
  solveRelaxed();
}

bool CplexLatticeSolver::checkFeasibilityStep(const VectorXd& current_feas_scores, const unordered_map<int,int>& bound_map){
  if (bound_map.size() > 0) return false;
  for (int j = 0; j < A.innerSize(); j ++){
    if (current_feas_scores(j) < 0) return false;
  }
  return true;
}

bool CplexLatticeSolver::latticeStep(int d, VectorXd& x, VectorXd& current_feas_scores, double& current_c_score, unordered_map<int,int>& bound_map){
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
  for (int j = 0; j < A.innerSize(); j ++) current_feas_scores(j) -= s*A(j, i);
  current_c_score += s*c(i);
  return true;
}

void CplexLatticeSolver::solve(double time_budget) {
  chrono::high_resolution_clock::time_point ep = chrono::high_resolution_clock::now();
  try_count = 0;
  avg_step_count = 0;
  #pragma omp parallel num_threads(CORE_COUNT)
  {
    // Local declaration
    int local_try_count = 0;
    int local_step_count = 0;
    double local_best_c_score = -DBL_MAX;
    VectorXd local_best_x;
    default_random_engine gen {static_cast<long unsigned int>(time(0)*omp_get_thread_num())};
    uniform_real_distribution u_dist(0.0, 1.0);

    while (true) {
      if (chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - ep).count() / 1000000000.0 > time_budget)
        break;
      local_try_count ++;
      vector<double> splits(n-1);
      for (int i = 0; i < n-1; i++) splits[i] = u_dist(gen);
      VectorXd coefs = bucketSort(splits);
      for (int i = n-1; i >= 1; i--) coefs(i) -= coefs(i-1);
      if (!kIgnoreDegenerate && degenerate_count < n){
        for (int i = 0; i < n; i ++){
          // Weird bug when coefs(i) = 0
          if (degenerate_violation[i].size() > 0) coefs(i) = kFloatEps;
        }
        double sum_coefs = coefs.sum();
        for (int i = 0; i < n; i ++) coefs(i) /= sum_coefs;
      }
      VectorXd dir = lattice_dirs.vectorProduct(coefs);
      AbstractWalker* walker;
      if (n < kLatticeUB) walker = new Walker(dir);
      else if (n < kRetargetUB) walker = new PseudoWalker(dir, true);
      else walker = new PseudoWalker(dir, false);
      VectorXd x = floor_r0;
      VectorXd current_feas_scores = feas_scores;
      double current_c_score = c_score;
      unordered_map<int, int> bound_map;
      for (int i = 0; i < n; i ++){
        if (u_dist(gen) < r0(i) - floor_r0(i)){
          latticeStep(i+1, x, current_feas_scores, current_c_score, bound_map);
        }
      }
      if (checkFeasibilityStep(current_feas_scores, bound_map)){
        if (local_best_c_score < current_c_score){
          local_best_c_score = current_c_score;
          local_best_x = x;
        }
      }
      int step_count = 0;
      double stubborn_scale = 1.0;
      while (true){
        int d = walker->step();
        if (!latticeStep(d, x, current_feas_scores, current_c_score, bound_map)) break;
        if (checkFeasibilityStep(current_feas_scores, bound_map)){
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
      local_step_count += step_count;
    }
    #pragma omp critical 
    {
      if (best_c_score < local_best_c_score){
        best_c_score = local_best_c_score;
        best_x = local_best_x;
      }
      try_count += local_try_count;
      avg_step_count += local_step_count;
    }
  }
  avg_step_count /= try_count;
  if (best_c_score != -DBL_MAX){
    best_x = uncrushX(best_x);
    best_c_score = getObjValue(best_x);
    if (sol_status == LS_NOT_FOUND) sol_status = LS_FOUND;
  }
}

void CplexLatticeSolver::milpSolve(){
  exe_milp = exeTime([](CPXENVptr& env, CPXLPptr& model){
    assert(!CPXXmipopt(env, model));
  }, env, model);
  x0.resize(numcols);
  assert(!CPXXgetx(env, model, &x0(0), 0, numcols-1));
  milp_score = getObjValue(x0);
}