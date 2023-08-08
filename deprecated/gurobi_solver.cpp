#define FMT_HEADER_ONLY

#include <iostream>
#include <random>
#include <float.h>

#include "fmt/core.h"
#include "gurobi_solver.h"
#include "utility.h"
#include "omp.h"

using namespace Eigen;
using namespace std;

GurobiSolver::GurobiSolver(const MatrixXd& A, const VectorXd& bl, const VectorXd& bu, const VectorXd& c, const VectorXd& l, const VectorXd& u, bool is_presolve): c(c){
  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
  relaxed_status = LS_NOT_FOUND;
  ilp_status = LS_NOT_FOUND;
  n = c.size();
  m = bl.size();
  env = NULL;
  model = NULL;
  assert(!GRBemptyenv(&env));
  assert(!GRBsetparam(env, "OutputFlag", "0"));
  assert(!GRBstartenv(env));
  assert(!GRBnewmodel(env, &model, "origin", 0, NULL, NULL, NULL, NULL, NULL));
  for (int i = 0; i < n; i++) {
    assert(!GRBaddvar(model, 0, NULL, NULL, c(i), l(i), u(i), GRB_INTEGER, NULL));
  }
  assert(!GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE));
  for (int i = 0; i < A.innerSize(); i++) {
    int* ind = new int[n];
    double* val = new double[n];
    for (int j = 0; j < n; j++) {
      ind[j] = j;
      val[j] = A(i, j);
    }
    if (bu(i) != DBL_MAX) assert(!GRBaddconstr(model, n, ind, val, GRB_LESS_EQUAL, bu(i), NULL));
    if (bl(i) != -DBL_MAX) assert(!GRBaddconstr(model, n, ind, val, GRB_GREATER_EQUAL, bl(i), NULL));
    delete ind;
    delete val;
  }
  if (!is_presolve) assert(!GRBsetintparam(env, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF));
  assert(!GRBupdatemodel(model));
  exe_init = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count() / 1000000.0;
}

GurobiSolver::GurobiSolver(const MatrixXd& A, const VectorXd& b, const VectorXd& c, const VectorXd& u, bool is_presolve): c(c) {
  chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
  relaxed_status = LS_NOT_FOUND;
  ilp_status = LS_NOT_FOUND;
  n = c.size();
  m = b.size();
  env = NULL;
  model = NULL;
  assert(!GRBemptyenv(&env));
  assert(!GRBsetparam(env, "OutputFlag", "0"));
  assert(!GRBstartenv(env));
  assert(!GRBnewmodel(env, &model, "origin", 0, NULL, NULL, NULL, NULL, NULL));
  for (int i = 0; i < n; i++) {
    assert(!GRBaddvar(model, 0, NULL, NULL, c(i), 0, u(i), GRB_INTEGER, NULL));
  }
  assert(!GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE));
  for (int i = 0; i < A.innerSize(); i++) {
    int* ind = new int[n];
    double* val = new double[n];
    for (int j = 0; j < n; j++) {
      ind[j] = j;
      val[j] = A(i, j);
    }
    assert(!GRBaddconstr(model, n, ind, val, GRB_LESS_EQUAL, b(i), NULL));
    delete ind;
    delete val;
  }
  if (!is_presolve) assert(!GRBsetintparam(env, GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF));
  assert(!GRBupdatemodel(model));
  exe_init = chrono::duration_cast<chrono::nanoseconds>(chrono::high_resolution_clock::now() - start).count() / 1000000.0;
}

void GurobiSolver::solveIlp(bool is_concurrent){
  if (!is_concurrent) assert(!GRBsetintparam(env, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL));
  else assert(!GRBsetintparam(env, GRB_INT_PAR_METHOD, GRB_METHOD_CONCURRENT));
  exe_ilp = exeTime([](GRBmodel* model){
    assert(!GRBoptimize(model));
  }, model);
  assert(!GRBgetintattr(model, GRB_INT_ATTR_STATUS, &ilp_status));
  if (ilp_status == GRB_OPTIMAL){
    ilp_status = LS_FOUND;
    x0.resize(n);
    assert(!GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, n, &x0(0)));
    ilp_cscore = c.dot(x0);
  } else if (ilp_status == GRB_INFEASIBLE) ilp_status = LS_INFEASIBLE;
  else if (ilp_status == GRB_UNBOUNDED) ilp_status = LS_UNBOUNDED;
}

void GurobiSolver::solveRelaxed(bool is_concurrent){
  if (!is_concurrent) assert(!GRBsetintparam(env, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL));
  else assert(!GRBsetintparam(env, GRB_INT_PAR_METHOD, GRB_METHOD_CONCURRENT));
  GRBmodel* relaxed;
  assert(!GRBrelaxmodel(model, &relaxed));
  exe_relaxed = exeTime([](GRBmodel* relaxed) {
    assert(!GRBoptimize(relaxed));
  }, relaxed);
  assert(!GRBgetdblattr(relaxed, GRB_DBL_ATTR_ITERCOUNT, &iteration_count));
  assert(!GRBgetintattr(relaxed, GRB_INT_ATTR_STATUS, &relaxed_status));
  if (relaxed_status == GRB_OPTIMAL){
    relaxed_status = LS_FOUND;
    r0.resize(n);
    assert(!GRBgetdblattrarray(relaxed, GRB_DBL_ATTR_X, 0, n, &r0(0)));
    relaxed_cscore = c.dot(r0);
  } else if (relaxed_status == GRB_INFEASIBLE) relaxed_status = LS_INFEASIBLE;
  else if (relaxed_status == GRB_UNBOUNDED) relaxed_status = LS_UNBOUNDED;
}