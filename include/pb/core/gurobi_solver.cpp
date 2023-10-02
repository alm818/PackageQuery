#include "pb/util/umisc.h"
#include "pb/util/unumeric.h"
#include "pb/core/gurobi_solver.h"

GurobiSolver::~GurobiSolver(){
  assert(!GRBfreemodel(model));
  GRBfreeenv(env);
}

GurobiSolver::GurobiSolver(const DetProb &prob, bool with_objective): Checker(prob){
  std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
  lp_status = NotFound;
  ilp_status = NotFound;
  env = NULL;
  model = NULL;
  assert(!GRBemptyenv(&env));
  assert(!GRBsetparam(env, "OutputFlag", "0"));
  assert(!GRBstartenv(env));
  assert(!GRBnewmodel(env, &model, "origin", 0, NULL, NULL, NULL, NULL, NULL));
  if (with_objective){
    for (int i = 0; i < n; i++) assert(!GRBaddvar(model, 0, NULL, NULL, prob.c(i), prob.l(i), prob.u(i), GRB_INTEGER, NULL));
  } else{
    assert(!GRBsetintparam(env, GRB_INT_PAR_SOLUTIONLIMIT, 1));
    for (int i = 0; i < n; i++) assert(!GRBaddvar(model, 0, NULL, NULL, 0.0, prob.l(i), prob.u(i), GRB_INTEGER, NULL));
  }

  assert(!GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE));
  for (int i = 0; i < m; i++) {
    int *ind = new int[n];
    double *val = new double[n];
    for (int j = 0; j < n; j++) {
      ind[j] = j;
      val[j] = prob.A(i, j);
    }
    if (prob.bu(i) != DBL_MAX) assert(!GRBaddconstr(model, n, ind, val, GRB_LESS_EQUAL, prob.bu(i), NULL));
    if (prob.bl(i) != -DBL_MAX) assert(!GRBaddconstr(model, n, ind, val, GRB_GREATER_EQUAL, prob.bl(i), NULL));
    delete ind;
    delete val;
  }
  assert(!GRBupdatemodel(model));
  exe_init = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000000.0;
  exe_lp = 0;
  exe_ilp = 0;
}

void GurobiSolver::writeModel(string file_name){
  assert(!GRBwrite(model, file_name.c_str()));
}

void GurobiSolver::solveIlp(double mipGap, double time_limit){
  // cout << "START GS " << n << " " << time_limit << endl;
  if (time_limit > 0){
    assert(!GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_TIMELIMIT, time_limit));
  }
  assert(!GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_MIPGAP, mipGap));
  exe_ilp = exeTime([](GRBmodel *model){
    assert(!GRBoptimize(model));
  }, model) + exe_init;
  // cout << "END GS " << n << " " << time_limit << endl;
  assert(!GRBgetintattr(model, GRB_INT_ATTR_STATUS, &ilp_status));
  if (ilp_status == GRB_OPTIMAL){
    ilp_status = Found;
    ilp_sol.resize(n);
    assert(!GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, n, &ilp_sol(0)));
    ilp_score = getScore(ilp_sol);
  } else if (ilp_status == GRB_INFEASIBLE) ilp_status = Infeasible;
  else if (ilp_status == GRB_UNBOUNDED) ilp_status = Unbounded;
  else if (ilp_status == GRB_TIME_LIMIT){
    int sol_count;
    assert(!GRBgetintattr(model, GRB_INT_ATTR_SOLCOUNT, &sol_count));
    if (sol_count){
      ilp_status = Found;
      ilp_sol.resize(n);
      assert(!GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, n, &ilp_sol(0)));
      ilp_score = getScore(ilp_sol);
    } else ilp_status = Timeout;
  }
}

void GurobiSolver::solveLp(){
  GRBmodel *relaxed;
  assert(!GRBrelaxmodel(model, &relaxed));
  GRBenv *local_env = GRBgetenv(relaxed);
  // Best configuration according to Gurobi
  assert(!GRBsetintparam(local_env, GRB_INT_PAR_PRESOLVE, 0));
  // assert(!GRBsetintparam(local_env, GRB_INT_PAR_SIFTING, 1));
  assert(!GRBsetintparam(local_env, GRB_INT_PAR_METHOD, 0));
  exe_lp = exeTime([](GRBmodel *relaxed) {
    assert(!GRBoptimize(relaxed));
  }, relaxed) + exe_init;
  assert(!GRBgetintattr(relaxed, GRB_INT_ATTR_STATUS, &lp_status));
  if (lp_status == GRB_OPTIMAL){
    lp_status = Found;
    double v;
    assert(!GRBgetdblattr(relaxed, GRB_DBL_ATTR_ITERCOUNT, &v));
    iteration_count = (int) v;
    lp_sol.resize(n);
    assert(!GRBgetdblattrarray(relaxed, GRB_DBL_ATTR_X, 0, n, &lp_sol(0)));
    lp_score = getScore(lp_sol);
  } else if (lp_status == GRB_INFEASIBLE) lp_status = Infeasible;
  else if (lp_status == GRB_UNBOUNDED) lp_status = Unbounded;
}

double GurobiSolver::getScore(const VectorXd &sol){
  if (sol.size() != n) return 0;
  VectorXd c (n);
  assert(!GRBgetdblattrarray(model, GRB_DBL_ATTR_OBJ, 0, n, &c(0)));
  return c.dot(sol);
}

int GurobiSolver::checkLpFeasibility(const VectorXd &sol){
  if (sol.size() != n) return Unsolved;
  double* sol_ptr = const_cast<double*>(&sol(0));
  GRBmodel *relaxed;
  assert(!GRBrelaxmodel(model, &relaxed));
  assert(!GRBsetdblattrarray(relaxed, GRB_DBL_ATTR_LB, 0, n, sol_ptr));
  assert(!GRBsetdblattrarray(relaxed, GRB_DBL_ATTR_UB, 0, n, sol_ptr));
  assert(!GRBoptimize(relaxed));
  int status;
  assert(!GRBgetintattr(relaxed, GRB_INT_ATTR_STATUS, &status));
  if (status == GRB_OPTIMAL) return Feasibility;
  else return Infeasibility;
}

int GurobiSolver::checkIlpFeasibility(const VectorXd &sol){
  if (sol.size() != n) return Unsolved;
  for (int i = 0; i < n; i ++){
    if (!isInteger(sol(i), epsilon)) return Integrality;
  }
  return checkLpFeasibility(sol);
}

bool GurobiSolver::hasIlpSolution(double time_limit){
  GurobiSolver gs = GurobiSolver(prob, false);
  gs.solveIlp(1e-4, time_limit);
  return gs.ilp_status == Found;
}