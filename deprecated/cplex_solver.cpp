#include "pb/util/umisc.h"
#include "pb/util/unumeric.h"
#include "pb/core/cplex_solver.h"

CplexSolver::~CplexSolver(){
  assert(!CPXXfreeprob(env, &model));
  assert(!CPXXcloseCPLEX(&env));
}

CplexSolver::CplexSolver(const DetProb &prob): Checker(prob){
  std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
  lp_status = NotFound;
  ilp_status = NotFound;
  int status = 0;
  env = CPXXopenCPLEX(&status); assert(!status);
  assert(!CPXXsetintparam(env, CPXPARAM_ScreenOutput, CPX_OFF));
  model = CPXXcreateprob(env, &status, "prob"); assert(!status);
  assert(!CPXXnewcols(env, model, n, &prob.c(0), &prob.l(0), &prob.u(0), NULL, NULL));
  assert(!CPXXchgobjsen(env, model, CPX_MAX));
  VectorXi ind (n); iota(ind.begin(), ind.end(), 0);
  long long zero = 0;
  for (int i = 0; i < m; i ++){
    if (prob.bl(i) != -DBL_MAX && prob.bu(i) != DBL_MAX){
      // Ranged constraint
      assert(!CPXXaddrows(env, model, 0, 1, n, &prob.bl(i), "R", &zero, &ind(0), &prob.A(i, 0), NULL, NULL));
      double rng = prob.bu(i) - prob.bl(i);
      assert(!CPXXchgrngval(env, model, 1, &i, &rng));
    } else if (prob.bl(i) != -DBL_MAX){
      // Lower bound constraint
      assert(!CPXXaddrows(env, model, 0, 1, n, &prob.bl(i), "G", &zero, &ind(0), &prob.A(i, 0), NULL, NULL));
    } else if (prob.bu(i) != DBL_MAX){
      // Upper bound constraint
      assert(!CPXXaddrows(env, model, 0, 1, n, &prob.bu(i), "L", &zero, &ind(0), &prob.A(i, 0), NULL, NULL));
    }
  }
  exe_init = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000000.0;
}

void CplexSolver::solveIlp(){
  VectorXi ind (n); iota(ind.begin(), ind.end(), 0);
  char *xtype = new char[n]; memset(xtype, CPX_INTEGER, n);
  assert(!CPXXchgctype(env, model, n, &ind(0), xtype));
  delete[] xtype;
  assert(!CPXXchgprobtype(env, model, CPXPROB_MILP));
  exe_ilp = exeTime([](CPXENVptr env, CPXLPptr model){
    assert(!CPXXmipopt(env, model));
  }, env, model);
  ilp_status = CPXXgetstat(env, model);
  if (ilp_status == CPXMIP_OPTIMAL || ilp_status == CPXMIP_OPTIMAL_TOL){
    ilp_sol.resize(n);
    assert(!CPXXsolution(env, model, &ilp_status, &ilp_score, &ilp_sol(0), NULL, NULL, NULL));
    ilp_status = Found;
  } else if (ilp_status == CPXMIP_INFEASIBLE) ilp_status = Infeasible;
  else if (ilp_status == CPXMIP_UNBOUNDED) ilp_status = Unbounded;
}

void CplexSolver::solveLp(){
  VectorXi ind (n); iota(ind.begin(), ind.end(), 0);
  char *xtype = new char[n]; memset(xtype, CPX_CONTINUOUS, n);
  assert(!CPXXchgctype(env, model, n, &ind(0), xtype));
  delete[] xtype;
  assert(!CPXXchgprobtype(env, model, CPXPROB_LP));
  exe_ilp = exeTime([](CPXENVptr env, CPXLPptr model){
    assert(!CPXXlpopt(env, model));
  }, env, model);
  lp_status = CPXXgetstat(env, model);
  if (lp_status == CPX_STAT_OPTIMAL){
    lp_sol.resize(n);
    assert(!CPXXsolution(env, model, &lp_status, &lp_score, &lp_sol(0), NULL, NULL, NULL));
    lp_status = Found;
  } else if (lp_status == CPX_STAT_INFEASIBLE) lp_status = Infeasible;
  else if (lp_status == CPX_STAT_UNBOUNDED) lp_status = Unbounded;
}

double CplexSolver::getLpTime(){
  return exe_init + exe_lp;
}

double CplexSolver::getIlpTime(){
  return exe_init + exe_ilp;
}

double CplexSolver::getScore(const VectorXd &sol){
  if (sol.size() != n) return 0;
  VectorXd c (n);
  assert(!CPXXgetobj(env, model, &c(0), 0, n-1));
  return c.dot(sol);
}

int CplexSolver::checkLpFeasibility(const VectorXd &sol){
  if (sol.size() != n) return Unsolved;
  VectorXd constraint_infeas (m);
  assert(!CPXXgetrowinfeas(env, model, &sol(0), &constraint_infeas(0), 0, m-1));
  for (int i = 0; i < m; i ++){
    if (isLess(constraint_infeas(i), 0, kEpsilon)) return LbConstraint;
    if (isGreater(constraint_infeas(i), 0, kEpsilon)) return UbConstraint;
  }
  VectorXd variable_infeas (n);
  assert(!CPXXgetcolinfeas(env, model, &sol(0), &variable_infeas(0), 0, n-1));
  for (int i = 0; i < n; i ++){
    if (isLess(variable_infeas(i), 0, kEpsilon)) return LbVariable;
    if (isGreater(variable_infeas(i), 0, kEpsilon)) return UbVariable;
  }
  return Feasibility;
}

int CplexSolver::checkIlpFeasibility(const VectorXd &sol){
  if (sol.size() != n) return Unsolved;
  for (int i = 0; i < n; i ++){
    if (!isInteger(sol(i), kEpsilon)) return Integrality;
  }
  return checkLpFeasibility(sol);
}