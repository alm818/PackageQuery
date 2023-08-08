#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <random>
#include <numeric>

#include "utility.h"
#include "lattice_solver.h"
#include "gurobi_lattice_solver.h"
#include "gurobi_solver.h"
#include "fmt/core.h"
#include "fitsio.h"
#include "omp.h"

using namespace std;
using namespace Eigen;

void run(){
  GalaxyDB gdb = GalaxyDB(10000000);
  MatrixXd A;
  VectorXd b, c, u;
  gdb.generateQuery(0.3, 3, GalaxyDB::Q2, 50, A, b, c, u);
  // for (int i = 0; i < gdb.column_names.size(); i ++){
  //   cout << i << " " << gdb.column_names[i] << " " << gdb.mv.mean(i) << " " << gdb.mv.var(i) << endl;
  // }
  LatticeSolver ls = LatticeSolver(80, A, b, c, u);
  ls.solve(2);
  cout << "LS FINISHED" << endl;
  GurobiSolver gs = GurobiSolver(A, b, c, u);
  gs.solveIlp();
  ls.compareReport(gs.x0, gs.exe_ilp+gs.exe_init);
}

int main(){
  run();
}