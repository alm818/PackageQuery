#define _USE_MATH_DEFINES
#define FMT_HEADER_ONLY

#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include <float.h>
#include <queue>
#include <random>
#include <chrono>
#include <numeric>

#include "utility.h"
#include "parallel_pq.h"
#include "fmt/core.h"
#include "fitsio.h"
#include "omp.h"
#include "lattice_solver.h"
#include "pseudo_walker.h"
#include "simplex.h"
#include "gurobi_lattice_solver.h"
#include "gurobi_solver.h"
#include "reducer.h"
#include "fmt/core.h"
#include "matplotlibcpp.h"
#include <cmath>

using namespace std;
using namespace Eigen;

namespace plt = matplotlibcpp;

int main(){
  // Prepare data.
  int n = 5000;
  std::vector<double> x(n), y(n), z(n), w(n,2);
  for(int i=0; i<n; ++i) {
    x.at(i) = i*i;
    y.at(i) = sin(2*M_PI*i/360.0);
    z.at(i) = log(i);
  }

  // Set the size of output image to 1200x780 pixels
  // plt::figure_size(1200, 780);
  // // Plot line from given x and y data. Color is selected automatically.
  // plt::plot(x, y);
  // // Plot a red dashed line from given x and y data.
  // plt::plot(x, w,"r--");
  // // Plot a line whose name will show up as "log(x)" in the legend.
  // plt::named_plot("log(x)", x, z);
  // // Set x-axis to interval [0,1000000]
  // plt::xlim(0, 1000*1000);
  // // Add graph title
  // plt::title("Sample figure");
  // // Enable legend.
  // plt::legend();
  // // Save the image (file format is determined by the extension)
  // plt::save(fmt::format("{}/resource/matplotlib/basic.png", kProjectHome));
  cout << "OK" << endl;
}
