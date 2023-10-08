#define FMT_HEADER_ONLY

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <string>
#include <numeric>
#include <algorithm>
#include <map>
#include <set>
#include <random>
#include <float.h>
#include <queue>
#include <omp.h>
#include <chrono>
#include <thread>
#include <ctime>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <functional>

#include "Eigen/Dense"
#include "fmt/core.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RMatrixXd;

enum DetSense { LowerBounded, UpperBounded, Bounded};
enum SolStatus { NotFound, Found, Feasible, Infeasible, Unbounded, DualUnbounded, Timeout, NoPartitionFound, IncompatiblePartition, NumericalUnstability, HalfFeasible};
enum FeasStatus { Unsolved, Feasibility, Infeasibility, LbConstraint, UbConstraint, LbVariable, UbVariable, Integrality, BadFilter};
enum GroupStatus { Unitialized, Unlocked, Locked};

enum StoSense { VaR1, VaR2, VaR3, VaR4, CVaR1, CVaR2, CVaR3, CVaR4};
/*
VaR1: Pr(G >= v) >= p
VaR2: Pr(G <= v) >= p
VaR3: Pr(G >= v) <= p
VaR4: Pr(G <= v) <= p
CVaR1: E[G | G <= Q_G(p)] >= v
CVaR2: E[G | G >= Q_G(p)] <= v
CVaR3: E[G | G <= Q_G(p)] <= v
CVaR4: E[G | G >= Q_G(p)] >= v
*/

namespace pb{
  using std::map;
  using std::set;
  using std::unordered_map;
  using std::unordered_set;
  using std::priority_queue;
  using std::vector;
  using std::string;
  using std::pair;
  using std::tuple;
  using std::tie;
  using std::shared_ptr;
  
  using std::to_string;
  using std::iota;
  using std::max_element;
  using std::min_element;
  using std::sort;
  using std::min;
  using std::max;
  using std::nth_element;
  using std::distance;
  using std::hash;
  using std::isnan;

  using std::make_heap;
  using std::push_heap;
  using std::pop_heap;
  using std::make_shared;
  using std::fill;
  using std::shuffle;
  using std::random_shuffle;

  using Eigen::VectorXd;
  using Eigen::VectorXi;
  using Eigen::MatrixXd;

  using std::random_device;
  using std::default_random_engine;
  using std::uniform_real_distribution;
  using std::uniform_int_distribution;
  using std::normal_distribution;
  using std::lognormal_distribution;
  using std::seed_seq;

  using std::fstream;
  using std::ios;
  using std::cout;
  using std::cin;
  using std::endl;

  using std::invalid_argument;
}