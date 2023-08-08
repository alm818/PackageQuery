#define FMT_HEADER_ONLY

#include <float.h>
#include <iostream>
#include "utility.h"
#include "simplex.h"
#include "fmt/core.h"

using namespace std;

#define at(row, col) (numcols*(row)+(col))
#define slack(index) (n+(index))
#define artif(index) (npm+1+(index))

void draw(double* tableau, int m, int n){
  cout << "--------------------------------" << endl;
  for (int i = 0; i < m; i ++){
    for (int j = 0; j < n; j ++){
      fmt::print("{: 8.3Lf} ", tableau[i*n+j]);
    }
    cout << endl;
  }
  cout << "--------------------------------" << endl;
}

Simplex::~Simplex(){
  delete[] bhead;
  delete[] tableau;
}

void Simplex::pivot(int horizon, const VectorXd& u, int ent_col){
  #pragma omp master
  {
    lea_row = -1;
    min_ratio = DBL_MAX;
    double col_obj = tableau[at(0, ent_col)];
    if (col_obj < 0){
      // Increase entering variable
      // Check upper bound on entering variable
      if (ent_col < n && u(ent_col) >= 0){
        double ratio = u(ent_col) - tableau[at(1, ent_col)];
        if (min_ratio > ratio){
          min_ratio = ratio;
          lea_row = 1;
        }
      }
    } else{
      // Decrease entering variable
      // Check lower bound on entering variable
      double ratio = tableau[at(1, ent_col)];
      if (min_ratio > ratio){
        min_ratio = ratio;
        lea_row = 1;
      }
    }
    for (int j = 2; j <= m+1; j ++){
      int basic_col = bhead[j-2];
      double step = sign(col_obj)*tableau[at(j, ent_col)];
      double cur = tableau[at(1, basic_col)];
      double ratio = -1;
      if (basic_col < n){
        if (isGreater(step, 0) && u(basic_col) >= 0) ratio = (u(basic_col) - cur) / step;
        else if (isLess(step,0)) ratio = -cur / step;
      } else{
        if (isLess(step,0)) ratio = -cur / step;
      }
      if (isGreaterEqual(ratio,0) && min_ratio > ratio){
        min_ratio = ratio;
        lea_row = j;
      }
    }

  }
  #pragma omp barrier
  if (lea_row == 1){
    #pragma omp master
    {
      double col_obj = tableau[at(0, ent_col)];
      double step = min_ratio * sign(col_obj);
      tableau[at(1, ent_col)] -= step;
      // No pivot, just change non-basic from lb to ub or ub to lb
      for (int j = 2; j <= m+1; j ++){
        // Update cur
        int basic_col = bhead[j-2];
        tableau[at(1, basic_col)] += step*tableau[at(j, ent_col)];
        if (isEqual(tableau[at(1, basic_col)], 0)) tableau[at(1, basic_col)] = 0;
      }
    }
  } else if (lea_row == -1){
    #pragma omp master
    {
      status = LS_UNBOUNDED;
    }
  } else{
    // Pivot happening, entering variable ent_col, leaving variable lea_row
    #pragma omp master
    {
      iteration_count ++;
      double col_obj = tableau[at(0, ent_col)];
      double step = min_ratio * sign(col_obj);
      tableau[at(1, ent_col)] -= step;
      for (int j = 2; j <= m+1; j ++){
        // Update cur
        int basic_col = bhead[j-2];
        tableau[at(1, basic_col)] += step*tableau[at(j, ent_col)];
        if (isEqual(tableau[at(1, basic_col)], 0)) tableau[at(1, basic_col)] = 0;
      }
      bhead[lea_row-2] = ent_col;
      d_cache = tableau[at(lea_row, ent_col)];
      v_cache.resize(m+2);
      for (int j = 0; j <= m+1; j ++) v_cache(j) = tableau[at(j, ent_col)];
    }
    #pragma omp barrier
    #pragma omp for
    for (int i = 0; i < horizon; i ++){
      tableau[at(lea_row, i)] /= d_cache;
    }
    #pragma omp for nowait
    for (int i = 0; i < horizon; i ++){
      for (int j = m+1; j >= 0; j --){
        if (j != 1 && j != lea_row)
          tableau[at(j, i)] -= tableau[at(lea_row, i)] * v_cache(j);
      }
    }
  }
  #pragma omp barrier
}

void Simplex::selectEnteringColumn(int horizon, const VectorXd& c, const VectorXd& u, double& ent_value, int& ent_col){
  // Choosing entering variable
  #pragma omp master
  {
    ent_value = 0;
    ent_col = -1;
  }
  #pragma omp barrier
  double local_ent_value = 0;
  int local_ent_col = -1;
  #pragma omp for
  for (int i = 0; i < horizon; i ++){
    if (i == npm) continue;
    double cur = tableau[at(1, i)];
    // Positive col_obj wants to decrease, Negative col_obj wants to increase
    double col_obj = tableau[at(0, i)];
    // If the entering variable already lower bounded and obj is non negative
    if (isEqual(cur, 0) && isGreaterEqual(col_obj, 0)) continue;
    // If the entering variable already upper bounded and obj is non positive
    if (i < n && isEqual(cur, u(i)) && isLessEqual(col_obj, 0)) continue;

    // "Dantzig's rule"
    // if (local_ent_value < fabs(col_obj)){
    //   local_ent_value = fabs(col_obj);
    //   local_ent_col = i;
    // }

    // "Steepest-edge rule"
    // TO-DO: IMPLEMENT RECURSIVE UPDATE FOR STEEPEST-EDGE RULE
    double s = sign(col_obj);
    double dot = 0;
    double sq_norm = 0;
    if (i < n){
      dot -= s*c(i);
      sq_norm += 1;
    }
    for (int j = 2; j <= m+1; j ++){
      int basic_col = bhead[j-2];
      if (basic_col < n){
        double val = tableau[at(j, i)];
        dot += s*val*c(basic_col);
        sq_norm += val*val;
      }
    }
    if (sq_norm > 0){
      double cmp = dot / sqrt(sq_norm);
      if (local_ent_value < cmp){
        local_ent_value = cmp;
        local_ent_col = i;
      }
    }

    // "Another"
    // double sq_norm = 1;
    // for (int j = 2; j <= m+1; j ++){
    //   double val = tableau[at(j, i)];
    //   sq_norm += val*val;
    // }
    // double cmp = abs(tableau[at(0, i)]) / sqrt(sq_norm);
    // if (local_ent_value < cmp){
    //   local_ent_value = cmp;
    //   local_ent_col = i;
    // }
  }
  #pragma omp critical
  {
    if (ent_value < local_ent_value){
      ent_value = local_ent_value;
      ent_col = local_ent_col;
    }
  }
  #pragma omp barrier
}

Simplex::Simplex(int core, const MatrixXd& A, const VectorXd& b, const VectorXd& c, const VectorXd& u){
  status = LS_NOT_FOUND;
  iteration_count = 0;
  n = c.size();
  m = b.size();
  npm = n+m;
  bhead = new double[m];
  vector<int> ma_map;
  ma = 0;
  for (int i = 0; i < m; i ++){
    if (b(i) >= 0){
      bhead[i] = slack(i);
    } else{
      bhead[i] = artif(ma);
      ma_map.push_back(i);
      ma ++;
    }
  }
  numcols = npm+1+ma;
  tableau = new double[(m+2)*numcols];
  double ent_value = 0;
  int ent_col = -1;

  #pragma omp parallel num_threads(core)
  {
    // Phase-I: Initialization
    #pragma omp for
    for (int i = 0; i < numcols; i ++){
      for (int j = m+1; j >= 0; j --){
        if (j > 1){
          int index = j - 2;
          double s = nonNegativeSign(b(index));
          if (i < n) tableau[at(j, i)] = A(index, i) * s;
          else if (i < npm){
            if (i-n == index) tableau[at(j, i)] = s;
            else tableau[at(j, i)] = 0;
          } else if (i == npm) tableau[at(j, i)] = s * b(index);
          else{
            int ma_index = ma_map[i-npm-1];
            if (ma_index == index) tableau[at(j, i)] = 1;
            else tableau[at(j, i)] = 0;
          }
        } else{
          if (j == 1){
            if (i >= n && i < npm && b(i-n) >= 0) tableau[at(j, i)] = b(i-n);
            else if (i > npm) tableau[at(j, i)] = -b(ma_map[i-npm-1]);
            else tableau[at(j, i)] = 0;
          } else{
            if (i <= npm){
              double val = 0;
              for (auto ma_index : ma_map) val -= tableau[at(ma_index+2, i)];
              tableau[at(j, i)] = val;
            } else tableau[at(j, i)] = 0;
          }
        }
      }
    }
    // Phase-I
    while (true){
      selectEnteringColumn(numcols, c, u, ent_value, ent_col);
      #pragma omp master
      {
        if (ent_value == 0){
          if (isGreaterEqual(tableau[at(0, npm)], 0)) status = LS_FEASIBLE;
          else status = LS_INFEASIBLE;
        }
      }
      #pragma omp barrier
      if (status != LS_NOT_FOUND) break;
      pivot(numcols, u, ent_col);
      if (status == LS_UNBOUNDED) break;
    }
    if (status == LS_FEASIBLE){
      // Phase-II: Initialization
      #pragma omp master
      {
        sum = 0;
      }
      #pragma omp barrier
      double local_sum = 0;
      #pragma omp for nowait
      for (int i = 0; i <= npm; i ++){
        if (i < n) tableau[at(0, i)] = -c(i);
        else tableau[at(0, i)] = 0;
        local_sum += tableau[at(0, i)] * tableau[at(1, i)];
      }
      #pragma omp atomic
      sum += local_sum;
      #pragma omp barrier
      #pragma omp master
      {
        tableau[at(0, npm)] = sum;
        v_cache.resize(m);
        for (int j = 2; j <= m+1; j ++){
          int basic_col = bhead[j-2];
          v_cache(j-2) = tableau[at(0, basic_col)];
        }
      }
      #pragma omp barrier
      #pragma omp for
      for (int i = 0; i <= npm; i ++){
        double sub_obj = 0;
        for (int j = 2; j <= m+1; j ++){
          int basic_col = bhead[j-2];
          if (basic_col < n) sub_obj += tableau[at(j, i)] * v_cache(j-2);
        }
        tableau[at(0, i)] -= sub_obj;
      }

      // Phase-II
      while (true){
        selectEnteringColumn(npm+1, c, u, ent_value, ent_col);
        #pragma omp master
        {
          if (ent_value == 0){
            status = LS_FOUND;
          }
        }
        #pragma omp barrier
        if (status == LS_FOUND) break;
        pivot(npm+1, u, ent_col);
        if (status == LS_UNBOUNDED) break;
      }
      // #pragma omp single
      // draw(tableau, m+2, numcols);
      // #pragma omp single
      // {
      //   for (int i = 0; i < m; i ++){
      //     cout << (int)bhead[i] << " ";
      //   }
      //   cout << endl;
      // }
    }
  }
}
