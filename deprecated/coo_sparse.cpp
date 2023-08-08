#include <iostream>
#include "coo_sparse.h"

CooSparse::CooSparse(){
}

CooSparse::CooSparse(int n, int m){
  this->n = n;
  this->m = m;
  rows.resize(n);
  for (int i = 0; i < n; i ++) rows[i].reserve(m+1);
}

void CooSparse::addEntry(int r, int c, double v){
  rows[r].emplace_back(c, v);
}

int CooSparse::getSize(){
  int size = 0;
  for (int i = 0; i < n; i ++) size += rows[i].size();
  return size;
}

// Product, rows as direction, cols as component
void CooSparse::inplaceVectorProduct(const VectorXd& x, VectorXd& res){
  for (int i = 0; i < n; i ++){
    for (const auto& pi : rows[i]){
      res(pi.first) += x(i) * pi.second;
    }
  }
}

// Transpose product, rows as component, cols as direction
VectorXd CooSparse::vectorProduct(const VectorXd& x){
  VectorXd res (n); res.fill(0);
  for (int i = 0; i < n; i ++){
    for (const auto& pi : rows[i]){
      res(i) += x(pi.first) * pi.second;
    }
  }
  return res;
}

void CooSparse::print(){
  for (int i = 0; i < n; i ++){
    for (const auto& pi : rows[i]){
      cout << "(" << pi.first << "," << pi.second << ") ";
    }
    cout << endl;
  }
}