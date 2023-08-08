#include <iostream>

#include "parallel_pq.h"
#include "utility.h"
#include "omp.h"

using namespace std;

ParallelPQ::~ParallelPQ(){
  if (q) delete q;
}

inline void ParallelPQ::heapify(int i){
  while (true){
    int lc = (i << 1) + 1;
    int rc = (i+1) << 1;
    if (lc < n){
      if (rc < n){
        if ((*q)[lc] < (*q)[rc]){
          if ((*q)[rc] > (*q)[i]){
            swap((*q)[rc], (*q)[i]);
            i = rc;
          } else break;
        } else{
          if ((*q)[lc] > (*q)[i]){
            swap((*q)[lc], (*q)[i]);
            i = lc;
          } else break;
        }
      } else{
        if ((*q)[lc] > (*q)[i]){
          swap((*q)[lc], (*q)[i]);
          i = lc;
        } else break;
      }
    } else break;
  }
}

ParallelPQ::ParallelPQ(int core, vector<pair<double, int>>* arr): q(arr){
  n = q->size();
  if (n <= 2 && n > 0){
    if (n == 2){
      if ((*q)[0] < (*q)[1]) swap((*q)[0], (*q)[1]);
    }
  } else if (n > 2){
    int layer = 1;
    while (n > 1){
      n >>= 1;
      layer <<= 1;
    }
    n = q->size();
    int left = -1; 
    int right = -1;
    #pragma omp parallel num_threads(core)
    {
      while (layer > 1){
        #pragma omp barrier
        #pragma omp master
        {
          left = (layer >> 1) - 1;
          right = layer - 1;
          layer >>= 1;
        }
        #pragma omp barrier
        #pragma omp for
        for (int index = left; index < right; index ++){
          heapify(index);
        }
      }
    }
    lazy = (*q)[0].first - max((*q)[1].first, (*q)[2].first);
  }
}

pair<double,int> ParallelPQ::peak(){
  if (n > 0) return (*q)[0];
  return {0,-1};
}

void ParallelPQ::pop(){
  if (n > 0){
    swap((*q)[0], (*q)[n-1]);
    n --;
    heapify(0);
  }
}

void ParallelPQ::subtractPeak(double pos_val){
  // cout << "SUBTRACTING: " << pos_val << " TO " << lazy << endl;
  // cout << "3 BIGS:" << endl;
  // cout << (*q)[0].first << ":" << (*q)[0].second << endl;
  // cout << (*q)[1].first << ":" << (*q)[1].second << endl;
  // cout << (*q)[2].first << ":" << (*q)[2].second << endl;
  lazy -= pos_val;
  (*q)[0].first -= pos_val;
  if (lazy < 0){
    heapify(0);
    lazy = (*q)[0].first - max((*q)[1].first, (*q)[2].first);
  }
}

int ParallelPQ::size(){
  return n;
}

