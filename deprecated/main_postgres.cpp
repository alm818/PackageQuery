#include "pb/util/udeclare.h"
#include "pb/det/synthetic.h"
#include "pb/det/dlv.h"

using namespace pb;

void create_table(){
  Synthetic syn = Synthetic();
  // long long N = 1000000000;
  // vector<string> cols = {"tmass_prox", "j", "h", "k"};
  // syn.createSubtable("ssds", 7, cols, 1);
  vector<string> cols = {"price", "quantity", "discount", "tax"};
  syn.createSubtable("tpch", 7, cols, 1);
  syn.pro.print();
} 

void test_dlv(){
  DynamicLowVariance dlv = DynamicLowVariance();
  dlv.partition(Synthetic::table_name, "P0");
  dlv.pro.print();
}

int main(){
  create_table();
  // test_dlv();
}