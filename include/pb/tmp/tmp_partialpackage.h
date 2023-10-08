#pragma once

#include "pb/det/det_prob.h"
#include "pb/det/lsr_prob.h"
#include "pb/det/kd_tree.h"


#include "pb/util/udebug.h"
#include "pb/util/udeclare.h"
#include "pb/util/upostgres.h"
#include "pb/util/uconfig.h"
#include "libpq-fe.h"
#include "float.h"
#include "Eigen/Dense"

using namespace pb;

class TmpPartialPackage {
    public:
        PGconn *_conn;
        PGresult *_res;
        VectorXd sum;
        map<long long, long long> sketch_sol; // mapping gid to its corresponding occurrences
        vector<long long> initial_ids; // the initial ids of sketch_det_prob, indicates the index to gid mapping
        RMatrixXd temp_A;
        VectorXd refine_sol;
        LsrProb lsr_prob;
        DetProb sketch_det_prob, refine_det_prob;
        long long num_total_group;
        //vector<int> group_indices;
        vector<string> g_cols; // att_cols aliased as g
        unordered_set<long long> sketch_gids; // used for indexing the RMatrix
        string filter_conds;

    public:
        ~TmpPartialPackage();
        TmpPartialPackage(LsrProb &lsr_rob);
        void init(PGconn *conn, DetProb &sketch_det_prob, map<long long, long long> &sketch_sol, unordered_set<long long> &sketch_gids);
        bool refine(map<long long, long long> &sol, int core);
};
