#pragma once

#include "pb/det/det_prob.h"
#include "pb/det/lsr_prob.h"
#include "pb/det/kd_tree.h"

#include "pb/util/udebug.h"
#include "pb/util/udeclare.h"
#include "pb/util/upostgres.h"
#include "pb/util/uconfig.h"
#include "libpq-fe.h"

using namespace pb;

class TmpSketchRefine {
    private:
        PgManager *pg;
        PGconn *_conn;
        PGresult *_res;
        string _sql, group_table_name, partition_table_name;
        int core;
        bool checkTupleFiltered(VectorXd &tuple, VectorXd &bl, VectorXd &bu);
    public:
        DetProb det_prob;
        LsrProb prob;
        double exec_sr, exec_refine, exec_sketch;
        ~TmpSketchRefine();
        TmpSketchRefine(LsrProb &lsr_prob, int core=kPCore);
        void init();
        bool sketchAndRefine(map<long long, long long> &sol);
        void formulateDetProb(LsrProb &prob, string sql);
};
