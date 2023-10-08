#include "tmp_partialpackage.h"
#include "pb/core/gurobi_solver.h"
#include "pb/core/checker.h"

#define VERBOSE 0

TmpPartialPackage::~TmpPartialPackage() {}
// #TODO: does the order of column is addressed here? i.e. if A has different column order than bl & bu, could we still solve it?
void TmpPartialPackage::init(PGconn *conn, DetProb &sketch_det_prob, map<long long, long long> &sketch_sol, unordered_set<long long> &sketch_gids)
{
    this->_conn = conn;
    this->sketch_det_prob = sketch_det_prob;
    //this->sketch_det_prob.copyBounds(sketch_det_prob.bl, sketch_det_prob.bu, sketch_det_prob)
    this->refine_det_prob = sketch_det_prob;
    this->sketch_sol = sketch_sol;
    this->sketch_gids = sketch_gids;
    this->temp_A = sketch_det_prob.A;
    this->initial_ids = sketch_det_prob.ids;
}

// TmpPartialPackage::TmpPartialPackage(DetProb &det_prob, map<long long, long long> &sketch_sol, vector<int> &group_indices, vector<string> &g_cols, string filter_conds): sketch_sol(sketch_sol), group_indices(group_indices), g_cols(g_cols), filter_conds(filter_conds){
//     num_total_group = det_prob.ids.size();
//     sketch_det_prob = det_prob;
//     refine_det_prob = det_prob;
//     sketch_gids = det_prob.ids;
// }

TmpPartialPackage::TmpPartialPackage(LsrProb &lsr_prob) : lsr_prob(lsr_prob)
{
    g_cols = vector<string>(lsr_prob.det_sql.att_cols.size());
    for (size_t i = 0; i < lsr_prob.det_sql.att_cols.size(); i++)
    {
        g_cols[i] = "g." + lsr_prob.det_sql.att_cols[i];
    }
    filter_conds = getFilterConds(lsr_prob.det_sql.filter_cols, lsr_prob.det_sql.filter_intervals, kPrecision);
}

bool TmpPartialPackage::refine(map<long long, long long> &sol, int core)
{
    if (sketch_gids.empty())
        return true;

    string g_cols_name = join(g_cols, ",");
    string group_table_name = fmt::format("[1G]_{}_{}", lsr_prob.det_sql.table_name, lsr_prob.partition_name);
    string partition_table_name = fmt::format("[1P]_{}_{}", lsr_prob.det_sql.table_name, lsr_prob.partition_name);

    /* Procedure for refining the groups:
            1. Query Postgres for the actual tuples for a specific group id -> det_prob.A
            2. Remove the group repr tuples from sketch
            3. recalculate the bounds based on other sketch tuples -> det_prob.bl & bu
            4. Copy bounds from lsr_prob for cl, cu, u, l
            5. Solve it using Gurobi

            100 million
            ~10-200 million operations

            1 million = 100 millisecond
    */

    // RMatrixXd effective_A = refine_det_prob.A(Eigen::all, group_indices);
    // auto mean = effective_A.colwise().sum() / effective_A.cols();
    long E = sketch_det_prob.ids.size(); // Package size
    
    #if VERBOSE
    fmt::print("Total package size: {}\n", E);
    #endif

    int m = g_cols.size() + 1; // Account for cl & cu
    int num_group_refined_phase1 = 0;
    int num_iter = 1;

    // Phase 1 start here

    while (sketch_gids.size() > 0)
    {
        // cout << "Sketch size " << sketch_gids.size() << endl;
        // As long as there are unrefined groups, start refining from the first one
        long long refine_group_start_idx = 0, refine_group_end_idx;
        int num_group_refined_iter = 0;
    
        while (refine_group_start_idx < E)
        {
            // cout << "OKc " << refine_group_start_idx << " " << E << endl; 
            long long refine_gid = sketch_det_prob.ids[refine_group_start_idx];
            // If this group is already refined, skip
            if (sketch_gids.find(refine_gid) == sketch_gids.end())
            {
                refine_group_start_idx ++;
                continue;
            }
            long long num_repr_tuple = sketch_sol.at(refine_gid);
            refine_group_end_idx = refine_group_start_idx + num_repr_tuple;

            // fmt::print("refine gid: {}\tnum_repr_tuple:{}\trefine_group_start_idx:{}\trefine_group_end_idx:{}\n", refine_gid, num_repr_tuple, refine_group_start_idx, refine_group_end_idx);
            // vector<int> group_indices;
            long seq_len = E - (long)num_repr_tuple;
            int *group_idx_seq = new int[seq_len];
            iota(group_idx_seq, group_idx_seq + refine_group_start_idx, 0);
            iota(group_idx_seq + refine_group_start_idx, group_idx_seq + E - num_repr_tuple, refine_group_end_idx);
            vector<int> seq_v (group_idx_seq, group_idx_seq+seq_len);
            delete[] group_idx_seq;
            // for(int i=0; i<E; i++) {
            //     if(i>=refine_group_start_idx && i<refine_group_end_idx) {
            //         continue;
            //     }
            //     group_indices.push_back(i);
            // }
            // TODO: Handle cl & cu
            RMatrixXd sketch_tmp_A = sketch_det_prob.A(Eigen::all, seq_v);
            // cout << "sketch_det_prob_A size: " << sketch_det_prob.A.rows() << " by " << sketch_det_prob.A.cols() << endl;
            // cout << "tmp_A size: " << sketch_tmp_A.rows() << " by " << sketch_tmp_A.cols() << endl;
            // cout << "OKd\n"; 
            auto sum = sketch_tmp_A.rowwise().sum();

            //cout << "SUM: \n" << sum << endl;

            string sql = fmt::format("SELECT {}, {}, {} FROM \"{}\" p INNER JOIN \"{}\" g ON p.tid=g.id WHERE p.gid={} AND {};", kId, lsr_prob.det_sql.obj_col, g_cols_name, partition_table_name, lsr_prob.det_sql.table_name, refine_gid, filter_conds);
            _res = PQexec(_conn, sql.c_str());
            assert(PQresultStatus(_res) == PGRES_TUPLES_OK);

            int n = PQntuples(_res);
            refine_det_prob.resize(m, n);

            refine_det_prob.bl = sketch_det_prob.bl - sum;
            refine_det_prob.bl(m-1) = 0; // cl is always 0
            refine_det_prob.bu = sketch_det_prob.bu - sum;
            refine_det_prob.bu(m-1) = num_repr_tuple; // cu is the group occurrences in sketch solution
            
            // fmt::print("Phase 1 - iteration {}:\n", num_group_refined_phase1);
            // cout << "\tbl:\n" << refine_det_prob.bl << endl;
            // cout << "\tbu:\n" << refine_det_prob.bu << endl;

            for (int i = 0; i < n; i++)
            {
                refine_det_prob.ids[i] = atoll(PQgetvalue(_res, i, 0));
                refine_det_prob.c[i] = atof(PQgetvalue(_res, i, 1));
                for (int j = 0; j < m - 1; j++)
                {
                    refine_det_prob.A(j, i) = atof(PQgetvalue(_res, i, 2 + j));
                }
                if (m > 0) refine_det_prob.A(m - 1, i) = 1.0; // Only change for cl & cu when there is bl & bu
            }
            PQclear(_res);
            // cout << "OKe\n"; 

            refine_det_prob.u.fill((double)lsr_prob.det_sql.u); // Each tuple cannot exceed u times
            refine_det_prob.truncate();
            n = refine_det_prob.A.cols();
            m = refine_det_prob.A.rows();

            // fmt::print("Refining for gid: {}, size {}\n", refine_gid, num_repr_tuple);
            // cout << "bl:\n" << refine_det_prob.bl << "\nbu:\n" << refine_det_prob.bu << endl;
            // cout << "u: " << refine_det_prob.u << endl << endl;
            // cout << "A:\n" << refine_det_prob.A << endl;

            // Solve refine for each group
            GurobiSolver gs = GurobiSolver(refine_det_prob, true);
            // cout << "Start refine" << endl;
            gs.solveIlp(1e-4, kTimeLimit/10);
            // cout << "Finish refine" << endl;
            // cout << "Refine status:" << gs.ilp_status << ": " << solMessage(gs.ilp_status) << endl;

            Checker ch = Checker(refine_det_prob);
            int feasStatus = ch.checkIlpFeasibility(gs.ilp_sol);
            if (feasStatus != 1)
            {
                // fmt::print("Refining for gid: {}, size {}\n", refine_gid, num_repr_tuple);
                // cout << "bl:\n" << refine_det_prob.bl << "\nbu:\n" << refine_det_prob.bu << endl;
                // cout << "u: " << refine_det_prob.u << endl << endl;
                // cout << "A:\n" << refine_det_prob.A << endl;
                #if VERBOSE
                fmt::print("Refine failed for group {} in iteration {}\n", refine_gid, num_iter);
                #endif

                refine_group_start_idx = refine_group_end_idx;
                continue;
            }
            else
            {
                num_group_refined_iter++;
                assert(sketch_gids.erase(refine_gid) > 0);
            }

            // cout << "OKa\n"; 

            // Replace the repr tuples with the actual tuples
            int tmp = refine_group_start_idx;
            for (int i = 0; i < gs.ilp_sol.size() && tmp < refine_group_end_idx; i++)
            {
                if (gs.ilp_sol(i) != 0)
                {
                    // int occurrence = (int)gs.ilp_sol(i);
                    for (int j = 0; j < m; j++)
                    {
                        this->temp_A(j, tmp) = refine_det_prob.A(j, i);
                        // sketch_det_prob.A(j, tmp) = refine_det_prob.A(j, i);
                    }
                    sketch_det_prob.ids[tmp] = refine_det_prob.ids[i]; // Assigning the actual tuples' IDs to sketch_det_prob does not affect the calculation of ILP constraints
                    tmp += 1;
                }
            }
            // cout << "OKb\n"; 
            num_group_refined_phase1 += 1;
            refine_group_start_idx = refine_group_end_idx;
        }

        if (num_group_refined_iter == 0)
        {
            #if VERBOSE
            fmt::print("[Failed] No group could be refined in iteration {}\n", num_iter);
            #endif

            return false;
        }else {
            #if VERBOSE
            fmt::print("Refined {} groups during iteration {}\n", num_group_refined_iter, num_iter);
            #endif
        }
        // Update the sketch_det_prob to actual tuples refined in this round
        sketch_det_prob.A = temp_A;
        num_iter++;
    }

    // Phase 1 end here

    // cout << "OK3\n"; 
    for (auto id : sketch_det_prob.ids)
    {
        auto pair = sol.emplace(id, 0);
        pair.first->second += 1;
    }

    LsrChecker lsr_ch = LsrChecker(lsr_prob);
    int feasStatus = lsr_ch.checkIlpFeasibility(sol);
    if (feasStatus == 1)
    {
        #if VERBOSE
        fmt::print("[Suceed] Phase 1 feasible, finished refining for {} groups, with {} tuples in the package\n", num_group_refined_phase1, sketch_det_prob.ids.size());
        #endif

        return true;
    }

    #if VERBOSE
    fmt::print("Phase 1 failed, feasibility status: {}, trying to re-refine\n", feasMessage(feasStatus));
    #endif

    // Phase 2: Going though each group again for replacement
    long long refine_group_start_idx = 0, refine_group_end_idx;
    while (refine_group_start_idx < E)
    {
        // cout << "RGSI " << refine_group_start_idx << " E " << E << "\n"; 
        // long long refine_gid = sketch_det_prob.ids[refine_group_start_idx];
        long long re_refine_gid = initial_ids[refine_group_start_idx];
        long long num_repr_tuple = sketch_sol[re_refine_gid];
        assert(num_repr_tuple > 0);
        refine_group_end_idx = refine_group_start_idx + num_repr_tuple;

        // fmt::print("refine gid: {}\tnum_repr_tuple:{}\trefine_group_start_idx:{}\trefine_group_end_idx:{}\n", refine_gid, num_repr_tuple, refine_group_start_idx, refine_group_end_idx);
        // vector<int> group_indices;
        long seq_len = E - (long)num_repr_tuple;
        int *group_idx_seq = new int[seq_len];
        iota(group_idx_seq, group_idx_seq + refine_group_start_idx, 0);
        iota(group_idx_seq + refine_group_start_idx, group_idx_seq + E - num_repr_tuple, refine_group_end_idx);
        vector<int> seq_v (group_idx_seq, group_idx_seq+seq_len);
        delete[] group_idx_seq;

        RMatrixXd sketch_tmp_A = sketch_det_prob.A(Eigen::all, seq_v);
        auto sum = sketch_tmp_A.rowwise().sum();

        string sql = fmt::format("SELECT {}, {}, {} FROM \"{}\" p INNER JOIN \"{}\" g ON p.tid=g.id WHERE p.gid={} AND {};", kId, lsr_prob.det_sql.obj_col, g_cols_name, partition_table_name, lsr_prob.det_sql.table_name, re_refine_gid, filter_conds);
        _res = PQexec(_conn, sql.c_str());
        assert(PQresultStatus(_res) == PGRES_TUPLES_OK);

        int n = PQntuples(_res);
        refine_det_prob.resize(m, n);
        refine_det_prob.bl = sketch_det_prob.bl - sum;
        refine_det_prob.bu = sketch_det_prob.bu - sum;
        // cout << "Sketch\n";
        // cout << "bl:\n" << sketch_det_prob.bl << "\nbu:\n" << sketch_det_prob.bu << endl;
        // cout << refine_det_prob.A << endl;
        for (int i = 0; i < n; i++)
        {
            refine_det_prob.ids[i] = atoll(PQgetvalue(_res, i, 0));
            refine_det_prob.c[i] = atof(PQgetvalue(_res, i, 1));
            for (int j = 0; j < m - 1; j++)
            {
                refine_det_prob.A(j, i) = atof(PQgetvalue(_res, i, 2 + j));
            }
            if (m > 0)
                refine_det_prob.A(m - 1, i) = 1.0; // Only change for cl & cu when there is bl & bu
        }
        PQclear(_res);

        refine_det_prob.u.fill((double)lsr_prob.det_sql.u);
        refine_det_prob.truncate();
        n = refine_det_prob.A.cols();
        m = refine_det_prob.A.rows();

        // Solve refine for each group
        GurobiSolver gs = GurobiSolver(refine_det_prob);
        gs.solveIlp(1e-4, kTimeLimit/10);
        // cout << "Refine status:" << gs.ilp_status << ": " << solMessage(gs.ilp_status) << endl;
        Checker ch = Checker(refine_det_prob);
        int feasStatus = ch.checkIlpFeasibility(gs.ilp_sol);
        if (feasStatus != 1)
        {
            // fmt::print("Phase2: Re-refine failed {} at group {}\n", feasMessage(feasStatus), re_refine_gid);
            
            refine_group_start_idx = refine_group_end_idx;
            continue;
        }

        //cout << "u: " << refine_det_prob.u << endl << endl;
        //cout << "A:\n" << refine_det_prob.A << endl;

        // cout << "Before re-refine\n" ;
        // print(sol);

        // Replace the repr tuples with the actual tuples
        int tmp = refine_group_start_idx;
        for (int i = 0; i < gs.ilp_sol.size() && tmp < refine_group_end_idx; i++)
        {
            if (gs.ilp_sol(i) != 0)
            {
                int occurrence = (int)gs.ilp_sol(i);
                for (int j = 0; j < m; j++)
                {
                    sketch_det_prob.A(j, tmp) = refine_det_prob.A(j, i);
                    // cout << refine_det_prob.A(j, i) << " ";
                }
                // long long old_tid = sketch_det_prob.ids[tmp];
                // long long new_tid = refine_det_prob.ids[i];
                assert(sol.erase(sketch_det_prob.ids[tmp]) > 0 );
                sketch_det_prob.ids[tmp] = refine_det_prob.ids[i];
                auto pair = sol.emplace(refine_det_prob.ids[i], 0);
                pair.first->second += occurrence;
                //fmt::print("Substituted gid {}: tid {} -> tid {}\n", re_refine_gid, old_tid, new_tid );
                tmp += 1;
            }
        }
        // Solution updated, return
        #if VERBOSE
        fmt::print("[Suceed] Phase 2 suceeded, re-refined group {} with size {}, feasMessage: {}\n", re_refine_gid, num_repr_tuple,feasMessage(feasStatus));
        #endif
        //fmt::print("Refining for gid: {}, size {}\n", refine_gid, num_repr_tuple);

        return true;

    }

    // Exited the loop means that no re-refining can be done on any group
    #if VERBOSE
    fmt::print("[Failed] Phase 2 failed, no re-refining can be done.\n");
    #endif
    
    return false;
}
