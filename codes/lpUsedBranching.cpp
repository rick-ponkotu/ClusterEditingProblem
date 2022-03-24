#include <vector>
#include <list>
#include <iostream>
#include <algorithm>

#include "main.h"
#include "graph.h"
#include "random.hpp"
#include "randomPivot.h"


using namespace std;

int rec_cnt = 0;
void ret_rec_cnt_v3() {
    cout << rec_cnt << endl;
}
int bcd_branching(Graph& G, const Graph& G_orig,int max_obj, vector <edge>& sol, int rec_depth, vector<vector<int>>& uv, vector<vector<double>> lp_sol){
    rec_cnt ++;

    if(G.num_nodes <= 1)  return 0;

    if(G.num_nodes == 2){
        if(G.Weight(0,1) > 0 && G.Flag(0,1) == 0){
            G.permanent(0, 1, sol, G_orig);
            G.reset_flag(0,1);
        }
        if(G.Weight(0,1) <= 0 && G.Flag(0,1) == 0){
            G.forbid(0, 1, sol, G_orig);
            G.reset_flag(0,1);
        }
        return 0;
    }

    vector <int> triple;

    if(!G.conflict_triple(triple)){
        FOR(u, 0, G.num_nodes-1) FOR(v, u+1, G.num_nodes){
            if(G.Weight(u,v) > 0 && G.Flag(u,v) == 0){
                G.permanent(u, v, sol, G_orig);
                G.reset_flag(u,v);
            }
            if(G.Weight(u,v) <= 0 && G.Flag(u,v) == 0){
                G.forbid(u, v, sol, G_orig);
                G.reset_flag(u,v);
            }
        }
        return 0;
    }

    if(max_obj <= 0)  return -1;

    int best = -1;
    vector <edge> best_sol;
    if (rec_depth % 10 == 0){
        uv.clear();
        FOR(u, 0, G.num_nodes-1) FOR(v, u+1, G.num_nodes){
            if(G.Weight(u,v) <= 0 || G.Flag(u,v) != 0) continue;
            vector<int> s2_tmp;
            int cost = 0;
            //bool flag_no_edge = false;
            FOR(z, 0, G.num_nodes){
                if(u == z || v == z) continue;
                if((G.Weight(u,z) > 0 && G.Weight(v,z) <= 0) || (G.Weight(u,z) <= 0 && G.Weight(v,z) > 0)){
                    //flag_no_edge = true;
                    cost += min(abs(G.Weight(u,z)), abs(G.Weight(v,z)));
                }
            }
            //if(flag_no_edge == false) continue;
            s2_tmp.push_back(cost);
            s2_tmp.push_back(G.node_names[u]);
            s2_tmp.push_back(G.node_names[v]);
            uv.push_back(s2_tmp);
            s2_tmp.clear();
        }
        sort(uv.begin(), uv.end(), [](const vector<int> &alpha, const vector<int > &beta){return alpha[0] > beta[0];});
    }

    bool flag = false;
    if (uv.size() > 0) {
        FOR (i,0,uv.size()) {
            int u_name = uv[i][1];
            int v_name = uv[i][2];
            int u = G.name_to_ind[u_name];
            int v = G.name_to_ind[v_name];
            if (u == -1 || v == -1) continue;
            if(G.weight[u_name][v_name] <= 0 || G.flag[u_name][v_name] != 0) continue;
            if(G.Weight(u,v) <= 0 || G.Flag(u,v) != 0) continue;
            flag = true;
            bool merge_flag = true;
            FOR(z, 0, G.num_nodes){
                if(u_name == z || v_name == z) continue;
                if((G.flag[u_name][z] == -1 && G.flag[v_name][z] == 1) || (G.flag[u_name][z] == 1 && G.flag[v_name][z] == -1)) merge_flag = false;
            }

            if (lp_sol[u][v] == 1) {
                if(u_name > v_name) SWAP(int, u_name, v_name);
                int w = G.Weight(u,v);
                G.forbid(u, v, best_sol, G_orig);
                G.flip_edge(u,v);
                best = bcd_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_sol);
                if(best != -1) {
                    best += w;
                    if(best < max_obj) max_obj = best;
                }
                G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);

                if(merge_flag){
                    vector <edge> tmpsol;
                    G.permanent(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, G_orig);
                    MergeData mg(G_orig.num_nodes);
                    int mergecost = G.merge_nodes(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, mg, G_orig);
                    int cost = bcd_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_sol);
                    if((cost != -1) && (best == -1 || cost + mergecost < best)){
                        best = cost + mergecost;
                        best_sol = tmpsol;
                    }

                    G.expand_nodes(u_name, v_name, mg);
                    G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
                }
            }
            else if (lp_sol[u][v] == 0) {
                //merge forbid 入れ替え(Graph copyで対応)
                if(merge_flag){
                    Graph G_copy = G;
                    vector <edge> tmpsol;
                    G_copy.permanent(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, G_orig);
                    MergeData mg(G_orig.num_nodes);
                    int mergecost = G_copy.merge_nodes(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, mg, G_orig);
                    int cost = bcd_branching(G_copy, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_sol);
                    if((cost != -1) && (best == -1 || cost + mergecost < best)){
                        best = cost + mergecost;
                        best_sol = tmpsol;
                    }

                    G_copy.expand_nodes(u_name, v_name, mg);
                    G_copy.reset_flag(G_copy.name_to_ind[u_name],G_copy.num_nodes-1);
                }

                if(G.Flag(u,v) == 0) {
                    if(u_name > v_name) SWAP(int, u_name, v_name);
                    int w = G.Weight(u,v);
                    G.forbid(u, v, best_sol, G_orig);
                    G.flip_edge(u,v);
                    best = bcd_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_sol);
                    if(best != -1) {
                        best += w;
                        if(best < max_obj) max_obj = best;
                    }
                    G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                    G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                }
            }
            else {
                if(u_name > v_name) SWAP(int, u_name, v_name);
                int w = G.Weight(u,v);
                G.forbid(u, v, best_sol, G_orig);
                G.flip_edge(u,v);
                best = bcd_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_sol);
                if(best != -1) {
                    best += w;
                    if(best < max_obj) max_obj = best;
                }
                G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);

                if(merge_flag){
                    vector <edge> tmpsol;
                    G.permanent(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, G_orig);
                    MergeData mg(G_orig.num_nodes);
                    int mergecost = G.merge_nodes(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, mg, G_orig);
                    int cost = bcd_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_sol);
                    if((cost != -1) && (best == -1 || cost + mergecost < best)){
                        best = cost + mergecost;
                        best_sol = tmpsol;
                    }

                    G.expand_nodes(u_name, v_name, mg);
                    G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
                }
            }

            if(best != -1) sol.insert(sol.end(), best_sol.begin(), best_sol.end());
            return best;
        }
    }

    if(flag == false){
        uv.clear();
        FOR(u, 0, G.num_nodes-1) FOR(v, u+1, G.num_nodes){
            if(G.Weight(u,v) <= 0 || G.Flag(u,v) != 0) continue;
            vector<int> s2_tmp;
            int cost = 0;
            //bool flag_no_edge = false;
            FOR(z, 0, G.num_nodes){
                if(u == z || v == z) continue;
                if((G.Weight(u,z) > 0 && G.Weight(v,z) <= 0) || (G.Weight(u,z) <= 0 && G.Weight(v,z) > 0)){
                    //flag_no_edge = true;
                    cost += min(abs(G.Weight(u,z)), abs(G.Weight(v,z)));
                }
            }
            s2_tmp.push_back(cost);
            s2_tmp.push_back(G.node_names[u]);
            s2_tmp.push_back(G.node_names[v]);
            uv.push_back(s2_tmp);
            s2_tmp.clear();
        }
        sort(uv.begin(), uv.end(), [](const vector<int> &alpha, const vector<int > &beta){return alpha[0] > beta[0];});
    }

    if (uv.size() > 0) {
        FOR (i,0,uv.size()) {
            int u_name = uv[i][1];
            int v_name = uv[i][2];
            int u = G.name_to_ind[u_name];
            int v = G.name_to_ind[v_name];

            if (u == -1 || v == -1) continue;
            if(G.weight[u_name][v_name] <= 0 || G.flag[u_name][v_name] != 0) continue;
            if(G.Weight(u,v) <= 0 || G.Flag(u,v) != 0) continue;
            bool merge_flag = true;
            FOR(z, 0, G.num_nodes){
                if(u_name == z || v_name == z) continue;
                if((G.flag[u_name][z] == -1 && G.flag[v_name][z] == 1) || (G.flag[u_name][z] == 1 && G.flag[v_name][z] == -1)) merge_flag = false;
            }

            if(lp_sol[u][v] == 0) {
            //merge forbid 入れ替え(Graph copyで対応)
            if(merge_flag){
                Graph G_copy = G;
                vector <edge> tmpsol;
                G_copy.permanent(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, G_orig);
                MergeData mg(G_orig.num_nodes);
                int mergecost = G_copy.merge_nodes(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, mg, G_orig);
                int cost = bcd_branching(G_copy, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_sol);
                if((cost != -1) && (best == -1 || cost + mergecost < best)){
                    best = cost + mergecost;
                    best_sol = tmpsol;
                }

                G_copy.expand_nodes(u_name, v_name, mg);
                G_copy.reset_flag(G_copy.name_to_ind[u_name],G_copy.num_nodes-1);
            }

            if(G.Flag(u,v) == 0) {
                if(u_name > v_name) SWAP(int, u_name, v_name);
                int w = G.Weight(u,v);
                G.forbid(u, v, best_sol, G_orig);
                G.flip_edge(u,v);
                best = bcd_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_sol);
                if(best != -1) {
                    best += w;
                    if(best < max_obj) max_obj = best;
                }
                G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);
            }
            } else if (lp_sol[u][v] == 1) {
                //merge forbid 入れ替え
                if(u_name > v_name) SWAP(int, u_name, v_name);
                int w = G.Weight(u,v);
                G.forbid(u, v, best_sol, G_orig);
                G.flip_edge(u,v);

                best = bcd_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_sol);
                if(best != -1) {
                    best += w;
                    if(best < max_obj) max_obj = best;
                }
                G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);

                if(merge_flag){
                    vector <edge> tmpsol;
                    G.permanent(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, G_orig);
                    MergeData mg(G_orig.num_nodes);
                    int mergecost = G.merge_nodes(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, mg, G_orig);
                    int cost = bcd_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_sol);
                    if((cost != -1) && (best == -1 || cost + mergecost < best)){
                        best = cost + mergecost;
                        best_sol = tmpsol;
                    }

                    G.expand_nodes(u_name, v_name, mg);
                    G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
                }
            }
            else {
                if(u_name > v_name) SWAP(int, u_name, v_name);
                int w = G.Weight(u,v);
                G.forbid(u, v, best_sol, G_orig);
                G.flip_edge(u,v);

                best = bcd_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_sol);
                if(best != -1) {
                    best += w;
                    if(best < max_obj) max_obj = best;
                }
                G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);

                if(merge_flag){
                    vector <edge> tmpsol;
                    G.permanent(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, G_orig);
                    MergeData mg(G_orig.num_nodes);
                    int mergecost = G.merge_nodes(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, mg, G_orig);
                    int cost = bcd_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_sol);
                    if((cost != -1) && (best == -1 || cost + mergecost < best)){
                        best = cost + mergecost;
                        best_sol = tmpsol;
                    }

                    G.expand_nodes(u_name, v_name, mg);
                    G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
                }
            }
        if(best != -1) sol.insert(sol.end(), best_sol.begin(), best_sol.end());
        return best;
        }
    }
    if(best != -1) sol.insert(sol.end(), best_sol.begin(), best_sol.end());
    return best;
}

int bd_branching(Graph& G, const Graph& G_orig,int max_obj, vector <edge>& sol, vector<vector<int>>& uv, vector<vector<double>> lp_sol) {
    rec_cnt ++;

    if(G.num_nodes <= 1)  return 0;

    if(G.num_nodes == 2){
        if(G.Weight(0,1) > 0 && G.Flag(0,1) == 0){
            G.permanent(0, 1, sol, G_orig);
            G.reset_flag(0,1);
        }
        if(G.Weight(0,1) <= 0 && G.Flag(0,1) == 0){
            G.forbid(0, 1, sol, G_orig);
            G.reset_flag(0,1);
        }
        return 0;
    }

    vector <int> triple;

    if(!G.conflict_triple(triple)){
        FOR(u, 0, G.num_nodes-1) FOR(v, u+1, G.num_nodes){
            if(G.Weight(u,v) > 0 && G.Flag(u,v) == 0){
                G.permanent(u, v, sol, G_orig);
                G.reset_flag(u,v);
            }
            if(G.Weight(u,v) <= 0 && G.Flag(u,v) == 0){
                G.forbid(u, v, sol, G_orig);
                G.reset_flag(u,v);
            }
        }
        return 0;
    }

    if(max_obj <= 0)  return -1;

    int best = -1;
    vector <edge> best_sol;

    uv.clear();
    FOR(u, 0, G.num_nodes-1) FOR(v, u+1, G.num_nodes){
        if(G.Weight(u,v) <= 0 || G.Flag(u,v) != 0) continue;
        vector<int> s2_tmp;
        //bool flag_no_edge = false;
        int cost = 0;
        FOR(z, 0, G.num_nodes){
            if(u == z || v == z) continue;
            if((G.Weight(u,z) > 0 && G.Weight(v,z) <= 0) || (G.Weight(u,z) <= 0 && G.Weight(v,z) > 0)){
                //flag_no_edge = true;
                cost += min(abs(G.Weight(u,z)), abs(G.Weight(v,z)));
            }
        }
        //if(flag_no_edge == false) continue;
        s2_tmp.push_back(cost);
        s2_tmp.push_back(G.node_names[u]);
        s2_tmp.push_back(G.node_names[v]);
        uv.push_back(s2_tmp);
        s2_tmp.clear();
    }
    sort(uv.begin(), uv.end(), [](const vector<int> &alpha, const vector<int > &beta){return alpha[0] > beta[0];});

    FOR (i,0,uv.size()) {
        int u_name = uv[i][1];
        int v_name = uv[i][2];
        int u = G.name_to_ind[u_name];
        int v = G.name_to_ind[v_name];
        if (u == -1 || v == -1) continue;
        if(G.weight[u_name][v_name] <= 0 || G.flag[u_name][v_name] != 0) continue;
        if(G.Weight(u,v) <= 0 || G.Flag(u,v) != 0) continue;
        bool merge_flag = true;
        FOR(z, 0, G.num_nodes){
            if(u_name == z || v_name == z) continue;
            if((G.flag[u_name][z] == -1 && G.flag[v_name][z] == 1) || (G.flag[u_name][z] == 1 && G.flag[v_name][z] == -1)) merge_flag = false;
        }

        if (lp_sol[u][v] == 1) {
            if(u_name > v_name) SWAP(int, u_name, v_name);
            int w = G.Weight(u,v);
            G.forbid(u, v, best_sol, G_orig);
            G.flip_edge(u,v);
            best = bd_branching(G, G_orig, max_obj-w, best_sol, uv, lp_sol);
            if(best != -1) {
                best += w;
                if(best < max_obj) max_obj = best;
            }
            G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
            G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);

            if(merge_flag){
                vector <edge> tmpsol;
                G.permanent(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, G_orig);
                MergeData mg(G_orig.num_nodes);
                int mergecost = G.merge_nodes(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, mg, G_orig);
                int cost = bd_branching(G, G_orig, max_obj-mergecost, tmpsol, uv, lp_sol);
                if((cost != -1) && (best == -1 || cost + mergecost < best)){
                    best = cost + mergecost;
                    best_sol = tmpsol;
                }

                G.expand_nodes(u_name, v_name, mg);
                G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
            }
        }
        else if (lp_sol[u][v] == 0) {
            //merge forbid 入れ替え(Graph copyで対応)
            if(merge_flag){
                Graph G_copy = G;
                vector <edge> tmpsol;
                G_copy.permanent(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, G_orig);
                MergeData mg(G_orig.num_nodes);
                int mergecost = G_copy.merge_nodes(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, mg, G_orig);
                int cost = bd_branching(G_copy, G_orig, max_obj-mergecost, tmpsol, uv, lp_sol);
                if((cost != -1) && (best == -1 || cost + mergecost < best)){
                    best = cost + mergecost;
                    best_sol = tmpsol;
                }

                G_copy.expand_nodes(u_name, v_name, mg);
                G_copy.reset_flag(G_copy.name_to_ind[u_name],G_copy.num_nodes-1);
            }

            if(G.Flag(u,v) == 0) {
                if(u_name > v_name) SWAP(int, u_name, v_name);
                int w = G.Weight(u,v);
                G.forbid(u, v, best_sol, G_orig);
                G.flip_edge(u,v);
                best = bd_branching(G, G_orig, max_obj-w, best_sol, uv, lp_sol);
                if(best != -1) {
                    best += w;
                    if(best < max_obj) max_obj = best;
                }
                G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);
            }
        }
        else {
            if(u_name > v_name) SWAP(int, u_name, v_name);
            int w = G.Weight(u,v);
            G.forbid(u, v, best_sol, G_orig);
            G.flip_edge(u,v);
            best = bd_branching(G, G_orig, max_obj-w, best_sol, uv, lp_sol);
            if(best != -1) {
                best += w;
                if(best < max_obj) max_obj = best;
            }
            G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
            G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);

            if(merge_flag){
                vector <edge> tmpsol;
                G.permanent(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, G_orig);
                MergeData mg(G_orig.num_nodes);
                int mergecost = G.merge_nodes(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, mg, G_orig);
                int cost = bd_branching(G, G_orig, max_obj-mergecost, tmpsol, uv, lp_sol);
                if((cost != -1) && (best == -1 || cost + mergecost < best)){
                    best = cost + mergecost;
                    best_sol = tmpsol;
                }

                G.expand_nodes(u_name, v_name, mg);
                G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
            }
        }

        if(best != -1) sol.insert(sol.end(), best_sol.begin(), best_sol.end());
        return best;
    }
    if(best != -1) sol.insert(sol.end(), best_sol.begin(), best_sol.end());
    return best;
}

int d_branching(Graph& G, const Graph& G_orig,int max_obj, vector <edge>& sol, vector<vector<int>>& uv, vector<vector<double>> lp_sol){
    rec_cnt ++;

    if(G.num_nodes <= 1)  return 0;

    if(G.num_nodes == 2){
        if(G.Weight(0,1) > 0 && G.Flag(0,1) == 0){
            G.permanent(0, 1, sol, G_orig);
            G.reset_flag(0,1);
        }
        if(G.Weight(0,1) <= 0 && G.Flag(0,1) == 0){
            G.forbid(0, 1, sol, G_orig);
            G.reset_flag(0,1);
        }
        return 0;
    }

    vector <int> triple;

    if(!G.conflict_triple(triple)){
        FOR(u, 0, G.num_nodes-1) FOR(v, u+1, G.num_nodes){
            if(G.Weight(u,v) > 0 && G.Flag(u,v) == 0){
                G.permanent(u, v, sol, G_orig);
                G.reset_flag(u,v);
            }
            if(G.Weight(u,v) <= 0 && G.Flag(u,v) == 0){
                G.forbid(u, v, sol, G_orig);
                G.reset_flag(u,v);
            }
        }
        return 0;
    }

    if(max_obj <= 0)  return -1;

    int best = -1;
    vector <edge> best_sol;

    uv.clear();
    FOR(u, 0, G.num_nodes-1) FOR(v, u+1, G.num_nodes){
        if(G.Weight(u,v) <= 0 || G.Flag(u,v) != 0) continue;
        vector<int> s2_tmp;
        //bool flag = false;
        FOR(z, 0, G.num_nodes){
            if(u == z || v == z) continue;
            //if((G.Weight(u,z) > 0 && G.Weight(v,z) <= 0) || (G.Weight(u,z) <= 0 && G.Weight(v,z) > 0)) flag = true;
        }
        //if(flag == false) continue;
        s2_tmp.push_back(G.node_names[u]);
        s2_tmp.push_back(G.node_names[v]);
        uv.push_back(s2_tmp);
        s2_tmp.clear();
    }

    FOR (i,0,uv.size()) {
        int u_name = uv[i][0];
        int v_name = uv[i][1];
        int u = G.name_to_ind[u_name];
        int v = G.name_to_ind[v_name];

        if (u == -1 || v == -1) continue;
        if(G.weight[u_name][v_name] <= 0 || G.flag[u_name][v_name] != 0) continue;
        if(G.Weight(u,v) <= 0 || G.Flag(u,v) != 0) continue;

        bool merge_flag = true;
        FOR(z, 0, G.num_nodes){
            if(u_name == z || v_name == z) continue;
            if((G.flag[u_name][z] == -1 && G.flag[v_name][z] == 1) || (G.flag[u_name][z] == 1 && G.flag[v_name][z] == -1)) merge_flag = false;
        }

        if(lp_sol[u][v] == -1){
            if(u_name > v_name) SWAP(int, u_name, v_name);
            int w = G.Weight(u,v);
            G.forbid(u, v, best_sol, G_orig);
            G.flip_edge(u,v);
            best = d_branching(G, G_orig, max_obj-w, best_sol, uv, lp_sol);
            if(best != -1) {
                best += w;
                if(best < max_obj) max_obj = best;
            }
            G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
            G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);

            if(merge_flag){
                vector <edge> tmpsol;
                G.permanent(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, G_orig);
                MergeData mg(G_orig.num_nodes);
                int mergecost = G.merge_nodes(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, mg, G_orig);
                int cost = d_branching(G, G_orig, max_obj-mergecost, tmpsol, uv, lp_sol);
                if((cost != -1) && (best == -1 || cost + mergecost < best)){
                    best = cost + mergecost;
                    best_sol = tmpsol;
                }
                G.expand_nodes(u_name, v_name, mg);
                G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
            }
        }
        else if (lp_sol[u][v] == 0) {
            //merge forbid 入れ替え(Graph copyで対応)
            if(merge_flag){
                Graph G_copy = G;
                vector <edge> tmpsol;
                G_copy.permanent(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, G_orig);
                MergeData mg(G_orig.num_nodes);
                int mergecost = G_copy.merge_nodes(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, mg, G_orig);
                int cost = d_branching(G_copy, G_orig, max_obj-mergecost, tmpsol, uv, lp_sol);
                if((cost != -1) && (best == -1 || cost + mergecost < best)){
                    best = cost + mergecost;
                    best_sol = tmpsol;
                }

                G_copy.expand_nodes(u_name, v_name, mg);
                G_copy.reset_flag(G_copy.name_to_ind[u_name],G_copy.num_nodes-1);
            }

            if(G.Flag(u,v) == 0) {
                if(u_name > v_name) SWAP(int, u_name, v_name);
                int w = G.Weight(u,v);
                G.forbid(u, v, best_sol, G_orig);
                G.flip_edge(u,v);
                best = d_branching(G, G_orig, max_obj-w, best_sol, uv, lp_sol);
                if(best != -1) {
                    best += w;
                    if(best < max_obj) max_obj = best;
                }
                G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);
            }
        }
        else {
            if(u_name > v_name) SWAP(int, u_name, v_name);
            int w = G.Weight(u,v);
            G.forbid(u, v, best_sol, G_orig);
            G.flip_edge(u,v);
            best = d_branching(G, G_orig, max_obj-w, best_sol, uv, lp_sol);
            if(best != -1) {
                best += w;
                if(best < max_obj) max_obj = best;
            }
            G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
            G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);

            if(merge_flag){
                vector <edge> tmpsol;
                G.permanent(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, G_orig);
                MergeData mg(G_orig.num_nodes);
                int mergecost = G.merge_nodes(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, mg, G_orig);
                int cost = d_branching(G, G_orig, max_obj-mergecost, tmpsol, uv, lp_sol);
                if((cost != -1) && (best == -1 || cost + mergecost < best)){
                    best = cost + mergecost;
                    best_sol = tmpsol;
                }

                G.expand_nodes(u_name, v_name, mg);
                G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
            }
        }
        if(best != -1) sol.insert(sol.end(), best_sol.begin(), best_sol.end());
        return best;
    }
    if(best != -1) sol.insert(sol.end(), best_sol.begin(), best_sol.end());
    return best;
}

int cd_branching(Graph& G, const Graph& G_orig,int max_obj, vector <edge>& sol, int rec_depth, vector<vector<int>>& uv, vector<vector<double>> lp_sol){
    rec_cnt ++;

    if(G.num_nodes <= 1)  return 0;

    if(G.num_nodes == 2){
        if(G.Weight(0,1) > 0 && G.Flag(0,1) == 0){
            G.permanent(0, 1, sol, G_orig);
            G.reset_flag(0,1);
        }
        if(G.Weight(0,1) <= 0 && G.Flag(0,1) == 0){
            G.forbid(0, 1, sol, G_orig);
            G.reset_flag(0,1);
        }
        return 0;
    }

    vector <int> triple;

    if(!G.conflict_triple(triple)){
        FOR(u, 0, G.num_nodes-1) FOR(v, u+1, G.num_nodes){
            if(G.Weight(u,v) > 0 && G.Flag(u,v) == 0){
                G.permanent(u, v, sol, G_orig);
                G.reset_flag(u,v);
            }
            if(G.Weight(u,v) <= 0 && G.Flag(u,v) == 0){
                G.forbid(u, v, sol, G_orig);
                G.reset_flag(u,v);
            }
        }
        return 0;
    }

    if(max_obj <= 0)  return -1;

    int best = -1;
    vector <edge> best_sol;
    if (rec_depth % 10 == 0){
        uv.clear();
        FOR(u, 0, G.num_nodes-1) FOR(v, u+1, G.num_nodes){
            if(G.Weight(u,v) <= 0 || G.Flag(u,v) != 0) continue;
            vector<int> s2_tmp;
            //bool flag_no_edge = false;
            FOR(z, 0, G.num_nodes){
                if(u == z || v == z) continue;
                //if((G.Weight(u,z) > 0 && G.Weight(v,z) <= 0) || (G.Weight(u,z) <= 0 && G.Weight(v,z) > 0)) flag_no_edge = true;
            }
            //if(flag_no_edge == false) continue;
            s2_tmp.push_back(G.node_names[u]);
            s2_tmp.push_back(G.node_names[v]);
            uv.push_back(s2_tmp);
            s2_tmp.clear();
        }
    }

    bool flag = false;
    if (uv.size() > 0) {
        FOR (i,0,uv.size()) {
            int u_name = uv[i][0];
            int v_name = uv[i][1];
            int u = G.name_to_ind[u_name];
            int v = G.name_to_ind[v_name];
            if (u == -1 || v == -1) continue;
            if(G.weight[u_name][v_name] <= 0 || G.flag[u_name][v_name] != 0) continue;
            if(G.Weight(u,v) <= 0 || G.Flag(u,v) != 0) continue;
            flag = true;
            bool merge_flag = true;
            FOR(z, 0, G.num_nodes){
                if(u_name == z || v_name == z) continue;
                if((G.flag[u_name][z] == -1 && G.flag[v_name][z] == 1) || (G.flag[u_name][z] == 1 && G.flag[v_name][z] == -1)) merge_flag = false;
            }

            if (lp_sol[u][v] == 1) {
                if(u_name > v_name) SWAP(int, u_name, v_name);
                int w = G.Weight(u,v);
                G.forbid(u, v, best_sol, G_orig);
                G.flip_edge(u,v);
                best = cd_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_sol);
                if(best != -1) {
                    best += w;
                    if(best < max_obj) max_obj = best;
                }
                G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);

                if(merge_flag){
                    vector <edge> tmpsol;
                    G.permanent(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, G_orig);
                    MergeData mg(G_orig.num_nodes);
                    int mergecost = G.merge_nodes(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, mg, G_orig);
                    int cost = cd_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_sol);
                    if((cost != -1) && (best == -1 || cost + mergecost < best)){
                        best = cost + mergecost;
                        best_sol = tmpsol;
                    }

                    G.expand_nodes(u_name, v_name, mg);
                    G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
                }
            }
            else if (lp_sol[u][v] == 0) {
                //merge forbid 入れ替え(Graph copyで対応)
                if(merge_flag){
                    Graph G_copy = G;
                    vector <edge> tmpsol;
                    G_copy.permanent(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, G_orig);
                    MergeData mg(G_orig.num_nodes);
                    int mergecost = G_copy.merge_nodes(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, mg, G_orig);
                    int cost = cd_branching(G_copy, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_sol);
                    if((cost != -1) && (best == -1 || cost + mergecost < best)){
                        best = cost + mergecost;
                        best_sol = tmpsol;
                    }

                    G_copy.expand_nodes(u_name, v_name, mg);
                    G_copy.reset_flag(G_copy.name_to_ind[u_name],G_copy.num_nodes-1);
                }

                if(G.Flag(u,v) == 0) {
                    if(u_name > v_name) SWAP(int, u_name, v_name);
                    int w = G.Weight(u,v);
                    G.forbid(u, v, best_sol, G_orig);
                    G.flip_edge(u,v);
                    best = cd_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_sol);
                    if(best != -1) {
                        best += w;
                        if(best < max_obj) max_obj = best;
                    }
                    G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                    G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                }
            }
            else {
                if(u_name > v_name) SWAP(int, u_name, v_name);
                int w = G.Weight(u,v);
                G.forbid(u, v, best_sol, G_orig);
                G.flip_edge(u,v);
                best = cd_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_sol);
                if(best != -1) {
                    best += w;
                    if(best < max_obj) max_obj = best;
                }
                G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);

                if(merge_flag){
                    vector <edge> tmpsol;
                    G.permanent(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, G_orig);
                    MergeData mg(G_orig.num_nodes);
                    int mergecost = G.merge_nodes(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, mg, G_orig);
                    int cost = cd_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_sol);
                    if((cost != -1) && (best == -1 || cost + mergecost < best)){
                        best = cost + mergecost;
                        best_sol = tmpsol;
                    }

                    G.expand_nodes(u_name, v_name, mg);
                    G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
                }
            }

            if(best != -1) sol.insert(sol.end(), best_sol.begin(), best_sol.end());
            return best;
        }
    }

    if(flag == false){
        uv.clear();
        FOR(u, 0, G.num_nodes-1) FOR(v, u+1, G.num_nodes){
            if(G.Weight(u,v) <= 0 || G.Flag(u,v) != 0) continue;
            vector<int> s2_tmp;
            //bool flag_no_edge = false;
            FOR(z, 0, G.num_nodes){
                if(u == z || v == z) continue;
                //if((G.Weight(u,z) > 0 && G.Weight(v,z) <= 0) || (G.Weight(u,z) <= 0 && G.Weight(v,z) > 0)) flag_no_edge = true;
            }
            //if(flag_no_edge == false) continue;
            s2_tmp.push_back(G.node_names[u]);
            s2_tmp.push_back(G.node_names[v]);
            uv.push_back(s2_tmp);
            s2_tmp.clear();
        }
    }

    if (uv.size() > 0) {
        FOR (i,0,uv.size()) {
            int u_name = uv[i][0];
            int v_name = uv[i][1];
            int u = G.name_to_ind[u_name];
            int v = G.name_to_ind[v_name];

            if (u == -1 || v == -1) continue;
            if(G.weight[u_name][v_name] <= 0 || G.flag[u_name][v_name] != 0) continue;
            if(G.Weight(u,v) <= 0 || G.Flag(u,v) != 0) continue;
            bool merge_flag = true;
            FOR(z, 0, G.num_nodes){
                if(u_name == z || v_name == z) continue;
                if((G.flag[u_name][z] == -1 && G.flag[v_name][z] == 1) || (G.flag[u_name][z] == 1 && G.flag[v_name][z] == -1)) merge_flag = false;
            }

            if(lp_sol[u][v] == 0) {
            //merge forbid 入れ替え(Graph copyで対応)
                if(merge_flag){
                    Graph G_copy = G;
                    vector <edge> tmpsol;
                    G_copy.permanent(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, G_orig);
                    MergeData mg(G_orig.num_nodes);
                    int mergecost = G_copy.merge_nodes(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, mg, G_orig);
                    int cost = cd_branching(G_copy, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_sol);
                    if((cost != -1) && (best == -1 || cost + mergecost < best)){
                        best = cost + mergecost;
                        best_sol = tmpsol;
                    }

                    G_copy.expand_nodes(u_name, v_name, mg);
                    G_copy.reset_flag(G_copy.name_to_ind[u_name],G_copy.num_nodes-1);
                }

                if(G.Flag(u,v) == 0) {
                    if(u_name > v_name) SWAP(int, u_name, v_name);
                    int w = G.Weight(u,v);
                    G.forbid(u, v, best_sol, G_orig);
                    G.flip_edge(u,v);
                    best = cd_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_sol);
                    if(best != -1) {
                        best += w;
                        if(best < max_obj) max_obj = best;
                    }
                    G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                    G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                }
            }
            else if (lp_sol[u][v] == 1) {
                //merge forbid 入れ替え
                if(u_name > v_name) SWAP(int, u_name, v_name);
                int w = G.Weight(u,v);
                G.forbid(u, v, best_sol, G_orig);
                G.flip_edge(u,v);

                best = cd_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_sol);
                if(best != -1) {
                    best += w;
                    if(best < max_obj) max_obj = best;
                }
                G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);

                if(merge_flag){
                    vector <edge> tmpsol;
                    G.permanent(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, G_orig);
                    MergeData mg(G_orig.num_nodes);
                    int mergecost = G.merge_nodes(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, mg, G_orig);
                    int cost = cd_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_sol);
                    if((cost != -1) && (best == -1 || cost + mergecost < best)){
                        best = cost + mergecost;
                        best_sol = tmpsol;
                    }

                    G.expand_nodes(u_name, v_name, mg);
                    G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
                }
            }
            else {
                if(u_name > v_name) SWAP(int, u_name, v_name);
                int w = G.Weight(u,v);
                G.forbid(u, v, best_sol, G_orig);
                G.flip_edge(u,v);

                best = cd_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_sol);
                if(best != -1) {
                    best += w;
                    if(best < max_obj) max_obj = best;
                }
                G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
                G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);

                if(merge_flag){
                    vector <edge> tmpsol;
                    G.permanent(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, G_orig);
                    MergeData mg(G_orig.num_nodes);
                    int mergecost = G.merge_nodes(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, mg, G_orig);
                    int cost = cd_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_sol);
                    if((cost != -1) && (best == -1 || cost + mergecost < best)){
                        best = cost + mergecost;
                        best_sol = tmpsol;
                    }

                    G.expand_nodes(u_name, v_name, mg);
                    G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
                }
            if(best != -1) sol.insert(sol.end(), best_sol.begin(), best_sol.end());
            return best;
            }
        }
    }
    if(best != -1) sol.insert(sol.end(), best_sol.begin(), best_sol.end());
    return best;
}