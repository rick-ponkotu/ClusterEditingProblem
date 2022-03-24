#include <vector>
#include <list>
#include <iostream>
#include <algorithm>
#include <chrono>

#include "main.h"
#include "graph.h"
#include "heuristicBranching.h"
#include "random.hpp"
#include "randomPivot.h"
#include "solver.h"

#include <ilcplex/ilocplex.h>

using namespace std;
using namespace chrono;

static int rec_cnt = 0;
static double create_verctor_time = 0;
int idx;
vector<vector<int>> s2_uv;
vector <int> tmp_uv;
system_clock::time_point start, end;


void ret_cnt() {
  cout << rec_cnt << endl;
}

double ret_time() {
  return create_verctor_time / 1000000000;
}

void show(vector<vector<int>> arr) {
  FOR(i,0,arr.size()) {
    FOR(j,0,arr[i].size()) {
      cout << arr[i][j] << ' ';
    }
    cout <<endl;
  }
}

int naive_branching(Graph& G, const Graph& G_orig,int max_obj, vector <edge>& sol, int rec_depth, vector<vector<int>>& uv, vector<vector<int>>& lp_zero, vector<vector<int>>& lp_one, vector<vector<double>> lp_sol){
  rec_cnt++;

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


  //create vector for branching edges
  //lpの値が1の頂点ペアを選ぶ（forbid)が先に動くため
  //選んだペアのlpの値によってforbid,mergeの順番を変える

  if (rec_depth % 5 == 0){
    //cout << rec_depth << endl;
    uv.clear();
    //auto start = system_clock::now();
    FOR(u, 0, G.num_nodes-1) FOR(v, u+1, G.num_nodes){
      if(G.Weight(u,v) <= 0 || G.Flag(u,v) != 0) continue;
      vector<int> s2_tmp;
      int cost = 0;
      FOR(z, 0, G.num_nodes){
        if(u == z || v == z) continue;
        if((G.Weight(u,z) > 0 && G.Weight(v,z) <= 0) || (G.Weight(u,z) <= 0 && G.Weight(v,z) > 0)){
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
      //通常ver
      lp_sol[u][v] = -1;
      //heuristic
      if (lp_sol[u][v] == 1) {
        /*
        //heuristic
        if(u_name > v_name) SWAP(int, u_name, v_name);
        int w = G.Weight(u,v);
        G.forbid(u, v, best_sol, G_orig);
        G.flip_edge(u,v);
        best = naive_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
        if(best != -1) {
          best += w;
          if(best < max_obj) max_obj = best;
        }
        G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
        G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);
        */

        if(u_name > v_name) SWAP(int, u_name, v_name);
        int w = G.Weight(u,v);
        G.forbid(u, v, best_sol, G_orig);
        G.flip_edge(u,v);
        best = naive_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
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
          int cost = naive_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
          if((cost != -1) && (best == -1 || cost + mergecost < best)){
            best = cost + mergecost;
            best_sol = tmpsol;
          }

          G.expand_nodes(u_name, v_name, mg);
          G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
        }
      }
      else if (lp_sol[u][v] == 0) {
        /*
        //heuristic
        if(merge_flag){
          vector <edge> tmpsol;
          G.permanent(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, G_orig);
          MergeData mg(G_orig.num_nodes);
          int mergecost = G.merge_nodes(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, mg, G_orig);
          int cost = naive_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
          if((cost != -1) && (best == -1 || cost + mergecost < best)){
            best = cost + mergecost;
            best_sol = tmpsol;
          }

          G.expand_nodes(u_name, v_name, mg);
          G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
        }
        */

        //merge forbid 入れ替え(Graph copyで対応)
        if(merge_flag){
          Graph G_copy = G;
          vector <edge> tmpsol;
          G_copy.permanent(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, G_orig);
          MergeData mg(G_orig.num_nodes);
          int mergecost = G_copy.merge_nodes(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, mg, G_orig);
          int cost = naive_branching(G_copy, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
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
          best = naive_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
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
        best = naive_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
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
          int cost = naive_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
          if((cost != -1) && (best == -1 || cost + mergecost < best)){
            best = cost + mergecost;
            best_sol = tmpsol;
          }

          G.expand_nodes(u_name, v_name, mg);
          G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
        }
      }

      /* normal
      if(u_name > v_name) SWAP(int, u_name, v_name);
      int w = G.Weight(u,v);
      G.forbid(u, v, best_sol, G_orig);
      G.flip_edge(u,v);
      best = naive_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
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
        int cost = naive_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
        if((cost != -1) && (best == -1 || cost + mergecost < best)){
          best = cost + mergecost;
          best_sol = tmpsol;
        }

        G.expand_nodes(u_name, v_name, mg);
        G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
      }
      */

      if(best != -1) sol.insert(sol.end(), best_sol.begin(), best_sol.end());
      return best;
    }
  }

  if(flag == false){
    uv.clear();
    //auto start = system_clock::now();
    FOR(u, 0, G.num_nodes-1) FOR(v, u+1, G.num_nodes){
      if(G.Weight(u,v) <= 0 || G.Flag(u,v) != 0) continue;
      vector<int> s2_tmp;
      int cost = 0;
      FOR(z, 0, G.num_nodes){
        if(u == z || v == z) continue;
        if((G.Weight(u,z) > 0 && G.Weight(v,z) <= 0) || (G.Weight(u,z) <= 0 && G.Weight(v,z) > 0)){
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
      flag = true;
      bool merge_flag = true;
      FOR(z, 0, G.num_nodes){
        if(u_name == z || v_name == z) continue;
        if((G.flag[u_name][z] == -1 && G.flag[v_name][z] == 1) || (G.flag[u_name][z] == 1 && G.flag[v_name][z] == -1)) merge_flag = false;
      }

      //通常ver
      lp_sol[u][v] = -1;
      //heur
      if(lp_sol[u][v] == 0) {
        /*
        if(merge_flag){
          vector <edge> tmpsol;
          G.permanent(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, G_orig);
          MergeData mg(G_orig.num_nodes);
          int mergecost = G.merge_nodes(G.name_to_ind[u_name], G.name_to_ind[v_name], tmpsol, mg, G_orig);
          int cost = naive_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
          if((cost != -1) && (best == -1 || cost + mergecost < best)){
            best = cost + mergecost;
            best_sol = tmpsol;
          }

          G.expand_nodes(u_name, v_name, mg);
          G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
        }
        */

        //merge forbid 入れ替え(Graph copyで対応)
        if(merge_flag){
          Graph G_copy = G;
          vector <edge> tmpsol;
          G_copy.permanent(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, G_orig);
          MergeData mg(G_orig.num_nodes);
          int mergecost = G_copy.merge_nodes(G_copy.name_to_ind[u_name], G_copy.name_to_ind[v_name], tmpsol, mg, G_orig);
          int cost = naive_branching(G_copy, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
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
          best = naive_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
          if(best != -1) {
            best += w;
            if(best < max_obj) max_obj = best;
          }
          G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
          G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);
        }
      } else if (lp_sol[u][v] == 1) {
        /*
        if(u_name > v_name) SWAP(int, u_name, v_name);
        int w = G.Weight(u,v);
        G.forbid(u, v, best_sol, G_orig);
        G.flip_edge(u,v);

        best = naive_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
        if(best != -1) {
          best += w;
          if(best < max_obj) max_obj = best;
        }
        G.reset_flag(G.name_to_ind[u_name],G.name_to_ind[v_name]);
        G.flip_edge(G.name_to_ind[u_name],G.name_to_ind[v_name]);
        */

        //merge forbid 入れ替え
        if(u_name > v_name) SWAP(int, u_name, v_name);
        int w = G.Weight(u,v);
        G.forbid(u, v, best_sol, G_orig);
        G.flip_edge(u,v);

        best = naive_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
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
          int cost = naive_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
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

        best = naive_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
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
          int cost = naive_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
          if((cost != -1) && (best == -1 || cost + mergecost < best)){
            best = cost + mergecost;
            best_sol = tmpsol;
          }

          G.expand_nodes(u_name, v_name, mg);
          G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
        }
      }

      /* normal
      if(u_name > v_name) SWAP(int, u_name, v_name);
      int w = G.Weight(u,v);
      G.forbid(u, v, best_sol, G_orig);
      G.flip_edge(u,v);
      best = naive_branching(G, G_orig, max_obj-w, best_sol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
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
        int cost = naive_branching(G, G_orig, max_obj-mergecost, tmpsol, rec_depth+1, uv, lp_zero, lp_one, lp_sol);
        if((cost != -1) && (best == -1 || cost + mergecost < best)){
          best = cost + mergecost;
          best_sol = tmpsol;
        }

        G.expand_nodes(u_name, v_name, mg);
        G.reset_flag(G.name_to_ind[u_name],G.num_nodes-1);
      }

      */
      if(best != -1) sol.insert(sol.end(), best_sol.begin(), best_sol.end());
      return best;
    }
  }
}

/*
int random_pivot(Graph& G, const Graph& G_orig, std::vector <edge>& sol){
  rnd.seed(RANDOM_SEED);
  int best_cost = -1;
  vector <edge> best_sol;

  REP(t, NUM_PIVOT){
    int n = G.num_nodes;
    vector <int> nodes;
    REP(i, n) nodes.push_back(i);
    vector <edge> tmpsol;
    int cost = 0;

    while(n > 1){
      int k = rnd(0, n-1);
      int pivot = nodes[k];
      vector <int> neighbor;
      nodes.erase(nodes.begin()+k);
      n--;

      for(int i=n-1; i >= 0; i--){
        if(G.Weight(nodes[i],pivot) > 0) {
          neighbor.push_back(nodes[i]);
          if(G.Flag(nodes[i],pivot) == 0){
            G.permanent(nodes[i],pivot,sol,G_orig);
            G.reset_flag(nodes[i], pivot);
          }
          nodes.erase(nodes.begin()+i);
        }
      }
      n -= neighbor.size();

      if(neighbor.size() > 1){
        FOR(i,0, (int) neighbor.size()-1) FOR(j, i+1, (int) neighbor.size()){
          if(G.Flag(neighbor[i],neighbor[j]) == 0){
            G.permanent(neighbor[i],neighbor[j],sol,G_orig);
            G.reset_flag(neighbor[i], neighbor[j]);
          }

          if(G.Weight(neighbor[i],neighbor[j]) < 0){
            cost-= G.Weight(neighbor[i],neighbor[j]);
          }
        }
      }
      if(neighbor.size() > 0 && n > 0){
        REP(i, (int) neighbor.size()) REP(j, n){
          if(G.Flag(neighbor[i],nodes[j]) == 0){
            G.forbid(neighbor[i],nodes[j],sol,G_orig);
            G.reset_flag(neighbor[i],nodes[j]);
          }

          if(G.Weight(neighbor[i],nodes[j]) > 0){
            cost += G.Weight(neighbor[i],nodes[j]);
          }
        }
      }
    }
    if(best_cost == -1 || cost < best_cost){
      best_cost = cost;
      best_sol = tmpsol;
    }
  }

  if(best_cost != -1) sol.insert(sol.end(), best_sol.begin(), best_sol.end());
  return best_cost;
}
*/

/*
double lp_solve(const Graph &G, vector <vector <double> >& lp_sol){
  int n = G.num_nodes;
  lp_sol = vector <vector <double> > (n, vector <double>(n, 0));


  IloEnv env;	//enviromental object
  double obj_val = 0;

  env = IloEnv();
  IloModel model = IloModel(env);
  IloNumVarArray x (env);
  IloExprArray coef (env);
  IloExpr f (env);
  IloRangeArray constraint (env);
  //if(DEBUG) env.setOut(cerr);
  //else env.setOut(env.getNullStream());
  env.setOut(env.getNullStream());

  FOR (i, 0, n*(n-1)/2) x.add(IloNumVar(env , 0 , 1 , ILOFLOAT));

  FOR(i, 0, n-1) FOR(j, i+1, n){
    // ij corresponds to index i(2n-i-1)/2 + j - i - 1
    x.add(IloNumVar(env , 0 , 1 , ILOFLOAT) );
    if(G.Weight(i,j) > 0) f += x[i*(2*n-i-1)/2+j-i-1];
    else{
      f -= x[i*(2*n-i-1)/2+j-i-1];
      obj_val += 1;
    }
  }

  IloObjective obj(env, f, IloObjective::Minimize);
  model.add(obj);

  FOR(i, 0, n-1) FOR(j, i+1, n) FOR(k, 0, n){
    if(k == i || k == j) continue;
    if(i < k && j < k) constraint.add(x[i*(2*n-i-1)/2+j-i-1] - x[i*(2*n-i-1)/2+k-i-1] - x[j*(2*n-j-1)/2+k-j-1] <= 0);
    else if(i < k && j > k) constraint.add(x[i*(2*n-i-1)/2+j-i-1] - x[i*(2*n-i-1)/2+k-i-1] - x[k*(2*n-k-1)/2+j-k-1] <= 0);
    else if(i > k && j > k) constraint.add(x[i*(2*n-i-1)/2+j-i-1] - x[k*(2*n-k-1)/2+i-k-1] - x[k*(2*n-k-1)/2+j-k-1] <= 0);
  }
  model.add(constraint);

  IloCplex cplex = IloCplex(model);
  cplex.solve();

  obj_val += cplex.getObjValue();

  FOR(i, 0, n-1) FOR(j, i+1, n){
    lp_sol[i][j] = cplex.getValue(x[i*(2*n-i-1)/2+j-i-1]);
    lp_sol[j][i] = cplex.getValue(x[i*(2*n-i-1)/2+j-i-1]);
    // cout << i << "," << j << ": " << lp_sol[i][j] << endl;
  }

  return obj_val;
}

int lp_pivot(Graph &G,  const Graph& G_orig, vector <edge>& sol, const vector <vector <double> >& lp_sol){
  rnd.seed(RANDOM_SEED);
  int best_cost = -1;
  vector <edge> best_sol;

  REP(t, NUM_PIVOT){
    int n = G.num_nodes;
    vector <int> nodes;
    REP(i, n) nodes.push_back(i);
    vector <edge> tmpsol;
    int cost = 0;

    while(n > 1){
      int k = rnd(0, n-1);
      int pivot = nodes[k];
      vector <int> neighbor;
      nodes.erase(nodes.begin()+k);
      n--;
      neighbor.push_back(pivot);

      for(int i=n-1; i >= 0; i--){
        if(rnd(0,0xffffff) > lp_sol[nodes[i]][pivot] * 0xffffff) {
          neighbor.push_back(nodes[i]);
          nodes.erase(nodes.begin()+i);
          n--;
        }
      }

      if(neighbor.size() > 1){
        FOR(i,0, (int) neighbor.size()-1) FOR(j, i+1, (int) neighbor.size()){
          if(G.Flag(neighbor[i],neighbor[j]) == 0){
            G.permanent(neighbor[i],neighbor[j],sol,G_orig);
            G.reset_flag(neighbor[i], neighbor[j]);
          }
          if(G.Weight(neighbor[i],neighbor[j]) < 0){
            cost-= G.Weight(neighbor[i],neighbor[j]);
          }
        }
      }
      if(neighbor.size() > 0 && n > 0){
        REP(i, (int) neighbor.size()) REP(j, n){
          if(G.Flag(neighbor[i],nodes[j]) == 0){
            G.forbid(neighbor[i],nodes[j],sol,G_orig);
            G.reset_flag(neighbor[i],nodes[j]);
          }
          if(G.Weight(neighbor[i],nodes[j]) > 0){
            cost += G.Weight(neighbor[i],nodes[j]);
          }
        }
      }
    }
    if(best_cost == -1 || cost < best_cost){
      best_cost = cost;
      best_sol = tmpsol;
    }
  }


  if(best_cost != -1) sol.insert(sol.end(), best_sol.begin(), best_sol.end());
  return best_cost;
}
*/