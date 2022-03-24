#include <vector>
#include <list>
#include <iostream>
#include <algorithm>
#include "random.hpp"
#include "main.h"
#include "graph.h"
#include "randomPivot.h"

using namespace std;

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