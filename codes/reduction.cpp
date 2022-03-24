#include <vector>
#include <iostream>
#include <algorithm>
#include <list>

#include "main.h"
#include "graph.h"

using namespace std;

void search_clique(const Graph& G, const int u, vector <int>& clique);
int check_unaffordable(const Graph&G, const int u, const int v, const int obj);
int reduction(Graph& G, const Graph& G_orig, const int obj, vector <edge>& sol);
void show_vec(vector <int> vec);

//search clique in G
void search_clique(const Graph& G, const int u, vector<int>& clique) {
    int n = G.num_nodes;

    clique.clear();
    if(n <= 2) return;

    clique.push_back(u);

    //search neighbor
    for(int i=0; i < n; i++) {
        if(u != i && G.Weight(u,i) > 0) clique.push_back(i);
    }

    if(clique.size() == 1) return;

    for(int i=0; i < clique.size(); i++) {
        for(int j=0; j < n; j++) {
            if(clique[i] == j || u == j)  continue;
            if(G.Weight(u,j) > 0) {
                if(G.Weight(clique[i],j) <= 0){ 
                    clique.clear();
                    return;
                }
            }else {
                if(G.Weight(clique[i],j) > 0){ 
                    clique.clear();
                    return;
                }
            }
        }
    }
    return;
}

/*
set permanent or forbidden for each u, v
*/
int check_unaffordable(const Graph& G, const int u, const int v, const int obj) {
    int N = G.num_nodes;
    int icf;
    int icp;
    int sum_icf = 0;
    int sum_icp = 0;

    for(int w=0; w < N; w++) {
        if(u == w || v == w) continue;

        if(G.Weight(u,w) > 0 && G.Weight(v,w) > 0) {
            sum_icf += min(G.Weight(u,w), G.Weight(v,w));
        }
        else if(G.Weight(u,w) > 0 || G.Weight(v,w) > 0) {
            sum_icp += min(abs(G.Weight(u,w)), abs(G.Weight(v,w)));
        }
    }

    icf = max(0, G.Weight(u,v)) + sum_icf;
    icp = max(0, -G.Weight(u,v)) + sum_icp;

    if(icf > obj && icp > obj) {
        return 0;
    }else if(icf > obj) {
        return icp;
    }else if(icp > obj) {
        return -1;
    }else {
        return 0;
    }
    return 0;
}

/*
graph reduction
*/
int reduction(Graph& G, const Graph& G_orig, const int obj, vector <edge>& sol) {
    int N = G.num_nodes;
    if(obj < 0) return -1;

    for(int u=0; u < N; u++) {
        vector<int> clique;
        search_clique(G, u, clique);
        if(clique.size() > 1) {
            G.delete_nodes(clique, sol, G_orig);
            return 0;
        }
        else if(clique.size() == 1){
            G.delete_node(u, sol ,G_orig);
            return 0;
        }
    }

    for(int u=0; u < N-1; u++) {
      for(int v=u+1; v < N; v++) {
          int cost = check_unaffordable(G, u, v, obj);
          if(cost > 0) {
              if (G.Flag(u,v) == 0) G.permanent(u, v, sol, G_orig);
              MergeData mg(G_orig.num_nodes);
              G.merge_nodes(u, v, sol, mg, G_orig);
              return cost;
          }
      }
    }

    return -1;
}

int cal_reduction(Graph& G, const Graph& G_orig, const int obj, vector <edge>& sol) {
    int cost = 0;
    int flag = 0;
    do {
        cost += flag;
        flag = reduction(G, G_orig, obj-cost, sol);
    } while(flag != -1 && G.num_nodes >= 3);

    return cost;
}


//kernelization based on edge cuts
int kernelization_EdgeCuts(Graph& G, const Graph& G_orig, vector <edge>& sol) {
    vector<int> adj_set;
    int cost = 0;
    int delta;
    int gamma;

    for(auto v: G.node_names) {
        if(G.neighbors[v].size() == 0){
            G.delete_node(G.name_to_ind[v], sol, G_orig);
            return 0;
        }

        delta = 0;
        gamma = 0;
        vector<int> adj_Nv;
        vector<int> gamma_nn(G_orig.num_nodes,0);
        vector <vector <int>> nn(G_orig.num_nodes);

        for(auto i : G.node_names){
            if(i == v || G.weight[i][v] > 0) continue;
            for(auto j: G.neighbors[i]){
                if(G.weight[v][j] > 0){
                    gamma += G.weight[i][j];
                    gamma_nn[i] +=G.weight[i][j];
                    if(nn[i].size() == 0) adj_Nv.push_back(i);
                    nn[i].push_back(j);
                }
            }
        }
        bool flag_neighbor = false;
        if(G.neighbors[v].size() >= 2){
            for(auto i = G.neighbors[v].begin(); next(i) != G.neighbors[v].end(); i++) {
                for(auto j = next(i); j != G.neighbors[v].end(); j++) {
                    if(G.weight[*i][*j] <= 0){
                        if (G.flag[*i][*j] == -1){
                            flag_neighbor = true;
                            break;
                        }
                        delta -= G.weight[*i][*j];
                    }
                }   
                if(flag_neighbor) break;
            }
        }
        if(flag_neighbor) continue;

        if(2*delta + gamma < G.neighbors[v].size()+1) {
            //Step 1
            if(G.neighbors[v].size() >= 2){
                for(auto i = G.neighbors[v].begin(); next(i) != G.neighbors[v].end(); i++) {
                    for(auto j = next(i); j != G.neighbors[v].end(); j++) {
                        if(G.flag[*i][*j] == 0) G.permanent(G.name_to_ind[*i], G.name_to_ind[*j], sol, G_orig);
                        if(G.weight[*i][*j] <= 0)G.flip_edge(G.name_to_ind[*i], G.name_to_ind[*j]);
                    }
                }
            }
            cost += delta;

            //Step 2
            int x_remain = -1;
            for(auto& x: adj_Nv) {
                if(gamma_nn[x] <= (G.neighbors[v].size()+1) / 2) {
                    for(auto j: nn[x]) {
                        G.forbid(G.name_to_ind[x], G.name_to_ind[j], sol, G_orig);
                        G.flip_edge(G.name_to_ind[x], G.name_to_ind[j]);
                    }
                    cost += gamma_nn[x];
                }
                else x_remain = x;
            }

            //Step 3
            if(x_remain != -1) {
                // merge N[v] into a single vertex v'
                if(G.neighbors[v].size() > 0){
                    vector <int> neighb;
                    for(auto i : G.neighbors[v]) neighb.push_back(i);
                    sort(neighb.begin(), neighb.end(), greater<int>());
                    REP(i, neighb.size()-1) {
                        if(G.flag[neighb[i]][neighb[i+1]] == 0) G.permanent(G.name_to_ind[neighb[i]], G.name_to_ind[neighb[i+1]], sol, G_orig);
                        MergeData mg(G_orig.num_nodes);
                        cost += G.merge_nodes(G.name_to_ind[neighb[i]], G.name_to_ind[neighb[i+1]], sol, mg, G_orig);
                    }
                }
                
                if(G.flag[v][*G.neighbors[v].begin()] == 0) G.permanent(G.name_to_ind[v], G.name_to_ind[*G.neighbors[v].begin()], sol, G_orig);
                MergeData mg(G_orig.num_nodes);
                cost += G.merge_nodes(G.name_to_ind[v], G.name_to_ind[*G.neighbors[v].begin()], sol, mg, G_orig);

                return cost;

            }else {
                if(G.neighbors[v].size() > 0)  {
                    vector <int> adj_set;
                    adj_set.push_back(G.name_to_ind[v]);
                    for(int i : G.neighbors[v]) adj_set.push_back(G.name_to_ind[i]);
                    G.delete_nodes(adj_set, sol, G_orig);
                    return cost;
                }else{
                    G.delete_node(G.name_to_ind[v], sol, G_orig);
                    return cost;
                }
            }
        }
    }
    return cost;
}


int cal_ker(Graph& G, const Graph& G_orig, vector <edge>& sol) {
    int cost = 0;
    int flag = 0;
    do {
        cost += flag;
        flag = kernelization_EdgeCuts(G, G_orig, sol);
    } while(flag > 0);
    return cost;
}


void show_vec(vector <int> vec) {
    for(auto& v: vec) {
        cout << v << ' ';
    }
    cout << endl;
}
