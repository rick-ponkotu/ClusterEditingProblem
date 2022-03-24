#include <vector>
#include <iostream>
#include <algorithm>

#include "main.h"
#include "graph.h"

using namespace std;

void search_clique(const Graph& G, const int u, vector <int>& clique);
int check_unaffordable(const Graph&G, const int u, const int v, const int obj);
int reduction(Graph& G, const Graph& G_orig, const int obj, vector <edge>& sol);

//search clique in G
void search_clique(const Graph& G, const int u, vector<int>& clique) {
    int n = G.num_nodes;

    clique.clear();
    if(n <= 2) return;

    clique.push_back(u);

    //search neighbor
    for(int i=0; i < n; i++) {
        if(u != i && G.weight[u][i] > 0) clique.push_back(i);
    }

    if(clique.size() == 1) return;

    for(int i=0; i < clique.size(); i++) {
        for(int j=0; j < n; j++) {
            if(clique[i] == j || u == j)  continue;

            if(G.weight[u][j] > 0) {
                if(G.weight[clique[i]][j] <= 0){
                    clique.clear();
                    return;
                }
            }else {
                if(G.weight[clique[i]][j] > 0){
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

        if(G.weight[u][w] > 0 && G.weight[v][w] > 0) {
            sum_icf += min(G.weight[u][w], G.weight[v][w]);
        }
        else if(G.weight[u][w] > 0 || G.weight[v][w] > 0) {
            sum_icp += min(abs(G.weight[u][w]), abs(G.weight[v][w]));
        }
    }

    icf = max(0, G.weight[u][v]) + sum_icf;
    icp = max(0, -G.weight[u][v]) + sum_icp;

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
              if (G.flag[u][v] == 0) G.permanent(u, v, sol, G_orig);
              G.merge_nodes(u, v, sol, G_orig);
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

//output N(v)
vector<int> open_adjacency_set(Graph& G, int v) {
    int N = G.num_nodes;
    vector<int> N_v;

    for(int i=0; i < N; i++) {
        if(G.weight[v][i] > 1) N_v.push_back(G.weight[v][i]);
    }
    return N_v;
}

//delta(v) = \Pi({{x,y} | x,y \in N(v), [x,y] \notin E}). delta is positive.
int delta_v(Graph& G, vector<int> adj_v) {
    int delta_cost = 0;

    for(auto& s: adj_v) {
        for(auto& t: adj_v) {
            if(s != t) {
                if(G.weight[s][t] < 0) delta_cost += (-1) *G.weight[s][t];
            }
        }
    }
    return delta_cost;
}

//gamma(v) is total weight of edges to make N[v] disjoint.
int gamma_v(Graph& G, vector<int> closed_adj_set) {
    int gamma_cost = 0;
    vector<int> diff;
    set_difference(G.node_names.begin(), G.node_names.end(), closed_adj_set.begin(), closed_adj_set.end(), inserter(diff, diff.end()));

    for(auto& s: diff) {
        for(auto& t: closed_adj_set) {
            if(G.weight[s][t] > 0) gamma_cost += G.weight[s][t];
        }
    }
    return gamma_cost;
}

//output E(x, N[v]). E(x,N[v]) is the set of edges that has one end in v and the other in N[v].
vector<int> edge_from_x_to_openNv(Graph& G, int x, int v) {
    vector<int> N_v;
    vector<int> end_points;
    N_v = open_adjacency_set(G, v);

    for(auto& y: N_v) {
        if(G.weight[x][y] > 0)  end_points.push_back(y);
    }
    return end_points;
}

//add edges to make G[N[v]] a clique
int make_clique(Graph& G, vector<int> closed_N_v, vector <edge>& sol) {
    int clique_cost = 0;

    for(auto& s: closed_N_v) {
        for(auto& t: closed_N_v) {
            if(s != t) {
                if(G.weight[s][t] < 0) {
                    G.weight[s][t] = 1;
                    clique_cost += 1;
                }
            }
        }
    }
    return clique_cost;
}

//output N(N[v])
vector<int> adj_adj_set(Graph& G, vector<int> closed_adj_set) {
    vector<int> adj_Nv;
    vector<int> diff;
    set_difference(G.node_names.begin(), G.node_names.end(), closed_adj_set.begin(), closed_adj_set.end(), inserter(diff, diff.end()));

    for(auto& s: diff) {
        for(auto& t: closed_adj_set) {
            if(G.weight[s][t] > 0) adj_Nv.push_back(s);
        }
    }
    return adj_Nv;
}

int kernelization_EdgeCuts(Graph& G, const Graph& G_orig, vector <edge>& sol) {
    vector<int> adj_set;
    int cost = 0;
    int delta;
    int gamma;

    for(auto& v: G.node_names) {
        adj_set = open_adjacency_set(G, v); //open adjacency set
        delta = delta_v(G, adj_set);
        adj_set.push_back(v); //closed adjacency set
        gamma = gamma_v(G, adj_set);

        //Step 1
        if(2*delta + gamma < adj_set.size()) {
            cost += make_clique(G, adj_set, sol);


            //Step 2
            adj_set = open_adjacency_set(G, v);
            delta = delta_v(G, adj_set);
            adj_set.push_back(v);
            gamma = gamma_v(G, adj_set);

            if(2*delta + gamma < adj_set.size()) {
                vector<int> adj_Nv;
                adj_Nv = adj_adj_set(G, adj_set);
                for(auto& x: adj_Nv) {
                    vector<int> edge;
                    edge = edge_from_x_to_openNv(G,x,v);
                    int edge_cost = 0;
                    for(auto& e: edge) {
                        edge_cost += G.weight[x][e];
                    }
                    if(edge_cost <= adj_set.size() / 2) {
                        for(auto& e: edge) {
                            G.forbid(x, e, sol, G_orig);
                            G.delete_edge(x, e);
                        }
                    }
                }


                //Step 3
                adj_set = open_adjacency_set(G, v);
                delta = delta_v(G, adj_set);
                adj_set.push_back(v);
                gamma = gamma_v(G, adj_set);

                if(2*delta + gamma < adj_set.size()) {
                    adj_Nv = adj_adj_set(G, adj_set);
                    if(adj_Nv.size() > 0) {
                        // merge N[v] into a single vertex v'
                        sort(adj_set.begin(), adj_set.end(), greater<int>());
                        for(int i=0; i < adj_set.size()-1; i++) {
                            G.permanent(adj_set[i], adj_set[i+1], sol, G_orig);
                            cost += G.merge_nodes(adj_set[i], adj_set[i+1], sol, G_orig);
                            return cost;
                        }
                    }
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
