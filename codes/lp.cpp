#include <vector>
#include <list>
#include <iostream>
#include <algorithm>
#include "random.hpp"
#include "main.h"
#include "graph.h"
#include "lp.h"

#include <ilcplex/ilocplex.h>
using namespace std;

double lp_solve(const Graph &G, vector <vector <double> >& lp_sol){
    int n = G.num_nodes;
    lp_sol = vector <vector <double> > (n, vector <double>(n, 0));

    /***/
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