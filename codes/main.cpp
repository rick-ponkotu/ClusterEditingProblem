#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <time.h>

#include "cmdline.h"

#include "main.h"
#include "graph.h"
#include "solver.h"
#include "reduction.h"
#include "exactBranching.h"
#include "heuristicBranching.h"
#include "lpUsedBranching.h"


using namespace std;

int main(int argc, char *argv[]){
  // parser
  cmdline::parser op;
  op.add<string>("input", 'i', "input file name", false);
  op.add<string>("output", 'o', "output fine name", false);
  op.add("reduction", 'r', "reduction-only mode");
  op.add("help", 0, "print help");

  if (!op.parse(argc, argv)||op.exist("help")){
    std::cout<<op.error_full()<<op.usage();
    return 0;
  }

  // input
  Graph Gin;

  if(op.exist("input")){
    const char *fn = op.get<string>("input").c_str();
    Gin = Graph(fn);
  }
  else Gin = Graph("../instances/exact/exact007.gr");

  vector <edge> sol;

  //////////////////// begin
  clock_t c_start = clock();

  Graph G = Gin;

  vector <edge> random_pivot_sol;
  //clock_t h1_start = clock();
  int random_pivot_obj = random_pivot(G, Gin, random_pivot_sol);
  //clock_t h1_end = clock();

  vector <vector <double> > lp_sol;
  double lp_obj = lp_solve(Gin, lp_sol);

  vector <edge> sol1;
  //clock_t h2_start = clock();
  int obj1 = lp_pivot(G, Gin, sol1, lp_sol);
  //clock_t h2_end = clock();

/*
  int rec_depth_heuristic = 0;
  Graph G_heuristic = Gin;
  vector <edge> sol_heuristic;
  clock_t h3_start = clock();
  vector<vector<int> > uv_heuristic;
  uv_heuristic.clear();
  int obj_heuristic = heuristic_branching(G_heuristic, Gin, random_pivot_obj, sol_heuristic, rec_depth_heuristic, uv_heuristic, lp_sol);
  clock_t h3_end = clock();
*/

  //cout << "time: " << (double) (h1_end - h1_start) / CLOCKS_PER_SEC << endl;
  //cout << "time: " << (double) (h2_end - h2_start) / CLOCKS_PER_SEC << endl;
  //cout << "time: " << (double) (h3_end - h3_start) / CLOCKS_PER_SEC << endl;
  //cout << op.get<string>("input") << " " <<lp_obj << " " << obj1 << " " << random_pivot_obj << " " << obj_heuristic << endl;
  //cout << op.get<string>("input") << " " <<lp_obj << " " << obj1 << " " << random_pivot_obj << endl;

/*
  if (obj_heuristic != -1 && obj1 > obj_heuristic) {
    obj1 = obj_heuristic;
    sol1 = sol_heuristic;
  }
*/

  if(obj1 > random_pivot_obj){
    obj1 = random_pivot_obj;
    sol1 = random_pivot_sol;
  }

  vector <edge> sol2;
  int red_cost = cal_reduction(G, Gin, obj1, sol2);
  red_cost += cal_ker(G, Gin, sol2);
  //cout << red_cost << endl;
  if(op.exist("reduction")){
    clock_t c_end = clock();
    cout << op.get<string>("input") << " " << obj1 << " " << red_cost << " " << Gin.num_nodes << " " << G.num_nodes << " " << (double) (c_end - c_start)/ CLOCKS_PER_SEC << endl;
    return 0;
  }




/*
  int obj2;
  int obj_default = 0;
  clock_t c1_start = clock();
  while(true) {
    Graph G_default = Gin;
    vector<edge> sol_default;
    uv.clear();
    obj2 = d_branching(G_default, Gin, obj_default, sol_default, uv, lp_sol);
    //obj2 = c_branching(G_default, Gin, obj_default, sol_default, rec_depth, uv);
    //obj2 = default_branching(G_default, Gin, obj_default, sol_default, uv);
    //obj2 = v0_branching(G_default, Gin, obj_default, sol_default, uv);
    //int obj2 = v1_branching(G, Gin, obj1-red_cost, sol2, rec_depth, uv);
    //obj2 = d_branching(G, Gin, obj1-red_cost, sol2, uv, lp_sol);
    //int obj2 = v3_branching(G, Gin, obj1-red_cost, sol2, rec_depth, uv, lp_sol);
    if(obj2 != -1) break;
    obj_default ++;
  }
*/

/*
  int rec_depth = 0;
  vector<vector<int> > uv;
  uv.clear();
  clock_t c1_start = clock();
  //int obj2 = default_branching(G, Gin, obj1-red_cost, sol2, uv);
  //int obj2 = b_branching(G, Gin, obj1-red_cost, sol2, uv);
  //int obj2 = c_branching(G, Gin, obj1-red_cost, sol2, rec_depth, uv);
  //int obj2 = d_branching(G, Gin, obj1-red_cost, sol2, uv, lp_sol);
  //int obj2 = bc_branching(G, Gin, obj1-red_cost, sol2, rec_depth, uv);
  //int obj2 = bd_branching(G, Gin, obj1-red_cost, sol2, uv, lp_sol);
  //int obj2 = cd_branching(G, Gin, obj1-red_cost, sol2, rec_depth, uv, lp_sol);
  //int obj2 = bcd_branching(G, Gin, obj1-red_cost, sol2, rec_depth, uv, lp_sol);
  clock_t c1_end = clock();

  if (obj1 > obj2 + red_cost && obj2 != -1){
    obj1 = obj2 + red_cost;
    sol1 = sol2;
  }

  clock_t c_end = clock();
*/

  if(DEBUG){
    //double elapsed = (double) (c_end - c_start) / CLOCKS_PER_SEC;
    //cout << op.get<string>("input") << " "<< obj1 << " "<< (double) (c_end - c_start)/ CLOCKS_PER_SEC << endl;
    //cout << op.get<string>("input") << " "<< lp_obj << " " << random_pivot_obj << " " << obj1 << endl;
    //ret_cnt();
    //ret_rec_cnt();
    //ret_rec_cnt_v3();
    //write_sol(sol2, "../log/exact003sol.txt");
    //cout << (ret_time() / elapsed) * 100 << endl;
    cout << op.get<string>("input") << endl;
    cout << "reduction cost " << red_cost << endl;
    //cout << "branching time: " << (double) (c1_end - c1_start) / CLOCKS_PER_SEC << endl;
    //ret_rec_cnt();
  }

  if(op.exist("output")){
    const char *fo = op.get<string>("output").c_str();
    write_sol(sol1,fo);
  }
}


void show_sol(const vector <edge>& sol){
  for (auto u: sol){
    cout << u.first+1 << " " << u.second+1 << endl;
  }

  return;
}

void write_sol(const vector <edge>& sol, const char *fo){
  ofstream ot(fo);

  for (auto u: sol){
    ot << u.first+1 << " " << u.second+1 << endl;
  }
  ot.close();
}
