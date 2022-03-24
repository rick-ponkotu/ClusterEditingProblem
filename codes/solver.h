
// return:
//  - the opt val if it is <= max_obj
//  - -1 if opt val > max_obj
// assumption: max_obj >= 0

//int naive_branching(Graph& G, const Graph& G_orig, int max_obj, std::vector <edge>& sol, int rec_depth);
int naive_branching(Graph& G, const Graph& G_orig, int max_obj, std::vector <edge>& sol, int rec_depth, std::vector<std::vector<int>>& uv, std::vector<std::vector<int>>& lp_zero, std::vector<std::vector<int>>& lp_one, std::vector<std::vector<double>> lp_sol);
int edge_branching(Graph& G, const Graph& G_orig, int max_obj, std::vector <edge>& sol);

int random_pivot(Graph& G, const Graph& G_orig, std::vector <edge>& sol);

void ret_cnt();

double ret_time();

double lp_solve(const Graph &G, std::vector <std::vector <double> >& lp_sol);

int lp_pivot(Graph &G,  const Graph& G_orig, std::vector <edge>& sol, const std::vector <std::vector <double> >& lp_sol);
int lp_pivot(Graph &G,  const Graph& G_orig, std::vector <edge>& sol, const std::vector <std::vector <double> >& lp_sol);
