void ret_rec_cnt_v3();
int bd_branching(Graph& G, const Graph& G_orig,int max_obj, std::vector <edge>& sol, std::vector<std::vector<int>>& uv, std::vector<std::vector<double>> lp_sol);
int bcd_branching(Graph& G, const Graph& G_orig, int max_obj, std::vector <edge>& sol, int rec_depth, std::vector<std::vector<int>>& uv, std::vector<std::vector<double>> lp_sol);

int d_branching(Graph& G, const Graph& G_orig,int max_obj, std::vector <edge>& sol, std::vector<std::vector<int>>& uv, std::vector<std::vector<double>> lp_sol);
int cd_branching(Graph& G, const Graph& G_orig,int max_obj, std::vector <edge>& sol, int rec_depth, std::vector<std::vector<int>>& uv, std::vector<std::vector<double>> lp_sol);