void ret_rec_cnt();
int default_branching(Graph& G, const Graph& G_orig,int max_obj, std::vector <edge>& sol, std::vector<std::vector<int>>& uv);
int b_branching(Graph& G, const Graph& G_orig,int max_obj, std::vector <edge>& sol, std::vector<std::vector<int>>& uv);
int bc_branching(Graph& G, const Graph& G_orig,int max_obj, std::vector <edge>& sol, int rec_depth, std::vector<std::vector<int>>& uv);
int c_branching(Graph& G, const Graph& G_orig,int max_obj, std::vector <edge>& sol, int rec_depth, std::vector<std::vector<int>>& uv);