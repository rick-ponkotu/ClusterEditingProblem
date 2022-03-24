double lp_solve(const Graph &G, std::vector <std::vector <double> >& lp_sol);
int lp_pivot(Graph &G,  const Graph& G_orig, std::vector <edge>& sol, const std::vector <std::vector <double> >& lp_sol);