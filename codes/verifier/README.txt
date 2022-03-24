To compile this verifier for Cluster Editing execute "make all" in the command line.
Afterwards, there is an executable "verifier". The verifier expects 
as first argument a path to graph file in the following format:

https://pacechallenge.org/2021/tracks/#input-format

The second argument is expected to be a path to a solution file 
containing a list of edges in the following format:

https://pacechallenge.org/2021/tracks/#output-format-for-exact-track-and-heuristic-track

The verifier checks whether the graph is a cluster graph after the edges in the solution file are edited. The verifier does not check for optimality.
