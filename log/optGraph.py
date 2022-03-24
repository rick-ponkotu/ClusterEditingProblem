import networkx as nx

G = nx.read_edgelist('input003.txt', nodetype=int)


f = open('exact003sol.txt', 'r')
d = f.readlines()
for j in d:
    aa, bb = j.split()
    a = int(aa)
    b = int(bb)

    if (a,b) in G.edges:
        G.remove_edge(a,b)
    else:
        G.add_edge(a,b)

nx.nx_agraph.view_pygraphviz(G, prog='fdp')