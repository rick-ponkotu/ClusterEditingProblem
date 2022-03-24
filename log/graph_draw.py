import matplotlib.pyplot as plt
import networkx as nx

G = nx.DiGraph()
G = nx.read_weighted_edgelist('lp_exact007.txt', nodetype=int)

pos = nx.spring_layout(G, k=0.7)

'''
#weight 0.5以下のedgeを消す
removeEdges = []

for (u,v,d) in G.edges(data=True):
    if d['weight'] == 0.5:
        removeEdges.append([u,v])

#グラフから削除
G.remove_edges_from(removeEdges)
'''
nx.nx_agraph.view_pygraphviz(G, prog='fdp')