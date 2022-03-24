import sys
import networkx as nx
import matplotlib.pyplot as plt
import random

A = [255,0,0,218,188,238,0,64,128,255,255,255,106,140]
B = [0,255,0,112,143,130,191,224,128,215,160,20,90,80]
D = [0,0,255,214,143,238,255,208,0,0,0,147,205,45]


G = nx.read_edgelist(sys.argv[2], comments="p", nodetype=int)
Gcopy = G.copy()

if len(sys.argv) >= 4:
    f = open(sys.argv[3], 'r')
    d = f.readlines()
    for j in d:
        aa, bb = j.split()
        a = int(aa)
        b = int(bb)

        if (a,b) in Gcopy.edges:
            Gcopy.remove_edge(a,b)
        else:
            Gcopy.add_edge(a,b)

C = nx.connected_components(Gcopy)


t = 0
for i in C:
    if t < 14:
        r = A[t]
        g = B[t]
        b = D[t]
    else:
        r = random.randint(0,255)
        g = random.randint(0,255)
        b = random.randint(0,255)
    for j in i:
        G.nodes[j]['viz'] = {}
        G.nodes[j]['viz']['color'] = {}
        G.nodes[j]['viz']['color']['r']= r
        G.nodes[j]['viz']['color']['g']= g
        G.nodes[j]['viz']['color']['b']= b
        G.nodes[j]['viz']['color']['a']= 1
        G.nodes[j]['id'] = j
    t += 1


nx.write_gexf(G, sys.argv[1])
