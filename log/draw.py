import matplotlib.pyplot as plt
import networkx as nx

G = nx.DiGraph()
G = nx.read_edgelist('pivot005.txt', nodetype=int)

pos = nx.spring_layout(G)
nx.draw_networkx(G, pos, with_labels=True, alpha=0.5)
plt.axis("off")
plt.show()