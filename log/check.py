import csv
with open('exact005.csv') as f:
    reader = csv.reader(f, delimiter=' ')
    l = [row for row in reader]
f.closed

n = 20
graphEdge = [[-1 for i in range(n)] for j in range(n)]

for e in l:
    u = int(e[0]) - 1
    v = int(e[1]) - 1
    graphEdge[u][v] = 1
    graphEdge[v][u] = 1
#print(graphEdge)

with open('exact005_sol.csv') as f2:
    reader = csv.reader(f2, delimiter=' ')
    l2 = [row for row in reader]

#print(l2)

for e2 in l2:
    u = int(e2[0])
    v = int(e2[1])
    if graphEdge[u][v] == 1:
        graphEdge[u][v] = -1
        graphEdge[v][u] == -1
    else:
        graphEdge[u][v] = 1
        graphEdge[v][u] = 1
f2.closed
#print(graphEdge)

with open('lp_exact05.txt') as f:
    reader = csv.reader(f, delimiter=' ')
    lp = [row for row in reader]
f.closed

one = 0
zero = 0
Edge = 0
noEdge = 0
all_one = 0
all_zero = 0

for i in range(n):
    for j in range(n):
        if graphEdge[i][j] == 1:
            Edge += 1
        if lp[i][j] == '1':
            all_one += 1
        if graphEdge[i][j] == -1:
            noEdge += 1
        if lp[i][j] == '0':
            all_zero += 1
        if graphEdge[i][j] == -1 and lp[i][j] == '1':
            one += 1
        if graphEdge[i][j] == 1 and lp[i][j] == '0':
            zero += 1




print(one, noEdge,one*100/noEdge)
print(zero,Edge,zero*100/Edge)
print(one*100/all_one)
print(zero*100/all_zero)