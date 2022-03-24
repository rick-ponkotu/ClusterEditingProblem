import csv

with open('lp_exact007.csv') as f:
    reader = csv.reader(f, delimiter=' ')
    l = [row for row in reader]

out = []
for i in range(len(l)):
    for j in range(i+1, len(l[i])):
        if l[i][j] != '1':
            out.append(f"{i+1} {j+1} {l[i][j]}\n")
            #F.write(f"{i+1} {j+1} {l[i][j]}\n")

with open('lp_exact007.txt', mode='w') as F:
    for e in out:
        F.write(e)