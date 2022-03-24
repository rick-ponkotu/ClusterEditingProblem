from graphviz import Graph
# 無向グラフ
g = Graph(format='png')
g.attr("node", style="filled", fontname='MS Gothic')
g.node('A', 'アーサー王', shape="circle", color="pink")
g.node('B', 'ベディヴィア卿', shape="square", color="cyan")
g.node('C', '勇敢なランスロット卿 ', color="lightgreen")

# edgeを追加
g.edge('A', 'B')
g.edge('B', 'C')
g.edge('C', 'A')
g.render("sample", "png", view=True)
