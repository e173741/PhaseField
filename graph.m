
file1=['node_',num2str(1),'.mat'] ;
file2=['input_',num2str(1),'.mat'] ;
file3=['output_',num2str(1),'.mat'] ;
file4=['phasefieldsimulator_',num2str(1),'.mat'] ;
A1x=load(file1) ;
A1=A1x.node2;
A2x=load(file2) ;
A2=A2x.nodeinputfeature ;
A3x=load(file3);
A3=A3x.nodeoutputfeature ;
A4x=load(file4) ;
cord=A4x.cord ;
t = A1(:,1);
h = A1(:,2);
g = graph(t,h);
p = plot(g);
g.Nodes.value = A2(:,3) ;
g.Nodes.NodeColors = g.Nodes.value;
p.NodeCData = g.Nodes.NodeColors;
x = cord(:,1);
y = cord(:,2) ;
q= plot(g,'XData',x,'YData',y,'NodeCData',g.Nodes.NodeColors);
colorbar
q.MarkerSize=10 ;
q.EdgeColor = 'k';