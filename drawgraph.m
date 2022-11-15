i=950 ;
file1=['node_',num2str(i),'.mat'] ;
file2=['input_',num2str(i),'.mat'] ;
file3=['output_',num2str(i),'.mat'] ;
file4=['phasefieldsimulator_',num2str(i),'.mat'] ;
A1x=load(file1) ;
A1=A1x.node2;
A2x=load(file2) ;
A2=A2x.nodeinputfeature1 ;
A3x=load(file3);
A3=A3x.nodeoutputfeature ;
A4x=load(file4) ;
cord=A4x.cord ;
sz1=size(A1,1)*0.25;
sz2=size(cord,1) ;
j=2; 
t = A1(1:sz1,1);
h = A1(1:sz1,2);
g = graph(t,h);
p = plot(g);
g.Nodes.value = 0.8*A3((j-1)*sz2+1:j*sz2,3) ;
g.Nodes.NodeColors = g.Nodes.value ;
p.NodeCData = g.Nodes.NodeColors ;
x = cord(:,1) ;
y = cord(:,2) ;
q= plot(g,'XData',x,'YData',y,'NodeCData',g.Nodes.NodeColors);
daspect([1 1 1]) ;
colormap(jet(256))
colorbar
set(gca,'Visible','off')
q.MarkerSize=2 ;
q.EdgeColor = 'k';
q.LineWidth=1 ;
