figure(1)
clf;
axis equal
ed_node=zeros(nNode,1);
for i = 1:nNode
    [row, col]=find(Element_information(:,2:4)==i);
    ed_node(i)=mean(ed(row));
end

for i =1:nEle
    nodeplot=Element_information(i,2:4);
    xplot=[Node_information(nodeplot,2)];yplot=[Node_information(nodeplot,3)];
    edplot=ed_node(nodeplot);
    patch(xplot,yplot,edplot);
    hold on
end
colormap jet;
colorbar
title("Deviatoric strain")
hold off;

