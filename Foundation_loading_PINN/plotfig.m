figure(1)
clf
axis equal
for i = 1:nEle
    nodeplot=Element_information(i,2:4);
    xplot=[Node_information(nodeplot,2)];yplot=[Node_information(nodeplot,3)];
    U2plot=[U(nodeplot+nNode)];
    patch(xplot,yplot,U2plot,"EdgeColor","none");
    hold on
end
colormap jet;
colorbar
title("Vertical displacement-contour")
hold off;

figure(2)
clf;
axis equal
q_node=zeros(nNode,1);
for i = 1:nNode
    [row, col]=find(Element_information(:,2:4)==i);
    q_node(i,1)=mean(q(row));
end

for i =1:nEle
    nodeplot=Element_information(i,2:4);
    xplot=[Node_information(nodeplot,2)];yplot=[Node_information(nodeplot,3)];
    qplot=q_node(nodeplot);
    patch(xplot,yplot,qplot,"EdgeColor","none");
    hold on
end
colormap jet;
colorbar
title("Deviatoric stress")
hold off;


figure(3)
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
    patch(xplot,yplot,edplot,"EdgeColor","none");
    hold on
end
colormap jet;
colorbar
title("Plastic deviatoric strain")
hold off;

figure(4)
clf;
axis equal
ev_node=zeros(nNode,1);
for i = 1:nNode
    [row, col]=find(Element_information(:,2:4)==i);
    ev_node(i)=mean(epv(row));
end

for i =1:nEle
    nodeplot=Element_information(i,2:4);
    xplot=[Node_information(nodeplot,2)];yplot=[Node_information(nodeplot,3)];
    evplot=ev_node(nodeplot);
    patch(xplot,yplot,evplot,"EdgeColor","none");
    hold on
end
colormap jet;
colorbar
title("Plastic volumetric strain")
hold off;
