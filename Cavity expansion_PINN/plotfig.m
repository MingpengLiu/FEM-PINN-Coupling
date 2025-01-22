figure(1)
clf
axis equal
for i = 1:nEle
    nodeplot=Element_information(i,2:4);
    xplot=[Node_information(nodeplot,2)];yplot=[Node_information(nodeplot,3)];
    Uplot=sqrt(U(nodeplot+nNode).^2+U(nodeplot).^2);
    patch(xplot,yplot,Uplot);
    hold on
end
colormap jet;
colorbar
title("Displacement-contour")
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
    patch(xplot,yplot,qplot);
    hold on
end
colormap jet;
colorbar
clim([1 73])
title("Deviatoric stress")
hold off;

figure(3)
clf;
axis equal
ed_node=zeros(nNode,1);
for i = 1:nNode
    [row, col]=find(Element_information(:,2:4)==i);
    ed_node(i)=mean(epd(row));
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
clim([0 0.04])
title("Plastic deviatoric strain")
hold off;


