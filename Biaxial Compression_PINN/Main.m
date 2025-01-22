clear;clc;
%% loading neural network
NNcell = load("netcell.mat");
NNcell = NNcell.netcell;
%% meshing, node and element
node1=[0 0];node2=[0 1];node3=[1 0];node4=[1 1];
Node_information = [1 node1;2 node2;3 node3;4 node4];
Element_information = [1 1 2 3; 2 4 2 3];
nNode = length(Node_information);
nEle = size(Element_information,1);
%% loading
% boundary
R=0.6;
dis1 = 0.15;dis2=R*dis1;
nload = 100;
F_ext = zeros(2*nNode,1);
F_ext_b = zeros(2*nNode,1);

% initial state
p0=100;e0=0.85;
e_pre = repelem(e0,nEle)';
plstrain_pre = zeros(nEle,4);
ini_stress = [-p0 -p0 -p0 0];        % in plane strain condition
stress_pre = repmat(ini_stress,nEle,1);

% save data
U_pre =zeros(2*nNode,1);U = zeros(2*nNode,1);save_stress=zeros(nEle,4);save_e=zeros(nEle,1);
save_strain=zeros(nEle,3);save_strainp=zeros(nEle,4);save_dstrain=zeros(nEle,3);
cell_stress=cell(nload,1);cell_e=cell(nload,1);
cell_strain=cell(nload,1);cell_strainp=cell(nload,1);cell_dstrain=cell(nload,1);
cell_strainp2=cell(nload,1);

for step = 1:nload
    U(1)=0;U(2)=0;U(3)=dis2*step/nload;U(4)=dis2*step/nload;
    U(5)=0;U(6)=-dis1*step/nload;U(7)=0;U(8)=-dis1*step/nload;

    for ie=1:nEle
        node=Element_information(ie,2:4);
        nx0_e = Node_information(node,2:3);
        ue = [U(node);U(node+nNode)];
        ue_pre = [U_pre(node);U_pre(node+nNode)];

        [k,f,Ja,stress,strain,strainpre,dstrain,plstrain,e]=Assemble...
        (nx0_e,ue,ue_pre,e_pre(ie),plstrain_pre(ie,:),stress_pre(ie,:)', ...
        NNcell);
        
        save_stress(ie,:)=stress;save_e(ie,:)=e;
        save_strain(ie,:)=strain;save_strainp(ie,:)=plstrain;save_dstrain(ie,:)=dstrain;  
    end

    e_pre=save_e;
    stress_pre=save_stress;
    plstrain_pre=save_strainp;

    cell_stress(step,1)={save_stress};cell_e(step,1)={save_e};
    cell_strain(step,1)={save_strain};cell_strainp(step,1)={save_strainp};cell_dstrain(step,1)={save_dstrain};
    U_pre = U;
end

%% extract results
plotstress=zeros(nload,4);p=zeros(nload,1);q=zeros(nload,1);plote=zeros(nload,1);
plotstrain=zeros(nload,3);ev=zeros(nload,1);ed=zeros(nload,1);
plotstrainp=zeros(nload,4);epv=zeros(nload,1);epd=zeros(nload,1);epv2=zeros(nload,1);epd2=zeros(nload,1);
for j=1:nload
    plote(j) = cell_e{j}(2,1);

    plotstress(j,:)=-cell_stress{j}(2,:);
    [p_j,q_j]=SINV([cell_stress{j}(2,:),0,0]);
    p(j)=-p_j;q(j)=q_j;

    plotstrain(j,:)=-cell_strain{j}(2,:);
    inpstrain = [cell_strain{j}(2,1:2),0,cell_strain{j}(2,3),0,0];
    [ev_j,ed_j]=SINV(inpstrain);
    ev(j)=-ev_j;ed(j)=ed_j;

    plotstrainp(j,:)=-cell_strainp{j}(2,:);
    [epv_j,epd_j]=SINV([cell_strainp{j}(2,:),0,0]);
    epv(j)=-epv_j;epd(j)=epd_j;
end
plotdataall = [plotstrain(:,2),q,p,epd,plote];
plotFig;