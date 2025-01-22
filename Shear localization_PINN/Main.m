clear;clc;
%% Material parameter
e0=0.85;p0=100;

NNcell=load("netcell.mat");
NNcell=NNcell.netcell;
net=NNcell{1};
nor_input=NNcell{2};
nor_output=NNcell{3};

tolc = 1e-6;
maxit = 35;
%% meshing, node and element

[Node_information,Element_information,Gauss]=Meshing2;
nNode = length(Node_information);
nEle = size(Element_information,1);

xmin = min(Node_information(:,2));xmax = max(Node_information(:,2));
ymin = min(Node_information(:,3));ymax = max(Node_information(:,3));
%% loading
% boundary and loading
node_left = find(abs(Node_information(:,2)-xmin)<tolc);
node_right = find(abs(Node_information(:,2)-xmax)<tolc);
node_top_dis = find(abs(Node_information(:,3)-ymax)<tolc);
node_bottom = find(abs(Node_information(:,3)-ymin)<tolc);

node_ndof= [node_bottom+nNode;node_top_dis+nNode];

disp1 = 0.006;
nload = 10;
F_ext = zeros(2*nNode,1);
F_ext_b = zeros(2*nNode,1);

ini_stress = repmat([-p0 -p0 -p0 0 0 0], nEle,1);ini_e = repmat(e0, nEle,1);
stress_pre=ini_stress(:,1:4);
e_pre=ini_e;
plstrain_pre = zeros(nEle,4);

U_pre =zeros(2*nNode,1);U = zeros(2*nNode,1);cell_U=cell(nload,1);
save_e=zeros(nEle,1);save_stress=zeros(nEle,4);
save_strain=zeros(nEle,3);save_strainp=zeros(nEle,4);save_dstrain=zeros(nEle,3);

cell_stress=cell(nload,1);cell_e=cell(nload,1);
cell_strain=cell(nload,1);cell_strainp=cell(nload,1);cell_dstrain=cell(nload,1);

for step = 1:nload
    U(node_bottom+nNode)=0;
    U(node_top_dis+nNode)=-disp1*step/nload;
    
    % Newton-Rasphon iteration
    for itn = 1:maxit
        KK=zeros(2*nNode); 
        for ie= 1:nEle
            node=Element_information(ie,2:4);
            nx0_e = Node_information(node,2:3);
            ue = [U(node);U(node+nNode)];
            ue_pre = [U_pre(node);U_pre(node+nNode)];

            [k,f,Ja,stress,strain,strainpre,dstrain,plstrain,e]=Assemble...
            (nx0_e,ue,ue_pre,e_pre(ie),plstrain_pre(ie,:),stress_pre(ie,:)', ...
            NNcell);

            list = [node,node+nNode];
            K=zeros(2*nNode);
            for i =1:6
                for j =1:6
                    K(list(i),list(j))=k(i,j);
                end
            end
            KK = KK+K;

            save_stress(ie,:)=stress;save_e(ie,:)=e;
            save_strain(ie,:)=strain;save_strainp(ie,:)=plstrain;save_dstrain(ie,:)=dstrain;  
        end
        F_int = KK*U;
        R = F_int-F_ext-F_ext_b;
        [du,R_r]=NewtonRaphson(KK,R,node_ndof);
        
        % convergence check
        if norm(R_r)<tolc
            break;
        end
        U = U+du;
    end
    e_pre=save_e;stress_pre=save_stress;plstrain_pre=save_strainp;
    cell_U(step,1)={U};
    cell_stress(step,1)={save_stress};cell_e(step,1)={save_e};
    cell_strain(step,1)={save_strain};cell_strainp(step,1)={save_strainp};cell_dstrain(step,1)={save_dstrain};
    U_pre = U;
end

%% extract results
plotstrain=zeros(nload,3);ev=zeros(nload,1);ed=zeros(nload,1);
loadstep = nload;
for j =1:nEle
    plotstrain=cell_strain{loadstep};
    [ev_j,ed_j]=SINV([plotstrain(j,1:2),0,plotstrain(j,3),0,0]);
    ev(j)=ev_j;ed(j)=ed_j;
end
plotfig;
