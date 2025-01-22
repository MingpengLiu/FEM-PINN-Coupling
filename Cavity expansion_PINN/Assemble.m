function [k,f,Ja,stress,strain,strain_pre,dstrain,plstrain,e]=Assemble...
            (nx0_e,ue,ue_pre,e_pre,plstrain_pre,stress_pre, ...
            net,nor_input,nor_output)
% elastic parameters
kapa = 0.05;poi=0.3;
% shape function B
A=0.5*abs(det([1 nx0_e(1,1) nx0_e(1,2);1 nx0_e(2,1) nx0_e(2,2); 1 nx0_e(3,1) nx0_e(3,2)]));
B=(1/(2*A))*[nx0_e(2,2)-nx0_e(3,2) nx0_e(3,2)-nx0_e(1,2) nx0_e(1,2)-nx0_e(2,2) 0 0 0;...
             0 0 0 nx0_e(3,1)-nx0_e(2,1) nx0_e(1,1)-nx0_e(3,1) nx0_e(2,1)-nx0_e(1,1);...
nx0_e(3,1)-nx0_e(2,1) nx0_e(1,1)-nx0_e(3,1) nx0_e(2,1)-nx0_e(1,1) nx0_e(2,2)-nx0_e(3,2) nx0_e(3,2)-nx0_e(1,2) nx0_e(1,2)-nx0_e(2,2)];
% element strain
strain = B*ue;         % strain at current step
strain_pre = B*ue_pre; % strain at last step
dstrain = strain-strain_pre; % strain increment

Xinput = [e_pre,stress_pre',strain_pre',dstrain',plstrain_pre]; % input for NN
xinput = mapminmax("apply",Xinput',nor_input);                  % normalized input
ysim=predict(net,xinput');                                      % prediction via NN
Ysim = mapminmax("reverse",ysim',nor_output);                   % revserse normalization for NN

e = Ysim(1);
stress=Ysim(2:5);
plstrain = Ysim(6:9);

% elastic stiffness matrix
[p, q] = SINV([stress',0,0]);
K = -(1+e)*p/kapa;
G = K*3*(1-2*poi)/2/(1+poi);
De = GETDE(K, G);

% delimension reduction
Ja = De;
Ja(5:6,:)=[];
Ja(:,5:6)=[];
Ja(3,:)=[];
Ja(:,3)=[];
k= A* B'*Ja*B;
f = k*ue;
end