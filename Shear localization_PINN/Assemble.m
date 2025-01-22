function [k,f,Ja,stress,strain,strain_pre,dstrain,plstrain,e]=Assemble...
(nx0_e,ue,ue_pre,e_pre,plstrain_pre,stress_pre, ...
NNcell)

kapa = 0.05;poi=0.3;
% shape function B
A=0.5*abs(det([1 nx0_e(1,1) nx0_e(1,2);1 nx0_e(2,1) nx0_e(2,2); 1 nx0_e(3,1) nx0_e(3,2)]));
B=(1/(2*A))*[nx0_e(2,2)-nx0_e(3,2) nx0_e(3,2)-nx0_e(1,2) nx0_e(1,2)-nx0_e(2,2) 0 0 0;...
             0 0 0 nx0_e(3,1)-nx0_e(2,1) nx0_e(1,1)-nx0_e(3,1) nx0_e(2,1)-nx0_e(1,1);...
nx0_e(3,1)-nx0_e(2,1) nx0_e(1,1)-nx0_e(3,1) nx0_e(2,1)-nx0_e(1,1) nx0_e(2,2)-nx0_e(3,2) nx0_e(3,2)-nx0_e(1,2) nx0_e(1,2)-nx0_e(2,2)];
% element strain
strain = B*ue;
strain_pre = B*ue_pre;
dstrain = strain-strain_pre;

net=NNcell{1};
nor_input=NNcell{2};
nor_output=NNcell{3};

Xinput = [e_pre,stress_pre',strain_pre',dstrain',plstrain_pre]; % symbol transformation
xinput = mapminmax("apply",Xinput',nor_input);
ysim=predict(net,xinput');
Ysim = mapminmax("reverse",ysim',nor_output);

e = Ysim(1);
stress=Ysim(2:5);
plstrain = Ysim(6:9);

p = mean(stress(1:3));
K = -(1+e)*p/kapa;
G = K*3*(1-2*poi)/2/(1+poi);
De = GETDE(K, G);

Ja = De;        
Ja(5:6,:)=[];
Ja(:,5:6)=[];
Ja(3,:)=[];
Ja(:,3)=[];

k=A*B'*Ja*B;
f = k*ue;
end