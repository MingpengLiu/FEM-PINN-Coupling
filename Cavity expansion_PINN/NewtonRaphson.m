function [d,R_r]=NewtonRaphson(KK,R,node_ndof)
nd = size(KK,1);
fdof=(1:nd)';

d=zeros(size(fdof));

pdof=node_ndof;
fdof(pdof)=[];

R_r=R(fdof);
s = -lsqminnorm(KK(fdof,fdof),R_r);

d(pdof)=0;
d(fdof)=s;

end