function [Dep, dstrain_sp, F, G] = GetDep(stress, q, Ds, dfdp, dfdq, De, Kp, ddstrain)
A = [1/3 1/3 1/3 0 0 0];

B(1) = 1/2/q*(2*stress(1)-stress(2)-stress(3));
B(2) = 1/2/q*(2*stress(2)-stress(1)-stress(3));
B(3) = 1/2/q*(2*stress(3)-stress(1)-stress(2));
B(4:6) = 3/2/q*2*stress(4:6);

C = Ds*A+B;
D = dfdp*A+dfdq*B;  % Differential of F to stress

E = D*De; % dF_dstress * De  (D is the elastic matrix) 

F =sum(C.*E);  %  dF/dstress * De * C
G=C*De;        % D * C

H = G'*E;

Dep=De-1/(Kp+F)*H;

dstrain_sp=E*ddstrain/(Kp+F);

unload=D*ddstrain;

if unload < 0
    Dep=De;
    dstrain_sp=0;
end

end





