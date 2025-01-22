function [p, q] = solvepq(stress)
STRESS = reshape(stress, [6, 1]);
I = [1 1 1 0 0 0]';
p = (STRESS(1)+STRESS(2)+STRESS(3))/3;
SS = STRESS - p*I;
q = sqrt(3/2*(SS')*SS);
end