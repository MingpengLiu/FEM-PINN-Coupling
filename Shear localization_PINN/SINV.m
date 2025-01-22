function [SINV1, SINV2] = SINV(STRESS)
STRESS = reshape(STRESS, [6, 1]);
I = [1 1 1 0 0 0]';
SINV1 = (STRESS(1)+STRESS(2)+STRESS(3))/3;
SS = STRESS - SINV1*I;
SINV2 = sqrt(3/2*(SS')*SS);
end