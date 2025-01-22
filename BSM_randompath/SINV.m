function [SINV1, SINV2] = SINV(stress)
stress = reshape(stress, [6, 1]);
I = [1 1 1 0 0 0]';
SINV1 = (stress(1)+stress(2)+stress(3))/3;
SS = stress - SINV1*I;
SINV2 = sqrt(3/2*(SS')*SS);
end