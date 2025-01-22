function [X, Y] = SetXY(data, datanum)
M = length(datanum);
X = cell(M,1);
Y = cell(M,1);
for i = 1:length(datanum)
    e = data{datanum(i),1};
    e_x = e(1:end-1,:);
    e_y = e(2:end,:);

    stress = data{i,2};
    stress_x = [stress(1:end-1,1:3),stress(1:end-1,4)];
    stress_y = [stress(2:end,1:3), stress(2:end,4)];

    strain = data{i,3};
    strain_x = [strain(1:150,1:2),strain(1:150,4)];

    strain_inc = data{i,4};
    strain_inc_x = [strain_inc(:,1:2),strain_inc(:,4)];

%     pstrain_inc = data{i,6};
%     pstrain_inc_y = pstrain_inc(:,1:3);

    pstrain = data{i,5};
    pstrain_x = pstrain(1:end-1,1:4);
    pstrain_y = pstrain(2:end,1:4);

    X(i) = {[e_x,stress_x,strain_x,strain_inc_x,pstrain_x]};
    Y(i) = {[e_y,stress_y,pstrain_y]};
end
end