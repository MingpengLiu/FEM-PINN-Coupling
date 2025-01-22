function [loss, gradients,Y_pre,y_pre,y_pre_r] = modelLoss(net, x, y, X_max,X_min,Y_max,Y_min)

% loss between real value and predicted value
y_pre = predict(net, x);  

loss_e = mse(y(1,:), y_pre(1,:));
loss_stress = mse(y(2:5,:), y_pre(2:5,:));
loss_data = loss_e+loss_stress;

% extract data in real value
y_pre_r = extractdata(y_pre);
y_r = extractdata(y);
x_r = extractdata(x);

% reverse normalization to normal values
Y_pre = y_pre_r.*(Y_max-Y_min)+Y_min;
Y = y_r.*(Y_max-Y_min)+Y_min;
X = x_r.*(X_max-X_min)+X_min;

% elastic parameters
kapa = 0.05;
poi = 0.3;

% real variables at time n
e_n = X(1,:);
stress_n = X(2:5,:);
straininc_n = X(9:11,:);
plstrain_n = X(12:15,:);

% real variables at time n+1
e_n1 = Y(1,:);
stress_n1 = Y(2:5,:);
plstrain_n1 = y(6:9,:); %  normalized dlarray format

% predicted variables at time n+1
stress_n1_pre = Y_pre(2:5,:);
plstrain_n1_pre = y_pre(6:9,:);  %  normalized dlarray format

% stress increment
dstress_pre = stress_n1_pre-stress_n;

% elastic stiffness matrix
p_n = mean(stress_n(1:3,:),1);
K = (1+e_n).*p_n/kapa;
G = K*3*(1-2*poi)/2/(1+poi);

% inverse the plastic strain
plstrain_n1_cal = zeros(4,length(y));
for i = 1:length(y)
    De = zeros(4,4);
    De(1,1) = K(i)+4/3*G(i);
    De(2,2) = K(i)+4/3*G(i);
    De(3,3) = K(i)+4/3*G(i);
    De(4,4) = G(i);
    De(1,2) = K(i)-2/3*G(i);
    De(1,3) = K(i)-2/3*G(i);
    De(2,1) = K(i)-2/3*G(i);
    De(2,3) = K(i)-2/3*G(i);
    De(3,1) = K(i)-2/3*G(i);
    De(3,2) = K(i)-2/3*G(i);

    % predicted variables
    dstrain_e_pre = De\dstress_pre(:,i);     % elastic strain increment
    dstrain_p_pre = [straininc_n(1:2,i);0;straininc_n(3,i)]-dstrain_e_pre;  % plastic strain increment
    plstrain_n1_cal(1:4,i) = plstrain_n(:,i)+dstrain_p_pre;   % plastic strain at time n+1
end
plstrain_n1_cal_norm = (plstrain_n1_cal-Y_min(6:9))./(Y_max(6:9)-Y_min(6:9)); % nomalized plastic strain
plstrain_n1_cal_dl = dlarray(plstrain_n1_cal_norm,'CB');  % transfer to dlarray format

loss_plastrain = mse(plstrain_n1_cal_dl,plstrain_n1_pre); % loss in plastic strain

% total_loss = w_e*loss_e + w_stress*loss_stress + w_plastrain*loss_plastrain;            % total loss
total_loss = loss_data + loss_plastrain;            % total loss

gradients = dlgradient(total_loss, net.Learnables);             % update gradients

loss = [total_loss,loss_e,loss_stress,loss_plastrain];
end