function validation = DataValidation(X)
poi = 0.3;
kapa = 0.05;
for i=1:length(X)
    Data = X{i};
    e= Data(:,1);
    stress = Data(:,2:5);
    strain = Data(:,6:8);
    straininc = Data(:,9:11);
    plastrain = Data(:,12:15);
    
    for j=1:length(Data)-1
        stress_n=[stress(j,1:4),0,0];
        [p_n, q_n] = solvepq(stress_n);
        stress_n1=[stress(j+1,1:4),0,0];
        d_stress = stress_n1-stress_n;

        e_n = e(j);
        K_n = (1+e_n)*p_n/kapa;        %% Bulk modolus 
        G_n = K_n*3*(1-2*poi)/2/(1+poi); %% Shear modolus
        De_n = GetDe(K_n, G_n);

        dstrain_e = De_n\d_stress'; 
        dstrain_p = [straininc(j,1:2),0,straininc(j,3),0,0]-dstrain_e';
        plastrain_n1 = [plastrain(j,:),0,0]+dstrain_p;
        error_plastrain(j,:) = [plastrain(j+1,:),0,0]-plastrain_n1;
    end    
    errori(i)=norm(error_plastrain);
end
if norm(errori)<10-4
    validation=1;
else
    validation=0;
end
end