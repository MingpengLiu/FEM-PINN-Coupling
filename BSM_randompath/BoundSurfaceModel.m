function [stress_n1, DDSDDE,statev_n1] = BoundSurfaceModel(stress_n, strain_n, dstrain_n, Pro, statev_n)

%% read the initial stress at this step
stress = stress_n;
strain = strain_n;
dstrain = dstrain_n;

SSTOL = 1e-3;  %% error tolerance

lama = Pro(1); % compression index
kapa = Pro(2); % swelling index
poi = Pro(3);  % poisson's ratio
M = Pro(4);    % slope of critical state line 
N = Pro(5);    % intersection of NSL
Gama = Pro(6); % intersection of CSL
n = Pro(7);    % n in bounding surface
r = exp((N-Gama)/(lama-kapa)); % r in bounding surface

%% stress integration along all strain increment
% initialize the time 
time = 0;  
dtime = 1;

while time<1
    e1 = statev_n(1);     % vodi ratio
    pc1 = statev_n(2);    % pre-consolidation pressure
    e0 = statev_n(3);     % initial void ratio 
    %% Get the state at the beginning of substep
    % stress constant, p and q have been positive
    [p1, q1] = SINV(stress);  
    q1 = max(1e-5,q1);

    rho1 = sqrt(p1^2+q1^2);  % distance between stress state to original point 

    K1 = -(1+e1)*p1/kapa;          % Bulk modolus 
    G1 = K1*3*(1-2*poi)/2/(1+poi); % Shear modolus
    De1 = GetDe(K1, G1);           % Get the elastic matrix at the begining of substep

    dfdp = -n*((-q1/M/p1)^(n-1))*(-q1/M)/p1/p1+1/p1/log(r);     % deriation of yield surface tp p 
    dfdq = -n*((-q1/M/p1))^(n-1)/M/p1;                          % deriation of yield surface tp q
    ata = q1/p1;                                                % stress ratio
    
    % stress state at bounding surface 
    p1_BS = -pc1*exp(-((-ata/M)^n)*log(r));
    rho1_BS = -sqrt(1+ata^2)*p1_BS;                  

    Kp1 = -(1+e1)/(lama-kapa)/log(r)*(M^2*(rho1_BS/rho1)-ata^2)/2/ata;  % plastic modulus
    Ds1 = (M^2*(rho1/rho1_BS)-ata^2)/2/ata;               % stress-dilation parameter
    
    %% Get the first trial state at the end of substep
    error =1;
    while error > SSTOL
        ddstrain = dstrain*dtime;      % strain increment in substep
        
        %%% first trial computation
        [Dep1, dstrain_sp1, F1, G1] = GetDep(stress, q1, Ds1, dfdp, dfdq, De1, Kp1, ddstrain);  % get elastoplastic matrix

        Devp1 = Ds1*dstrain_sp1;      % plastic volume strain increment

        pc2 = pc1*exp(-(1+e1)/(lama-kapa)*Devp1);    % yield stress at the end of substep

        dstress1 = GetDstress(Dep1, ddstrain);               % stress increment
        stress1 = reshape(stress,1,6)+dstress1;
                                                     
        Dev = (ddstrain(1)+ddstrain(2)+ddstrain(3));  % volumetric strain increment       
        e2 = exp(Dev)*(1+e1)-1;                       % void ratio at the end point

        %%% second trial computation
        [p2, q2] = SINV(stress1);
        q2=max(1e-5,q2);

        K2 = -(1+e2)*p2/kapa;      
        G2 = K2*3*(1-2*poi)/2/(1+poi); 
        De2 = GetDe(K2, G2); 

        dfdp = -n*((-q2/M/p2)^(n-1))*(-q2/M)/p2/p2+1/p2/log(r); 
        dfdq = -n*((-q2/M/p2)^(n-1))/M/p2;                          
        ata = q2/p2;             

        rho2 = sqrt(p2^2+q2^2);
        p2_BS = -pc2*exp(-((-ata/M)^n)*log(r));
        rho2_BS = -sqrt(1+ata^2)*p2_BS;                  

        Kp2 = -(1+e2)/(lama-kapa)/log(r)*(M^2*(rho2_BS/rho2)-ata^2)/2/ata; 
        Ds2 = (M^2*(rho2/rho2_BS)-ata^2)/2/ata;               

        [Dep2, dstrain_sp2, F2, G2] = GetDep(stress1, q2, Ds2, dfdp, dfdq, De2, Kp2, ddstrain);  %% 获取弹塑性刚度矩阵

        Devp2 = Ds2*dstrain_sp2;

        dstress2 = GetDstress(Dep2, ddstrain);               

        %%%  error control
        estress=0.5*(dstress2-dstress1);
        stress2 = stress'+0.5*(dstress1+dstress2);
        error = sqrt(sum(estress.^2)/sum(stress2.^2));
        error = max(1e-8,error);     % error of R
        
        % if not satisfy error, reduce time step
        beta = 0.8*sqrt(SSTOL/error);
        beta = max(0.1,beta);
        dtime=beta*dtime;
    end
    % update time and go on next substep
    time=time+dtime;
    beta= min(2,beta);
 
    dtime=beta*dtime;

    stress=stress+0.5*dstress1'+0.5*dstress2';

    Devp=0.5*(Devp1+Devp2);
    pc=pc1*exp(-(1+e1)/(lama-kapa)*Devp);

    statev_n1(1)=e2;
    statev_n1(2)=pc;
    statev_n1(3)=e0;

    if dtime>(1-time) 
        dtime = 1-time;
    end
end

% final state
[p, q] = SINV(stress);
e = e2;
statev_n1(1)=e2;
statev_n1(2)=pc;
statev_n1(3)=e0;

K = -(1+e)*p/kapa;
G = K*3*(1-2*poi)/2/(1+poi);

De = GetDe(K, G);
dfdp = -n*((-q/M/p)^(n-1))*(-q/M)/p/p+1/p/log(r);  
dfdq = -n*((-q/M/p)^(n-1))/M/p;                          
ata = q/p;              
rho = sqrt(p^2+q^2);
p_BS = -pc*exp(-((-ata/M)^n)*log(r));
rho_BS = -sqrt(1+ata^2)*p_BS;                
Kp = -(1+e)/(lama-kapa)/log(r)*(M^2*(rho_BS/rho)-ata^2)/2/ata;  
Ds = (M^2*(rho/rho_BS)-ata^2)/2/ata;              
[DDSDDE, dstrain_sp, F, G] = GetDep(stress, q, Ds, dfdp, dfdq, De, Kp, dstrain);  

stress_n1 = stress;
end