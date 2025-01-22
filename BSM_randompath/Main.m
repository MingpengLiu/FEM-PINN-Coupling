clear;
%% model parameter
lama = 0.15;  
kapa = 0.05;  
poi = 0.3;   
M = 1.2;      
N = 2.69;     
Gama = 2.62;  
n = 1.6;      
Pro = [lama, kapa, poi, M, N, Gama, n];

% initial state parameters
pc= 100;   % consolidation pressure
e0=0.85;   % initial void ratio
p0=100;    % initial stress state - confining pressure

ns = 150;  % total step in each curve
nc = 4;    % control points of strain path 
np = 20;   % number of generated strain paths
strainrange = [-0.15 0.15];  % strain range -0.15 to 0.15
%% generate strain/stress data
ini_stress = [-p0 -p0 -p0 0 0 0]';    % tensile is positive，compression is negative
STATEV = [e0 pc e0];                  % state variables
stress = ini_stress;

strainpath = strain_path(ns,nc,strainrange);  % random strain path 
strainincre = strainpath(2:end,:)-strainpath(1:end-1,:);  % strain increment

% saving data
save_stress = zeros(ns, 6);
save_void   = zeros(ns, 1);
save_strain_p = zeros(ns, 6);
for i = 1:ns
    strain = -strainpath(i,:)';    % change to the sign for constitutive model
    dstrain = -strainincre(i,:)';  % strain increment at current step

    % using bounding surface model to get stress and state variables at
    % next step
    [stress_n1, DDSDDE_n1,STATEV_n1] = BoundSurfaceModel(stress, strain, dstrain, Pro, STATEV);
        
    d_stress = stress_n1-stress;   % stress increment
    
    % assemble tangient elastic matrix
    e_n = STATEV(1);
    [p_n, q_n] = SINV(stress);       % get p and q
    K_n = -(1+e_n)*p_n/kapa;         % Bulk modolus 
    G_n = K_n*3*(1-2*poi)/2/(1+poi); % Shear modolus
    De_n = GetDe(K_n, G_n);          % elastic stiffiness matrix

    dstrain_e = De_n\d_stress;       % elastic strain increment
    dstrain_p = dstrain-dstrain_e;   % plastic strain increment

    e_n1 = STATEV_n1(1);             % void ratio
    stress = stress_n1;              % update stress
    STATEV = STATEV_n1;              % update state variables
    
    % save variables
    save_stress(i,1:6)=-stress_n1;
    save_strain_p(i+1,1:6) = save_strain_p(i,1:6)-dstrain_p';
    save_void(i,1) = e_n1;
end

save_stress = [-ini_stress';save_stress];  
save_void   = [e0;save_void];
save_dstrain_p = save_strain_p(2:end,:)-save_strain_p(1:end-1,:);
            
plotfig;