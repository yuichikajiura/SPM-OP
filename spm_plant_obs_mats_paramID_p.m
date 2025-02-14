%% Matrices for Li Diffusion in Single Particle Model, c_s_p(r,t) & c_s_n(r,t)


function sys_p = spm_plant_obs_mats_paramID_p(p)

% Electrochemical Model Parameters
alpha_p = p.D_s_p / (p.R_s_p * p.delta_r_p)^2;

% Block matrices
M1_p = zeros(p.Nr-1);

for idx = 1:(p.Nr-1)
    
    % Lower diagonal
    if(idx ~= 1)
        M1_p(idx,idx-1) = (idx-1)/idx * alpha_p;  
    end
    
    % Main diagonal
    M1_p(idx,idx) = -2*alpha_p;
    
    % Upper diagonal
    if(idx ~= (p.Nr-1))
        M1_p(idx,idx+1) = (idx+1)/idx * alpha_p;
    end
end

M2_p = zeros(p.Nr-1,2);
M2_p(end,end) = p.Nr/(p.Nr-1) * alpha_p;
N1 = zeros(2,p.Nr-1);

% 2nd Order BCs
N1(1,1) = 4;
N1(1,2) = -1;
N1(2,end) = -4;
N1(2,end-1) = 1;
N2 = diag([-3,3]);
N3_p = [0; (2*p.delta_r_p * p.R_s_p)/(p.D_s_p * p.Faraday * p.a_s_p * p.Area * p.L_p)];

% A,B matrices for each electrode
A_p = M1_p - M2_p*(N2\N1);
B_p = M2_p*(N2\N3_p);

% C,D matrices for each electrode
C_p = -(N2\N1);
D_p = (N2\N3_p);
sys_p = ss(A_p,B_p,C_p,D_p);