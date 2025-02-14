%% Matrices for Li Diffusion in Single Particle Model, c_s_p(r,t) & c_s_n(r,t)


function sys_n = spm_plant_obs_mats_paramID_n(p)

% Electrochemical Model Parameters
alpha_n = p.D_s_n / (p.R_s_n * p.delta_r_n)^2;

% Block matrices
M1_n = zeros(p.Nr-1);

for idx = 1:(p.Nr-1)
    
    % Lower diagonal
    if(idx ~= 1)
        M1_n(idx,idx-1) = (idx-1)/idx * alpha_n;
    end
    
    % Main diagonal
    M1_n(idx,idx) = -2*alpha_n;
    
    % Upper diagonal
    if(idx ~= (p.Nr-1))
        M1_n(idx,idx+1) = (idx+1)/idx * alpha_n;
    end
end

M2_n = zeros(p.Nr-1,2);
M2_n(end,end) = p.Nr/(p.Nr-1) * alpha_n;
N1 = zeros(2,p.Nr-1);

% 2nd Order BCs
N1(1,1) = 4;
N1(1,2) = -1;
N1(2,end) = -4;
N1(2,end-1) = 1;
N2 = diag([-3,3]);
N3_n = [0; -(2*p.delta_r_n * p.R_s_n)/(p.D_s_n * p.Faraday * p.a_s_n * p.Area * p.L_n)];

% A,B matrices for each electrode
A_n = M1_n - M2_n*(N2\N1);
B_n = M2_n*(N2\N3_n);

% C,D matrices for each electrode
C_n = -(N2\N1);
D_n = (N2\N3_n);

sys_n = ss(A_n,B_n,C_n,D_n);