%% Nonlinear output for voltage in Single Particle Model
%   Created July 21, 2011 by Scott Moura
%   Modified June 25 2024 by Yuichi Kajiura

function [V] = nonlinearSPMOutputVoltage_paramID(p,c_ss_n,c_ss_p,cen_bar,ces_bar,cep_bar, I)

V = zeros(height(I), 1);

for k = 1:height(I)
    % Stochiometric Concentration Ratio
    theta_n = c_ss_n(k) / p.c_s_n_max;
    theta_p = c_ss_p(k) / p.c_s_p_max;


    % Equilibrium Potential
    Unref = refPotentialAnode(p,theta_n);
    Upref = refPotentialCathode(p,theta_p);
    
    % Exchange Current Density
    c_e = zeros(p.Nx,1);     % Fixed electrolyte concentration [mol/m^3]
    
    c_e(1:p.Nxn) = cen_bar;
    c_e((p.Nxn+1):(p.Nxn+p.Nxs-1)) = ces_bar;
    c_e(p.Nxn+p.Nxs+1:end) = cep_bar;
    
    [i_0n,i_0p] = exch_cur_dens(p,c_ss_n,c_ss_p,c_e);
    RTaF=(p.R*p.T_amp)/(p.alph*p.Faraday);
    
    % Voltage
    V(k) = RTaF * asinh(-I(k) / (2*p.a_s_p*p.Area*p.L_p*i_0p(end))) ...
        -RTaF * asinh(I(k) / (2*p.a_s_n*p.Area*p.L_n*i_0n(1))) ...
        + Upref - Unref ...
        - (p.R_f_n/(p.a_s_n*p.L_n*p.Area) + p.R_f_p/(p.a_s_p*p.L_p*p.Area))*I(k);
end