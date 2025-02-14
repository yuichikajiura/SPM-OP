%% Compute Initial Solid Concentrations from Voltage and Params
%   Created Mar 12, 2014

function cs0_electrode = init_cs_paramID(p, electrode_id)

csp0 = p.csp0;
csn0 = (p.n_Li_s - p.epsilon_s_p * p.L_p * p.Area * csp0) / (p.epsilon_s_n * p.L_n * p.Area);
cs0 = [csp0, csn0];
cs0_electrode = cs0(electrode_id);