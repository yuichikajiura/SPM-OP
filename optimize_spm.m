function [x_est] = optimize_spm(cell_target, cell_known, bound, ...
    profile, method, p_targets, Nr_model, Nr_data, tell_csp0)

    
    %% Parameters
    
    %%% setting parameters
    if cell_known == 1
        run params_setting % load params
    else
        run params_setting_2
    end
    
    p_upper = [p.n_Li_s; p.R_f_n; p.k_n; p.k_p; p.D_s_n; p.D_s_p].* bound;
    p_lower = [p.n_Li_s; p.R_f_n; p.k_n; p.k_p; p.D_s_n; p.D_s_p]./ bound;
    
    if cell_target == 1
        run params_setting
    else
        run params_setting_2
    end
    
    p_true = [p.n_Li_s; p.R_f_n; p.k_n; p.k_p; p.D_s_n; p.D_s_p];
    
    for i = 1:6
        if ~ismember(i, p_targets)
            p_upper(i) = p_true(i);
            p_lower(i) = p_true(i);
        end
        p_upper(i) = p_upper(i)/p_true(i);
        p_lower(i) = p_lower(i)/p_true(i);
    end
    
    %%% importing data
    data = readtable("data/SPM_FDM_nr"+Nr_data+"_simulation_"+profile+"_cell"+cell_target+".csv");
    
    
    %%% Finite Difference Model for particle
    p.Nr = Nr_model;
    p.delta_r_n = 1/p.Nr;
    p.delta_r_p = 1/p.Nr;
    r_vec = (0:p.delta_r_n:1)';
    r_vecx = r_vec(2:end-1);
    
    %%% Finite Difference points along x-coordinate
    p.Nxn = 10;
    p.Nxs = 5;
    p.Nxp = 10;
    p.Nx = p.Nxn + p.Nxs + p.Nxp;
    Nx = p.Nx - 3;
    
    p.delta_x_n = 1 / p.Nxn;
    p.delta_x_s = 1 / p.Nxs;
    p.delta_x_p = 1 / p.Nxp;
    
    % Data setup
    data.Properties.VariableNames = {'t', 'I', 'V_spm', 'css_n', ...
        'cs_ave_n', 'css_p', 'cs_ave_p'};
    I = data.I;
    v_true = data.V_spm;
    t = data.t;
    
    noise = 0.03;  % the magnitude of white noise to standard deviation of data
    std_current = 38;  % standard deviation of current data in Ampere
    std_voltage = 0.16;  % standard deviation of voltage data in Volt
    I = I + normrnd(0, std_current * noise, size(I));
    v_true = v_true + normrnd(0, std_voltage * noise, size(v_true));
    
    % Area for param search
    csp_0_true = data.cs_ave_p(1);
    p_true = [csp_0_true; p_true];

    csp_lower = 0;
    csp_upper = p.c_s_p_max;
    if tell_csp0
        lb = [csp_0_true; p_lower];
        ub = [csp_0_true; p_upper];
    else        
        lb = [csp_lower; p_lower];
        ub = [csp_upper; p_upper];
    end
    
    
    params = @(x) param_update(p, p_true, x);
    c_n0 = @(x) init_cs_paramID(params(x), 2);  %x = [c_p0_1; p_targets]
    c_p0 = @(x) init_cs_paramID(params(x), 1);  %x = [c_p0_1; p_targets]
    c_n0_Nr = @(x) c_n0(x) * ones(p.Nr-1,1);
    c_p0_Nr = @(x) c_p0(x) * ones(p.Nr-1,1);
    sys_n = @(x) spm_plant_obs_mats_paramID_n(params(x));
    sys_p = @(x) spm_plant_obs_mats_paramID_p(params(x));
    css_n = @(x, I, t) calc_css(sys_n(x), I, t, c_n0_Nr(x));
    css_p = @(x, I, t) calc_css(sys_p(x), I, t, c_p0_Nr(x));
    v_pred = @(x, I, t) nonlinearSPMOutputVoltage_paramID(params(x),css_n(x, I, t),css_p(x, I, t),p.c_e,p.c_e,p.c_e,I);
    objFun = @(x) my_rmse(v_pred(x, I, t), v_true);
    
    if method==1
        pop_size = (numel(p_targets)+1)*20;
        max_gen = (numel(p_targets)+1)*100;
        options = optimoptions("ga","Display","iter", ...
            "FunctionTolerance",1e-3, ...
            'PlotFcn','gaplotbestindiv', ...
            'UseParallel', true, ...
            'UseVectorized', false, ...
            'PopulationSize', pop_size, ...
            'OutputFcn', @my_output, ...            
            'MaxGenerations', max_gen);
        x_est = ga(objFun,7,[],[],[],[],lb,ub, [], options);
    else
        options = optimoptions("particleswarm","Display","iter", ...
            'UseParallel', true, ...
            'UseVectorized', false, ...
            'OutputFcn', @my_output, ...
            "FunctionTolerance",1e-3);
        x_est = particleswarm(objFun,7,lb,ub,options);
    end
    disp(x_est);
end

% Function used in optimization algorithm 
function css = calc_css(sys, I, t, c_0)
    [csz, t, c_x] = lsim(sys, I, t, c_0);
    css = csz(:,2);
end

function params = param_update(params, params_true, x)
    params.csp0 = x(1) * params_true(1);
    params.n_Li_s = x(2) * params_true(2);
    params.R_f_n = x(3) * params_true(3);
    params.k_n = x(4) * params_true(4);
    params.k_p = x(5) * params_true(5);
    params.D_s_n = x(6) * params_true(6);
    params.D_s_p = x(7) * params_true(7);
end

function rmse = my_rmse(actual, predicted)
    rmse = sqrt(mean((actual - predicted).^2));
end

function stop = my_output(optimValues, state) 
    stop = false; 
    if state == 'iter' 
        disp([' x = ',num2str(optimValues.bestx)]); 
    end
end