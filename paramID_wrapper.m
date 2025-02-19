clc;
clear;
close all;
% parpool


%% Configuration
cell_target = 2 ;% target cell to identify params. 1:fresh, 2:degraded
cell_known = 1; % cell of known parameter. 1:fresh, 2:degraded
bound = 1.5; %ratio to the known params of other cell for the boundary of param search range
profiles = ["UDDSx2", "FUDS", "Charge", "Charge2"]; %option:UDDSx2, FUDS, Charge, Charge2, US06_Extended
tell_csp0 = false; % true if assume we know initial concentration 
p_targets = {[1, 3, 4, 5, 6]}; %%1:nLi, 2:R_f_n, 3:k_n, 4:k_p, 5:D_s_n, 6:D_s_p (0 to estimate csp0)
Nr_model = 20;
Nr_data = 100;
methods = [1, 2]; %1:GA, 2:PSO

%% Variables for recording purpose 
method_names = ["GA", "PSO"];
results = cell2table(cell(0, 11), 'VariableNames', {'experiment_ID', 'p_target', 'profile', 'method', 'c_n_0', 'nLi', 'R_f_n', 'k_n', 'k_p', 'D_s_n', 'D_s_p'});

%% Run experiments 
l = 1;
for i = 1:length(p_targets)
    results = cell2table(cell(0, 11), 'VariableNames', {'experiment_ID', 'p_target', 'profile', 'method', 'c_n_0', 'nLi', 'R_f_n', 'k_n', 'k_p', 'D_s_n', 'D_s_p'});
    for j = 1:length(profiles)
        for k = 1:length(methods)
            X = sprintf('Run experiment for p_target %d on  %s profile with %s', cell2mat(p_targets(i)), profiles(j), method_names(methods(k)));
            disp(X);
            x_est = optimize_spm(cell_target, cell_known, bound, ...
                profiles(j), methods(k), cell2mat(p_targets(i)), Nr_model, Nr_data, tell_csp0);
            result = {l, cell2mat(p_targets(i)), profiles(j), method_names(k), ...
                x_est(1), x_est(2), x_est(3), x_est(4), x_est(5), x_est(6), x_est(7)};
            results = [results; result];
            l = l+1;
        end
    end
    writetable(results, "result_"+strjoin(string(cell2mat(p_targets(i))), "-")+".txt")
end
