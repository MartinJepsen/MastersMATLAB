clear; clc; close all

%% Load relevant variables
nsr = 0.05;
load(sprintf("simulation/SYSID/model_error/00_000_%03d", nsr*100))

s_facs = sort([1:0.05:1.4, 1.001, 1.002, 1.003], 'ascend');
for s_num = 1:numel(s_facs)
    for damel = [1:14]
        dam = [damel, 0.80];
        set_up
        %% load simulation results
        filename = sprintf('%02d_%03d_%03d',dam(1,1), dam(1,2)*100, nsr*100);
        disp("Loading " + filename)
        load("simulation/SYSID/model_error/"+filename)
        
        %% load gains
        
        load(sprintf("gaindesign/01_strain_cond/pole_3/gains_%0.3f.mat", s_facs(s_num)))
    %     load(sprintf("gaindesign/02_sens/constrained/gains_%02d", damel))
    %     load(sprintf("gaindesign/03_strain_norm/gains_%02d", damel))
        
        H_ref = (Mg*s^2 + Cg*s + Kg)^-1;                % reference OL transfer matrix
        H_CL_ref = (Mg*s^2 + Cg*s + Kg + B2*K*cdis)^-1; % reference CL transfer matrix
        A_CL_ex = SS_exact.A + SS_exact.B * B2 * K * cdis * SS_exact.C;
        Lambda_CL_ex = eig(A_CL_ex);                    % exact CL poles
    
        %% Obtain characteristic strains for every run
        tot_runs = 0;
    
        for run = 1:numel(SS_est)
            % OL
            SS = SS_est{run};
            H = SS.transfer_matrix(s); 
            SS_d = SS_est_d{run};
            H_d = SS_d.transfer_matrix(s);
            
            l_CL_est = eig(SS.A + SS.B*K*SS.C);
            [~, sorting] = sort(abs(l_CL_est), 'ascend');
            Lambda_CL_est(:, run) = l_CL_est(sorting);          % estimated CL poles
    
            DeltaH = H_d - H;                               % damage-induced transfer matrix shift (estimated)
            [~, ~, V] = svd(DeltaH);                        % DDLVs
            d_OL = zeros(n_dof, 1);
            d_OL(idx) = H_ref * B2 * V(:, end);             % full OL displacement vector
            eps_OL = B_strain * d_OL;                       % full OL strain vector
            
            % CL
            H_CL = (eye(size(H)) + H * K)^-1 * H;           % estimated CL transfer matrix, reference
            H_CL_d = (eye(size(H_d)) + H_d * K)^-1 * H_d;   % estimated CL transfer matrix, damaged
            DeltaH_CL = H_CL_d - H_CL;                      % CL damage-induced transfer matrix change
            [~, ~, V] = svd(DeltaH_CL);                     % CLDDLVs
            d_CL = zeros(n_dof, 1);
            d_CL(idx) = H_CL_ref * B2 * V(:, end);          % CL displacement vector
            eps_CL = B_strain * d_CL;                       % cCL strain vector
            
            % Compute strains
            strains(:, tot_runs+1, 1) = abs(eps_OL);        % array of characteristic OL strain vectors
            strains(:, tot_runs+1, 2) = abs(eps_CL);        % array of characteristic CL strain vectors
            tot_runs = tot_runs + 1;
            
            eps_norm_OL = abs(eps_OL) / max(abs(eps_OL));
            eps_norm_CL = abs(eps_CL) / max(abs(eps_CL));
            
            min_strain_OL(tot_runs, 1) = find(eps_norm_OL == min(eps_norm_OL));   % index of smallest OL strain
            min_strain_CL(tot_runs, 1) = find(eps_norm_CL == min(eps_norm_CL));   % index of smallest CL strain
        end
    
        %% Results post-processing
        n_el = 14;
        clearvars success_rates
        for i = 1:n_el
            success_rates(i, 1:2) = [sum(min_strain_OL == i), sum(min_strain_CL == i)]/tot_runs*100;
        end
        success_rates = array2table([[1:n_el]', round(success_rates)], 'VariableNames', {'element_number', 'detection_rate_OL', 'detection_rate_CL'});  
        results(damel, :) = success_rates(damel, :);
    end
    OL(:, s_num) = table2array(results(:,2));
    CL(:, s_num) = table2array(results(:,3));
end

save("s_values/results_pole3.mat", "OL", "CL", "s_facs")