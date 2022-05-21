clear; close all

%% Set simulation variables
nsr = 0.05;
err = 0.00;
dam_ = 0.80;
sensor = "dis";

show_plots = false;

%% Compute results
base_dir = sprintf("simulation/SYSID/model_error_%03d_%s", err*100, sensor);
load(sprintf("%s/00_000_%03d", base_dir, nsr*100))

% load gains
load('gaindesign/01_strain_cond/gains_1_1.120.mat')


ObjectiveFunction = @cost_function;
options = optimoptions('ga', 'Generations', 10,...
                        'PopulationSize', 10,...
                        'CrossoverFraction', 0.6,...
                        'FunctionTolerance',1e-7,...
                        'PlotFcn', @gaplotbestf);
    
[res, fval] = ga(ObjectiveFunction, 2, [], [], [], [], [], [], [], options);




function [J] = cost_function(X) 
    dam_ = evalin('base', 'dam_');
    nsr = evalin('base', 'nsr');
    err = evalin('base', 'err');
    sensor = evalin("base", "sensor");
    base_dir = evalin("base", "base_dir");
    K = evalin("base", "K");
    SS_est = evalin("base", "SS_est");
    n_el = 14;
    results = zeros(n_el, 2);
    for damel = [1:n_el]
        dam = [damel, dam_];
        set_up
        tot_runs = 1;
        % load simulation results
        filename = sprintf('%02d_%03d_%03d', dam(1,1), dam(1,2)*100, nsr*100);
        load(fullfile(base_dir, filename))

        s = complex(X(1), X(2));

        % account for output type
        if sensor == "dis"
            s_fac = 1;
        elseif sensor == "vel"
            s_fac = 1/s;
        elseif sensor == "acc"
            s_fac = 1 / (s^2);
        end
    
        % model transfer matrices
        H_ref = (Mg*s^2 + Cg*s + Kg)^-1;                % reference OL transfer matrix
        H_CL_ref = (Mg*s^2 + Cg*s + Kg + B2*K*cdis)^-1; % reference CL transfer matrix

        % Compute all estimated transfer matrices
        n_sim = numel(SS_est);
        n_sim_d = numel(SS_est_d);
        H_arr = cell(n_sim, 1);
        H_CL_arr = cell(n_sim, 1);
        H_d_arr = cell(n_sim_d, 1);
        H_CL_d_arr = cell(n_sim_d, 1);
        if tot_runs == 1
            for run_u = 1:n_sim
                H = s_fac * SS_est{run_u}.transfer_matrix(s);
                H_arr{run_u, 1} = H;
                H_CL_arr{run_u, 1} = (eye(size(H)) + H * K)^-1 * H;
            end
        end
        for run_d = 1:n_sim_d
            H_d = s_fac * SS_est_d{run_d}.transfer_matrix(s);
            H_d_arr{run_d, 1} = H_d;
            H_CL_d_arr{run_d, 1} = (eye(size(H_d)) + H_d * K)^-1 * H_d;   % estimated CL transfer matrix, damaged
        end
    
        strains = zeros(size(B_strain, 1), n_sim * n_sim_d, 2);
        min_strain_OL = zeros(n_sim * n_sim_d, 1);
        min_strain_CL = zeros(n_sim * n_sim_d, 1);
        for run_u = 1:n_sim
            H = H_arr{run_u};
            H_CL = H_CL_arr{run_u};
            
            for run_d = 1:n_sim_d
                H_d = H_d_arr{run_d};
                H_CL_d = H_CL_d_arr{run_d};
    
                DeltaH = H_d - H;                               % damage-induced transfer matrix shift (estimated)
                [~, ~, V] = svd(DeltaH);                        % DDLVs
                d_OL = zeros(n_dof, 1);
                d_OL(idx) = H_ref * B2 * V(:, end);             % full OL displacement vector
                eps_OL = B_strain * d_OL;                       % full OL strain vector
            
                DeltaH_CL = H_CL_d - H_CL;                      % CL damage-induced transfer matrix change
                [~, ~, V] = svd(DeltaH_CL);                     % CLDDLVs
                d_CL = zeros(n_dof, 1);
                d_CL(idx) = H_CL_ref * B2 * V(:, end);          % CL displacement vector
                eps_CL = B_strain * d_CL;                       % cCL strain vector
            
                % Compute strains
                strains(:, tot_runs, 1) = abs(eps_OL);        % array of characteristic OL strain vectors
                strains(:, tot_runs, 2) = abs(eps_CL);        % array of characteristic CL strain vectors
    
                eps_norm_OL = abs(eps_OL) / max(abs(eps_OL));
                eps_norm_CL = abs(eps_CL) / max(abs(eps_CL));
                
                min_strain_OL(tot_runs, 1) = find(eps_norm_OL == min(eps_norm_OL));   % index of smallest OL strain
                min_strain_CL(tot_runs, 1) = find(eps_norm_CL == min(eps_norm_CL));   % index of smallest CL strain
                tot_runs = tot_runs + 1;
            end
        end
    
        %% Results post-processin
        
        success_rates = zeros(n_el, 2);
        for i = 1:n_el
            success_rates(i, 1:2) = [sum(min_strain_OL == i), sum(min_strain_CL == i)];
        end
        success_rates = success_rates/size(min_strain_OL, 1)*100;
        results(damel, :) = success_rates(damel, :);  
    end
    J = - sum(results(:, 1));
end