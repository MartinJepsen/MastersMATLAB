clear; close all

%% Set simulation variables
nsr = 0.05;
err = 0.02;
dam_ = 0.80;
sensor = "dis";




%% Compute results
base_dir = sprintf("simulation/SYSID/model_error_%03d_%s", err*100, sensor);
load(fullfile(base_dir, "SetUp.mat"))

r = GeneralParameters.r;
m = GeneralParameters.m;
n_el_K = r*m;
np = 2*n_el_K + 2;
% load('gaindesign\01_strain_cond\gains_5_0.120.mat')
load('gaindesign/01_strain_cond/gains_1_1.120.mat')
initpop = [reshape(real(K), 1, numel(K)), reshape(imag(K), 1, numel(K)), real(s), imag(s)];

ObjectiveFunction = @costfunction;
options = optimoptions('ga', 'Generations', 80,...
                        'PopulationSize', 25,...
                        'FunctionTolerance',1e-7,...
                        'InitialPopulation', initpop,...
                        'PlotFcn', @gaplotbestf);


[res, fval] = ga(ObjectiveFunction, np, [], [], [], [], [], [], [], options);
% Reshape result into complex-values r x m matrix
re = reshape(res(1:n_el_K), r, m);
im = reshape(res(n_el_K+1:2*n_el_K), r, m);
K = complex(re, im);
s = complex(res(end-1), res(end));

save(sprintf("K2_%03d_%03d_%03d_%s", err*100, dam_*100, nsr*100, sensor), "K", "s", "res")

%%
costfunction(res)
%%
function [J] = costfunction(X)
tic
    nsr = evalin('base', 'nsr');
    dam_ = evalin("base", "dam_");
    sensor = evalin("base", "sensor");

    base_dir = evalin('base', 'base_dir');
    load(fullfile(base_dir, "SetUp.mat"))
    load(sprintf("%s/00_000_%03d", base_dir, nsr*100))

    
    n_el_K = evalin('base', 'n_el_K');
    r = evalin('base', 'r');
    m = evalin('base', 'm');    

    re = reshape(X(1:n_el_K), r, m);
    im = reshape(X(n_el_K+1:2*n_el_K), r, m);
    K = complex(re, im);
    s = complex(X(end-1), X(end));

for damel = [1:14]
    damage = [damel, dam_];
%     set_up
    tot_runs = 1;
    
    % load simulation results
    filename = sprintf('%02d_%03d_%03d', damage(1,1), damage(1,2)*100, nsr*100);
%     disp("Loading " + filename)
    load(fullfile(base_dir, filename))
    
    
    Kg = ReferenceModels.Kg;
    Cg = ReferenceModels.Cg;
    Mg = ReferenceModels.Mg;
    B2 = GeneralParameters.B2;
    cdis = GeneralParameters.cdis;
    B_strain = GeneralParameters.B_strain;
    n_dof = GeneralParameters.n_dof;
    idx = GeneralParameters.idx;

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
%     A_CL_ex = SS_exact.A + SS_exact.B * B2 * K * cdis * SS_exact.C;
%     Lambda_CL = eig(A_CL_ex);                       % exact CL poles

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

    %% Results post-processing
    n_el = size(B_strain, 1);
    
    clearvars success_rates
    for i = 1:n_el
        success_rates(i, 1:2) = [sum(min_strain_OL == i), sum(min_strain_CL == i)];
    end
    success_rates = array2table([[1:n_el]', round(success_rates/size(min_strain_OL, 1)*100)], 'VariableNames',...
                    {'el', 'OL', 'CL'});  
    results(damel, :) = success_rates(damel, :);

end
   J = -sum(results.CL);
   toc
end