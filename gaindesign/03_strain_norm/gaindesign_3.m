clc; clear; close all;

load("gaindesign/gain_pars")          % load system matrices
idx = setdiff(1:n_dof, bc);


% 1) set up reference and undamaged SS models
% 2) get transfer matrices

H = SS_exact.H;
H_d = SS_exact.H;


% H_ = zeros(n_dof, free_dof);
% H_(idx, :) = H;

%% Genetic algorithm
run = 0;
good = 0;
for run = 0;

np = r*m; % No. of parameters

ObjectiveFunction = @main_gain_design;
options = optimoptions('ga', 'Generations', 5000,...
                        'PopulationSize', 100,...
                        'FunctionTolerance',1e-20,...
                        'PlotFcn', @gaplotbestf);

% [res, fval] = ga(ObjectiveFunction, nvars, [], [], [], [], lb, ub, [], options);
[res, fval] = ga(ObjectiveFunction, np, [], [], [], [], [], [], [], options);

K = reshape(res, r, m);
    
    results{run+1,1} = K;
    results{run+1,2} = fval;

end


%% Sort results
[sens_norms, ind] = sort([results{:,2}].^-1, 'descend'); % sort the inverse fitness function values (norm of sensitivity) in descending order

for i = 1:numel(sens_norms)
    gains{i,2} = sens_norms(i);
    gains(i,1) = results(ind(i),1);
end

beep

%%
save(sprintf("optimised_DDLV/03_strain_norm/unit_perturbations/gains_%d",damel),"gains")

function [J] = main_gain_design(X)

    free_dof = evalin('base', 'free_dof');
    n_dof = evalin('base', 'n_dof');
    m = evalin('base', 'm');
    r = evalin('base', 'r');
    B2 = evalin('base', 'B2');
    B = evalin('base', 'B_strain');
    idx = evalin('base', 'idx');
    H_d = evalin('base', 'H_d');
    H = evalin('base','H');
    cdis = evalin('base','cdis');
    FE = evalin('base', 'FE')
    
    K = reshape(X, r, m);              % gain matrix
    
    
    % OL transfer matrices
    deltaH = H_d - H;
    
    % CL transfer matrices
    H_CL = (eye(size(H)) + H * b2*K*cdis) \ H;
    H_CL_d = (eye(size(H_d)) + H_d * b2*K*cdis) \ H_d;
    deltaH_CL = H_CL_d - H_CL;

    % Strain fields
    [~, ~, V] = svd(deltaH);
    d_OL = [0; G*V(:, end)];
    
    [~, ~, V] = svd(dG_CL);
    d_CL = [0; G_CL * V(:, end)]; 
    eps = zeros(n_dof, 1, 2);
    for el = 1:n_dof
        eps(el, 1, 1) = [1 -1] * d_OL(el:el+1);
        eps(el, 1, 2) = [1 -1] * d_CL(el:el+1);
    end
    
    eps = abs(eps) ./ max(abs(eps));
    
    % Variable to be minimised
    J = norm(eps(setdiff(1:end, damel),1,1))/norm(eps(setdiff(1:end, damel),1,2));
end
