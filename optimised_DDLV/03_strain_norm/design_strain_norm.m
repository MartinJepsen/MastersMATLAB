clc; clear; close all;

load("optimised_DDLV/gain_pars")          % load system matrices
idx = setdiff(1:n_dof, bc);

H = SS_exact.H;
H_ = zeros(n_dof, free_dof);
H_(idx, :) = H;

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
    H_ = evalin('base', 'H_');
    H = evalin('base','H');
    cdis = evalin('base','cdis');
    
    K = reshape(X, r, m);              % gain matrix
    
    % Obtain closed-loop stiffness (damaged and undamaged)
    dKg = Kg_d - Kg;
    dKg = dKg(out_dof, in_dof);
    
    % OL transfer matrices
    G = inv(Mg * s^2 + Cg * s + Kg);        
    G_d = inv(Mg * s^2 + Cg * s + Kg_d);
    dG = G_d - G;
    
    % CL transfer matrices
    G_CL = (eye(size(G)) + G * b2*K*cdis) \ G;
    G_CL_d = (eye(size(G_d)) + G_d * b2*K*cdis) \ G_d;
    dG_CL = G_CL_d - G_CL;

    % Strain fields
    [~, ~, V] = svd(dG);
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
