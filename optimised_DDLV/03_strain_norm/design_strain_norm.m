clc; clear; close all;

damel = 6;

%% Exact system
n_dof = 6;
out_dof = [1, 3, 5];
m = numel(out_dof);
in_dof = [1, 3, 5];
r = numel(in_dof);

%% Genetic algorithm
run = 0;
good = 0;

for run = 0:5

    np = r*m; % No. of parameters
 
    upper = 1000;
    lower = -1000;
%     
    lb = complex(lower, lower) * ones(1, np);
    ub = complex(upper,upper) * ones(1, np); 
    nvars = numel(lb);
    ObjectiveFunction = @main_gain_design;

    options = optimoptions('ga', 'Generations', 1000, 'PopulationSize', 100, 'ConstraintTolerance', 1e-20);
    [res, fval] = ga(ObjectiveFunction, nvars, [], [], [], [], [], [], [], options);
    
    results{run+1,1} = reshape(res, r, m);
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

    load("simulation/system_matrices/unit_perturbation/model_6")          % load system matrices
    n_dof = size(Kg,1);
    in_dof = [1, 3, 5];                                         % input DOF
    out_dof = [1, 3, 5];                                        % output DOF
    m = numel(out_dof);
    r = numel(in_dof);

    b2 = zeros(n_dof, r);
    for ii=1:r
        b2(in_dof(ii), ii) = 1;
    end
    
    cdis = zeros(m, n_dof);
    for ii=1:m
        cdis(ii, out_dof(ii)) = 1;
    end

    K = reshape(X, r, m);                                       % gain matrix
    s = complex(real(Lambda(1)), 1.1*imag(Lambda(1)));
    
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
