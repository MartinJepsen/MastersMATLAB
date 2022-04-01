clc; clear; close all;

damel = 1;

%% Exact system
n_dof = 6;
out_dof = [1, 3, 5];
m = numel(out_dof);
in_dof = [1, 3, 5];
r = numel(in_dof);

%% Genetic algorithm
run = 0;
good = 0;
for run = 0:5;

np = r*m; % No. of parameters

upper = 1000;
lower = -1000;

lb = complex(lower,lower) * ones(1, np);
ub = complex(upper,upper) * ones(1, np); 
nvars = numel(lb);
ObjectiveFunction = @main_gain_design;
options = optimoptions('ga', 'Generations', 1000, 'PopulationSize', 100, 'ConstraintTolerance', 1e-20);

% [res, fval] = ga(ObjectiveFunction, nvars, [], [], [], [], lb, ub, [], options);
[res, fval] = ga(ObjectiveFunction, nvars, [], [], [], [], [], [], [], options);

K = reshape(res, r, m);  

results{run+1,1} = K;
results{run+1,2} = fval;

end
[fvals, ind] = sort([results{:,2}], 'ascend');

%% Sort results
for i = 1:numel(fvals)
    gains{i,2} = fvals(i);
    gains(i,1) = results(ind(i),1);
end

beep


%%
save(sprintf("optimised_DDLV/01_strain_cond/unconstrained/gains_%d",damel),"gains")

function [J] = main_gain_design(X)

    load("simulation/system_matrices/unit_perturbation/model_1")          % load system matrices
    n_dof = size(Kg,1);
    in_dof = [1, 3, 5];                                         % input DOF
    out_dof = [1, 3, 5];                                        % output DOF
    m = numel(out_dof);
    r = numel(in_dof);

    K = reshape(X, r, m);                                       % gain matrix
    s = complex(real(Lambda(1)), 1.1*imag(Lambda(1)));         % pole
    c_dis = eye(n_dof);
    
    b2 = zeros(n_dof, r);
    for ii=1:r
        b2(in_dof(ii), ii) = 1;
    end
    
    cdis = zeros(m, n_dof);
    for ii=1:m
        cdis(ii, out_dof(ii)) = 1;
    end
    
    % Transfer matrix shift
    G = (Mg *s^2 + Cg*s + Kg)^-1;        % OL transfer matrix (reference)
    G_CL = (Mg*s^2 + Cg*s + Kg + b2*K*cdis)^-1;            % CL transfer matrix (reference)

    G_d = (Mg*s^2 + Cg*s + Kg_d)^-1;    % OL transfer matrix (damaged)
    G_CL_d = (Mg*s^2 + Cg*s + Kg_d + b2*K*cdis)^-1;    % CL transfer matrix (damaged)
    
    DeltaG = G_d - G;                       % OL transfer matrix shift
    DeltaG_CL = G_CL_d - G_CL;              % CL transfer matrix shift

    % Get rank of OL transfer matrix shift
    s_OL = svd(DeltaG);
    s_OL = s_OL/max(s_OL);
    nullity_DG = nnz(~round(s_OL, 8));
    rank_DG = min(size(DeltaG)) - nullity_DG;

    [U, S, V] = svd(DeltaG_CL);
    s_CL = diag(S)/max(diag(S));
    nullity_DG_CL = nnz(~round(s_CL, 8));

    B = zeros(n_dof, n_dof+1);
    B(1:size(B)+1:numel(B)) = 1;
    B = B - circshift(B,1,2);
    
    G = [zeros(1, n_dof); G];
    G_CL = [zeros(1, n_dof); G_CL];

    E = B*G*b2;
    E_CL = B*G_CL*b2;
    
    % reject members of the population that violate the nullity constraint
%     if nullity_DG_CL ~= nullity_DG
%         J = 1e40;
%         disp("Rank changed, population member rejected")
%         return
%     end

%     if cond(E_CL) > cond(E)
%         J = 1e40;
%         disp("Conditions number incompatible")
%         return
%     else
        J = cond(E_CL) / cond(E);
%         disp("Good")
%     end
end
