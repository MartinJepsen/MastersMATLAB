clc; clear; close all;

damel = 1;
load(sprintf("simulation/system_matrices/unit_perturbation/model_%d", damel))
s = complex(real(Lambda(1)), 0.1 + imag(Lambda(1)));         % pole




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
options = optimoptions('ga',...
                       'Generations', 1000,...
                       'PopulationSize', 100,...
                       'PlotFcn', 'gaplotbestf');

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
save(sprintf("optimised_DDLV/01_strain_cond/reduced/gains_%d",damel),"gains")

function [J] = main_gain_design(X)

    Kg = evalin('base', 'Kg');
    Cg = evalin('base', 'Cg');
    Mg = evalin('base', 'Mg');
    Kg_d = evalin('base', 'Kg_d');
    s = evalin('base', 's');
    
    n_dof = size(Kg,1);
    in_dof = [1, 3, 5];                                         % input DOF
    out_dof = [1, 3, 5];                                        % output DOF
    m = numel(out_dof);
    r = numel(in_dof);

    K = reshape(X, r, m);                                       % gain matrix
    
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

    B = zeros(n_dof, n_dof+1);
    B(1:size(B)+1:numel(B)) = -1;
    B = B - circshift(B,1,2);
    
    G = [zeros(1, n_dof); G];
    G_CL = [zeros(1, n_dof); G_CL];

    E = B*G*b2;
    E_CL = B*G_CL*b2;
    
    % reject members of the population that violate the nullity constraint

    J = cond(E_CL) / cond(E);
end
