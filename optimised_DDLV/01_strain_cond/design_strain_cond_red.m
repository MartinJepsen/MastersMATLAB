clc; clear; close all;

load('testing/legacy_test/gain_pars.mat')
damel=dam(1)
%% Exact system
G_ref = (Mg*s^2 + Cg*s + Kg)^-1;
% G = G_ref(out_dof, in_dof);
G = G_ref;

G_d = (Mg*s^2 + Cg*s + Kg_d)^-1;
% G_d = G_d(out_dof, in_dof);

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

    free_dof = evalin('base', 'free_dof');
    n_dof = evalin('base', 'n_dof');
    G = evalin('base', 'G');

    in_dof = evalin('base', 'in_dof');
    out_dof = evalin('base', 'out_dof');
    m = numel(out_dof);
    r = numel(in_dof);
    cdis = evalin('base', 'cdis');
    B2 = evalin('base', 'B2');
    idx = evalin('base', 'idx');
    B_strain = evalin('base', 'B_strain');
    s = evalin('base', 's');

    Kg = evalin('base', 'Kg');
    Cg = evalin('base', 'Cg');
    Mg = evalin('base', 'Mg');
   
    K = reshape(X, r, m);                                       % gain matrix
   
    % Transfer matrix shift
    G_CL = (Mg*s^2 + Cg*s + Kg + B2*K*cdis)^-1;        % CL transfer matrix (reference)
    G_CL_ = zeros(n_dof, free_dof);
    G_CL_(idx, :) = G_CL;

    G_ = zeros(n_dof, free_dof);
    G_(idx, :) = G;

    E = B_strain * G_ * B2;
    E_CL = B_strain * G_CL_ * B2;
    

    J = cond(E_CL) / cond(E);
end
