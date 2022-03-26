clc; clear; close all;

load("gaindesign/gain_pars")          % load system matrices
idx = setdiff(1:n_dof, bc);

H = (Mg*s^2 + Cg*s + Kg)^-1;
H_ = zeros(n_dof, free_dof);
H_(idx, :) = H;

H_d = (Mg*s^2 + Cg*s + Kg_d)^-1;

%% Genetic algorithm
run = 0;
good = 0;
tic
for run = 0:2
    seed = ceil(abs(randn * randn) * 10)
    rng(seed);
    np = r*m; % No. of parameters
    
    ObjectiveFunction = @main_gain_design;
    options = optimoptions('ga', 'Generations', 20000,...
                            'PopulationSize', 100,...
                            'FunctionTolerance',1e-5,...
                            'PlotFcn', @gaplotbestf);
    
    % [res, fval] = ga(ObjectiveFunction, nvars, [], [], [], [], lb, ub, [], options);
    [res, fval] = ga(ObjectiveFunction, np*2, [], [], [], [], [], [], [], options);
    
    % K = reshape(res, r, m);  
    
        re = reshape(res(1:np), r, m);
        im = reshape(res(np+1:end), r, m);
        K = complex(re, im);
    
    results{run+1,1} = K;
    results{run+1,2} = fval;

end
[fvals, ind] = sort([results{:,2}], 'ascend');
toc
%% Sort results
for i = 1:numel(fvals)
    gains{i,2} = fvals(i);
    gains(i,1) = results(ind(i),1);
end

beep


%%
savenum = 3;
K = gains{1,1};
save(sprintf("optimised_DDLV/01_strain_cond/gains_%d", savenum),"K", 'gains')


function [J] = main_gain_design(X)

    free_dof = evalin('base', 'free_dof');
    n_dof = evalin('base', 'n_dof');
    m = evalin('base', 'm');
    r = evalin('base', 'r');
    B2 = evalin('base', 'B2');
    B = evalin('base', 'B_strain');
    idx = evalin('base', 'idx');
    H_ = evalin('base', 'H_');
    H = evalin('base', 'H');
    H_d = evalin('base', 'H_d');
    cdis = evalin('base','cdis');
    FE = evalin('base', 'FE');
    s = evalin('base', 's');
    Kg = evalin('base', 'Kg');
    Cg = evalin('base', 'Cg');
    Mg = evalin('base', 'Mg');
    in_dof = evalin('base', 'in_dof');
    out_dof = evalin('base', 'out_dof');
    
    np = r*m;
    re = reshape(X(1:np), r, m);
    im = reshape(X(np+1:end), r, m);
    K = complex(re, im);

    % OLDDLV
    DeltaH = H - H_d;
    [~, ~, V] = svd(DeltaH);

    eps = B * H_ * B2 * V(:, end);
    eps = abs(eps) / max(abs(eps));

    % CL transfer matrices 
    H_CL_ref = (Mg*s^2 + Cg*s + Kg + (B2*K*cdis))^-1;   % full transfer matrix
    H_CL = H_CL_ref(out_dof, in_dof);                   % reduced transfer matrix

    H_CL_d = (Mg*s^2 + Cg*s + Kg_d + (B2*K*cdis))^-1;
    H_CL_d = H_CL_d(out_dof, in_dof);
    
    % CLDDLV
    DeltaH_CL = H_CL_d - H_CL;
    [~, ~, V] = svd(DeltaH_CL);
    H_CL_ = zeros(n_dof, free_dof);
    H_CL_(idx, :) = H_CL_ref;
    eps_CL = B * H_CL_ * B2 * V(:, end);
    eps_CL = abs(eps_CL) / max(abs(eps_CL));

    J = norm(eps(setdiff(1:end, damel),1,1))/norm(eps(setdiff(1:end, damel),1,2));

end