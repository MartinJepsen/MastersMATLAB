clc; clear; close all;

set_up

load("gaindesign/gain_pars")          % load system matrices
damel = dam(:, 1);

% reduced
SS1 = StateSpaceModel();
SS1.set_io(in_dof, out_dof);
SS1.dt_from_FE(Kg, Cg, Mg, dt);
SS1.to_ct();
SS1.transfer_matrix(s);

SS1_d = StateSpaceModel();
SS1_d.set_io(in_dof, out_dof);
SS1_d.dt_from_FE(Kg_d, Cg, Mg, dt);
SS1_d.to_ct();
SS1_d.transfer_matrix(s);

H = SS1.H;
H_d = SS1_d.H;

% full
H_ref = (Mg*s^2 + Cg*s + Kg)^-1;
H_ = zeros(n_dof, free_dof);
H_(idx, :) = H_ref;

% Apply normalisation of stiffness perturbation
DeltaKg = Kg_d - Kg;
DeltaKg(DeltaKg ~= 0) = DeltaKg(DeltaKg ~= 0) ./ abs(DeltaKg(DeltaKg ~= 0));
Kg_d = DeltaKg - Kg;

%% Genetic algorithm
run = 0;
good = 0;
tic
for run = 0:2
    seed = ceil(abs(randn * randn) * 10)
    rng(seed);
    np = r*m; % No. of parameters
    
    ObjectiveFunction = @main_gain_design;
    options = optimoptions('ga', 'Generations', 300,...
                            'PopulationSize', 100,...
                            'FunctionTolerance',0,...
                            'MaxStallGenerations', 1000,...
                            'CrossoverFraction', 0.50,...
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

%%
K = gains{1,1};
save(sprintf("gaindesign/03_strain_norm/gains_%02d", damel),"K", 'gains')
beep

function [J] = main_gain_design(X)
    in_dof = evalin('base', 'in_dof');
    out_dof = evalin('base', 'out_dof');
    free_dof = evalin('base', 'free_dof');
    n_dof = evalin('base', 'n_dof');
    m = evalin('base', 'm');
    r = evalin('base', 'r');
    B2 = evalin('base', 'B2');
    cdis = evalin('base','cdis');
    B = evalin('base', 'B_strain');
    idx = evalin('base', 'idx');

    H_ref = evalin('base', 'H_ref');
    H_ = evalin('base', 'H_');
    H = evalin('base', 'H');
    H_d = evalin('base', 'H_d');

    FE = evalin('base', 'FE');
    s = evalin('base', 's');
    Kg = evalin('base', 'Kg');
    Kg_d = evalin('base', 'Kg_d');
    Cg = evalin('base', 'Cg');
    Mg = evalin('base', 'Mg');
    damel = evalin('base', 'damel');

    np = r*m;
    re = reshape(X(1:np), r, m);
    im = reshape(X(np+1:end), r, m);
    K = complex(re, im);

    % OLDDLV
    DeltaH = H - H_d;
    [~, ~, V] = svd(DeltaH);

    d = H_ref * B2 * V(:, end);
    d_ = zeros(n_dof, 1);
    d_(idx) = d;
    eps = B * d_;
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
    
    dofs = [1:size(B,1)];
    dofs(damel) = [];
    J = norm(eps(dofs)) / norm(eps_CL(dofs));

end