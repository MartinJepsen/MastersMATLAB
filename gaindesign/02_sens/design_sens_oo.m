clc; clear; close all;

load("optimised_DDLV/gain_pars")          % load system matrices
idx = setdiff(1:n_dof, bc);
H = SS_exact.H;

dKg = Kg_d - Kg;

%% Genetic algorithm
run = 0;
good = 0;

for run = 0:5

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
save(sprintf("optimised_DDLV/02_sens/unit_perturbations/gains_%d",damel),"gains")

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
    dKG = evalin('base','dKg');
    
    K = reshape(X, r, m);                                       % gain matrix    
    
    H_CL = (eye(size(H)) + H * (B2*K*cdis))^-1 * H;

    dH = -H_CL * dKg * H_CL;                % Change in OL transfer matrix

    if cond(H_CL) > cond(HG)
        J = 1e40;
        return
    end

    % Variable to be minimised
    J = 1/norm(H_CL);
end
