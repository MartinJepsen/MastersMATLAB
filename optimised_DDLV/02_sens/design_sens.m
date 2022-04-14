clc; clear; close all;

damel = 5;
load(sprintf("simulation/system_matrices/unit_perturbation/model_%d", damel))

n_dof = 6;
out_dof = [1, 3, 5];
m = numel(out_dof);
in_dof = [1, 3, 5];
r = numel(in_dof);

s = complex(real(Lambda(1)), 0.1 + imag(Lambda(1)));         % pole
dKg = Kg_d - Kg;
dKg = dKg(out_dof, in_dof);
G = inv(Mg * s^2 + Cg * s + Kg);        % OL transfer matrix
G = G(out_dof, in_dof);


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

    options = optimoptions('ga',...
                           'Generations', 1000,...
                           'PopulationSize', 100,...
                           'PlotFcn', @gaplotbestf);

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
save(sprintf("optimised_DDLV/02_sens/unit_perturbations/gains_%d",damel),"gains")

function [J] = main_gain_design(X)

    in_dof = [1, 3, 5];                                         % input DOF
    out_dof = [1, 3, 5];                                        % output DOF
    m = numel(out_dof);
    r = numel(in_dof);

    K = reshape(X, r, m);                                       % gain matrix
    s = evalin('base', 's');
    
    % Obtain closed-loop stiffness (damaged and undamaged)
    G = evalin('base', 'G');
    dKg = evalin('base', 'dKg');
    
    G_CL = (eye(size(G))+G*K)\G;            % CL transfer matrix
    dG = -G_CL * dKg * G_CL;                % Change in OL transfer matrix

    if cond(G_CL) > cond(G)
        J = 1e40;
        return
    end

    % Variable to be minimised
    J = 1/norm(dG);
end
