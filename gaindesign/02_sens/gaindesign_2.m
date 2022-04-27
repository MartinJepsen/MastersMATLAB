set_up;

load("gaindesign/gain_pars")          % load system matrices
damel = dam(:, 1);

% Apply normalisation of stiffness perturbation
DeltaKg = Kg_d - Kg;
DeltaKg(DeltaKg ~= 0) = DeltaKg(DeltaKg ~= 0) ./ abs(DeltaKg(DeltaKg ~= 0));
Kg_d = DeltaKg - Kg;
DeltaKg = DeltaKg(out_dof, in_dof);

%% Genetic algorithm
run = 0;
good = 0;
tic
for run = 0
    seed = ceil(abs(randn * randn) * 10)
    rng(seed);
    np = r*m; % No. of parameters
    
    ObjectiveFunction = @main_gain_design;
    options = optimoptions('ga', 'Generations', 5000,...
                            'PopulationSize', 200,...
                            'FunctionTolerance',1e20,...
                            'MaxStallGenerations', 200,...
                            'PlotFcn', @gaplotbestf);
    
    [res, fval] = ga(ObjectiveFunction, np*2, [], [], [], [], [], [], [], options);
       
    re = reshape(res(1:np), r, m);
    im = reshape(res(np+1:end), r, m);
    K = complex(re, im);
    
    results{run+1,1} = K;
    results{run+1,2} = fval;

end

[fvals, ind] = sort([results{:,2}], 'ascend');
toc

%% Sort and save results
for i = 1:numel(fvals)
    gains{i,2} = fvals(i);
    gains(i,1) = results(ind(i),1);
end

K = gains{1,1};
save(sprintf("gaindesign/02_sens/constrained/gains_%02d", damel),"K", 'gains')
beep

function [J] = main_gain_design(X)

    in_dof = evalin('base', 'in_dof');
    out_dof = evalin('base', 'out_dof');
    m = evalin('base', 'm');
    r = evalin('base', 'r');
    B2 = evalin('base', 'B2');
    cdis = evalin('base','cdis');
    DeltaKg = evalin('base','DeltaKg');
    s = evalin('base', 's');
    Kg = evalin('base', 'Kg');
    Cg = evalin('base', 'Cg');
    Mg = evalin('base', 'Mg');
    
    % rearrange function argument to complex-valued gain matrix
    np = r*m;
    re = reshape(X(1:np), r, m);
    im = reshape(X(np+1:end), r, m);
    K = complex(re, im);
    
    % OL
    H = (Mg*s^2 + Cg*s + Kg)^-1;
    H = H(out_dof, in_dof);
    dH = -H * DeltaKg * H;

    % CL
    H_CL = (Mg*s^2 + Cg*s + Kg + B2*K*cdis)^-1;
    H_CL = H_CL(out_dof, in_dof);
    dH_CL = -H_CL * DeltaKg * H_CL;                 % CL transfer matrix sensitivity towards unit stiffness perturbation
    
    % reject population member if it increases the transfer matrix condition number
    if cond(H_CL) > cond(H)
        J = 1;
        return
    end

    J = norm(dH) / norm(dH_CL);
end
