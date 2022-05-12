set_up

s = complex(real(Lambda(1)), imag(1.12*Lambda(1)));
H = (Mg*s^2 + Cg*s + Kg)^-1;
H_ = zeros(n_dof, free_dof);
H_(idx, :) = H;

vars.n_dof = n_dof;
vars.free_dof = free_dof;
vars.m = m;
vars.r = r;
vars.B2 = B2;
vars.B = B_strain;
vars.idx  = idx ;
vars.H_ = H_;
vars.cdis = cdis;
vars.s = s;
vars.Kg = Kg;
vars.Cg = Cg;
vars.Mg = Mg;

%% Genetic algorithm
run = 0;
good = 0;
tic
for run = 0
    % change seed on every iteration
    seed = ceil(abs(randn * randn) * 10) 
    rng(seed);

    
    ObjectiveFunction = @main_gain_design;
    options = optimoptions('ga', 'Generations', 20000,...
                            'PopulationSize', 100,...
                            'CrossoverFraction', 0.6,...
                            'FunctionTolerance',1e-7);
    
    np = r*m; % No. of entries in gain matrix
    [res, fval] = ga(ObjectiveFunction, np*2, [], [], [], [], [], [], [], options);
    
    % Reshape result into complex-values r x m matrix
    re = reshape(res(1:np), r, m);
    im = reshape(res(np+1:end), r, m);
    K = complex(re, im);

    results{run+1, 1} = K;
    results{run+1, 2} = fval;
end

[fvals, ind] = sort([results{:,2}], 'ascend');
toc
%% Sort gain matrices by cost function value
for i = 1:numel(fvals)
    gains{i,2} = fvals(i);
    gains(i,1) = results(ind(i),1);
end
beep

%% Store results
savenum = 1;
K = gains{1,1};
save(sprintf("gaindesign/01_strain_cond/gains_%d", savenum),"K", "gains", "s")

function [J] = main_gain_design(X)
    % Load pre-defined variables from base workspace
    vars = evalin('base', 'vars');
    n_dof = vars.n_dof;
    free_dof = vars.free_dof;
    m = vars.m;
    r = vars.r;
    B2 = vars.B2;
    B = vars.B;
    idx  = vars.idx ;
    H_ = vars.H_;
    cdis = vars.cdis;
    s = vars.s;
    Kg = vars.Kg;
    Cg = vars.Cg;
    Mg = vars.Mg;
    
    % Reshape the function argument into a complex-valued r x m matrix.
    np = r*m;
    re = reshape(X(1:np), r, m);
    im = reshape(X(np+1:end), r, m);
    K = complex(re, im);

    % Transfer matrix shift    
    H_CL = (Mg*s^2 + Cg*s + Kg + (B2*K*cdis))^-1;
    H_CL_ = zeros(n_dof, free_dof);
    H_CL_(idx, :) = H_CL;
    
    % Load-to-strain map
    E =  B * H_ * B2;
    E_CL = B * H_CL_ * B2;
    
    % Cost function value
    J = cond(E_CL) / cond(E);

end
