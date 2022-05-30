clear
load("gaindesign/01_strain_cond/SetUp.mat")

for polenum = 1:2:19

im_fac = 1.12;
Lambda = ReferenceModels.Lambda;
r = GeneralParameters.r;
m = GeneralParameters.m;

try
    load('gaindesign/01_strain_cond/gains_1_1.120.mat',"K")
    initpop = [reshape(real(K), 1, numel(K)), reshape(imag(K), 1, numel(K))];
catch
    initpop = ones(1, 2*m*r);
end

s = complex(real(Lambda(polenum)), im_fac*imag(Lambda(polenum)))
GeneralParameters.s = s;

%%

Kg = ReferenceModels.Kg;
Mg = ReferenceModels.Mg;
Cg = ReferenceModels.Cg;

H = (Mg*s^2 + Cg*s + Kg)^-1;
H_ = zeros(GeneralParameters.n_dof, GeneralParameters.free_dof);
H_(GeneralParameters.idx, :) = H;

ga_vars.free_dof =  GeneralParameters.free_dof;
ga_vars.m = GeneralParameters.m;
ga_vars.r = GeneralParameters.r;
ga_vars.B2 = GeneralParameters.B2;
ga_vars.B = GeneralParameters.B_strain;
ga_vars.idx = GeneralParameters.idx;
ga_vars.H_ = H_;
ga_vars.cdis = GeneralParameters.cdis;
ga_vars.s = s;
ga_vars.Kg = Kg;
ga_vars.Cg = Cg;
ga_vars.Mg = Mg;
ga_vars.n_dof = GeneralParameters.n_dof;
ga_vars.B_strain = GeneralParameters.B_strain;

%% Genetic algorithm
tic
for run = 0
    options = optimoptions('ga', 'Generations', 20000,...
                            'PopulationSize', 500,...
                            'CrossoverFraction', 0.5,...
                            'FunctionTolerance',1e-5,...
                            'InitialPopulation', initpop,...
                            'PlotFcn', @gaplotbestf);
    
    np = r*m; % No. of entries in gain matrix
    [res, fval] = ga({@scheme1, ga_vars}, np*2, [], [], [], [], [], [], [], options);
    
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
K = gains{1,1};
save(sprintf("gaindesign/01_strain_cond/gains_%d_%0.3f.mat", polenum, im_fac),"K", "gains", "s")

end
function [J] = scheme1(X, ga_pars)
    % Load pre-defined variables from base workspace


    n_dof = ga_pars.n_dof;
    free_dof =  ga_pars.free_dof;
    m = ga_pars.m;
    r = ga_pars.r;
    B2 = ga_pars.B2;
    B = ga_pars.B_strain;
    idx = ga_pars.idx;
    cdis = ga_pars.cdis;
    Kg = ga_pars.Kg;
    Cg = ga_pars.Cg;
    Mg = ga_pars.Mg;
    H_ = ga_pars.H_;
    s = ga_pars.s;
    
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
