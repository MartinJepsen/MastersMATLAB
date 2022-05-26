set_up
polenum = 5;
im_fac = 0.12;
Lambda = ReferenceModels.Lambda;
s = complex(real(Lambda(polenum)), imag(Lambda(polenum)) + imag(im_fac*Lambda(1)));

GeneralParameters.s = s;

%%

Kg = ReferenceModels.Kg;
Mg = ReferenceModels.Mg;
Cg = ReferenceModels.Cg;

H = (Mg*s^2 + Cg*s + Kg)^-1;
H_ = zeros(GeneralParameters.n_dof, GeneralParameters.free_dof);
H_(GeneralParameters.idx, :) = H;

ReferenceModels.H = H;
ReferenceModels.H_ = H_;

% H = (Mg*s^2 + Cg*s + Kg)^-1;
% n_dof = GeneralParameters.n_dof;
% free_dof =  GeneralParameters.free_dof;
% idx  = GeneralParameters.idx;
% H_ = zeros(n_dof, free_dof);
% H_(idx, :) = H;

% vars.n_dof = n_dof;
% vars.free_dof =  free_dof;
% vars.m = GeneralParameters.m;
% vars.r = GeneralParameters.r;
% vars.B2 = GeneralParameters.B2;
% vars.B = GeneralParameters.B_strain;
% vars.idx = idx;
% vars.H_ = H_;
% vars.cdis = GeneralParameters.cdis;
% vars.s = s;
% vars.Kg = Kg;
% vars.Cg = Cg;
% vars.Mg = Mg;

r = GeneralParameters.r;
m = GeneralParameters.m;
%%


%% Genetic algorithm
run = 0;
good = 0;
tic
for run = 0
    % change seed on every iteration
    seed = ceil(abs(randn * randn) * 10) 
    rng(seed);

    ObjectiveFunction = @main_gain_design;
    options = optimoptions('ga', 'Generations', 5000,...
                            'PopulationSize', 100,...
                            'CrossoverFraction', 0.6,...
                            'FunctionTolerance',1e-7);%,...
%                             'PlotFcn', @gaplotbestf);
    
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
clearvars -except GeneralParameters ReferenceModels DamagedModels K

%% Store results
K = gains{1,1};
save(sprintf("gaindesign/01_strain_cond/gains_%d_%0.3f.mat", polenum, im_fac),"K", "gains", "s")

function [J] = main_gain_design(X)
    % Load pre-defined variables from base workspace
    GeneralParameters = evalin('base', 'GeneralParameters');
    ReferenceModels = evalin('base', 'ReferenceModels');

    n_dof = GeneralParameters.n_dof;
    free_dof =  GeneralParameters.free_dof;
    m = GeneralParameters.m;
    r = GeneralParameters.r;
    B2 = GeneralParameters.B2;
    B = GeneralParameters.B_strain;
    idx = GeneralParameters.idx;
    cdis = GeneralParameters.cdis;
    Kg = ReferenceModels.Kg;
    Cg = ReferenceModels.Cg;
    Mg = ReferenceModels.Mg;
    H_ = ReferenceModels.H_;
    s = GeneralParameters.s;
    
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
