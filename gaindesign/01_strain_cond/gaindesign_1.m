clear
load("gaindesign/01_strain_cond/SetUp.mat")
for polenum = 3:2:21
im_fac = 1.12;
Lambda = ReferenceModels.Lambda;
s = complex(real(Lambda(polenum)), im_fac*imag(Lambda(polenum)));

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

options = optimoptions('ga', 'Generations', 20000,...
                        'PopulationSize', 100,...
                        'CrossoverFraction', 0.5,...
                        'FunctionTolerance',1e-6,...
                        'PlotFcn', @gaplotbestf);

np = r*m; % No. of entries in gain matrix
[res, fval] = ga({@scheme1, GeneralParameters, ReferenceModels}, np*2, [], [], [], [], [], [], [], options);

% Reshape result into complex-values r x m matrix
re = reshape(res(1:np), r, m);
im = reshape(res(np+1:end), r, m);
K = complex(re, im);

results{run+1, 1} = K;
results{run+1, 2} = fval;

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

function [J] = scheme1(X, GeneralParameters, ReferenceModels)
    % Load pre-defined variables from base workspace

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
