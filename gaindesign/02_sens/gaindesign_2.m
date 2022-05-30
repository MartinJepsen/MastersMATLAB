% load('gaindesign/02_sens/SetUp.mat')
[DamagedModels] = generate_damaged_models(ReferenceModels.FE, ReferenceModels.FE_e, damage);

for polenum = 1:2:15

im_fac = 1.2;
Lambda = ReferenceModels.Lambda;
s = complex(real(Lambda(polenum)), im_fac*imag(Lambda(polenum)));

GeneralParameters.s = s;




% Apply normalisation of stiffness perturbation
Kg = ReferenceModels.Kg;
Kg_d = DamagedModels.Kg_d;
DeltaKg = Kg_d - Kg;
DeltaKg(DeltaKg ~= 0) = DeltaKg(DeltaKg ~= 0) ./ abs(DeltaKg(DeltaKg ~= 0));
Kg_d = DeltaKg - Kg;
DeltaKg = DeltaKg(out_dof, in_dof);
DamagedModels.DeltaKg = DeltaKg;

Mg = ReferenceModels.Mg;
Cg = ReferenceModels.Cg;
Kg = ReferenceModels.Kg;
out_dof = GeneralParameters.out_dof;
in_dof = GeneralParameters.in_dof;

H = (Mg*s^2 + Cg*s + Kg)^-1;
H = H(out_dof, in_dof);
dH = -H * DeltaKg * H;
ReferenceModels.dH = dH;
ReferenceModels.H = H;

%% Genetic algorithm
tic

for run = 0
r = GeneralParameters.r;
m = GeneralParameters.m;
np = r*m;

ObjectiveFunction = @main_gain_design;
options = optimoptions('ga', 'Generations', 5000,...
                        'PopulationSize', 200,...
                        'FunctionTolerance',1e-5,...
                        'CrossoverFraction',0.5,...
                        'MaxStallGenerations', 200);%,...
%                         'PlotFcn', @gaplotbestf);

[res, fval] = ga(ObjectiveFunction, np*2, [], [], [], [], [], [], [], options);
   
re = reshape(res(1:np), r, m);
im = reshape(res(np+1:end), r, m);
K = complex(re, im);

results{run+1,1} = K;
results{run+1,2} = fval;
end

[fvals, ind] = sort([results{:,2}], "ascend");
toc

%% Sort and save results
for i = 1:numel(fvals)
    gains{i,2} = fvals(i);
    gains(i,1) = results(ind(i),1);
end

K = gains{1,1};
save(sprintf("gaindesign/02_sens/gain%d_%d_%0.3f.mat", damage(1,1), polenum, im_fac),"K", "gains", "s")
beep
end


function [J] = main_gain_design(X)

    GeneralParameters = evalin('base', 'GeneralParameters');
    ReferenceModels = evalin('base', 'ReferenceModels');
    DeltaKg = evalin('base', 'DeltaKg');
    r = GeneralParameters.r;
    m = GeneralParameters.m;
    s = GeneralParameters.s;
    in_dof = GeneralParameters.in_dof;
    out_dof = GeneralParameters.out_dof;
    B2 = GeneralParameters.B2;
    cdis = GeneralParameters.cdis;
    Mg = ReferenceModels.Mg;
    Cg = ReferenceModels.Cg;
    Kg = ReferenceModels.Kg;

    
    % rearrange function argument to complex-valued gain matrix
    np = r*m;
    re = reshape(X(1:np), r, m);
    im = reshape(X(np+1:end), r, m);
    K = complex(re, im);
    
    % OL
    dH = ReferenceModels.dH;
    H = ReferenceModels.H;

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
