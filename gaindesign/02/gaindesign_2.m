% load('gaindesign/02_sens/SetUp.mat')
[DamagedModels] = generate_damaged_models(ReferenceModels(1).FE, ReferenceModels(1).FE_e, damage);

expand = false;
if expand
    GeneralParameters = expand_coordinates(GeneralParameters);
end

Mg = ReferenceModels(1).Mg;
Cg = ReferenceModels(1).Cg;
Kg = ReferenceModels(1).Kg;
Cg_d = DamagedModels(1).Cg_d;
Kg_d = DamagedModels(1).Kg_d;
DeltaKg = Kg_d - Kg;
DeltaKg(DeltaKg ~= 0) = DeltaKg(DeltaKg ~= 0) ./ abs(DeltaKg(DeltaKg ~= 0));
Kg_d = DeltaKg - Kg;

out_dof = GeneralParameters(1).out_dof;
in_dof = GeneralParameters(1).in_dof;

% collect ga variables
ga_vars.idx = GeneralParameters(1).idx;
ga_vars.n_dof = GeneralParameters(1).n_dof;
ga_vars.free_dof = GeneralParameters(1).free_dof;
ga_vars.B = GeneralParameters(1).B_strain;
ga_vars.damel = DamagedModels(1).damage(1,1);
ga_vars.r = GeneralParameters(1).r;
ga_vars.m = GeneralParameters(1).m;
ga_vars.in_dof = in_dof;
ga_vars.out_dof = out_dof;
ga_vars.B2 = GeneralParameters(1).B2;
ga_vars.cdis = GeneralParameters(1).cdis;
ga_vars.Mg = Mg;
ga_vars.Cg = Cg;
ga_vars.Kg = Kg;
ga_vars.Cg_d = Cg_d;
ga_vars.Kg_d = Kg_d;

r = GeneralParameters(1).r;
m = GeneralParameters(1).m;
np = r*m;

% ga options
options = optimoptions('ga', 'Generations', 5000,...
                        'PopulationSize', 200,...
                        'FunctionTolerance',1e-5,...
                        'CrossoverFraction',0.5,...
                        'MaxStallGenerations', 500);%,...
%                         'PlotFcn', @gaplotbestf);
poles = 1:2:15;
Lambda = ReferenceModels(1).Lambda;
im_fac = 1;

parfor polenum = 1:numel(poles) 
pole = poles(polenum);
s = complex(0*real(Lambda(pole)), im_fac*imag(Lambda(pole)));
H_ref = (Mg*s^2 + Cg*s + Kg)^-1;
H = H_ref(out_dof, in_dof);
dH = -H * DeltaKg * H;

%% Genetic algorithm
[res, fval] = ga({@scheme2, ga_vars, s, H, dH, DeltaKg}, np*2, [], [], [], [], [], [], [], options);
   
re = reshape(res(1:np), r, m);
im = reshape(res(np+1:end), r, m);
K = complex(re, im);

parsave(damage, pole, im_fac, K, s, fval, expand)
toc
end

function [J] = scheme2(X, ga_vars, s, H, dH, DeltaKg)

    r          = ga_vars.r;
    m          = ga_vars.m;
    in_dof     = ga_vars.in_dof;
    out_dof    = ga_vars.out_dof;
    B2         = ga_vars.B2;
    cdis       = ga_vars.cdis;
    Mg         = ga_vars.Mg;
    Cg         = ga_vars.Cg;
    Kg         = ga_vars.Kg;
    B          = ga_vars.B;
    
    % rearrange function argument to complex-valued gain matrix
    np = r*m;
    re = reshape(X(1:np), r, m);
    im = reshape(X(np+1:end), r, m);
    K = complex(re, im);
    
    % OL

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

function parsave(damage, pole, im_fac, K, s, fval, expand)
    fprintf("Saving gain%d_%d_%0.3f.mat\n", damage(1,1), pole, im_fac)
    if expand
        save(sprintf("gaindesign/02/exp_gains/gain%d_%d_%0.3f.mat", damage(1,1), pole, im_fac),"K", "fval", "s")
    else
        save(sprintf("gaindesign/02/gains/gain%d_%d_%0.3f.mat", damage(1,1), pole, im_fac),"K", "fval", "s")
    end
    
end
