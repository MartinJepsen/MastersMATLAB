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

Lambda = ReferenceModels(1).Lambda;

% collect ga variables

ga_vars.free_dof = GeneralParameters(1).free_dof;
ga_vars.B = GeneralParameters(1).B_strain;
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

poles = 1:2:13;
im_fac = 1.12;

for polenum = 1:numel(poles) 
tic
pole = poles(polenum);
s = complex(real(Lambda(pole)), im_fac*imag(Lambda(pole)));

H_ref = (Mg*s^2 + Cg*s + Kg)^-1;
H = H_ref(out_dof, in_dof);

H_d = (Mg*s^2 + Cg_d*s + Kg_d)^-1;
H_d = H_d(out_dof, in_dof);

DeltaH = H - H_d;



%% Genetic algorithm
r = GeneralParameters(1).r;
m = GeneralParameters(1).m;
np = r*m;

rng default
options = optimoptions('ga', 'Generations', 100000,...
                        'PopulationSize', 10,...
                        'FunctionTolerance',1e-7,...
                        'CrossoverFraction',0.5,...
                        'MaxStallGenerations', 500);%,'PlotFcn', @gaplotbestf);

[res, fval] = ga({@scheme4, ga_vars, s, np, H_ref, DeltaH}, np*2, [], [], [], [], [], [], [], options);
   
re = reshape(res(1:np), r, m);
im = reshape(res(np+1:end), r, m);
K = complex(re, im);

parsave(damage, pole, im_fac, K, s, fval, expand)

toc
end


function [J] = scheme4(X, ga_vars, s, np, H_ref, DeltaH)

    r          = ga_vars.r;
    m          = ga_vars.m;
    in_dof     = ga_vars.in_dof;
    out_dof    = ga_vars.out_dof;
    B2         = ga_vars.B2;
    cdis       = ga_vars.cdis;
    Mg         = ga_vars.Mg;
    Cg         = ga_vars.Cg;
    Kg         = ga_vars.Kg;
    Cg_d       = ga_vars.Cg_d;
    Kg_d       = ga_vars.Kg_d;

    re = reshape(X(1:np), r, m);
    im = reshape(X(np+1:end), r, m);
    K = complex(re, im);
    

    % CL transfer matrices 
    H_CL_ref = (Mg*s^2 + Cg*s + Kg + (B2*K*cdis))^-1;   % full transfer matrix
    H_CL = H_CL_ref(out_dof, in_dof);                   % reduced transfer matrix
    H_CL_d = (Mg*s^2 + Cg_d*s + Kg_d + (B2*K*cdis))^-1;
    H_CL_d = H_CL_d(out_dof, in_dof);
    
    % CLDDLV
    DeltaH_CL = H_CL_d - H_CL;

    if cond(DeltaH_CL) > cond(DeltaH)
        J = 1;
        return
    end
    J = norm(DeltaH) / norm(DeltaH_CL);

end

function parsave(damage, pole, im_fac, K, s, fval, expand)
    fprintf("Saving gain%d_%d_%0.3f.mat", damage(1,1), pole, im_fac)
    if expand
        save(sprintf("gaindesign/04/exp_gains/gain%d_%d_%0.3f.mat", damage(1,1), pole, im_fac),"K", "fval", "s")
    else
        save(sprintf("gaindesign/04/gains/gain%d_%d_%0.3f.mat", damage(1,1), pole, im_fac),"K", "fval", "s")
    end
end