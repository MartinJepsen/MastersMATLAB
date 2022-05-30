% load('gaindesign/02_sens/SetUp.mat')
[DamagedModels] = generate_damaged_models(ReferenceModels(1).FE, ReferenceModels(1).FE_e, damage);


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

poles = 1:2:15;
im_fac = 0;

parfor polenum = 1:numel(poles) 
tic
pole = poles(polenum);
s = complex(real(Lambda(pole)), im_fac*imag(Lambda(pole)));

H_ref = (Mg*s^2 + Cg*s + Kg)^-1;
H = H_ref(out_dof, in_dof);

H_d = (Mg*s^2 + Cg_d*s + Kg_d)^-1;

DeltaH = H - H_d;
[~, ~, V] = svd(DeltaH);

% ga_vars.V = V;
% ga_vars.H_ref = H_ref;

%% Genetic algorithm
r = GeneralParameters(1).r;
m = GeneralParameters(1).m;
np = r*m;

options = optimoptions('ga', 'Generations', 5000,...
                        'PopulationSize', 200,...
                        'FunctionTolerance',1e-5,...
                        'CrossoverFraction',0.5,...
                        'MaxStallGenerations', 500);%,...
%                         'PlotFcn', @gaplotbestf);

[res, fval] = ga({@scheme3, ga_vars, s, V, H_ref, np}, np*2, [], [], [], [], [], [], [], options);
   
re = reshape(res(1:np), r, m);
im = reshape(res(np+1:end), r, m);
K = complex(re, im);

parsave(damage, pole, im_fac, K, s, fval)

% save(sprintf("gaindesign/03_strain_norm/gain%d_%d_%0.3f.mat", damage(1,1), polenum, im_fac),"K", "gains", "s")
toc
end


function [J] = scheme3(X, ga_vars, s, V, H_ref, np)

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
    idx        = ga_vars.idx;
    n_dof      = ga_vars.n_dof;
    free_dof   = ga_vars.free_dof;
    B          = ga_vars.B;

    re = reshape(X(1:np), r, m);
    im = reshape(X(np+1:end), r, m);
    K = complex(re, im);
    
    % DDLV
    d = H_ref * B2 * V(:, end);
    d_ = zeros(n_dof, 1);
    d_(idx) = d;
    eps = B * d_;
    eps = abs(eps) / max(abs(eps));

    % CL transfer matrices 
    H_CL_ref = (Mg*s^2 + Cg*s + Kg + (B2*K*cdis))^-1;   % full transfer matrix
    H_CL = H_CL_ref(out_dof, in_dof);                   % reduced transfer matrix
    H_CL_d = (Mg*s^2 + Cg_d*s + Kg_d + (B2*K*cdis))^-1;
    H_CL_d = H_CL_d(out_dof, in_dof);
    
    % CLDDLV
    DeltaH_CL = H_CL_d - H_CL;
    [~, ~, V] = svd(DeltaH_CL);
    H_CL_ = zeros(n_dof, free_dof);
    H_CL_(idx, :) = H_CL_ref;
    eps_CL = B * H_CL_ * B2 * V(:, end);
    eps_CL = abs(eps_CL) / max(abs(eps_CL));
    
    dofs = [1:size(B,1)];
    dofs(ga_vars.damel) = [];
    J = norm(eps(dofs)) / norm(eps_CL(dofs));

end

function parsave(damage, pole, im_fac, K, s, fval)
    fprintf("Saving gain%d_%d_%0.3f.mat", damage(1,1), pole, im_fac)
    save(sprintf("gaindesign/03_strain_norm/gain%d_%d_%0.3f.mat", damage(1,1), pole, im_fac),"K", "fval", "s")
end