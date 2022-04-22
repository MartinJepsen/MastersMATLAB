free_dof = size(FE.Kg,1);
n_dof = FE.n_dof;
m = numel(out_dof);
r = numel(in_dof);
Kg = FE.Kg;
Cg = FE.Cg;
Mg = FE.Mg;
Kg_d = FE.Kg_d;
Lambda = SS_exact.modal_parameters.Lambda;
s = complex(real(Lambda(1)), 0.1 + imag(Lambda(1)));         % pole

B2 = StateSpaceModel().set_io(in_dof, out_dof, 24);
B_strain = FE.B;
bc = FE.mesh.bc_dof;

cdis = zeros(m, free_dof);
for ii=1:m
    cdis(ii, out_dof(ii)) = 1;
end

idx = setdiff(1:n_dof, bc);

save('gaindesign/gain_pars', 'Kg', 'Mg', 'Cg', 'Kg_d',...
    'free_dof', 'n_dof', 'm', 'r','out_dof', 'in_dof',...
    's', 'B2', 'B_strain', 'bc', 'SS_exact', 'SS_exact_d', 'cdis', 'dt', 'u', 't','FE', 'idx', 'dam')
