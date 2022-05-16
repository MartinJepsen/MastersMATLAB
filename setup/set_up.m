in_dof = [1, 2, 5, 6, 7, 8, 11, 12];
out_dof = in_dof;
dt = 0.0001;                            % time increment size
Nsamples = 50000;
t = 0:dt:(Nsamples * dt - dt);          % time sequence
u = randn(numel(in_dof), numel(t));
blockrows = 60;

if ~exist('dam', 'var')
    disp('No damage defined. Setting damaged configuration equal to the reference.')
    dam = [1,1];
end

FE = FiniteElementModel();
FE.from_xlsx('structures/paper_truss.xlsx')
FE.assembly('bar', dam);
FE.apply_bc([1, 2, 9, 10]);
FE.modal_damping(0.02);
FE.strains_from_disp([])

free_dof = size(FE.Kg, 1);
n_dof = FE.n_dof;
m = numel(out_dof);
r = numel(in_dof);
n_el = size(FE.mesh.topology, 1);
B_strain = FE.B;
bc = FE.mesh.bc_dof;

Kg = FE.Kg;
Kg_d = FE.Kg_d;
Cg = FE.Cg;
Cg_d = FE.Cg_d;
Mg = FE.Mg;

SS_exact = StateSpaceModel();
<<<<<<< HEAD
SS_exact.set_io(1:12, 1:12, 24);
SS_exact.dt_from_FE(Kg, Cg, Mg, dt);
=======
SS_exact.set_io(1:n_dof, 1:n_dof, 2 * n_dof);
SS_exact.dt_from_FE(Kg, Cg, Mg, dt, "acc");
>>>>>>> fb20e35 (implement acceleration and velocity output)
SS_exact.get_modal_parameters();
SS_exact.to_ct();
Lambda = SS_exact.modal_parameters.Lambda;

% Exact, damaged
SS_exact_d = StateSpaceModel();
<<<<<<< HEAD
SS_exact_d.set_io(1:12, 1:12, 24);
SS_exact_d.dt_from_FE(FE.Kg_d, FE.Cg_d, FE.Mg, dt);
=======
SS_exact_d.set_io(1:n_dof, 1:n_dof, 2 * n_dof);
SS_exact_d.dt_from_FE(Kg_d, Cg_d, Mg, dt, "acc");
>>>>>>> fb20e35 (implement acceleration and velocity output)
SS_exact_d.get_modal_parameters();
SS_exact_d.to_ct();
Lambda_d = SS_exact_d.modal_parameters.Lambda;

B2 = StateSpaceModel().set_io(in_dof, out_dof, 24);
cdis = zeros(m, free_dof);
for ii=1:m
    cdis(ii, out_dof(ii)) = 1;
end
idx = setdiff(1:n_dof, bc);
