in_dof = [1, 2, 5, 6, 7, 8, 11, 12];
out_dof = in_dof;
dt = 0.0001;                        % time increment size
Nsamples = 50000;
t = 0:dt:(Nsamples * dt - dt);          % time sequence
u = randn(numel(in_dof), numel(t));
blockrows = 60;

FE = FiniteElementModel();
FE.from_xlsx('structures/paper_truss.xlsx')
FE.assembly('bar', dam);
FE.apply_bc([1, 2, 9, 10]);
FE.modal_damping(0.02);
FE.strains_from_disp([])

Kg = FE.Kg;
Kg_d = FE.Kg_d;
Cg = FE.Cg;
Mg = FE.Mg;

SS_exact = StateSpaceModel();
SS_exact.set_io(1:12, 1:12, 24);
SS_exact.dt_from_FE(Kg, Cg, Mg, dt);
SS_exact.get_modal_parameters();
SS_exact.to_ct();
Lambda = SS_exact.modal_parameters.Lambda;

% Exact, damaged
SS_exact_d = StateSpaceModel();
SS_exact_d.set_io(1:12, 1:12, 24);
SS_exact_d.dt_from_FE(FE.Kg, FE.Cg, FE.Mg, dt);
SS_exact_d.to_ct();
Lambda_d = SS_exact_d.modal_parameters.Lambda;

export_gain_pars