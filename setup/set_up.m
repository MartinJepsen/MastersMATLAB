% time = clock;
% rng(time(end)*100*rand);         % uses current time as random number seed
% dams = [96, 95, 97, 95, 97, 97, 95] / 100
% i = 1;
% for damel = [1,2,9,11,12,13,14]

rng(1)
in_dof = [1, 2, 5, 6, 7, 8, 11, 12];
% in_dof = 1:12;
out_dof = in_dof;
dt = 0.0001;                        % time increment size
Nsamples = 50000;
t = 0:dt:(Nsamples * dt - dt);          % time sequence
u = randn(numel(in_dof), numel(t));
nsr = 0.05;
blockrows = 60;
% dam = [2, 0.95];

FE = FiniteElementModel('structures/paper_truss.xlsx');
FE.assembly('bar', dam);
FE.apply_bc([1, 2, 9, 10]);
FE.modal_damping(0.02);
FE.strains_from_disp([])
cdis = zeros(numel(out_dof), 12);
for ii=1:numel(out_dof);
    cdis(ii, out_dof(ii)) = 1;
end

Kg = FE.Kg;
Kg_d = FE.Kg_d;
Cg = FE.Cg;
Mg = FE.Mg;

SS_exact = StateSpaceModel();
SS_exact.set_io(1:12, 1:12, 24);
SS_exact.dt_from_FE(Kg, Cg, Mg, dt);
SS_exact.get_modal_parameters();
% SS_exact.time_response(u, t, 0, true);
SS_exact.to_ct();
Lambda = SS_exact.modal_parameters.Lambda;
s = complex(real(Lambda(1)), 1.1*imag(Lambda(1)));         % pole
SS_exact.transfer_matrix(s);

% Exact, damaged
SS_exact_d = StateSpaceModel();
SS_exact_d.set_io(1:12, 1:12, 24);
SS_exact_d.dt_from_FE(FE.Kg, FE.Cg, FE.Mg, dt);
SS_exact_d.to_ct();
SS_exact_d.transfer_matrix(s);

export_gain_pars
% i = i+1;
% end