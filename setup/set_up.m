clear

rng(1)

dam = [1, 0.85];

chain_system

in_dof = [1, 3, 5];
out_dof = in_dof;
dt = 0.02;                        % time increment size
Nsamples = 50000;
t = 0:dt:(Nsamples * dt - dt);          % time sequence
u = randn(numel(in_dof), numel(t));
blockrows = 40;

nsr = 0.05;

cdis = zeros(numel(out_dof), free_dof);
for ii=1:numel(out_dof);
    cdis(ii, out_dof(ii)) = 1;
end

Kg = FE.Kg;
Kg_d = FE.Kg_d;
Cg = FE.Cg;
Mg = FE.Mg;

SS_exact = StateSpaceModel();
SS_exact.set_io(1:free_dof, 1:free_dof, 2*free_dof);
SS_exact.dt_from_FE(Kg, Cg, Mg, dt);
SS_exact.get_modal_parameters();
% SS_exact.time_response(u, t, 0, true);
SS_exact.to_ct();
Lambda = SS_exact.modal_parameters.Lambda;
s = complex(real(Lambda(1)), 1.1*imag(Lambda(1)));         % pole
SS_exact.transfer_matrix(s);

% Exact, damaged
SS_exact_d = StateSpaceModel();
SS_exact_d.set_io(1:free_dof, 1:free_dof, 2*free_dof);
SS_exact_d.dt_from_FE(FE.Kg, FE.Cg, FE.Mg, dt);
SS_exact_d.to_ct();
SS_exact_d.transfer_matrix(s);

export_gain_pars
