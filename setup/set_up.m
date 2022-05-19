in_dof = [1, 2, 5, 6, 7, 8, 11, 12];
out_dof = in_dof;
dt = 0.0001;                            % time increment size
Nsamples = 25000;
t = 0:dt:(Nsamples * dt - dt);          % time sequence
u = randn(numel(in_dof), numel(t));
blockrows = 60;

if ~exist('dam', 'var')
    disp('No damage defined. Setting damaged configuration equal to the reference.')
    dam = [1,1];
end

generate_FE_models

m = numel(out_dof);
r = numel(in_dof);


SS_exact = StateSpaceModel();
SS_exact.set_io(1:12, 1:12, 24);
SS_exact.dt_from_FE(Kg, Cg, Mg, dt, sensor);
SS_exact.get_modal_parameters();
SS_exact.to_ct();
Lambda = SS_exact.modal_parameters.Lambda;

% Exact, damaged
SS_exact_d = StateSpaceModel();
SS_exact_d.set_io(1:12, 1:12, 24);
SS_exact_d.dt_from_FE(Kg_d, Cg_d, Mg, dt, sensor);
SS_exact_d.get_modal_parameters();
SS_exact_d.to_ct();
Lambda_d = SS_exact_d.modal_parameters.Lambda;

B2 = StateSpaceModel().set_io(in_dof, out_dof, 24);
cdis = zeros(m, free_dof);
for ii=1:m
    cdis(ii, out_dof(ii)) = 1;
end
idx = setdiff(1:n_dof, bc);
