%% FE model
FE = FiniteElementModel('paper_truss.xlsx');
dam = [4, 0.999];
FE.assembly('bar', dam);
FE.apply_bc([1, 2, 9, 10]);
FE.modal_damping(0.05);

%% SS model
% Exact
SS_exact = StateSpaceModel();
SS_exact.set_io(1:12, 1:12, 24)
SS_exact.dt_from_FE(FE.Kg, FE.Cg, FE.Mg, dt);
SS_exact.get_modal_parameters();
SS_exact.time_response(u, t, 0, true);
SS_exact.transfer_matrix(1);

SS_exact_d = StateSpaceModel();
SS_exact_d.set_io(1:12, 1:12, 24)
SS_exact_d.dt_from_FE(FE.Kg_d, FE.Cg, FE.Mg, dt);
SS_exact_d.get_modal_parameters();
SS_exact_d.time_response(u, t, 0, true);
SS_exact_d.transfer_matrix(1);

DeltaG = SS_exact_d.H - SS_exact.H;
[~, ~, V] = svd(DeltaG);

d = SS_exact.H * V(:, end);
FE.strains_from_disp(d);