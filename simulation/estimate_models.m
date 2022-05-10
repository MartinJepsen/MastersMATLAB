set_up;

omega_ref = SS_exact_d.modal_parameters.omega;
zeta_ref = SS_exact_d.modal_parameters.zeta;
omega_dev = [];
zeta_dev = [];

FE_e = FiniteElementModel();
FE_e.from_xlsx('structures/paper_truss.xlsx')
FE_e.mesh.element_properties.E = FE.mesh.element_properties.E .* unifrnd(0.98, 1.02, [n_el, 1]);
FE_e.assembly('bar', dam);
FE_e.apply_bc([1, 2, 9, 10]);
FE_e.modal_damping(0.02);

Kg_e = FE_e.Kg;
Kg_de = FE_e.Kg_d;
Cg_e = FE_e.Cg;
Cg_de = FE_e.Cg_d;
Mg_e = FE_e.Mg;
filename_u = sprintf("simulation/SYSID/model_error/00_000_%03d.mat", nsr*100);
tic
parfor run = 1:100
    % check if simulations of undamaged config exist:
    if exist(filename_u, "file") == 0
        disp('Generating reference model')
        SS = StateSpaceModel();
        SS.set_io(in_dof, out_dof);
        SS.dt_from_FE(Kg_e, Cg_e, Mg_e, dt);
        [u_n, y] = SS.time_response(u, t, nsr, false);
        SS.estimate(u_n, y, blockrows);
        SS.get_modal_parameters();
        SS.to_ct();
        SS_est{run} = SS;
        lambda_est(:, run) = SS.modal_parameters.Lambda;
    end

    SS_d = StateSpaceModel();
    SS_d.set_io(in_dof, out_dof);
    SS_d.dt_from_FE(Kg_de, Cg_de, Mg_e, dt)
    [u_n, y] = SS_d.time_response(u, t, nsr, false);
    SS_d.estimate(u_n, y, blockrows);
    SS_d.get_modal_parameters()
    SS_d.to_ct();
    SS_est_d{run} = SS_d;

    omega_dev = [omega_dev, dev(omega_ref, SS_d.modal_parameters.omega)];
    zeta_dev = [zeta_dev, dev(zeta_ref, SS_d.modal_parameters.zeta)];
    lambda_est_d(:, run) = SS_d.modal_parameters.Lambda;
end

disp(sprintf("Finished estimation in %0.2f s", toc))
%%
filename = sprintf("simulation/SYSID/model_error/%02d_%03d_%03d", dam(1,1), dam(1,2)*100, nsr*100)
save(filename, 'SS_est_d', 'lambda_est_d', 'omega_dev', 'zeta_dev')

try
    save(filename_u, 'SS_est', 'lambda_est')
catch
    disp('Reference models already exist')
end

beep