clear
set_up;

omega_ref = SS_exact.modal_parameters.omega;
zeta_ref = SS_exact.modal_parameters.zeta;

omega_dev = [];
zeta_dev = [];

for run = 1:100
    run
    tic
    SS = StateSpaceModel();
    SS.set_io(in_dof, out_dof);
    SS.dt_from_FE(Kg, Cg, Mg, dt);
    [u, y] = SS.time_response(u, t, nsr, false);
    SS.estimate(u, y, blockrows);
    SS.get_modal_parameters();
    SS.to_ct();

    SS_d = StateSpaceModel();
    SS_d.set_io(in_dof, out_dof);
    SS_d.dt_from_FE(Kg_d, Cg, Mg, dt)
    [u, y] = SS_d.time_response(u, t, nsr, false);
    SS_d.estimate(u, y, blockrows);
    SS_d.to_ct();

    omega_dev = [omega_dev, dev(omega_ref, SS.modal_parameters.omega)];
    zeta_dev = [zeta_dev, dev(zeta_ref, SS.modal_parameters.zeta)];
    lambda_est(:, run) = SS.modal_parameters.Lambda;

    SS_est{run} = SS;
    SS_est_d{run} = SS_d;
    toc
end

%%
filename = sprintf("simulation/SYSID/%02d_%03d_%03d",dam(1,1), dam(1,2)*100, nsr*100)
save(filename, 'SS_est', 'SS_est_d', 'lambda_est', 'omega_dev', 'zeta_dev')
beep