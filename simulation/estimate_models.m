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
    SS.time_response(u, t, nsr, true);
    SS.estimate([],[], blockrows);
    SS.get_modal_parameters();
    SS.to_ct();
    SS.transfer_matrix(s);
    H = SS.H;

    SS_d = StateSpaceModel();
    SS_d.set_io(in_dof, out_dof);
    SS_d.dt_from_FE(Kg_d, Cg, Mg, dt)
    SS_d.time_response(u, t, nsr, true);
    SS_d.estimate([],[], blockrows);
    SS_d.to_ct();
    SS_d.transfer_matrix(s);
    H_d = SS_d.H;

    omega_dev = [omega_dev, dev(omega_ref, SS.modal_parameters.omega)];
    zeta_dev = [zeta_dev, dev(zeta_ref, SS.modal_parameters.zeta)];
    lambda_est(:, run) = SS.modal_parameters.Lambda;

    H_est{run, 1} = H;
    H_est_d{run, 1} = H_d;

    SS_est{run} = SS;
    SS_est_d{run} = SS_d;
    toc
end

%%
filename = sprintf("simulation/SYSID/%02d_%03d_%03d",dam(1,1), dam(1,2)*100, nsr*100)
% save(filename, 'H_est', 'H_est_d', 'lambda_est', 'omega_dev', 'zeta_dev')
beep