set_up;

omega_ref = SS_exact_d.modal_parameters.omega;
zeta_ref = SS_exact_d.modal_parameters.zeta;
omega_dev = [];
zeta_dev = [];

base_dir = sprintf("simulation/SYSID/model_error_%03d_%s", err*100, sensor);

if ~isfolder(base_dir)
    mkdir(base_dir)
    addpath(base_dir)
end

filename_u = sprintf("00_000_%03d.mat", nsr*100);

t_0 = tic;
parfor run = 1:n_runs
    t_0_run = tic;
    % check if simulations of undamaged config exist:
    if exist(fullfile(base_dir, filename_u), "file") == 0
        disp('Generating reference model\n')
        SS = StateSpaceModel();
        SS.set_io(in_dof, out_dof);
        SS.dt_from_FE(Kg_e, Cg_e, Mg_e, dt, sensor);
        [u_n, y] = SS.time_response(u, t, nsr, false);
        SS.estimate(u_n, y, blockrows);
        SS.get_modal_parameters();
        SS.to_ct();
        SS_est{run} = SS;
        lambda_est(:, run) = SS.modal_parameters.Lambda;
    end

    SS_d = StateSpaceModel();
    SS_d.set_io(in_dof, out_dof);
    SS_d.dt_from_FE(Kg_de, Cg_de, Mg_e, dt, sensor)
    [u_n, y] = SS_d.time_response(u, t, nsr, false);
    SS_d.estimate(u_n, y, blockrows);
    SS_d.get_modal_parameters()
    SS_d.to_ct();
    SS_est_d{run} = SS_d;
    lambda_est_d(:, run) = SS_d.modal_parameters.Lambda;
    fprintf("Finished realisation no. %03d in %0.2f s\n", run, toc(t_0_run))
end
fprintf("Finished all runs in %0.2f s", toc(t_0))

%%
filename = sprintf("%02d_%03d_%03d", dam(1,1), dam(1,2)*100, nsr*100)

save(fullfile(base_dir, filename), 'SS_est_d', 'lambda_est_d', 'FE_e')

try
    save(fullfile(base_dir, filename_u), 'SS_est', 'lambda_est')
catch
    disp('Reference models already exist')
end

beep