set_up;
Kg = ReferenceModels.Kg;
Cg = ReferenceModels.Cg;
Mg = ReferenceModels.Mg;
Kg_e = ReferenceModels.Kg_e;
Cg_e = ReferenceModels.Cg_e;

[DamagedModels] = generate_damaged_models(ReferenceModels.FE, ReferenceModels.FE_e, damage);
Kg_d = DamagedModels.Kg_d;
Cg_d = DamagedModels.Cg_d;
Kg_de = DamagedModels.Kg_de;
Cg_de = DamagedModels.Cg_de;

in_dof = GeneralParameters.in_dof;
out_dof = GeneralParameters.out_dof;
dt = GeneralParameters.dt;
u = GeneralParameters.u;
t = GeneralParameters.t;
blockrows = GeneralParameters.blockrows;

base_dir = GeneralParameters.base_dir;
filename_u = sprintf("00_000_%03d.mat", round(nsr*100,0));

t_0 = tic;
if truncated_mode == 0
    order = 24;
else
    order = truncated_mode*2;
end

parfor run = 1:n_runs
    % check if simulations of undamaged config exist:
    if exist(fullfile(base_dir, filename_u), "file") == 0
        disp('Generating reference model')
        SS = StateSpaceModel();
        SS.set_io(in_dof, out_dof);
        SS.dt_from_FE(Kg_e, Cg_e, Mg, dt, sensor);
        [u_n, y] = SS.time_response(u, t, nsr, false);
        SS.estimate(u_n, y, blockrows, order);
        SS.get_modal_parameters();
        SS.to_ct();
        SS_est{run} = SS;
        lambda_est(:, run) = SS.modal_parameters.Lambda;
    end

    SS_d = StateSpaceModel();
    SS_d.set_io(in_dof, out_dof);
    SS_d.dt_from_FE(Kg_de, Cg_de, Mg, dt, sensor)
    [u_n, y] = SS_d.time_response(u, t, nsr, false);
    SS_d.estimate(u_n, y, blockrows, order);
    SS_d.get_modal_parameters()
    SS_d.to_ct();
    SS_est_d{run} = SS_d;
    lambda_est_d(:, run) = SS_d.modal_parameters.Lambda;
%     fprintf("Finished realisation no. %03d in %0.2f s\n", run, toc(t_0_run))
end
fprintf("Finished all runs in %0.2f s", toc(t_0))

%%
filename = sprintf("%02d_%03d_%03d.mat", damage(1,1), damage(1,2)*100, round(nsr*100,0));
fprintf("Saving %s\n", filename)
save(fullfile(base_dir, filename), 'SS_est_d', 'lambda_est_d', 'DamagedModels')

try
    save(fullfile(base_dir, filename_u), 'SS_est', 'lambda_est', "ReferenceModels")
    fprintf("Saving %s\n", filename_u)
catch
end

beep