clear
blabla = [0.5, 0.6, 0.75, 0.80, 0.90, 0.99, 1];
for abcd = 1:numel(blabla)
damages = [1:14]';
rng(1)
damages(:, 2) = blabla(abcd);

nsr = 0.05;
err = 0.00;
sensor = "dis";
in_dof = [1:12];
out_dof = in_dof;
dt = 0.0001;                            % time increment size
n_samples = 20000;
blockrows = 48;
n_runs = 50;
truncate = false;

base_dir = sprintf("simulation/SYSID/model_error_%03d_%s", err*100, sensor);
set_up

%% Run estimation
for foo = 1:size(damages, 1)
    damage = damages(foo, :);
    estimate_models
end
end