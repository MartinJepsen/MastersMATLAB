clear
damages = [3:14]';
rng(1)
damages(:, 2) = 0.70;

nsr = 0.05;
err = 0.02;
sensor = "dis";
in_dof = [1:12];
out_dof = in_dof;
dt = 0.0001;                            % time increment size
n_samples = 20000;
blockrows = 60;
n_runs = 50;
set_up

%% Run estimation
for foo = 1:size(damages, 1)
    damage = damages(foo, :);
    estimate_models
end