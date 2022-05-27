clear
damages = [1:10]';
rng(1)
damages(:, 2) = 0.75;
nsr = 0.05;
err = 0.00;

sensor = "dis";
in_dof = [1:10];
out_dof = in_dof;
dt = 0.005;                            % time increment size
n_samples = 5000;
blockrows = 40;
n_runs = 25;
set_up

%% Run estimation
for foo = 1:size(damages, 1)
    damage = damages(foo, :);
    estimate_models
end