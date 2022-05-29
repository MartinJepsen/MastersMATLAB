clear
damages = [1:8]';
rng(1)
damages(:, 2) = 0.8;
nsr = 0.05;
err = 0.35;

sensor = "dis";
in_dof = [1:8];
out_dof = in_dof;
dt = 0.01;                            % time increment size
n_samples = 5000;
blockrows = 48;
n_runs = 25;
set_up

%% Run estimation
for foo = 1:size(damages, 1)
    damage = damages(foo, :);
    estimate_models
end

beep
pause(0.2)
beep
pause(0.2)
beep