clear

damages = [1:14]';
rng(1)
damages(:, 2) = 0.8;

nsr = 0.02;
sensor = "dis";
in_dof = [2,6];
out_dof = [1:12];
dt = 0.0001;                            % time increment size
n_samples = 15000;
blockrows = 48;
n_runs = 50;
truncated_mode = 4;
err = [5]/100;
set_up


t_start = tic;
for foo = 1:size(damages, 1)
    damage = damages(foo, :);
%     estimate_models
%             gaindesign_1
%             gaindesign_2
            gaindesign_3
end
toc(t_start)

beep
pause(1)
beep
pause(1)
beep
