% clear

damages = [1:14]';
rng(1)
damages(:, 2) = 0.40;

nsr = 0.02;
err = 0.02;
sensor = "dis";
in_dof = [1:12];
out_dof = in_dof;
dt = 0.0001;                            % time increment size
n_samples = 15000;
blockrows = 48;
n_runs = 50;

for truncated_mode = [0]
    set_up
    t1 = tic;
    for foo = 1:size(damages, 1)
        damage = damages(foo, :);
        estimate_models
%         gaindesign_2;
%         gaindesign_3;
    end
    toc(t1)
end

beep
pause(1)
beep
pause(1)
beep
