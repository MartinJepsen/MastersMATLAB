clear
damages = [1:8]';
rng(1)
damages(:, 2) = 0.8;
nsr = 0.05;
err = 0.00;

sensor = "dis";
in_dof = [1:8];
out_dof = in_dof;
dt = 0.01;                            % time increment size
n_samples = 5000;
blockrows = 48;
n_runs = 25;

for truncated_mode = [0]
    set_up
    tic
    for foo = 1:size(damages, 1)
        damage = damages(foo, :);
%         estimate_models
%         gaindesign_2;
    %     gaindesign_3;
    gaindesign_4;
    end
    toc
end

beep
pause(1)
beep
pause(1)
beep