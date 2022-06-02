clear

    damages = [1:14]';
    rng(1)
    damages(:, 2) = 0.80;
    
    nsr = 0.05;
    err = 0.00;
    sensor = "dis";
    in_dof = [1:12];
    out_dof = in_dof;
    dt = 0.0001;                            % time increment size
    n_samples = 20000;
    blockrows = 48;
    n_runs = 50;

    for truncated_mode = [12, 10]
        set_up
        t_start = tic;
        for foo = 1:size(damages, 1)
            damage = damages(foo, :);
            estimate_models
        end
        toc(t_start)
    end

beep
pause(1)
beep
pause(1)
beep
