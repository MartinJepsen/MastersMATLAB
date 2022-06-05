clear

    damages = [1]';
    rng(1)
    damages(:, 2) = 0.4;
    
    nsr = 0.02;
    sensor = "dis";
    in_dof = [1,2,5,6,7,8,11,12];
    out_dof = in_dof;
    dt = 0.0001;                            % time increment size
    n_samples = 15000;
    blockrows = 48;
    n_runs = 50;
    truncated_mode = 0;

    for err = [20,30]/100
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
