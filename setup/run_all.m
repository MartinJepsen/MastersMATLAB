clear

    damages = [1:14]';
    rng(1)
    damages(:, 2) = 0.4;
    
    sensor = "dis";
    in_dof = [1,2,5,6,7,8,11,12];
    out_dof = in_dof;
    dt = 0.0001;                            % time increment size
    n_samples = 15000;
    blockrows = 48;
    n_runs = 50;
    truncated_mode = 0;

errs = [0, 10, 20, 30, 40]/100; 
nsrs = [2, 5]/100;
dams = [70, 60, 50]/100;


for i_d = 1:numel(dams)
    damages(:, 2) = dams(i_d);
    for i_n = 1:numel(nsrs)
        nsr = nsrs(i_n);
        for i_e = 1:numel(errs)
            err = errs(i_e);

            set_up
            t_start = tic;
            for foo = 1:size(damages, 1)
                damage = damages(foo, :);
                estimate_models
            end
            toc(t_start)
        end
    end
end

beep
pause(1)
beep
pause(1)
beep
