clear
blabla = [10]/100;

for abcd = 1:numel(blabla)
    damages = [6:14]';
    rng(1)
    damages(:, 2) = 0.80;
    
    nsr = blabla(abcd);
    err = 0.00;
    sensor = "dis";
    in_dof = [1:12];
    out_dof = in_dof;
    dt = 0.0001;                            % time increment size
    n_samples = 20000;
    blockrows = 48;
    n_runs = 50;

    for truncated_mode = [0]
        set_up
        tic
        for foo = 1:size(damages, 1)
            damage = damages(foo, :);
            estimate_models
        %     gaindesign_2;
        %     gaindesign_3;
        end
        toc
    end
end

beep
pause(1)
beep
pause(1)
beep
