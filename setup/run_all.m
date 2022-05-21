clear
damages = [1:14]';
rng(1)
damages(:, 2) = 0.80;
nsr = 0.05;
err = 0.02;
sensor = "acc";

n_runs = 50;
%% Run estimation

for foo = 1:size(damages, 1)
    dam = damages(foo, :);
    estimate_models

%     gaindesign_2
%     gaindesign_3
end

%% Run gain design procedure 1
gaindesign_1