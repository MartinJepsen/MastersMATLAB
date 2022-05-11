clear
damages = [4, 6]';
rng(1)
damages(:, 2) = 0.95;
nsr = 0.05;
for foo = 1:size(damages, 1)
    dam = damages(foo, :);
    estimate_models
%     gaindesign_1
%     gaindesign_2
%     gaindesign_3
end