clear
damages = [1:14]';
rng(1)
damages(:, 2) = 0.80;
nsr = 0.05;
for foo = 1:size(damages, 1)
    dam = damages(foo, :);
    estimate_models
%     gaindesign_1
%     gaindesign_2
%     gaindesign_3
end