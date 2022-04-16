clear
damages = [1, 2, 12, 13, 14]';
rng(1)
damages(:, 2) = 0.95;


for foo = 1:size(damages, 1)
    dam = damages(foo, :);
    estimate_models
    gaindesign_2
    gaindesign_3
end