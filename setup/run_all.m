clear
damages = [3, 5, 7, 9]';
rng(1)
damages(:, 2) = 0.80;


for foo = 1:size(damages, 1)
    dam = damages(foo, :);
    estimate_models
end