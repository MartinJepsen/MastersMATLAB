clear
damages = [3]';
rng(1)
damages(:, 2) = 0.80;
nsr = 0.05;
for foo = 1:size(damages, 1)
    dam = damages(foo, :);
    estimate_models
end