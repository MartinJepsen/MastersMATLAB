clear
damages = [3, 5, 7, 9, 12]';
rng(1)
damages(:, 2) = 0.95;


for i = 1:numel(damages)
    dam = damages(i, :);
    estimate_models
end