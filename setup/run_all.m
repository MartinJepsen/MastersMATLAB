clear
% damages = [3, 5, 7, 9]';
rng(1)
% damages(:, 2) = 0.50;
nsr = 0.01;
damages =[7, 0.01]
for foo = 1:size(damages, 1)
    dam = damages(foo, :);
    estimate_models
end