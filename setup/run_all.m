clear
damages = [3, 5, 7, 9, 12]';
rng(1)
damages(:, 2) = 0.95;


for foo = 1:size(damages, 1)
    dam = damages(foo, :);
%     estimate_models
gaindesign_2
gaindesign_3
end