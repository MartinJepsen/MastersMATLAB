clear
damages = [2]';
rng(1)
damages(:, 2) = 0.80;
nsr = 0.00;
err = 0.00;

% gaindesign_1
% return
for foo = 1:size(damages, 1)
    dam = damages(foo, :);
    estimate_models

%     gaindesign_2
%     gaindesign_3
end