clear
damages = [4, 6]';
rng(1)
damages(:, 2) = 0.95;
nsr = 0.05;

s_factors = sort([1:0.05:1.4, 1.001, 1.002, 1.003], 'ascend');
dam = damages(1, :);

for foo = 1:numel(s_factors)
    s_fac = s_factors(foo);
    
%     estimate_models
    gaindesign_1
%     gaindesign_2
%     gaindesign_3
end