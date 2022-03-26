function x = TRAC(ref,obs)
% TRAC(ref,obs) computes the correlation between the reference signal (ref)
% and the observed signal (obs). Signals must be equally long.

% transpose if a row vectors were supplied
if size(ref,1) < size(ref,2); ref = ref'; end
if size(obs,1) < size(obs,2); obs = obs'; end

if size(obs, 2) ~= size(ref,2)
    disp('ERROR (TRAC): Unequal signal lengts')
    return
end

N = size(obs, 2);
x = nan(N);
for row = 1 : size(ref, 2)
    r = ref(:, row);
    for col = 1 : size(obs, 2)
        o = obs(:, col);
        x(row, col) = (r' * o)^2 / ((r' * r)*(o' * o));
    end
end



end