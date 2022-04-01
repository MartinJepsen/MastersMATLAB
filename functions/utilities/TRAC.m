function x = TRAC(ref,obs)
% TRAC(ref,obs) computes the correlation between the reference signal (ref)
% and the observed signal (obs). Signals must be equally long.

% transpose if a row vectors were supplied
if size(ref,1) < size(ref,2); ref = ref'; end
if size(obs,1) < size(obs,2); obs = obs'; end

x = (ref'*obs)^2/((ref'*ref)*(obs'*obs));


end