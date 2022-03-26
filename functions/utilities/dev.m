function err = dev(ref, obs)
% dev(ref, obs) computes percentage error
err = (obs-ref)./ref*100;