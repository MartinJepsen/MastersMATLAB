function y = noise(signal, nsr)
%     Y = NOISE(SIGNAL, NSR)
% time = clock;
% rng(time(end)*100*rand);         % uses current time as random number seed
% rng(1);
[n, N] = size(signal);

if N < n; signal = signal'; end

[n, N] = size(signal);

sig = std(signal, 0, 2);% standard deviation of each channel

snr = sig.*nsr;

v = randn(n,N).*snr;     % noise

y = signal+v;