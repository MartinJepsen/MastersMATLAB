% clear 
% close all
tic
% Updated 31_01_2022 11:00
%% Create new save folder
simnum = 0;
savedir = sprintf("simulation/SYSID/results_%03d", simnum);
while exist(fullfile(savedir), 'dir') == 7
    simnum = simnum + 1;
    savedir = sprintf("simulation/SYSID/results_%03d", simnum);
end
mkdir(savedir);
addpath(savedir);

%% Load systems

load('testing/legacy_test/gain_pars.mat')
damel = dam(1);

resultdir = savedir + sprintf("/%02d",[damel]);
mkdir(resultdir)
addpath(resultdir)



%% Frequency and time controls
omegas = sqrt(eig(Kg,Mg));
fmax = max(omegas)/(2*pi);          % highest frequency in the response
tmax = min(omegas/(2*pi))^-1;       % longest period in the response
dt = 0.0001;                          % time increment size
Nsamples = 50000;                   %
t = 0:dt:(Nsamples*dt-dt);                       % time sequence
% t = t(1:2^(nextpow2(numel(t))-1));  % Trim sample length to numel(t) = 2^n where n is an integer.
fn = (1/dt)/2;                      % Nyquist frequency

%% Input/output controls
freedof = size(Kg,1);
in_dof = [1, 3, 5];
out_dof = [1, 3, 5];

z0 = zeros(2*freedof,1);                  % initial conditions
rng(1);
u_raw = randn(numel(in_dof),numel(t));  % unfiltered input

%% Filter design
fc = 2*fmax;                          % Cutoff frequency

[z,p,k] = butter(3, fc/fn, 'low');      % transfer function den/num
[b,a] = zp2tf(z,p,k);
u = filter(b,a,u_raw,[],2);             % filtered input

filtering = false;

%% Simulation of time response
% modelorder = 2*size(Kg,1);
modelorder = 60;
nsr = 0.05;

[d, ~, ~, sys_ex_u] = SS_timeseries(Kg, Cg, Mg, z0, u_raw, t, in_dof);
[d_d, ~, ~, sys_ex_d] = SS_timeseries(Kg_d, Cg, Mg, z0, u_raw, t, in_dof);

d = d(out_dof,:);
d_d = d_d(out_dof,:);

model.z0 = z0;
model.t = t;
model.dt = dt;
model.in_dof = in_dof;
model.out_dof = out_dof;
model.nsr = nsr;
model.dam = dam;
model.damel = damel;
model.modelorder = modelorder;
model.filtering = filtering;

for i = 0:99
clear sys_d sys_u A B C D
% add noise
d_noise = noise(d, nsr);
d_noise_d = noise(d_d, nsr);
u_noise = noise(u_raw, nsr);

%% Estimate models
i

% with noise, damaged
[A, B, C, D] = n4sid_(d_noise_d, u_noise, modelorder);
sys_d.A = A;
sys_d.B = B;
sys_d.C = C;
sys_d.D = D;

% with noise, undamaged
clear A B C D
[A, B, C, D] = n4sid_(d_noise, u_noise, modelorder);
sys_u.A = A;
sys_u.B = B;
sys_u.C = C;
sys_u.D = D;

model.sys_u = sys_u;
model.sys_d = sys_d;

% [A, B, C, D] = n4sid_(d, u_raw, modelorder);
% sys_ref.A = A;
% sys_ref.B = B;
% sys_ref.C = C;
% sys_ref.D = D;

% model.sys_ref = sys_ref;
% model.sys_ex = sys_ex_u;
%% Save results

% resultdir = "Main/Chain/damage_patterns/statespace/dt_0_01/no_noise/" + sprintf("/run_%03d",[i]);

resultfile = resultdir + sprintf("/run_%03d",[i]);
save(resultfile, "model");

end

toc
beep
