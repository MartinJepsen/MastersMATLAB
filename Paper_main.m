close all
clear
tic

%% Setup
set_up

%% FE model
FE = FiniteElementModel('structures/paper_truss.xlsx');
FE.assembly('bar', dam);
FE.apply_bc([1, 2, 9, 10]);
FE.modal_damping(0.02);
FE.strains_from_disp([])
cdis = zeros(numel(out_dof), 12);
for ii=1:numel(out_dof);
    cdis(ii, out_dof(ii)) = 1;
end

%% Exact SS model
% Exact, reference
SS_exact = StateSpaceModel();
SS_exact.set_io(1:12, 1:12, 24);
SS_exact.dt_from_FE(FE.Kg, FE.Cg, FE.Mg, dt);
SS_exact.get_modal_parameters();
% SS_exact.time_response(u, t, 0, true);
SS_exact.to_ct();
Lambda = SS_exact.modal_parameters.Lambda;
s = complex(real(Lambda(1)), 1.1*imag(Lambda(1)));         % pole
SS_exact.transfer_matrix(s);

% Exact, damaged
SS_exact_d = StateSpaceModel();
SS_exact_d.set_io(1:12, 1:12, 24);
SS_exact_d.dt_from_FE(FE.Kg, FE.Cg, FE.Mg, dt);
SS_exact_d.to_ct();
SS_exact_d.transfer_matrix(s);

%% Estimated SS models
SS = StateSpaceModel();
B2 = SS.set_io(in_dof, out_dof, 24);
SS.dt_from_FE(FE.Kg, FE.Cg, FE.Mg, dt);
SS.time_response(u, t, nsr, true);
SS.estimate([],[], blockrows);
SS.get_modal_parameters();
SS.to_ct();
SS.transfer_matrix(s);

export_gain_pars
return

% omega_est = SS.modal_parameters.omega;
% omega_ref = SS_exact.modal_parameters.omega;
% zeta_est = SS.modal_parameters.zeta;
% zeta_ref = SS_exact.modal_parameters.zeta;
% [dev(omega_ref, omega_est), dev(zeta_ref, zeta_est)]

% Damaged
SS_d = StateSpaceModel();
SS_d.set_io(in_dof, out_dof, 24);
SS_d.dt_from_FE(FE.Kg_d, FE.Cg, FE.Mg, dt);
SS_d.time_response(u, t, nsr, true);
SS_d.estimate([],[], blockrows);
SS_d.get_modal_parameters();
SS_d.to_ct();
SS_d.transfer_matrix(s);

H = SS.H;
H_d = SS_d.H;

%% DDLV
gainpath = "optimised_DDLV\01_strain_cond\gains_1.mat";
load(gainpath)
H_ref = (FE.Mg * s^2 + FE.Cg * s + FE.Kg)^-1;
H_ref_CL = (FE.Mg * s^2 + FE.Cg * s + FE.Kg + B2 * K * cdis)^-1;

DeltaH = H_d - H;
[~, ~, V] = svd(DeltaH);

d = SS_exact.H * SS.B2 * V(:, end);
FE.strains_from_disp(d);
eps_OL = FE.results.eps;

eps_OL = abs(eps_OL) / max(abs(eps_OL));

% CL
SS.to_cl(K);
SS_d.to_cl(K);
H_CL = SS.H;
H_CL_d = SS_d.H;

DeltaH_CL = SS_d.H - SS.H;
[~, ~, V] = svd(DeltaH_CL);

% SS_exact.to_cl(K);
H_CL_full = SS_exact.H;
d_CL = H_ref_CL * SS.B2 * V(:, end);   % 
FE.strains_from_disp(d_CL);
eps_CL = FE.results.eps;

eps_CL = abs(eps_CL) / max(abs(eps_CL));

figure
hold on
x = [1:14];  % positions of the bars
b1 = bar(x, eps_OL, 'k');
% b2 = bar(x+x(end), eps_CL, 'w');
b2 = bar(x, eps_CL, 0.6, 'w');
set(gca, 'YScale', 'log')
ylim([1e-6, 2])
legend('OL', 'CL')
toc


