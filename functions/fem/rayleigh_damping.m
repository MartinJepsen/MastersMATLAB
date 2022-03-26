function [Cg] = rayleigh_damping(K,M,modes,zetas)
[Phi, Lambda] = eig(K,M);
% Computes damping
zeta1 = zetas(1);
zeta2 = zetas(2);
omega1 = sqrt(Lambda(modes(1),modes(1)));
omega2 = sqrt(Lambda(modes(2),modes(2)));

damp_coef = [1/(2*omega1) omega1/2; 1/(2*omega2) omega2/2]\[zeta1;zeta2];

Cg = damp_coef(1)*M+damp_coef(2)*K;
C_tilde = Phi'*Cg*Phi;

damping_ratios = diag(C_tilde)./(2*sqrt(diag(Lambda)));
criticalmode = find(damping_ratios>=1);

if isempty(criticalmode)==1
disp(['MESSAGE: Damping matrix computed.'...
    ,' There are no critically damped modes.'])
disp(['         zeta_max=',...
    num2str(max(damping_ratios))])
else
disp(['MESSAGE: Damping matrix computed.'...
    ,' Critical damping occurs at mode ', num2str(criticalmode(1)),...
    ', with zeta=',num2str(damping_ratios(criticalmode(1)))])
end
end