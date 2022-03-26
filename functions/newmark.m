function [x, v, a] = newmark(K,C,M,t,F,x0,v0,beta,gamma)
% [x, v, a] = NEWMARK(K,C,M,t,F,x0,v0,beta,gamma,snr) 
%   Computes the displacement (d), velocity (v) and acceleration (a) time 
%   series with white Gaussian noise for a structural system.
% INPUTS:
%   K: n x n stiffness matrix
%   C: n x n damping matrix  (set to "ZEROS(n)" for undamped system
%   M: n x n mass matrix
%   t: 1 x N time vector. must have uniform time step size
%   F: load vector of size n*N. Each row represents a DOF.
%   x0, v0: n x 1 initial displacement and velocity vectors
%   beta: integration parameter. 1/4 for constant acc., 1/6 for linear acc.
%   gamma: integration parameter. 1/2 for no algorithmic damping

%% Stability check

if size(F,1) ~= size(K,2)
    disp('ERROR: Inconsistent load vector dimensions')
    return
end

% Determine largest time increment
dt = t(2)-t(1);                       % time increment
N = numel(t);
% Determine damping ratios
[Phi, Lambda] = eig(K,M);
omega = sqrt(diag(Lambda));
zeta = diag(Phi'*C*Phi)./(2*omega);
Omega_cr = zeros(numel(zeta),1);
for i = 1:numel(zeta)
    Omega_cr(i) = (zeta(i)*(gamma-1/2)+sqrt(gamma/2-beta+zeta(i)^2*(gamma-1/2)^2))/...
                (gamma/2-beta);
end

if 2*beta >= gamma >= 1/2
    disp('MESSAGE: Unconditional stability.')
elseif beta <= 1/2*gamma && gamma >= 1/2 && dt <= min(Omega_cr./omega)
    disp(['MESSAGE: The Newmark algorithm is stable. Largest allowable dt = ',num2str(min(Omega_cr./omega))])    
else
    disp(['WARNING: The Newmark algorithm is unstable. Choose beta = 1/4 or dt <= ',num2str(min(Omega_cr./omega))])
end

% Newmark algorithm

    dof = size(K,1);                        % number of unconstrained DOF

    a0 = M\(F(:,1)-C*v0-K*x0);              % isolating a from eq. of motion
    
    % Helping parameters
    a1 = 1/(beta*dt^2); a2 = 1/(beta*dt); a3 = 1/(2*beta); a4 = 1/beta;
    
    K_eff = a1*M+a2*gamma*C+K;             % effective stiffness

    % Time response arrays
    x = zeros(dof,N); x(:,1) = x0;
    v = zeros(dof,N); v(:,1) = v0;
    a = zeros(dof,N); a(:,1) = a0;
    
for i = 1:N-1    
    % effective kinematic quantities
    a_eff = M*(a1*x(:,i)+a2*v(:,i)+(a3-1)*a(:,i));
    v_eff = C*(gamma*a2*x(:,i)+(gamma*a4-1)*v(:,i)+dt*(gamma*a3-1)*a(:,i));

    x(:,i+1) = K_eff\(F(:,i+1)+a_eff+v_eff);
    v(:,i+1) = gamma*a2*(x(:,i+1)-x(:,i))-(gamma*a4-1)*v(:,i)-dt*(gamma*a3-1)*a(:,i);
    a(:,i+1) = a1*(x(:,i+1)-x(:,i)-dt*v(:,i))-(a3-1)*a(:,i);
end

% Adds noise to the output
%  a = [a(:,1:floor(N/10)), awgn(a(:,ceil(N/10):end),snr,'measured')];
%  v = [v(:,1:floor(N/10)), awgn(v(:,ceil(N/10):end),snr,'measured')];
%  x = [x(:,1:floor(N/10)), awgn(x(:,ceil(N/10):end),snr,'measured')];

%  a = awgn(a,snr,'measured');
%  v = awgn(v,snr,'measured');
%  x = awgn(x,snr,'measured');

disp('MESSAGE: Time integration complete')

end