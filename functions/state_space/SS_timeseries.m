function [d, v, a, ABCD] = SS_timeseries(K,C,M,z0,u,t,in_dof)
% [d, v, a, ABCD] = SS_timeseries(K,C,M,z0,u,t,in_dof)
% Arguments:
%   K,C,M: stiffness, damping and mass matrix
%   u: vector of independent inputs.
%   in_dof: contains the DOF, where the components in u act. If in_dof(1)=3
%           , then the 1st component in u acts in the 3rd dof.
%   out_dof: the rows are [disp; vel; acc], with the columns representing 
%            dof. [0 1; 1 0; 0 1] will give velocity of the first dof and
%            displacement and acceleration for the second dof.

n = size(K,1);
dt = t(2)-t(1);
z = zeros(2*n,numel(t));
z(:,1) = z0;
y = zeros(3*n,numel(t)-1);

% Load distribution matrix
r = size(u,1);
B2 = zeros(n,r);
if isempty(in_dof) == 0
    for i = 1:r
        B2(in_dof(i),i)=1;
    end
end
% State space matrices


% Ouput distribution matrices
C_dis = zeros(3*n,n);
C_dis(1:n,1:n) = eye(n);
C_vel = zeros(3*n,n);
C_vel(n+1:2*n,1:end) = eye(n);
C_acc = zeros(3*n,n);
C_acc(2*n+1:end,1:end) = eye(n);


A_c = [zeros(n) eye(n); -M\K -M\C];
A_d = expm(A_c*dt);


B_c = [zeros(n,r); M\B2];
C_d = [C_dis C_vel]+[-C_acc*(M\K) -C_acc*(M\C)];
D_d = C_acc*(M\B2);
B_d = inv(A_c)*((A_d-eye(size(A_c))))*B_c;
for k = 1:numel(t)
    z(:,k+1) = A_d*z(:,k)+B_d*u(:,k);
    y(:,k) = C_d*z(:,k)+D_d*u(:,k);
end
d = y(1:n,:);
v = y(n+1:2*n,:);
a = y(2*n+1:end,:);

ABCD = struct('A',A_d,'B',B_d,'C',C_d,'D',D_d);
end