function ABCD = make_ss_ct(K, M, C, out_dof, in_dof)
% ABCD = MAKE_SS_CT(K, M, C, out_dof, in_dof)
%  Make struct containing state space matrices

n = size(K,1);          % system order
r = size(in_dof,1);     % number of inputs   
m = size(out_dof,1);    % number of sensors

% Load distribution matrix
B2 = zeros(n,r);
if isempty(in_dof) == 0
    for i = 1:r
        B2(in_dof(i),i)=1;
    end
end

% State space matrices
C_dis = zeros(3*n,n);
C_dis(1:n,1:n) = eye(n);
C_vel = zeros(3*n,n);
C_vel(n+1:2*n,1:end) = eye(n);
C_acc = zeros(3*n,n);
C_acc(2*n+1:end,1:end) = eye(n);

A = [zeros(n) eye(n); -M\K -M\C];
B = [zeros(n,r); M\B2];
C = [C_dis C_vel]+[-C_acc*(M\K) -C_acc*(M\C)];
C = C(out_dof,:);
D = C_acc*(M\B2);

ABCD = struct('A',A,'B',B,'C',C,'D',D);