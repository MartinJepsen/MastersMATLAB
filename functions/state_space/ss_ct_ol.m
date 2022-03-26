function ABCD = ss_ct_ol(K, C, M, out_dof, in_dof, type)
% ABCD = 	SS_CT_OL(K, C, M, out_dof, in_dof, type)
%  Make struct containing state space matrices and transfer matrix
%  K, C, M: system matrices
%  out_dof, in_dof: vector containing numbers of output and input dof
% type: 'dis', 'vel' or 'acc'. Selects the sensor type

n = size(K,1);         % system order
r = numel(in_dof);     % number of inputs   
m = numel(out_dof);    % number of sensors


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
D = C_acc*(M\B2);

%% Rewrite C and D to account for sensor locations
if type == 'vel'
    out_dof = out_dof+n
elseif type == 'acc'
    out_dof = out_dof+2*n
elseif type ~= 'dis'
    disp('ERROR: Invalid sensor type.')
    return
end

C = C(out_dof,:);
D = D(out_dof,:);

ABCD = struct('A',A,'B',B,'C',C,'D',D);