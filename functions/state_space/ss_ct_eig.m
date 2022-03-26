function [Lambda] = ss_ct_eig(K, C, M)
% [Lambda] = SS_CT_EIG(K, C, M] 
% Computes damped eigenfrequencies for the 
% CT state space model.
n = size(K,1);                                     % number of DOF
A_c = [zeros(n) eye(n); -M\K -M\C];                % CT state matrix
[Lambda] = eig(A_c);                            % poles
[~, ind] = sort(abs(Lambda),'ascend');
Lambda = Lambda(ind);
% omega_d = sortrows(imag(diag(Lambda)),'ascend');   % imaginary parts of poles
% omega_ct = omega_d(n+1:end);
