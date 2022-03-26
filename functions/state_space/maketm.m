function G = maketm(A, B, C, D, s)
% G = MAKETM(A, B, C, D, s)
% Creates transfer matrix.
%
% Arguments:
% A,B,C,D: matrices of continous-time state space model
% dt: time increment
% s: pole where G is evaluated. Can be symbolic but may take longer.

% [Phi, Lambda] = eig(Ac);
% 
% % Sort by eigenfrequency in ascending order
% [~, ind] = sortrows(abs(diag(Lambda)));
% L = diag(Lambda);
% Lambda = diag(L(ind));
% Phi = Phi(:,ind);

% transfer matrix
G = C*((eye(size(A))*s-A)\B)+D;