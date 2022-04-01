function G = maketf_colocated(A, B, C, D, dt, s, in, out, ord)
% G = MAKETF(A, B, C, D, dt)
% Creates transfer matrix for a system whose set of colocated inputs and 
% outputs is not empty.
%
% Arguments:
% A,B,C,D: matrices of discrete-time state space model
% dt: time increment
% s: pole where G is evaluated. Can be symbolic by may take longer.
% in: list of input DOF.
% out: list of output DOF.
% ord: order of DOF in the transfer matrix

Ac = logm(A)/dt;                % CT state matrix
Bc = (Ac\A-inv(Ac))\B;          % CT input matrix

n = size(Ac,1);                 % system order
r = size(Bc, 2);                % number of inputs
m = size(C, 1);                 % number of outputs

% Get eigenvectors and sort by eigenfrequency
[Phi, Lambda] = eig(Ac);
[~, ind] = sortrows(abs(diag(Lambda)));
L = diag(Lambda);
Lambda = diag(L(ind));
Phi = Phi(:,ind);

Psi = C*Phi;
Gamma = Phi\Bc;

% Select the colocated Gamma
GammaC = Gamma(:, v3);
% Select the colocated Psi
PsiC = Psi(v1,:);

% Compute the scaling constants (LS if there is more than one
% collocation)
for j = 1:N
    alpha(j) = PsiC(:,j).'*GammaC(j,:)/((PsiC(:,j)).'*PsiC(:,j));
end

% Obtain entries in the latent vectors at non-collocated inputs
tot = 1:r;
NC = setdiff(tot, v3);  % obtain DOF that are not in the set of inputs (?)

% Determine if there are any non-colocated inputs
if isempty(NC) == 0
    GammaNc = Gamma(:, NC);
    for j = 1:N
        PsiExtra(:,j) = (GammaNc(j, :)).'/alpha(j);
    end
else
    PsiExtra=[];   
end

% Organize the latent vectors such that it is the order of the outputs in
% Cc followed by any coordinates associated with non-collocated inputs (may
% not exist)
Fi = [Psi; PsiExtra];

% Solution
Cexp = Fi*inv(Phi);
Bexp = Phi*diag(alpha)*Phi.'*Cexp.';

% Reorder if desired
if isempty(ord) == 0
    Cexp = Cexp(ord, :);
    Bexp = Bexp(:, ord);
end




% transfer matrix
G = C*((eye(size(Ac))*s-logm(Ac)/dt)\Bc)+D;