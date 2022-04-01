function [Bexp,Cexp]=ExpandBandC(Ac,Bc,Cc,out,in,ord)
% [Bexp, Cexp] = ExpandBandC(Ac, Bc, Cc, out, in, ord)
% Expands the input and output matrices.
% Arguments:
% Ac: continuous-time state matrix
% Bc: continuous-time input matrix
% Cc: output matrix
% in: input DOF
% out: output DOF
% ord: ordering of DOF in the resulting B and C

% Get rows in C and column in B that are associated with DOF where both an
% input and output are located
[~, v1, v3] = intersect(out, in);

N = size(Ac, 1);        % system order
m = size(Cc, 1);        % number of outputs
r = size(Bc, 2);        % number of inputs

m = size(Cc,1);
r = size(Bc,2);


% [Bexp, Cexp] = EXPANDBANDC(Ac, Bc, Cc, v1, v3, ord)
% Expands the input and output matrices.
% Arguments:
% Ac: continuous-time state matrix
% Bc: continuous-time input matrix
% Cc: output matrix
% v1: Rows in C of DOF that have an input
% v3: Rows in B of DOF that have an output
% ord:


N = size(Ac, 1);        % system order
m = size(Cc, 1);        % number of outputs
r = size(Bc, 2);        % number of inputs

t = m + r - numel(v3);       % ?

if  isempty(ord) == 0 && numel(ord) - t ~= 0
    disp('ERROR: Invalid number of DOF in ord.')
    return
end

[Phi, Lambda] = eig(Ac);              % eigenvectors and eigenvalues
Lambda = diag(Lambda);

% sort by eigenfrequency in ascending order
[~, ind] = sort(abs(Lambda),'ascend');   
Lambda = Lambda(ind);
Phi = Phi(:, ind);


% for i = 1:size(Phi,1)
%         if imag(Phi(i,1)) < 1e-10
%             scale = complex( - real(Phi(i, :)) ./ imag(Phi(i, :)), 1);
%         break
%     end
% end
% 
% Phi = Phi*diag(scale);
% Phi(1:N/2,:) = real(Phi(1:N/2,:));


% *********************************************************************** %
Psi = (Cc*Phi);
Gamma = (Phi\Bc);

% Select the colocated Gamma
GammaC = Gamma(:, v3);
% GammaC = Gamma(v3, :);
% Select the colocated Psi
PsiC = Psi(v1,:);
% PsiC = Psi(:,v1);
% Compute the scaling constants (LS if there is more than one
% collocation)

for j = 1:N
    
    a1 = PsiC(: , j) * PsiC(: , j).';
    a2 = diag(a1);
    
    b1 = PsiC(: , j) * GammaC(j , :);
    b2 = diag(b1);
    
    alpha(j) = mean(b2 ./ sqrt(a2));
    
%     alpha(j) = PsiC(:,j)*GammaC(j,:)/((PsiC(:,j)).'*PsiC(:,j));
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