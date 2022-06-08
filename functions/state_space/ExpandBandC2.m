function [Bexp,Cexp]=ExpandBandC2(Ac,Bc,Cc,v1,v3,ord)

% v1: index of rows in C associated with colocated DOF
% v3: index of columns in B associated with colocated DOF

[mu,~]=size(Ac);
[m,~]=size(Cc);
[~,r]=size(Bc);
t=m+r-length(v3);

xx=isempty(ord);
if  xx~=1
    v4=length(ord);
    tes=v4-t;
       if tes~=0
       displ('error in specifying ord')
       return
       else
       end
else
end
[Phi,Lambda]=eig(Ac);

for i = 1:size(Phi,1)
%     if all(imag(Phi(i,:)) ~= 0) == 1
    if ~any(imag(Phi(i,:))==0)
        scale = complex( - real(Phi(i, :)) ./ imag(Phi(i, :)), 1);
        break
    else
        scale = complex( - real(Phi(1, :)) ./ imag(Phi(1, :)), 1);
    end
end

Phi = Phi * diag(scale);

Lambda=diag(Lambda);
[~,Y]=sort(abs(Lambda),'ascend');
Lambda=Lambda(Y);

Phi=Phi(:,Y);
Psi=Cc*Phi;
Gam=inv(Phi)*Bc;
% Select the collocated Gam
GamC=Gam(:,v3);
% Select the collocated Psi
PsiC=Psi(v1,:);
% Compute the scaling constants (least squares if there is more than one
% collocation)

for i=1:size(GamC,1)
    for j = 1:size(GamC,2)
%     alpha(j)=PsiC(:,j).'*GamC(j,:)/((PsiC(:,j)).'*PsiC(:,j));
        alpha(i,j) = GamC(i,j)/GamC(j,i);
    end
end

% Obtain entries in the latent vectors at non-collocated inputs
tot=1:r;
NC=setdiff(tot,v3);
% Determine if there are any non-collocated inputs

if ~isempty(NC)
    GamNC=Gam(:,NC);
    for j=1:mu
        d=(GamNC(j,:)).'/alpha(j);
        PsiExtra(:,j)=d;
    end
else
    PsiExtra=[];   
end

% Organize the latent vectors such that it is the order of the outputs in
% Cc followed by any coordinates associated with non-collocated inputs (may
% not exist)
Fi=[Psi;PsiExtra];
% Solution
Cexp=Fi*inv(Phi);
Bexp=Phi*diag(alpha)*Phi.'*Cexp.';
% Reorder if desired
if xx~=1
    Cexp=Cexp(ord,:);
    Bexp=Bexp(:,ord);
else

Bexp = real(Bexp);
Cexp = real(Cexp);
end
