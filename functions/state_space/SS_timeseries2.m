function y = SS_timeseries2(A,B,C,D,z0,u,t)
%  y = SS_timeseries2(A,B,C,D,z0,u,t)
% Arguments:
%   K,C,M: stiffness, damping and mass matrix
%   u: vector of independent inputs.
%   in_dof: contains the DOF, where the components in u act. If in_dof(1)=3
%           , then the 1st component in u acts in the 3rd dof.
%   out_dof: the rows are [disp; vel; acc], with the columns representing 
%            dof. [0 1; 1 0; 0 1] will give velocity of the first dof and
%            displacement and acceleration for the second dof.

n = size(A,1);
z = zeros(n,numel(t));
z(:,1) = z0;
y = zeros(size(C,1),numel(t)-1);


for k = 1:numel(t)
    z(:,k+1) = A*z(:,k)+B*u(:,k);
    y(:,k) = C*z(:,k)+D*u(:,k);
end
% d = y(1:n,:);
% v = y(n+1:2*n,:);
% a = y(2*n+1:end,:);

end