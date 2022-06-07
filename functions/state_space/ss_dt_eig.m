function [Phi, Lambda] = ss_dt_eig(A, C, dt)
% [Phi, Lambda] = SS_DT_EIG(A, C, dt)
% A: DT state matrix
% C: output matrix
% dt: time increment
% Lambda: poles, ordered by eigenfrequency
% Phi: eigenvectors, ordered by eigenfrequency

n = size(A,1);                  % system order
dof = floor(n/2);
m = size(C,1);                  % number of outputs

[Psi, Alpha] = eig(A);
% sort eigenvalues
l = log(diag(Alpha))/dt;
[~, ind] = sortrows(abs(l));
Lambda = l(ind);

% sort eigenvectors by magnitude of associated eigenvalue
Psi = Psi(:, ind);

% remove complex conjugates
Psi = Psi(:, 1:2:end);
Phi = C*Psi;
% ensure real-valued displacement partition
for i = 1:size(Phi,1)
%     if all(imag(Phi(i,:)) ~= 0) == 1
        if 0 < imag(Phi(i,1)) < 1e-10
            scale = complex( - real(Phi(i, :)) ./ imag(Phi(i, :)), 1);
            break
        else
            scale = complex( - real(Phi(1, :)) ./ imag(Phi(1, :)), 1);
        end
end

Phi = Phi*diag(scale);

% remove imaginary parts that result for computational error
Phi(1:m,:) = real(Phi(1:m,:));

% normalise displacement partition by largest displacement
Phi = Phi./max(abs(Phi(1:m,:)));

% ensure largest displacement is positive
for i = 1:size(Psi,2)
    if abs(max(real(Phi(1:m,i)))) < abs(min(real(Phi(1:m,i))))
        Phi(:,i) = -Phi(:,i);
    end
end