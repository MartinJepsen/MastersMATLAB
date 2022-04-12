

clear
% 31-01-2022 1438
k = [50, 40, 30, 20, 30, 60];
m = [0.48, 0.27, 0.0700, 0.0700, 0.090, 0.150];

n_dof = numel(m);

dam = 0.85;

save_dir = sprintf("/simulation/system_matrices/%05d", dam*10000)
if exist(fullfile(save_dir), 'dir') ~= 7
    mkdir(save_dir);
    addpath(save_dir);
end

for damel = 1:numel(k)
Kg = zeros(n_dof+1);
Cg = Kg;
Mg = diag(m);
Kg_d = Kg;

for i =1:n_dof
   kl = [1 -1; -1 1]*k(i);

%    ml = eye(2)*m(i);
 
   Kg(i:i+1,i:i+1) = Kg(i:i+1,i:i+1) + kl;
%    Mg(i:i+1,i:i+1) = Mg(i:i+1,i:i+1) + ml;
   if i == damel
       kl = dam*kl;
   end
   Kg_d(i:i+1,i:i+1) = Kg_d(i:i+1,i:i+1) + kl;
end

Kg(1,:)=[]; Kg(:,1)=[];
Kg_d(1,:)=[]; Kg_d(:,1)=[];
% Mg(1,:)=[]; Mg(:,1)=[];

% DeltaK = (Kg_d - Kg) / k(damel)
% Kg_d = Kg + DeltaK;

[Phi, Lambda] = eig(Kg, Mg);
C_tilde = diag([0.01, 0.015, 0.020, 0.025, 0.030, 0.035]) .* (2 * sqrt(diag(Lambda)));
Cg = inv(Phi)' * C_tilde * inv(Phi);

Lambda = ss_ct_eig(Kg, Cg, Mg);


savepath = sprintf(save_dir +"/model_%0d", damel);
save(savepath, 'Kg', 'Kg_d', 'Cg','Mg','damel','dam','Lambda')
end