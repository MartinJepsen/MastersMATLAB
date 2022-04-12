

k = [50, 40, 30, 20, 30, 60];
m = [0.48, 0.27, 0.0700, 0.0700, 0.090, 0.150];
free_dof = numel(m);

Kg = zeros(free_dof+1);
Cg = Kg;
Mg = diag(m);
Kg_d = Kg;

for i =1:free_dof
   kl = [1 -1; -1 1]*k(i);
   Kg(i:i+1,i:i+1) = Kg(i:i+1,i:i+1) + kl;

   if i == dam(1)
       kl = dam(2)*kl;
   end
   Kg_d(i:i+1,i:i+1) = Kg_d(i:i+1,i:i+1) + kl;
end

Kg(1,:)=[]; Kg(:,1)=[];
Kg_d(1,:)=[]; Kg_d(:,1)=[];

[Phi, Lambda] = eig(Kg, Mg);
C_tilde = diag([0.01, 0.015, 0.020, 0.025, 0.030, 0.035]) .* (2 * sqrt(diag(Lambda)));
Cg = inv(Phi)' * C_tilde * inv(Phi);


FE = FiniteElementModel();
FE.Kg = Kg;
FE.Kg_d = Kg_d;
FE.Cg = Cg;
FE.Mg = Mg;
FE.n_dof = free_dof + 1;
FE.mesh.topology = [1, 2;
                    2, 3;
                    3, 4;
                    4, 5;
                    5, 6;
                    6, 7];
FE.mesh.coordinates = [0, 0;
                       1, 0;
                       2, 0;
                       3, 0;
                       4, 0;
                       5, 0;
                       6, 0];
n_el = size(FE.mesh.topology, 1);
FE.mesh.n_el = n_el;
FE.mesh.n_node = size(FE.mesh.coordinates, 1);
FE.mesh.element_properties.L = ones(n_el,1);
FE.mesh.element_properties.rot = [ones(n_el, 1), zeros(n_el, 1)];

B = zeros(free_dof, FE.n_dof);
B(1:size(B)+1:numel(B)) = -1;
B = B - circshift(B,1,2);
FE.B = B;
FE.mesh.bc_dof = 1;

clearvars -except FE free_dof dam

