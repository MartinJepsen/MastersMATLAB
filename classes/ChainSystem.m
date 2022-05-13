classdef ChainSystem < handle
    properties
        Kg                  % stiffness matrix
        Kg_d                % stiffness matrix (damaged)
        Cg                  % damping matrix
        Cg_d                % damping matrix (damaged)
        Mg                  % mass matrix
        B                   % strain mapping matrix
        n_dof               % number of free DOF + constrained DOF
        damage              % damaged elements
        modal_parameters    % modal parameters
        results             % displacements/strains
        end
    methods
        function self = ChainSystem()
        end

        function assembly(self, k, m, damage)
            % ChainSystem.assembly(k, m)
            % generate stiffness and mass matrices for chain system
            % k: spring stiffnes(es)
            % m: point mass(es)
            % damage: [damaged_element, damaged_stiffness/reference_stiffness]
            
            n_dof = max(numel(k), numel(m));
            n_dof = n_dof(1);
            self.n_dof = n_dof;

            if exist('damage', 'var') 
                dam = [[1:n_dof]', ones(n_dof, 1)];
                dam(damage(:, 1), :) = damage;
            else
                disp('ChainSystem: no damage defined.')
                dam = ones(n_dof, 2);
            end
            self.damage = dam;
            
            % Formulate stiffness matrix
            if numel(k) ~= n_dof && numel(k) == 1
                k = ones(n_dof, 1) * k;
            elseif numel(k) ~= n_dof
                error('numel(k) ~= n_dof or 1')
            end
            Kg = zeros(n_dof + 1);
            Kg_d = Kg;
            for el = 1:n_dof
                k_el = [1, -1; -1,  1] * k(el);

                Kg([el:el + 1], [el:el + 1]) =...
                    Kg([el:el + 1], [el:el + 1]) + k_el;
                Kg_d([el:el + 1], [el:el + 1]) =...
                    Kg_d([el:el + 1], [el:el + 1]) + dam(el, 2) * k_el;
            end
            % apply BC's
            Kg(1, :) = []; Kg(:, 1) = [];
            Kg_d(1, :) = []; Kg_d(:, 1) = [];
            self.Kg = Kg;
            self.Kg_d = Kg_d;
            
            % Formulate mass matrix
            if numel(m) ~= n_dof
                if numel(m) == 1
                    m = ones(n_dof, 1) * m;
                else
                    error('numel(m) ~= n_dof or 1')
                end
            end
            Mg = diag(m);
            self.Mg = Mg;

            B = zeros(n_dof + 1);
            B(1:size(B)+1:numel(B)) = -1;
            B = B - circshift(B,1,2);
            self.B = B;
        end

        function modal_damping(self, zetas)
            % modal_damping(zetas)
            % applies damping ratio zeta to all modes
            if size(zetas, 1) < size(zetas, 2)
                zetas = zetas';
            end

            n_dof = self.n_dof;
            if 1 < numel(zetas) < n_dof
                a = numel(zetas);
                b = n_dof - numel(zetas);
                zetas = [zetas;
                        ones(b, 1) * zetas(end)];
            end
            
            % undamaged system
            Kg = self.Kg;
            Mg = self.Mg;
            [Phi, Lambda] = eig(Kg, Mg);
            omega = sqrt(diag(Lambda));
            % mass-normalise eigenvectors
            m_n = sqrt(diag(Phi' * Mg * Phi)); 
            Phi = Phi * diag(1 ./ m_n);
            
            % get damping matrix
            C_tilde = diag(2 * zetas .* omega);
            self.Cg = inv(Phi)' * C_tilde * inv(Phi);
            zetas = diag(C_tilde) ./ (2 * omega);
            
            % damaged system
            Kg_d = self.Kg_d;
            [Phi_d, Lambda_d] = eig(Kg_d, Mg);
            omega_d = sqrt(diag(Lambda_d));
            m_n = sqrt(diag(Phi_d' * Mg * Phi_d)); % mass-normalised eigenvectors
            Phi_d = Phi_d * diag(1./m_n);
            self.Cg_d = inv(Phi_d)' * C_tilde * inv(Phi_d);

            zetas_d = diag(Phi_d' * self.Cg * Phi) ./ (2 * omega_d);
            self.modal_parameters.zeta = zetas;
            self.modal_parameters.zeta_d = zetas_d;
            self.modal_parameters.Phi = Phi;
            self.modal_parameters.Phi_d = Phi_d;
            self.modal_parameters.Lambda = Lambda;
            self.modal_parameters.Lambda_d = Lambda_d;
            self.modal_parameters.omega = omega;
            self.modal_parameters.omega_d = omega_d;
            self.modal_parameters.omega_hz = omega / (2 * pi);
            self.modal_parameters.omega_hz_d = omega_d / (2 * pi);
        end
    end
end