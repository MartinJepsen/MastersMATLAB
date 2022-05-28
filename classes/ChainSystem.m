classdef ChainSystem < handle
    properties
        Kg                  % stiffness matrix
        Cg                  % damping matrix
        Mg                  % mass matrix
        B                   % strain mapping matrix
        n_dof               % number of free DOF + constrained DOF
        damage              % damaged elements
        modal_parameters    % modal parameters
        element_properties  % element properties
        results             % displacements/strains
        end
    methods
        function self = ChainSystem()
        end

        function assembly(self, k, m)
            % ChainSystem.assembly(k, m)
            % generate stiffness and mass matrices for chain system
            % k: spring stiffnes(es)
            % m: point mass(es)
            
            n_dof = max(numel(k), numel(m));
            n_dof = n_dof(1);
            self.n_dof = n_dof;
            self.element_properties.k = k;
            self.element_properties.m = m;

            
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
            end
            % apply BC's
            Kg(1, :) = []; Kg(:, 1) = [];
            self.Kg = Kg;
            
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

            B = zeros(n_dof, n_dof + 1);
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
            
            self.modal_parameters.zeta = zetas;
            self.modal_parameters.Phi = Phi;
            self.modal_parameters.Lambda = Lambda;
            self.modal_parameters.omega = omega;
            self.modal_parameters.omega_hz = omega / (2 * pi);
        end

        function rayleigh_damping(self, damping_ratios, modes)
            if numel(damping_ratios) ~= 2 || numel(modes) ~= 2
                throw('rayleigh_damping: number of damping ratios or specified modes is not equal to 2')
            end
            
            % if eigenfrequencies already exist, use them, if not, compute them
            if isfield(self.modal_parameters, 'Lambda')
                omega1 = sqrt(self.modal_parameters.Lambda(modes(1),modes(1)));
                omega2 = sqrt(self.modal_parameters.Lambda(modes(2),modes(2)));
            else
                [Phi, Lambda] = eig(self.Kg, self.Mg);
                self.modal_parameters.Phi = Phi;
                self.modal_parameters.Lambda = Lambda;
                omega1 = sqrt(Lambda(modes(1),modes(1)));
                omega2 = sqrt(Lambda(modes(2),modes(2)));
            end
            
            zeta1 = damping_ratios(1);
            zeta2 = damping_ratios(2);
            damp_coef = [1/(2*omega1), omega1/2; 1/(2*omega2), omega2/2] \ [zeta1;zeta2];
            
            fprintf("Generating Rayleigh damping matrix with alpha=%0.5f and beta=%0.5f\n", damp_coef)
            Cg = damp_coef(1) * self.Mg + damp_coef(2) * self.Kg;
            self.Cg = Cg;
            C_tilde = Phi' * Cg * Phi;
            
            zetas = diag(C_tilde) ./ (2 * sqrt(diag(Lambda)));
            self.modal_parameters.zeta = zetas;
            criticalmode = find(zetas >= 1);
            
            if isempty(criticalmode)
                disp(['MESSAGE: Damping matrix computed.'...
                    ,' There are no critically damped modes.'])
                disp(['         zeta_max=',...
                    num2str(max(zetas))])
            else
                disp(['MESSAGE: Damping matrix computed.'...
                    ,' Critical damping occurs at mode ', num2str(criticalmode(1)),...
                    ', with zeta=',num2str(zetas(criticalmode(1)))])
            end
            self.modal_parameters.omega = sqrt(diag(Lambda));
            self.modal_parameters.omega_hz = sqrt(diag(Lambda)) / (2 * pi);
            self.modal_parameters.alpha = damp_coef(1);
            self.modal_parameters.beta = damp_coef(2);
        end

        function get_modal_parameters(self)
            Kg = self.Kg;
            Cg = self.Cg;
            Mg = self.Mg;

            [Phi, Lambda] = eig(Kg, Mg);
            omega = sqrt(diag(Lambda));
            omega_hz = omega / (2*pi);
            % mass-normalise eigenvectors
            m_n = sqrt(diag(Phi' * Mg * Phi)); 
            Phi = Phi * diag(1 ./ m_n);
            
            % get damping matrix
            C_tilde = Phi' * Cg * Phi;
            zetas = diag(C_tilde) ./ (2 * omega);

            self.modal_parameters.zeta = zetas;
            self.modal_parameters.Phi = Phi;
            self.modal_parameters.Lambda = Lambda;
            self.modal_parameters.omega = omega;
            self.modal_parameters.omega_hz = omega_hz;
            
        end
    end
end