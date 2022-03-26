classdef StateSpaceModel < handle
    properties 
        A       % DT state matrix
        B       % DT input matrix
        C       % output matrix
        D       % feedthrough matrix
        B2     % input distribution matrix
        H       % transfer matrix
        s       % Laplace variable
        in_dof  % set of input DOF
        out_dof % set of output DOF
        sysid   % SYSID and simulation parameters
        modal_parameters 
    end

    methods
        function self = StateSpaceModel()
            % constructor
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function B2 = set_io(self, in_dof, out_dof, order)
            self.in_dof = in_dof;
            self.out_dof = out_dof;

            if exist('order', 'var') && isempty(self.B2)
                B2 = zeros(order/2, numel(in_dof));
                if isempty(in_dof) == 0
                    for i = 1:numel(in_dof)
                        B2(in_dof(i), i)=1;
                    end
                end
                self.B2 = B2;
            end
        end        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dt_from_FE(self, Kg, Cg, Mg, dt)
            % FROM_FE(K, C, M)
            % Make DT state space model from FE system matrices
            self.A = [];
            self.B = [];
            self.C = [];
            self.D = [];

            in_dof = self.in_dof;
            out_dof = self.out_dof;
            B2 = self.B2;

            n = size(Kg,1);         % system order
            r = numel(in_dof);     % number of inputs   
            m = numel(out_dof);    % number of sensors
            
            if isempty(B2)
                B2 = zeros(n,r);
                if ~isempty(in_dof)
                    for i = 1:r
                        B2(in_dof(i),i)=1;
                    end
                else
                    disp('ERROR (dt_from_FE): B2 matrix due to undefined input/output configuration')
                    return
                end
            end

            % State space matrices
            C_dis = zeros(3*n,n);
            C_dis(1:n,1:n) = eye(n);
            C_vel = zeros(3*n,n);
            C_vel(n+1:2*n,1:end) = eye(n);
            C_acc = zeros(3*n,n);
            C_acc(2*n+1:end,1:end) = eye(n);
            
            A_c = [zeros(n) eye(n); -Mg\Kg -Mg\Cg];
            A = expm(A_c * dt);
            B_c = [zeros(n,r); Mg\B2];
            B = inv(A_c) * ((A - eye(size(A_c)))) * B_c;
            C = [C_dis C_vel]+[-C_acc*(Mg\Kg) -C_acc*(Mg\Cg)];
            D = C_acc*(Mg\B2);

            C = C(out_dof, :);
            D = D(out_dof, :);

            self.A = A;
            self.B = B;
            self.B2 = B2;
            self.C = C;
            self.D = D;
            self.sysid.dt = dt;

        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function from_SS(self, SS)
            % inherit state space matrices from other state space model
            self.A = SS.A;
            self.B = SS.B;
            self.C = SS.C;
            self.D = SS.D;
            self.B2 = SS.B2;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [A, B, C, D] = estimate(self, y, u, blockrows)
            self.A = [];
            self.B = [];
            self.C = [];
            self.D = [];

            if isempty(y) && isempty(u)
                try
                    y = self.sysid.y;
                    u = self.sysid.u;
                catch
                    disp('ERROR: (estimate) no inputs or loads were found. Returning')
                    return
                end
            end
            [A, B, C, D] = n4sid_(y, u, blockrows);

            self.A = A;
            self.B = B;
            self.C = C;
            self.D = D;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [y, u] = time_response(self, u, t, nsr, store)
            A = self.A;
            B = self.B;
            C = self.C;
            D = self.D;
            dt = t(2) - t(1);
            self.sysid.dt = dt;

            n = size(A,1);
            z = zeros(n, numel(t));
            z(:, 1) = zeros(n, 1);
            y = zeros(size(C,1),numel(t)-1);
            
            for k = 1:numel(t)
                z(:, k+1) = A*z(:, k) + B*u(:, k);
                y(:, k) = C*z(:, k) + D*u(:, k);
            end
            
            y_n = noise(y, nsr);
            u_n = noise(u, nsr);
            self.sysid.nsr = nsr;
            self.sysid.N_s = numel(t);
            
            % store time response
            try
                if store
                    self.sysid.y = y_n;
                    self.sysid.u = u_n;
                end
            catch
            end
             
        end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function H = transfer_matrix(self, s)
            % Create transfer matrix from A, B, C, D, s

            % Get s-value from self if not supplied
            if ~exist('s', 'var') && ~isempty(self.s)
                s = self.s;
            elseif ~exist('s', 'var') && isempty(self.s)
                disp('WARNING (transfer_matrix): No s-value could be found. s = 0.')
                s = 0;
            end
            
            try
                H = self.C*((eye(size(self.A))*s-self.A)\self.B)+self.D;
                self.H = H;
                self.s = s;
            catch
                H = [];
                disp('ERROR (transfer_matrix): Missing parameter (A,B,C,D,s).')
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function expand_coords(self, ord)
            % expand_coords(ord)
            % Expands the input and output matrices.
            % ord: ordering of DOF in the resulting B and C
            
            if ~exist('ord','var')
                ord = [];
            end
            
            % Get rows in C and column in B that are associated with DOF where both an
            % input and output are located
            
            % Get sensor/input configuration
            in = self.in_dof;
            out = self.out_dof;
            [~, v1, v3] = intersect(out, in);
            
            if isempty(v1)
                disp('ERROR (expand_coords): No co-located inputs and outputs.')
                return
            end
            
            % v1 contains indices of rows in C for coordinates where an
            % input acts
            % v3 contains indices of columns in B for coordinates where a
            % sensor is located
            
            % Get state space matrices and convert to CT
            A_d = self.A;            
            A_c = logm(A_d) / self.sysid.dt;
            B_d = self.B;
            B_c = (A_c\A_d - inv(A_c)) \ B_d;
            C = self.C;
            
            N = size(A_d, 1);       % system order
            r = size(B_c, 2);       % number of inputs
            m = size(C, 1);         % number of outputs
           
            t = m + r - numel(v3);       % ?

            if  ~isempty(ord) && numel(ord) - t ~= 0
                disp('ERROR: Invalid number of DOF in ord.')
                return
            end
            
            [Phi, Lambda] = eig(A_c);
            Lambda = diag(Lambda);
            [~, ind] = sort(abs(Lambda), 'ascend');   
            Lambda = Lambda(ind);
            Phi = Phi(:, ind);

            Psi = C * Phi;
            Gamma = Phi \ B_c;
            % Select the colocated Gamma
            GammaC = Gamma(:, v3);
            % Select the colocated Psi
            PsiC = Psi(v1, :);
            
            % Compute the scaling constants (LS if there is more than one
            % collocation)
            for j = 1:N
                a1 = PsiC(: , j) * PsiC(: , j).';
                a2 = diag(a1);
                b1 = PsiC(: , j) * GammaC(j , :);
                b2 = diag(b1);
                alpha(j) = mean(b2 ./ sqrt(a2));
% %                 alpha(j)=fiC(:,j).'*TauC(j,:)/((fiC(:,j)).'*fiC(:,j));
%                 alpha(j) = PsiC(:, j).' * GammaC(j, :) / ((PsiC(:, j)).' * PsiC(:, j));
            end


            % Obtain entries in the latent vectors at non-collocated inputs
            tot = 1:r;
            NC = setdiff(tot, v3);  % obtain DOF that are not in the set of inputs

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
            Cexp = Fi * inv(Phi);
            Bexp = Phi * diag(alpha)*Phi.'*Cexp.';
            
            % Reorder if desired
            if ~isempty(ord)
                Cexp = Cexp(ord, :);
                Bexp = Bexp(:, ord); 
            end
            
            % Update class attributes
%             self.B = inv(A_c) * ((self.A - eye(size(A_c)))) * Bexp;
            B_c = A_c \ ((A_d - eye(size(A_c)))) * Bexp;
            
            Cexp = round(Cexp, 15);
            B_c = round(B_c, 15);
            self.B = real(B_c);
            self.C = real(Cexp);
            self.set_io(union(in, out), union(in, out));
            
            % TODO: implement update of D-matrix
            self.D = 0;
%             C_acc = zeros(3*N, N);
%             C_acc(2*N + 1:end,1:end) = eye(N);
%             D = C_acc * (self.Mg \ self.B2);
%             self.D = D(self.out_dof, :);
            self.transfer_matrix(self.s);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function get_modal_parameters(self, dt)
            
            if ~exist('dt', 'var')
                try
                    dt = self.sysid.dt;
                catch
                    disp('WARNING (modal_parameters): No dt is defined')
                end
            end

            [Phi, Lambda] = ss_dt_eig(self.A, self.C, dt);
            self.modal_parameters.Phi = Phi;
            self.modal_parameters.Lambda = Lambda;
            self.modal_parameters.omega = abs(Lambda);
            self.modal_parameters.zeta = - real(Lambda) ./ abs(Lambda);
        end

        function to_cl(self, K)

            self.A = self.A - self.B * K * self.C;
            try
                H = self.H;
                self.H = (eye(size(H)) + H * K)^-1*H;
            catch
                self.H = [];
            end
        end

        function to_ct(self)
            A_d = self.A;
            A_c = logm(A_d) / self.sysid.dt;
            B_d = self.B;
            B_c = (A_c\A_d-inv(A_c))\B_d;
            self.A = A_c;
            self.B = B_c;
        end

    end
end