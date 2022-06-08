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
        function dt_from_FE(self, Kg, Cg, Mg, dt, sensortype)

            if ~exist("sensortype", "var")
                disp("StateSpaceModel: No sensor type given. Outputting displacements.")
                sensortype = 'dis';
            end
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
            
            % acceleration output:
            if sensortype == "vel"
                out_dof = out_dof + n;
            elseif sensortype == "acc"
                out_dof = out_dof + 2*n;
            elseif sensortype == "dis"
                % do nothing
            else
                error("StateSpaceModel: No valid sensor type given")
            end
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
        function [A, B, C, D] = estimate(self, u, y, blockrows, order)
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
            [A, B, C, D] = n4sid_(y, u, blockrows, order);

            self.A = A;
            self.B = B;
            self.C = C;
            self.D = D;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [u, y] = time_response(self, u, t, nsr, store)
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
            
            if nsr ~= 0
                y = noise(y, nsr);
                u = noise(u, nsr);
            end

            self.sysid.nsr = nsr;
            self.sysid.N_s = numel(t);
            
            % store time response
            try
                if store
                    self.sysid.y = y;
                    self.sysid.u = u;
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
            catch
                try
                   H = self.C*((eye(size(self.A))*s-self.A)\self.B);
                catch
                    throw("Error (StateSpaceModel.transfer_matrix(): something went wrong")
                end
            end
            self.H = H;
            self.s = s;

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

        function expand(self)
            in_dof = self.in_dof;
            out_dof = self.out_dof;

            [~,v1,v3] = intersect(out_dof,in_dof);
            unsorted_dof = [out_dof, setdiff(in_dof,out_dof)];
            [~,order] = sort(unsorted_dof);

            [Bexp,Cexp]=ExpandBandC2(self.A,self.B,self.C,v1,v3,order);
            self.B = Bexp;
            self.C = Cexp;
        end

    end
end