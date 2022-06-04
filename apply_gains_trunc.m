clear; close all

%% Set simulation variables
nsr = 0.02;
err = 0.02;
dam_ = 0.5;
sensor = "dis";
poles = 1;
elements = 1:14;
mode = 6;

show_plots = false;

%% Compute results
if mode ~= 0
    base_dir = sprintf("simulation/SYSID/t%d_model_error_%03d_%s", mode, err*100, sensor);
elseif mode == 0
    base_dir = sprintf("simulation/SYSID/model_error_%03d_%s", err*100, sensor);
end

load(sprintf("%s/00_000_%03d", base_dir, round(nsr*100,0)))
load(fullfile(base_dir, "SetUp.mat"))

Kg = ReferenceModels.Kg;
Cg = ReferenceModels.Cg;
Mg = ReferenceModels.Mg;
B2 = GeneralParameters.B2;
cdis = GeneralParameters.cdis;
B_strain = GeneralParameters.B_strain;
n_dof = GeneralParameters.n_dof;
idx = GeneralParameters.idx;
SS_exact = ReferenceModels.SS_exact;

s_vals = [];

for damel = elements
    tot_runs = 1;

    dam = [damel, dam_];

    % load simulation results
    filename = sprintf('%02d_%03d_%03d', dam(1,1), dam(1,2)*100, round(nsr*100,0));
    disp("Loading " + filename)
    load(fullfile(base_dir, filename))

    % Compute all estimated transfer matrices
    n_sim = numel(SS_est);
    n_sim_d = numel(SS_est_d);
    n_runs = n_sim * n_sim_d * numel(poles);



    strains = zeros(size(B_strain, 1), n_runs, 2);
    min_strain_OL = zeros(n_runs, 1);
    min_strain_CL = zeros(n_runs, 1);

    for pole = poles
        
        % load gainss
        %     load('gaindesign/01_strain_cond/gains_5_0.120.mat')

            load(sprintf("gaindesign/01_strain_cond/gains_%d_1.120.mat", pole))
        %     load(sprintf("gaindesign/02_sens/constrained/gains_%02d", damel))
        %     load(sprintf("gaindesign/03_strain_norm/gains_%02d", damel))
        %     load(sprintf("Ks_%03d_%03d_%03d_%s", err*100, dam_*100, nsr*100, sensor))
        s_vals((pole+1)/2) = s;
        
        % account for output type
        if sensor == "dis"
            s_fac = 1;
        elseif sensor == "vel"
            s_fac = 1/s;
        elseif sensor == "acc"
            s_fac = 1 / (s^2);
        end
        
        % create new transfer matrices for each s-value
        H_arr = cell(n_sim, 1);
        H_CL_arr = cell(n_sim, 1);
        H_d_arr = cell(n_sim_d, 1);
        H_CL_d_arr = cell(n_sim_d, 1);
        for run_u = 1:n_sim
            H = s_fac * SS_est{run_u}.transfer_matrix(s);

            A = SS_est_d{run_u}.A;
            B = SS_est_d{run_u}.B;
            C = SS_est_d{run_u}.C;
% 
            [PsiUS,LamUS]=eig(A);
            lamUS=diag(LamUS);
            [~,D1]=sort(abs(imag(lamUS)));
            lam=lamUS(D1);
            Psi=PsiUS(:,D1);
            Phi=inv(Psi);
            Lam=diag(lam);
            
            G=C*(eye(size(A))*s-A)^-1*B;
            Gm=C*Psi*(eye(size(A))*s-Lam)^-1*Phi*B;
            
            Z=1:14;
            H=C*Psi(:,Z)*(eye(numel(Z))*s-Lam(Z,Z))^-1*Phi(Z,:)*B;

            H_arr{run_u, 1} = H;
            H_CL_arr{run_u, 1} = (eye(size(H)) + H * K)^-1 * H;
        end
        for run_d = 1:n_sim_d
%             H_d = s_fac * SS_est_d{run_d}.transfer_matrix(s);
            A = SS_est_d{run_d}.A;
            B = SS_est_d{run_d}.B;
            C = SS_est_d{run_d}.C;
% 
            [PsiUS,LamUS]=eig(A);
            lamUS=diag(LamUS);
            [~,D1]=sort(abs(imag(lamUS)));
            lam=lamUS(D1);
            Psi=PsiUS(:,D1);
            Phi=inv(Psi);
            Lam=diag(lam);
            
            G=C*(eye(size(A))*s-A)^-1*B;
            Gm=C*Psi*(eye(size(A))*s-Lam)^-1*Phi*B;
            
            H_d=C*Psi(:,Z)*(eye(numel(Z))*s-Lam(Z,Z))^-1*Phi(Z,:)*B;

            H_d_arr{run_d, 1} = H_d;
            H_CL_d_arr{run_d, 1} = (eye(size(H_d)) + H_d * K)^-1 * H_d;   % estimated CL transfer matrix, damaged
        end

        A_CL_ex = SS_exact.A + SS_exact.B * B2 * K * cdis * SS_exact.C;
        Lambda_CL = eig(A_CL_ex);                       % exact CL poles
    
        % model transfer matrices
        H_ref = (Mg*s^2 + Cg*s + Kg)^-1;                % reference OL transfer matrix
        SS = ReferenceModels.SS_exact;
        A = SS.A;
        B = SS.B;
        C = SS.C;   
        [PsiUS,LamUS]=eig(A);
        lamUS=diag(LamUS);
        [~,D1]=sort(abs(imag(lamUS)));
        lam=lamUS(D1);
        Psi=PsiUS(:,D1);
        Phi=inv(Psi);
        Lam=diag(lam);
        
        G=C*(eye(size(A))*s-A)^-1*B;
        Gm=C*Psi*(eye(size(A))*s-Lam)^-1*Phi*B;
        
        H_ref=C*Psi(:,Z)*(eye(numel(Z))*s-Lam(Z,Z))^-1*Phi(Z,:)*B;
        G = C*((eye(size(A))*s-A)\B);

%         H_CL_ref = (Mg*s^2 + Cg*s + Kg + B2*K*cdis)^-1; % reference CL transfer matrix
        H_CL_ref = (eye(size(H_ref)) + H_ref * B2*K*cdis)^-1 * H_ref

        for run_u = 1:n_sim
            H = H_arr{run_u};
            H_CL = H_CL_arr{run_u};
            
            for run_d = 1:n_sim_d
                H_d = H_d_arr{run_d};
                H_CL_d = H_CL_d_arr{run_d};
    
                DeltaH = H_d - H;                               % damage-induced transfer matrix shift (estimated)
                [~, ~, V] = svd(DeltaH);                        % DDLVs
                d_OL = zeros(n_dof, 1);
                d_OL(idx) = H_ref * B2 * V(:, end);             % full OL displacement vector
                eps_OL = B_strain * d_OL;                       % full OL strain vector
            
                DeltaH_CL = H_CL_d - H_CL;                      % CL damage-induced transfer matrix change
                [~, ~, V] = svd(DeltaH_CL);                     % CLDDLVs
                d_CL = zeros(n_dof, 1);
                d_CL(idx) = H_CL_ref * B2 * V(:, end);          % CL displacement vector
                eps_CL = B_strain * d_CL;                       % cCL strain vector
            
                % Compute strains
                strains(:, tot_runs, 1) = abs(eps_OL);        % array of characteristic OL strain vectors
                strains(:, tot_runs, 2) = abs(eps_CL);        % array of characteristic CL strain vectors
    
                eps_norm_OL = abs(eps_OL) / max(abs(eps_OL));
                eps_norm_CL = abs(eps_CL) / max(abs(eps_CL));
                
                min_strain_OL(tot_runs, 1) = find(eps_norm_OL == min(eps_norm_OL));   % index of smallest OL strain
                min_strain_CL(tot_runs, 1) = find(eps_norm_CL == min(eps_norm_CL));   % index of smallest CL strain
                tot_runs = tot_runs + 1;
            end
        end
    
        %% Results post-processing
        n_el = size(B_strain, 1);
    
        clearvars success_rates
        for i = 1:n_el
            success_rates(i, 1:2) = [sum(min_strain_OL == i), sum(min_strain_CL == i)];
        end
        success_rates = array2table([[1:n_el]', round(success_rates/size(min_strain_OL, 1)*100)], 'VariableNames',...
                        {'el', 'OL', 'CL'});  
        results(damel, :) = success_rates(damel, :);
    end
        if show_plots
            m = mean(strains,2);            % mean of each characteristic strain
            s_norm = strains ./ max(m);     % normalise rows by largest mean value
            s_s = std(s_norm,0,2);          % standard deviation of un-normalised strain array
            m_norm = mean(s_norm,2);        % mean of rows normalised by largest mean (the strain field to be plotted)
    
            % Plot results
            f2 = figure;
            hold on
    
            x = [1:n_el];  % positions of the bars
            b1 = bar(x, m_norm(:, 1), 'k');
            b2 = bar(x+x(end), m_norm(:,2), 'w');
    
            for i = 1:2*n_el
                end_pos = [m_norm(i), m_norm(i)]; % y-values of whiskers
                end_pos = end_pos + [s_s(i), -s_s(i)];
                end_pos(find(end_pos<0)) = 1e-15;
            
                line([i, i], end_pos, 'Marker', '_', 'Color','r')
            end
            
            % legend
            l = legend('OL', 'CL', 'Coeff. of variation');
    
            a2 = gca;
            a2.GridColor = 'k';
            grid on
            a2.XGrid = 'off';
    
            % x axis
            xticks([1:2*n_el])
            xticklabels([string([1:n_el, 1:n_el])])
            xlabel("Element number")
            a2.XTickLabelRotation = 0;
    
            % y axis
            set(gca, 'YScale', 'log')
            ylim([1e-3, 2.1])
            ylabel("Characteristic strain")
            yticks(logspace(-18, 0, 19))
        end
%         plot_poles(Lambda, lambda_est, s_vals, {'Exact', 'Estimated', 's'});
end
results.delta = results.CL - results.OL;
results
Lambda = ReferenceModels.Lambda;
plot_poles(Lambda, lambda_est, s_vals, {'Exact', 'Estimated', 's'});