clear; close all

%% Set simulation variables
% nsr = 0.02;
err = 0.00;
dam_ = 0.75;
sensor = "dis";
i = 1;
poles = 1:2:13;
pole_fac = 1.12;
% nsr = [0.01, 0.03, 0.05, 0.07];
nsr = 0.05;

% pole_fac = [1.01, 1.05, 1.12];
damages = [0.5, 0.6, 0.75, 0.8, 0.9, 0.95, 0.99, 1];
pole_fac = 1.12;
for dam_ = damages
    results = get_results(nsr, err, dam_, sensor, poles, pole_fac)
    OL(:, i) = results.OL;
    CL(:, i) = results.CL;
    DEL(:, i) = results.CL - results.OL;
    i = i + 1;
end

%%
close all
fig = figure;
DEL = DEL';
OL = OL';

del_neg = find(DEL<0);
del_zero = find(DEL==0);
OL_2 = OL;
OL_2(del_neg) = OL_2(del_neg) + DEL(del_neg);
DEL_2 = abs(DEL);
DEL_2(del_zero) = DEL_2(del_zero) + 0.01;

x = 1:numel(damages);
y = 1:8;
z(:, :, 1) = OL_2;
z(:, :, 2) = DEL_2;
[y1,x1]=meshgrid(y,x);
ngroups = 2;
z1=zeros(size(z,1),size(z,2));    % initial 'zldata'
for i1=1:ngroups
    z2=z1;
    z1=z1+squeeze(z(:,:,i1));
    h(i1)=CREATESTACKEDMULTIBAR3d(x1, y1, z2, z1, i1.*ones(numel(z1(:)),1), 1, ngroups);
    hold on
end

a = gca;
n_patches = prod(size(OL));

for ii = 1:floor(n_patches)
    a.Children(ii).FaceColor = [33,255,82]/255;
end
for ii = (n_patches+1):2*n_patches
    a.Children(ii).FaceColor = [200,200,200]/255;
    a.Children(ii).EdgeColor = 'k';
    a.Children(ii).FaceAlpha = 0.3;
end
for ii = 1:numel(del_neg)
    a.Children(n_patches+1-del_neg(ii)).FaceColor = 'r';
end
for ii = 1:numel(del_zero)
    a.Children(n_patches+1-del_zero(ii)).FaceColor = 'b';
end

axis tight
xlabel('Damage (%)')
xticks(x);
xticklabels(string((1-damages)*100));
a.XLabel.Rotation = -16;
a.XLabel.VerticalAlignment = 'middle';

ylabel('Damage pattern')
yticks(y);
a.YLabel.Rotation = 19;
a.YLabel.VerticalAlignment = 'middle';

zlim([0, 100])
zticks([0:10:100])
a.ZAxis.MinorTick = 'on';
a.GridColor = 'k';
a.GridAlpha = 1;
a.XGrid = 'off';
a.YGrid = 'off';
a.ZGrid = 'on';
zlabel('POL (%)')

view(45,45)
box on
l = legend([a.Children(end), a.Children(1), a.Children(n_patches+1-del_zero(1)), a.Children(n_patches+1-del_neg(1))], ...
    {'DDLV', 'CLDDLV (better)', 'CLDDLV (same)', 'CLDDLV (worse)'});
l.Position([1,2]) = [.04, .8];

% exportgraphics(fig, 'C:\Users\MAJP\MAJP_personal\Masters_thesis\MastersLaTeX\figures\damage_size.eps')


%%
function results = get_results(nsr, err, dam_, sensor, poles, im_fac)
    %% Compute results
    base_dir = sprintf("simulation/SYSID/model_error_%03d_%s", err*100, sensor);
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
    
    elements = 1:8;
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
    %             load(sprintf("gaindesign/01_strain_cond/gains_%d_1.010.mat", pole))
            load(sprintf("gaindesign/01_strain_cond/gains_%d_%0.3f.mat", pole, im_fac))
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
                H_arr{run_u, 1} = H;
                H_CL_arr{run_u, 1} = (eye(size(H)) + H * K)^-1 * H;
            end
            for run_d = 1:n_sim_d
                H_d = s_fac * SS_est_d{run_d}.transfer_matrix(s);
                H_d_arr{run_d, 1} = H_d;
                H_CL_d_arr{run_d, 1} = (eye(size(H_d)) + H_d * K)^-1 * H_d;   % estimated CL transfer matrix, damaged
            end
    
            A_CL_ex = SS_exact.A + SS_exact.B * B2 * K * cdis * SS_exact.C;
            Lambda_CL = eig(A_CL_ex);                       % exact CL poles
        
            % model transfer matrices
            H_ref = (Mg*s^2 + Cg*s + Kg)^-1;                % reference OL transfer matrix
            H_CL_ref = (Mg*s^2 + Cg*s + Kg + B2*K*cdis)^-1; % reference CL transfer matrix
    
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
    end
    results.delta = results.CL - results.OL;
end

% results
% Lambda = ReferenceModels.Lambda;
% plot_poles(Lambda, lambda_est, s_vals, {'Exact', 'Estimated', 's'});