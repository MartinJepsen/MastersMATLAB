function results = get_results(nsr, err, dam_, sensor, poles, im_fac, scheme, mode, elements, expand, deg)
    %% Compute results
if mode ~= 0
    base_dir = sprintf("simulation/SYSID/t%d_model_error_%03d_%s", mode, err*100, sensor);
elseif mode == 0
    base_dir = sprintf("simulation/SYSID/model_error_%03d_%s", err*100, sensor);
end

load(sprintf("%s/00_000_%03d", base_dir, round(nsr*100,0)))
load(fullfile(base_dir, "SetUp.mat"))


if expand
    for i_u = 1:numel(SS_est)
        SS_est{i_u}.expand();
    end
    GeneralParameters = expand_coordinates(GeneralParameters);
end


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

for i_e = elements
    tot_runs = 1;

    dam = [i_e, dam_];

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
    
    polenum = 1;
    for i_p = poles
        
        % load gainss
        if scheme == 1
            if expand
                load(sprintf("gaindesign/01/exp_gains/gains_%d_%0.3f.mat", i_p, im_fac))
            else
                load(sprintf("gaindesign/01/gains/K_%d_%d_%0.2f.mat", i_p,deg, im_fac))
            end     
        elseif scheme == 2
            if expand
                load(sprintf("gaindesign/02/exp_gains/gain%d_%d_%0.3f.mat", i_e, i_p, im_fac))
            else
                load(sprintf("gaindesign/02/gains/gain%d_%d_%0.3f.mat", i_e, i_p, im_fac))
            end                
        elseif scheme == 3
            if expand
                load(sprintf("gaindesign/03/exp_gains/gain%d_%d_%0.3f.mat", i_e, i_p, im_fac))
            else
                load(sprintf("gaindesign/03/gains/gain%d_%d_%0.3f.mat", i_e, i_p, im_fac))
            end
        end
        s_vals(polenum) = s;
        
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
        for i_u = 1:n_sim
            H = s_fac * SS_est{i_u}.transfer_matrix(s);
            H_arr{i_u, 1} = H;
            H_CL_arr{i_u, 1} = (eye(size(H,1)) + H * K)^-1 * H;
            A_CL_est = SS_est{i_u}.A +  SS_est{i_u}.B * K * SS_est{i_u}.C;
            Lambda_CL_est(:,i_u) = eig(A_CL_est);                       % exact CL poles
        end
        for i_d = 1:n_sim_d
            SS_d = SS_est_d{i_d};
            if expand && polenum==1
                SS_d.expand();
            end
            H_d = s_fac * SS_d.transfer_matrix(s);
            H_d_arr{i_d, 1} = H_d;
            H_CL_d_arr{i_d, 1} = (eye(size(H_d,1)) + H_d * K)^-1 * H_d;   % estimated CL transfer matrix, damaged
        end

        A_CL_ex = SS_exact.A + SS_exact.B * B2 * K * cdis * SS_exact.C;
%         A_CL_ex = SS_exact.A + SS_exact.B * K * SS_exact.C;
        Lambda_CL = eig(A_CL_ex);                       % exact CL poles
    
        % model transfer matrices
        H_ref = (Mg*s^2 + Cg*s + Kg)^-1;                % reference OL transfer matrix
        H_CL_ref = (Mg*s^2 + Cg*s + Kg + B2*K*cdis)^-1; % reference CL transfer matrix

        for i_u = 1:n_sim
            H = H_arr{i_u};
            H_CL = H_CL_arr{i_u};
            
            for i_d = 1:n_sim_d
                H_d = H_d_arr{i_d};
                H_CL_d = H_CL_d_arr{i_d};
    
                DeltaH = H_d - H;                               % damage-induced transfer matrix shift (estimated)
                [~, ~, V] = svd(DeltaH);                        % DDLVs
                d_OL = zeros(n_dof, 1);
                d_OL(idx) = H_ref * B2 * V(:, end);             % full OL displacement vector
                eps_OL = B_strain * d_OL;                       % full OL strain vector
            
                DeltaH_CL = H_CL_d - H_CL;                      % CL damage-induced transfer matrix change
                [~, ~, V] = svd(DeltaH_CL);                     % CLDDLVs
                d_CL = zeros(n_dof, 1);
                d_CL(idx) = H_CL_ref * B2 * V(:, end);          % CL displacement vector
                eps_CL = B_strain * d_CL;                       % CL strain vector
            
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
%         idx_s = ((polenum-1)*n_sim*n_sim_d+1):(tot_runs-1);
%         strains(:, idx_s, 1) = strains(:, idx_s,1) / max(strains(:, idx_s,1),[],'all');
%         strains(:, idx_s, 2) = strains(:, idx_s,2) / max(strains(:, idx_s,2),[],'all');
        polenum = polenum + 1;
    end

    
    %% Results post-processing
    clearvars success_rates
    n_el = size(B_strain, 1);
    for i = 1:n_el
        success_rates(i, 1:2) = [sum(min_strain_OL == i), sum(min_strain_CL == i)];
    end
    success_rates = array2table([[1:n_el]', round(success_rates/size(min_strain_OL, 1)*100)], 'VariableNames',...
                    {'el', 'OL', 'CL'});
    results(i_e, :) = success_rates(i_e, :);
    
end
results.delta = results.CL - results.OL;
end