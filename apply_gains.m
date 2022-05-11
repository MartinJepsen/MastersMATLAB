clear; clc; close all

%% Load relevant variables
nsr = 0.05;
load(sprintf("simulation/SYSID/model_error/00_000_%03d", nsr*100))

for damel = [1, 5, 14]
    dam = [damel, 0.95];
    set_up
    
    %% load simulation results
    filename = sprintf('%02d_%03d_%03d',dam(1,1), dam(1,2)*100, nsr*100);
    disp("Loading" + filename)
    load("simulation/SYSID/model_error/"+filename)
    
    %% load gains
    load('gaindesign\01_strain_cond\gains_4.mat')
%     load(sprintf("gaindesign/02_sens/constrained/gains_%02d", damel))
%     load(sprintf("gaindesign/03_strain_norm/gains_%02d", damel))
    
    H_ref = (Mg*s^2 + Cg*s + Kg)^-1;                % reference OL transfer matrix
    H_CL_ref = (Mg*s^2 + Cg*s + Kg + B2*K*cdis)^-1; % reference CL transfer matrix
    A_CL_ex = SS_exact.A + SS_exact.B * B2 * K * cdis * SS_exact.C;
    Lambda_CL_ex = eig(A_CL_ex);                    % exact CL poles

    %% Obtain characteristic strains for every run
    tot_runs = 0;

    for run = 1:numel(SS_est)
        % OL
        SS = SS_est{run};
        H = SS.transfer_matrix(s); 
        SS_d = SS_est_d{run};
        H_d = SS_d.transfer_matrix(s);
        
        l_CL_est = eig(SS.A + SS.B*K*SS.C);
        [~, sorting] = sort(abs(l_CL_est), 'ascend')
        Lambda_CL_est(:, run) = l_CL_est(sorting);          % estimated CL poles

        DeltaH = H_d - H;                               % damage-induced transfer matrix shift (estimated)
        [~, ~, V] = svd(DeltaH);                        % DDLVs
        d_OL = zeros(n_dof, 1);
        d_OL(idx) = H_ref * B2 * V(:, end);             % full OL displacement vector
        eps_OL = B_strain * d_OL;                       % full OL strain vector
        
        % CL
        H_CL = (eye(size(H)) + H * K)^-1 * H;           % estimated CL transfer matrix, reference
        H_CL_d = (eye(size(H_d)) + H_d * K)^-1 * H_d;   % estimated CL transfer matrix, damaged
        DeltaH_CL = H_CL_d - H_CL;                      % CL damage-induced transfer matrix change
        [~, ~, V] = svd(DeltaH_CL);                     % CLDDLVs
        d_CL = zeros(n_dof, 1);
        d_CL(idx) = H_CL_ref * B2 * V(:, end);          % CL displacement vector
        eps_CL = B_strain * d_CL;                       % cCL strain vector
        
        % Compute strains
        strains(:, tot_runs+1, 1) = abs(eps_OL);        % array of characteristic OL strain vectors
        strains(:, tot_runs+1, 2) = abs(eps_CL);        % array of characteristic CL strain vectors
        tot_runs = tot_runs + 1;
        
        eps_norm_OL = abs(eps_OL) / max(abs(eps_OL));
        eps_norm_CL = abs(eps_CL) / max(abs(eps_CL));
        
        min_strain_OL(tot_runs, 1) = find(eps_norm_OL == min(eps_norm_OL));   % index of smallest OL strain
        min_strain_CL(tot_runs, 1) = find(eps_norm_CL == min(eps_norm_CL));   % index of smallest CL strain
    end

    %% Results post-processing
    n_el = 14;
    clearvars success_rates
    for i = 1:n_el
        success_rates(i, 1:2) = [sum(min_strain_OL == i), sum(min_strain_CL == i)]/tot_runs*100;
    end
    success_rates = array2table([[1:n_el]', round(success_rates)], 'VariableNames', {'element_number', 'detection_rate_OL', 'detection_rate_CL'});  
    results(damel, :) = success_rates(damel, :)

    m = mean(strains,2);            % mean of each characteristic strain
    s_norm = strains ./ max(m);     % normalise rows by largest mean value
    s_s = std(s_norm,0,2);          % standard deviation of un-normalised strain array
    m_norm = mean(s_norm,2);        % mean of rows normalised by largest mean (the strain field to be plotted)
    
    [m_sorted, idx] = sort(m_norm(setdiff([1:n_el],damel),1,:));
    locatability = 1 ./ [m_norm(damel,1) / m_sorted(1, 1), m_norm(damel,2) / m_sorted(1, 2)];
    
    % Plot results
    close all
    f2 = figure;
    hold on
    f2.Position([3,4]) = [12, 5.5];
    
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
    percentile = 16;
    l = legend('OL', 'CL', 'Coeff. of variation', 'location', 'south west');
    
    fs_big = 10;
    fs_small = 7;
    a2 = gca;
    grid on
    a2.MinorGridLineStyle = '-';
    a2.GridColor = 'k';
    a2.GridAlpha = 1;
    a2.XGrid = 'off';
    
    % x axis
    xticks([1:2*n_el])
    xticklabels([string([1:n_el, 1:n_el])])
    xlabel("Element number", 'FontSize', fs_small)
    a2.XTickLabelRotation = 0;
    
    % y axis
    set(gca, 'YScale', 'log')
    ylim([1e-3, 2.1])
    ylabel("Normalized characteristic strain", 'FontSize', fs_small)
    yticks([1e-3, 1e-2, 1e-1, 1e0])
    a2.FontSize = 7;
    % exportgraphics(f2, sprintf('figures/sens_strains_%d.pdf', damel), 'Resolution', 200)
    % title(sprintf('Damaged element: %d', damel))
end