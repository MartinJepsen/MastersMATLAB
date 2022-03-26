clear; clc; close all
% Stand: 23_02_2022_1023
tot_runs = 0;
n_dof = 6;

results_number = 0;

for damel = 10
%     load("simulation/system_matrices/09500/model_" + num2str(damel))
load('D:\OneDrive - Aalborg Universitet\Speciale\MatLab\optimised_DDLV\gain_pars.mat')
    for gainset = 3
        load("optimised_DDLV/01_strain_cond/gains_" + num2str(gainset))
%         load("optimised_DDLV/02_sens/unit_perturbations/gains_" + num2str(gainset))  
%         load("optimised_DDLV/03_strain_norm/unit_perturbations/gains_" + num2str(gainset))  
%         K = gains{1,1};
%         norm_K = gains{1,2};

        for s_val = 1
%             s = [complex(real(Lambda(s_val)), 1.1*imag(Lambda(s_val)))];
            G_ref = (Mg*s^2+Cg*s+Kg)^-1;
            for run = 0:5
                
                % Load estimated models
%                 simpath = sprintf("simulation/SYSID/results_%03d/", results_number) ...
%                             + sprintf("%02d/", damel) + sprintf("run_%03d",run);
                simpath = sprintf("simulation/SYSID/results_%03d/", results_number) ...
                            + sprintf("%02d/", damel) + sprintf("run_%03d",run);
                load(simpath)
    
                in_dof = model.in_dof; r =  numel(in_dof);
                out_dof = model.out_dof; m = numel(out_dof);
                
                b2 = zeros(free_dof, r);
                for ii=1:r
                    b2(in_dof(ii), ii) = 1;
                end
                
                cdis = zeros(m, free_dof);
                for ii=1:m
                    cdis(ii, out_dof(ii)) = 1;
                end
    
        
                % Estimated, undamaged open-loop state space model
                [A, B] = dt2ct(model.sys_u.A, model.sys_u.B, model.dt);
                G_est = maketm(A, B, model.sys_u.C, 0, s);
                    
                % Open-loop, damaged, estimated
                [A, B] = dt2ct(model.sys_d.A, model.sys_d.B, model.dt);
                G_est_d = maketm(A, B, model.sys_d.C, 0, s);
    
                %DDLV
                DeltaG_est = G_est_d - G_est;
                [~, ~, V] = svd(DeltaG_est);
                F = zeros(free_dof, 1);
                v = V(:, end);
                F(in_dof) = v;
                d_OL_est = [0; G_ref * F];
        
                % Closed-loop transfer matrices
                G_CL_ref = (Mg*s^2 + Cg*s + Kg + b2*K*cdis)^-1;
                G_CL_est = (eye(size(G_est))+G_est*K)^-1*G_est;
                G_CL_est_d = (eye(size(G_est_d))+G_est_d*K)^-1*G_est_d;

                G_CL_ = zeros(n_dof, free_dof);
                G_CL_(idx, :) = G_CL;
                
                DeltaG_CL = G_CL_est_d-G_CL_est;                                % closed-loop transfer matrix change
    
                [~, S, V] = svd(DeltaG_CL);
                ss =  diag(S)./max(diag(S));
        
                F = zeros(n_dof,1);
                v = V(:, end);
                F(in_dof) = v;
                d_CL_est = [0; G_CL_ref*F];        % closed-loop displacement field
        
                for el = 1:free_dof
                    eps_OL(el,:) = [1 -1]*d_OL_est(el:el+1);
                    eps_CL(el,:) = [1 -1]*d_CL_est(el:el+1);
                end
        
%                 eps_OL = abs(eps_OL)./max(abs(eps_OL));
%                 eps_CL = abs(eps_CL)./max(abs(eps_CL));
        
                strains(:, tot_runs+1, 1) = abs(eps_OL);
                strains(:, tot_runs+1, 2) = abs(eps_CL);
    
                tot_runs = tot_runs + 1;
    
                min_strain_OL(tot_runs, 1) = find(eps_OL == min(eps_OL));
                min_strain_CL(tot_runs, 1) = find(eps_CL == min(eps_CL));
      
            end
        end
    end
end

%% sucess rates
for i = 1:free_dof
    success_rates(i, 1:2) = [sum(min_strain_OL ==i), sum(min_strain_CL ==i)]/tot_runs*100;
end
round(success_rates)

m = mean(strains,2);            % mean of each row
s_norm = strains ./ max(m);     % normalise rows by largest mean value
s_s = std(s_norm,0,2);         % standard deviation of un-normalised strain array
% cov = s_s ./ max(m);            % normalise standard deviation by largest mean value
m_norm = mean(s_norm,2);        % mean of rows normalised by largest mean (the strain field to be plotted)

[m_sorted, idx] = sort(m_norm(setdiff([1:free_dof],damel),1,:));
locatability = 1 ./ [m_norm(damel,1) / m_sorted(1, 1), m_norm(damel,2) / m_sorted(1, 2)]

% Plot results
close all

f2 = figure;
hold on
f2.Position(4) = 5;

x = [1:6];  % positions of the bars
b1 = bar(x, m_norm(:, 1), 'k');
b2 = bar(x+free_dof, m_norm(:,2), 'w');

for i = 1:2*free_dof
    end_pos = [m_norm(i), m_norm(i)]; % y-values of whiskers

    % 4) coefficients of variation (standard deviations normalised by largest normalised means)
%     if i<=n_dof
%         end_pos = end_pos + ([s_strains(i), -s_strains(i)] ./ [max(m_norm(1:n_dof))]);
%     else
%         end_pos = end_pos + ([s_strains(i), -s_strains(i)] ./ [max(m_norm(n_dof+1:end))]);
%     end
    % 5) coefficients of variation (standard deviations normalised by actual means)
    end_pos = end_pos + [s_s(i), -s_s(i)];
    end_pos(find(end_pos<0)) = 1e-15;

    line([i, i], end_pos, 'Marker', '_', 'Color','r')
end



% legend
percentile = 16;
% percentile_str = sprintf("%d$^{th}$", percentile) + " \& " + sprintf("%d$^{th}$ percentile", 100-percentile);
l = legend('OL', 'CL', 'Coeff. of variation', 'location', 'south west');

fs_big = 10;
fs_small = 9;
a2 = gca;
% title
% a2.Title.String = {sprintf("Strain field with damage in element %d", damel)};
% subtitle('Objective function 1')
% grid
grid on
a2.MinorGridLineStyle = '-';
a2.GridColor = 'k';
a2.GridAlpha = 1;
a2.XGrid = 'off';

% x axis
xticks([1:12])
xticklabels(["1", "2", "3", "4", "5", "6", "1", "2", "3", "4", "5", "6"])
xlabel("Element number", 'FontSize', fs_small)

% y axis
set(gca, 'YScale', 'log')
ylim([1e-3, 2.1])
ylabel("Normalized characteristic strain", 'FontSize', fs_small)
yticks([1e-3, 1e-2, 1e-1, 1e0])


% exportgraphics(f2, sprintf('figures/sens_strains_%d.pdf', damel), 'Resolution', 200)

return
