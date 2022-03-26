clear; clc; close all
% Stand: 23_02_2022_1023
set_up

filename = sprintf('%02d_%03d_%03d',dam(1,1), dam(1,2)*100, nsr*100)
load("simulation/SYSID/"+filename)
load('D:\OneDrive - Aalborg Universitet\Speciale\MatLab\gaindesign\gain_pars.mat')
load('D:\OneDrive - Aalborg Universitet\Speciale\MatLab\gaindesign\01_strain_cond\gains_3.mat')
H_ref = (Mg*s^2 + Cg*s + Kg)^-1;
H_CL_ref = (Mg*s^2 + Cg + Kg + B2*K*cdis)^-1;

SS = StateSpaceModel();
SS.set_io(in_dof, out_dof);

% export_gain_pars
SS.dt_from_FE(Kg, Cg, Mg, dt);
SS_d = StateSpaceModel();
SS_d.set_io(in_dof, out_dof);
SS_d.dt_from_FE(Kg_d, Cg, Mg, dt);

% idx = setdiff(1:n_dof, bc);

tot_runs = 0;

for damel = dam(:, 1)

    for gainset = 3
        K = gains{gainset}
        for run = 1:numel(H_est)
            tic
            % OL
            H = H_est{run};
            H_d = H_est_d{run};

            DeltaH = H_d - H;
            [~, ~, V] = svd(DeltaH);
            d_OL = zeros(n_dof, 1);
            d_OL(idx) = H_ref * B2 * V(:, end);
            eps_OL = B_strain * d_OL;
            
            H_CL = (eye(size(H)) + H * K)^-1 * H;
            H_CL_d = (eye(size(H_d))+H_d*K)^-1*H_d;
            DeltaH_CL = H_CL_d-H_CL;                                % closed-loop transfer matrix change
            [~, ~, V] = svd(DeltaH_CL);
            d_CL = zeros(n_dof, 1);
            d_CL(idx) = H_CL_ref * B2 * V(:, end);
            eps_CL = B_strain * d_CL;        % closed-loop displacement field
            
            

            strains(:, tot_runs+1, 1) = abs(eps_OL);
            strains(:, tot_runs+1, 2) = abs(eps_CL);
            tot_runs = tot_runs + 1;
            min_strain_OL(tot_runs, 1) = find(eps_OL == min(eps_OL));
            min_strain_CL(tot_runs, 1) = find(eps_CL == min(eps_CL));

%             n_el = 14;
%             figure
%             hold on
%             x = [1:n_el];  % positions of the bars
%             b1 = bar(x, abs(eps_OL)/max(abs(eps_OL)), 'k');
%             b2 = bar(x+x(end), abs(eps_CL)/max(abs(eps_CL)), 'w');
%             set(gca, 'YScale','log')
%             xticks([1:2*n_el])
%             xticklabels([string([1:n_el, 1:n_el])])
%             xlabel("Element number", 'FontSize', 10)
%             set(gca, 'XTickLabelRotation', 0);
%             ylim([1e-2, 2])
%             toc
        end
    end
end

%% sucess rates
n_el = 14;
for i = 1:n_el
    success_rates(i, 1:2) = [sum(min_strain_OL ==i), sum(min_strain_CL ==i)]/tot_runs*100;
end
success_rates = array2table([[1:n_el]', round(success_rates)], 'VariableNames', {'element_number', 'detection_rate_OL', 'detection_rate_CL'})

m = mean(strains,2);            % mean of each row
s_norm = strains ./ max(m);     % normalise rows by largest mean value
s_s = std(s_norm,0,2);         % standard deviation of un-normalised strain array
% cov = s_s ./ max(m);            % normalise standard deviation by largest mean value
m_norm = mean(s_norm,2);        % mean of rows normalised by largest mean (the strain field to be plotted)

[m_sorted, idx] = sort(m_norm(setdiff([1:n_el],damel),1,:));
locatability = 1 ./ [m_norm(damel,1) / m_sorted(1, 1), m_norm(damel,2) / m_sorted(1, 2)]

% Plot results
close all

f2 = figure;
hold on
f2.Position(4) = 5;

x = [1:n_el];  % positions of the bars
b1 = bar(x, m_norm(:, 1), 'k');
b2 = bar(x+x(end), m_norm(:,2), 'w');

for i = 1:2*n_el
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
xticks([1:2*n_el])
xticklabels([string([1:n_el, 1:n_el])])
xlabel("Element number", 'FontSize', fs_small)
a2.XTickLabelRotation = 0;

% y axis
set(gca, 'YScale', 'log')
ylim([1e-3, 2.1])
ylabel("Normalized characteristic strain", 'FontSize', fs_small)
yticks([1e-3, 1e-2, 1e-1, 1e0])


% exportgraphics(f2, sprintf('figures/sens_strains_%d.pdf', damel), 'Resolution', 200)

return
