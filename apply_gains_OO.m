clear; clc; close all
% Stand: 23_02_2022_1023
set_up

dam = [10, 0.8];
FE = FiniteElementModel('paper_truss.xlsx');
FE.assembly('bar', dam);
FE.apply_bc([1, 2, 9, 10]);
FE.modal_damping(0.02);
FE.strains_from_disp([])

SS_exact = StateSpaceModel();
SS_exact.set_io(1:12, 1:12, 24);
SS_exact.dt_from_FE(FE.Kg, FE.Cg, FE.Mg, dt);
SS_exact.get_modal_parameters();
Lambda = SS_exact.modal_parameters.Lambda;
s = complex(real(Lambda(1)), 1.1*imag(Lambda(1)));         % pole
z = exp(s * dt);
SS_exact.transfer_matrix(z);

H_ref = SS_exact.H;

SS = StateSpaceModel();
SS.set_io(in_dof, out_dof);

export_gain_pars
SS.dt_from_FE(Kg, Cg, Mg, dt);

SS_d = StateSpaceModel();
SS_d.set_io(in_dof, out_dof);
SS_d.dt_from_FE(Kg_d, Cg, Mg, dt);

idx = setdiff(1:n_dof, bc);

cdis = zeros(m, free_dof);
for ii=1:m
    cdis(ii, out_dof(ii)) = 1;
end

tot_runs = 0;

for damel = dam(:, 1)

    for gainset = 1
        load("optimised_DDLV/01_strain_cond/gains_" + num2str(gainset))
%         load("optimised_DDLV/02_sens/unit_perturbations/gains_" + num2str(gainset))  
%         load("optimised_DDLV/03_strain_norm/unit_perturbations/gains_" + num2str(gainset))  
        load('D:\OneDrive - Aalborg Universitet\Speciale\MatLab\optimised_DDLV\gain_pars.mat')
        
        for run = 0:5
            tic

            SS.time_response(u, t, nsr, true);
            SS.estimate([],[], blockrows);
            SS.get_modal_parameters();
            SS.transfer_matrix(z);

            SS_d.time_response(u, t, nsr, true);
            SS_d.estimate([],[], blockrows);
            SS_d.get_modal_parameters();
            SS_d.transfer_matrix(z);
    
            H_est = SS.H;
            H_est_d = SS_d.H;
            DeltaH_est = H_est_d - H_est;

            %DDLV
            [~, ~, V] = svd(DeltaH_est);
            F = zeros(free_dof, 1);
            F(in_dof) = V(:, end);
            d_OL = zeros(n_dof, 1);
            d_OL(idx) = H_ref * F;
            eps_OL = B_strain * d_OL;
            
%%             
            % Closed-loop transfer matrices
%             SS_exact_CL = StateSpaceModel();
%             SS_exact_CL.set_io(in_dof, out_dof);
%             SS_exact_CL.dt_from_FE(Kg, Cg, Mg, dt);
%             SS_exact_CL.transfer_matrix(z);
%             SS_exact_CL.to_cl(K);
%             A_CL = SS_exact.A - SS_exact.B * B2 * K * cdis;
            H_CL_ref = (eye(size(H_ref))+H_ref*B2*K*cdis)^-1*H_ref;
            
%             H_CL_ref = SS_exact_CL.H;
            H_CL_est = (eye(size(H_est))+H_est*K)^-1*H_est;
            H_CL_est_d = (eye(size(H_est_d))+H_est_d*K)^-1*H_est_d;
            
            DeltaH_CL = H_CL_est_d-H_CL_est;                                % closed-loop transfer matrix change

            [~, ~, V] = svd(DeltaH_CL);
    
            F = zeros(free_dof, 1);
            F(in_dof) = V(:, end);
            d_CL = zeros(n_dof, 1);
            d_CL(idx) = H_CL_ref * F;
            eps_CL = B_strain * d_CL;        % closed-loop displacement field
            

            strains(:, tot_runs+1, 1) = abs(eps_OL);
            strains(:, tot_runs+1, 2) = abs(eps_CL);
            tot_runs = tot_runs + 1;
            min_strain_OL(tot_runs, 1) = find(eps_OL == min(eps_OL));
            min_strain_CL(tot_runs, 1) = find(eps_CL == min(eps_CL));

            toc
        end
    end
end

%% sucess rates
n_el = 14;
for i = 1:n_el
    success_rates(i, 1:2) = [sum(min_strain_OL ==i), sum(min_strain_CL ==i)]/tot_runs*100;
end
round(success_rates)

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
