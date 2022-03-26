clear
clc
close all


tot_runs = 0;
n_dof = 6;
% s = complex(-0.032216900785070, 3.543681888973263);         % pole
% s = -0.6782 +23i
tic

results_number = 0;

for damel = 2
    load("simulation/system_matrices/08500/model_" + num2str(damel))
%     s = 0;
    

    for gainset = damel
        load("optimised_DDLV/gain_sensitivity/unit_perturbations/gains_" + num2str(gainset))  
        K = gains{1,1};
        norm_K = gains{1,2};

        for s_val = 1:2:numel(Lambda)
            s = [complex(real(Lambda(s_val)), 1.1*imag(Lambda(s_val)))];
            G_ref = (Mg*s^2+Cg*s+Kg)^-1;
            for run = 0:999
                
                % Load estimated models
                simpath = sprintf("simulation/SYSID/results_%03d/", results_number) ...
                            + sprintf("%02d/", damel) + sprintf("run_%03d",run);
                load(simpath)
    
                in_dof = model.in_dof; r =  numel(in_dof);
                out_dof = model.out_dof; m = numel(out_dof);
                
                b2 = zeros(n_dof, r);
                for ii=1:r
                    b2(in_dof(ii), ii) = 1;
                end
                
                cdis = zeros(m, n_dof);
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
                F = zeros(n_dof, 1);
                v = V(:, end);
                F(in_dof) = v;
                d_OL_est = [0; G_ref * F];
        
                % Closed-loop transfer matrices
                G_CL_ref = (Mg*s^2 + Cg*s + Kg + b2*K*cdis)^-1;
                G_CL_est = (eye(size(G_est))+G_est*K)^-1*G_est;
                G_CL_est_d = (eye(size(G_est_d))+G_est_d*K)^-1*G_est_d;
                DeltaG_CL = G_CL_est_d-G_CL_est;                                % closed-loop transfer matrix change
    
                [~, S, V] = svd(DeltaG_CL);
                ss =  diag(S)./max(diag(S));
        
                F = zeros(n_dof,1);
                v = V(:, end);
                F(in_dof) = v;
                d_CL_est = [0; G_CL_ref*F];        % closed-loop displacement field
        
                for el = 1:n_dof
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

for i = 1:n_dof
    success_rates(i, 1:2) = [sum(min_strain_OL ==i), sum(min_strain_CL ==i)]/tot_runs*100;
end

strains2 = [strains(:,:,1); strains(:,:,2)];
strains_mean = mean(strains2,2);
strains_mean2 = [strains_mean(1:n_dof) / max(strains_mean(1:n_dof));
                strains_mean(n_dof+1:end) / max(strains_mean(n_dof+1:end))];

s_strains = std(strains2,0,2);
cov = s_strains ./ strains_mean;
% cov = [cov(1:n_dof), cov(n_dof+1:end)] ./ [max(strains_mean(1:n_dof)), max(strains_mean(n_dof+1:end))];
% cov = [s_strains(1:n_dof), s_strains(n_dof+1:end)] ./ [max(strains_mean(1:n_dof)), max(strains_mean(n_dof+1:end))];

success_rates = round(success_rates,2)

% model
toc


% PLot results
close all

f2 = figure;
hold on
f2.Position(4) = 6;

x = [1:6];  % positions of the bars
bar_data = [strains_mean2(1:n_dof); strains_mean2(n_dof+1:end)];
b1 = bar(x, bar_data(1:n_dof,:), 'k');
b2 = bar(x+n_dof, bar_data(n_dof+1:end, :), 'w');

% arrows
% strains2 = [strains(:,:,1); strains(:,:,2)];
% for i = 1:2*n_dof
%     end_pos = [strains_mean(i) - cov(i), strains_mean(i) + cov(i)];
%     line([i, i], end_pos, 'Marker', '_', 'Color','r')
% end

% strains_mean = [strains_mean(:,1); strains_mean(:,2)];
strains3 = strains2;
for i = 1:2*n_dof
%     end_pos = [strains_mean2(i), strains_mean2(i)]  + ([prctile(strains3(i,:), 84),  - prctile(strains3(i,:), 16)] ./ [max(strains_mean(1:n_dof)), max(strains_mean(n_dof+1:end))]);
% end_pos = [strains_mean2(i), strains_mean2(i)]  + ([prctile(strains3(i,:), 84),  - prctile(strains3(i,:), 16)] ./ [strains_mean(i)]);
% end_pos = [strains_mean2(i), strains_mean2(i)]  + ([s_strains(i),  - s_strains(i)] ./ [max(strains_mean(1:n_dof)), max(strains_mean(n_dof+1:end))]);
end_pos = [strains_mean2(i), strains_mean2(i)]  + ([s_strains(i),  - s_strains(i)] ./ [strains_mean(i)]);
    line([i, i], end_pos, 'Marker', '_', 'Color','r')
end

% x axis
xticks([1:12])
xticklabels(["1", "2", "3", "4", "5", "6", "1", "2", "3", "4", "5", "6"])
xlabel("Element number")

% y axis
set(gca, 'YScale', 'log')
ylim([1e-2, 3])
ylabel("Normalised strain")

% legend
% percentile_str = sprintf("%d$^{th}$", percentile) + " \& " + sprintf("%d$^{th}$ percentile", 100-percentile);
l = legend('Mean strain (OL)', 'Mean strain (CL)', 'Coefficient of variation', 'location', 'south west');
% l.Title.String = 'Configuration';

% title
a2 = gca;
a2.Title.String = {sprintf("Damage in element %d", damel)};

% grid
grid on
a2.MinorGridLineStyle = '-';
a2.GridColor = 'k';
a2.GridAlpha = 1;
a2.XGrid = 'off';

return
