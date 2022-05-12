clear; close all;

load("s_values/results_pole3.mat")

f = figure;
hold on
i = 14;
    plot(s_facs, CL(i, :), 'x-')
    plot(s_facs, OL(i, :), 'x-')
    legend_str{i} = sprintf("%d", i);
    title(sprintf('Damage in element %d', i))

xlabel('$\frac{\Im(s)}{\Im(\lambda_1)}$','Interpreter','latex', 'FontSize', 14)
ylabel('Rate of sucessful localisation')
set(gca, 'XScale', 'log')
legend('CL', 'OL', 'location', 'south east')

grid on
% legend(legend_str)