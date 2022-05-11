clear; close all;
results = table2array(readtable("svalues.xlsx"));

OL = results(:, [2:2:end]);
CL = results(:, [3:2:end]);

im = [1.5, 1.1, 1.05, 1.01, 1];

f = figure;
hold on
i = 14;
    plot(im, CL(i, :), 'x-')
    plot(im, OL(i, :), 'x-')
    legend_str{i} = sprintf("%d", i);
    title(sprintf('Damage in elenent %d', i))
%     set(gca, 'YScale', 'log')
xlabel('$\frac{\Im(s)}{\Im(\lambda_1)}$','Interpreter','latex', 'FontSize', 14)
ylabel('Rate of sucessful localisation')
set(gca, 'XScale', 'log')
legend('CL', 'OL', 'location', 'south east')

grid on
% legend(legend_str)