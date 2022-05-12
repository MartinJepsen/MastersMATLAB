clear; close all;
results = table2array(readtable("svalues.xlsx", 'Sheet', 'Sheet2'));

OL = results(:, [2:2:end]);
CL = results(:, [3:2:end]);

im = sort([1:0.05:1.5, 1.001, 1.002, 1.003], 'ascend');

f = figure;
hold on
i = 7;
    plot(im, CL(i, :), 'x-')
    plot(im, OL(i, :), 'x-')
    legend_str{i} = sprintf("%d", i);
    title(sprintf('Damage in element %d', i))
%     set(gca, 'YScale', 'log')
xlabel('$\frac{\Im(s)}{\Im(\lambda_1)}$','Interpreter','latex', 'FontSize', 14)
ylabel('Rate of sucessful localisation')
set(gca, 'XScale', 'log')
legend('CL', 'OL', 'location', 'south east')

grid on
% legend(legend_str)