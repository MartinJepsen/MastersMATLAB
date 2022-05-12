
function f = plot_poles(Lambda_ref, Lambda_ex, s)
% PLOT_POLES(Lambda_ref, Lambda_ex)
% plot poles of two systems in the complex plane.
% Arguments:
% Lambda_ref: column vector containing poles (sorted by eigenfrequency)
% Lambda_ex: column vector containing poles  (sorted by eigenfrequency)

f = figure;
f.Units = 'centimeters';
f.Position([3, 4]) = [15, 8];
hold on; grid on
plot(real(Lambda_ref), imag(Lambda_ref), 'k.', 'MarkerSize', 15)
plot(real(Lambda_ex), imag(Lambda_ex), 'r.', 'MarkerSize', 12)

try
    plot(real(s), imag(s), 'bx', 'MarkerSize', 15);
catch
end
   

set(gca, 'FontName','Palatino Linotype')
a = gca;
a.XLabel.Interpreter = 'Latex';
a.XLabel.String = '$\Re(s)$';
% a.XLabel.HorizontalAlignment = 'left';
% a.XLabel.VerticalAlignment = 'middle';
% a.XLabel.Position = [1.05*a.XLim(2), 0, 0];
a.XLabel.Color = 'k';

a.YLabel.Interpreter = 'Latex';
a.YLabel.String = '$\Im(s)$';
a.YLabel.Color = 'k';
a.YLabel.HorizontalAlignment = 'right';
a.YLabel.VerticalAlignment = 'middle';
a.YLabel.Rotation = 0;

xline(0, 'k')
yline(0, 'k')

if exist('s', 'var')
    legend('Exact OL', 'Exact CL', 's', 'Location','southeast')
else
    legend('Exact OL', 'Exact CL', 'Location','southeast')
end