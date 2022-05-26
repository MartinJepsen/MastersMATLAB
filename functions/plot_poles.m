
function f = plot_poles(Lambda_ref, Lambda_ex, s, legend_str)
% PLOT_POLES(Lambda_ref, Lambda_ex)
% plot poles of two systems in the complex plane.
% Arguments:
% Lambda_ref: column vector containing poles (sorted by eigenfrequency)
% Lambda_ex: column vector containing poles  (sorted by eigenfrequency)

Lambda_ref = reshape(Lambda_ref, [numel(Lambda_ref), 1]);
Lambda_ex = reshape(Lambda_ex, [numel(Lambda_ex), 1]);

f = figure;
f.Units = 'centimeters';
f.Position([3, 4]) = [15, 8];
hold on; grid on
plot(real(Lambda_ref), imag(Lambda_ref), 'k.', 'MarkerSize', 15)
plot(real(Lambda_ex), imag(Lambda_ex), 'r.', 'MarkerSize', 10)

try
    plot(real(s), imag(s), 'bx', 'MarkerSize', 15);
catch
end
   

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


legend(legend_str)
end