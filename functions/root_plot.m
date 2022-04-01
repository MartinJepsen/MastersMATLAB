function root_plot(l_ex, l_n4sid)


% Plot roots
roots = [l_ex];
f = figure;
f.Units = 'centimeters';
f.Position = [2, 2, 15, 8];
hold on; grid on
% title(['$i=$',num2str(modelorder),', $\Delta t = $',num2str(dt),', $t_j=$',num2str(t(end))],'interpreter','latex','fontsize',12)
n = numel(roots);

for i = 1:2:n   
  q(i) = plot( [real(roots(i)) real(roots(i+1))],[imag(roots(i)) imag(roots(i+1))],...
      'k+-','LineWidth',1.5,'MarkerSize',8,'LineStyle',':') ;
end

roots = [l_n4sid];
n = numel(roots);
for i = 1:2:n
    e(i) = plot( [real(roots(i)) real(roots(i+1))],[imag(roots(i)) imag(roots(i+1))],'rx--','MarkerSize',8);
    
end

xlim([-5, 2])
xL = xlim;
yL = ylim;
% line([0 0], yL);  %x-axis
% line(xL, [0 0]);  %y-axis
set(gca, 'XAxisLocation','origin')
set(gca, 'YAxisLocation','origin')
set(gca, 'FontSize',11)
set(gca, 'FontName','Cambria Math')

a = gca;
a.XLabel.Interpreter = 'Latex';
a.XLabel.String = 'Re';
a.XLabel.FontSize = 12;
a.XLabel.HorizontalAlignment = 'left';
a.XLabel.VerticalAlignment = 'middle';
a.XLabel.Position = [a.XLim(2), 0, 0];
a.XLabel.Color = 'k';


a.YLabel.Interpreter = 'Latex';
a.YLabel.String = 'Im';
a.YLabel.FontSize = 12;
a.YLabel.HorizontalAlignment = 'center';
a.YLabel.VerticalAlignment = 'bottom';
a.YLabel.Position = [0, a.YLim(2), 0];
a.YLabel.Color = 'k';


legend([q(1), e(1)], {'Reference', 'Estimated'},'Location','North west','Interpreter','Latex','FontSize',11)

% plot( [real(roots(1)) real(roots(1+ceil(n/2)))],[imag(roots(1)) imag(roots(1+ceil(n/2)))] , '--*k')