
function plot_poles(Lambda_ref, Lambda_ex)
% PLOT_POLES(Lambda_ref, Lambda_ex)
% plot poles of two systems in the complex plane.
% Arguments:
% Lambda_ref: column vector containing poles (sorted by eigenfrequency)
% Lambda_ex: column vector containing poles  (sorted by eigenfrequency)

f = figure;
f.Position(4) = 6.5;
hold on; grid on
% fs = 12;

roots = [Lambda_ref];
n = numel(roots);

for i = 1:2:n   
  q(i) = plot( [real(roots(i)) real(roots(i+1))],[imag(roots(i)) imag(roots(i+1))],...
      'k.-','LineWidth',1,'MarkerSize',8,'LineStyle',':') ;
end

roots = [Lambda_ex];
n = size(roots,1);


%% Plots individually
for j = 1:size(roots,2)
    for i = 1:n       
        e(i) = plot(real(roots(i,j)), imag(roots(i,j)),'r.','MarkerSize',2);
%         string = strjoin({'\lambda_',num2str(i)});
    %     string = strjoin({'\lambda_{',num2str((i-1)/2+1)},'}');
%         text(real(roots(i,j)), imag(roots(i,j)),string,'Interpreter','tex','FontSize',fs,'FontName','palatino linotype')
        
    end
end

set(gca, 'XAxisLocation','origin')
set(gca, 'YAxisLocation','origin')
set(gca, 'FontSize',10)
set(gca, 'FontName','Cambria Math')

a = gca;
% a.XScale = 'log';
a.XLabel.Interpreter = 'Latex';
a.XLabel.String = 'Re';
a.XLabel.FontSize = 12;
a.XLabel.HorizontalAlignment = 'left';
a.XLabel.VerticalAlignment = 'middle';
a.XLabel.Position = [-0.11*a.XLim(1), 0, 0];
a.XLabel.Color = 'k';
% a.XLim(2) = -0.1*a.XLim(1);

a.YLabel.Interpreter = 'Latex';
a.YLabel.String = 'Im';
a.YLabel.FontSize = 12;
a.YLabel.HorizontalAlignment = 'center';
a.YLabel.VerticalAlignment = 'bottom';
a.YLabel.Position = [0, a.YLim(2), 0];
a.YLabel.Color = 'k';

 title('Pole map for exact and estimated models')

legend([q(1), e(1)], {'Exact OL poles', 'Estimated OL poles'},'Location','south west')
