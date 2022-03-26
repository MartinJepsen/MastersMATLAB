
function plot_poles(Lambda_ref, Lambda_ex)
% PLOT_POLES(Lambda_ref, Lambda_ex)
% plot poles of two systems in the complex plane.
% Arguments:
% Lambda_ref: column vector containing poles (sorted by eigenfrequency)
% Lambda_ex: column vector containing poles  (sorted by eigenfrequency)

f = figure;
f.Position(4) = 6;
hold on; grid on
% fs = 12;

roots = [Lambda_ref];
n = numel(roots);

for i = 1:2:n   
  q(i) = plot( [real(roots(i)) real(roots(i+1))],[imag(roots(i)) imag(roots(i+1))],...
      'k.-','LineWidth',1,'MarkerSize',15,'LineStyle',':') ;
end

roots = [Lambda_ex];
n = size(roots,1);

%% Plots in pairs
% for j = 1:size(roots,2)
% for i = 1:2:n
%     
%     e(i) = plot( [real(roots(i,j)) real(roots(i+1,j))],[imag(roots(i,j)) imag(roots(i+1,j))],'r.','MarkerSize',15);
%     string = strjoin({'\lambda_',num2str((i-1)/2+1),' snr_',num2str(j)},'');
%     string = strjoin({'\lambda_{',num2str((i-1)/2+1)},'}');
%     text(real(roots(i,j)),imag(roots(i,j)),string,'Interpreter','tex','FontSize',fs,'FontName','palatino linotype')
%     
% end
% end

%% Plots individually
for j = 1:size(roots,2)
    for i = 1:n       
        e(i) = plot(real(roots(i,j)), imag(roots(i,j)),'r.','MarkerSize',5);
        string = strjoin({'\lambda_',num2str(i)});
    %     string = strjoin({'\lambda_{',num2str((i-1)/2+1)},'}');
%         text(real(roots(i,j)), imag(roots(i,j)),string,'Interpreter','tex','FontSize',fs,'FontName','palatino linotype')
        
    end
end

set(gca, 'XAxisLocation','origin')
set(gca, 'YAxisLocation','origin')
% set(gca, 'FontSize',10)
set(gca, 'FontName','Palatino Linotype')

a = gca;
a.Position(1) = 0.07
% a.XScale = 'log';
a.XLabel.Interpreter = 'Latex';
a.XLabel.String = 'Re($s$)';
% a.XLabel.FontSize = 10;
a.XLabel.HorizontalAlignment = 'left';
a.XLabel.VerticalAlignment = 'middle';
a.XLabel.Position = [-0.05*a.XLim(1), 0, 0];
a.XLabel.Color = 'k';
% a.XLim(2) = -0.1*a.XLim(1);

a.YLabel.Interpreter = 'Latex';
a.YLabel.String = 'Im($s$)';
% a.YLabel.FontSize = 10;
a.YLabel.HorizontalAlignment = 'center';
a.YLabel.VerticalAlignment = 'bottom';
% a.YLabel.Position = [0.075*a.XLim(1), 0.9*a.YLim(2), 0];
a.YLabel.Position = [0, 1.08*a.YLim(2), 0]
a.YLabel.Color = 'k';
a.YLabel.Rotation = 0;


%%
% ylim([-45, 45])
%  title('Pole map for exact and estimated models')

legend([q(1), e(1)], {'Exact OL', 'Estimated OL'},'Location','south east')
