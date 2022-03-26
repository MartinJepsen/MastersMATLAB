function plotmesh(input)
coords = input.coords;
topology = input.topology;



sensors = [4, 6, 8, 9, 12];


f = figure;
hold on
axis equal
% f.Units = 'centimeters';
f.Position(4) = 15;
fs = 12;
xlabel('$x$','Interpreter','latex','FontSize',fs)
ylabel('$y$','Interpreter','latex','FontSize',fs)
set(gca, 'XMinorTick','off')
set(gca, 'YMinorTick','off')
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontName','Palatino Linotype','FontSize',10,'XColor','k')
% a = get(gca,'YTickLabel');
% set(gca,'YTickLabel',a,'FontName','Palatino Linotype','FontSize',10,'YColor','k')
set(gca, 'XMinorGrid', 'off')
set(gca, 'YMinorGrid', 'off')
grid on

fs2 = 9;

% ax = gca;
% ax.FontSize = 10; 
% ax.FontName = 'Palatino Linotype';
% ax.XColor = 'k';
% ax.YColor = 'k';



yscale = max(coords(:,2))-min(coords(:,2));
ymin = min(coords(:,2));
ymax = max(coords(:,2));
xscale = max(coords(:,1))-min(coords(:,1));
xmin = min(coords(:,1));
xmax = max(coords(:,1));
ylim([ymin-0.2*yscale, ymax+0.2*yscale])
xlim([xmin-0.05*xscale, xmax+0.05*xscale])



for i = 1:size(topology,1)

    
if ismember(i, sensors) == 1; nodecol = 'red';
else nodecol = 'black'; end
    
% makes x-y-z coordinates for the end and start node of the ith element as 
% [x1 y1 z1 x2 y2 z1]    
lines(i,:) = [coords(topology(i,1),:) coords(topology(i,2),:)];

% Average coordinates for placing element numbers
avcoord(i,:) = [(lines(i,1)+lines(i,4))/2, (lines(i,2)+lines(i,5))/2,...
                (lines(i,3)+lines(i,6))/2];

% Plots the lines and nodes
p = plot(lines(i,[1,4]),lines(i,[2,5]),'.-k','MarkerSize',40,'linewidth',1.5);

% Plots nodal numbers for all starting nodes
text(lines(i,1),lines(i,2), num2str(topology(i,1)),'color','white',...
    'HorizontalAlignment','center','FontSize',fs2,'Interpreter','Latex')
% Plots nodal numbers for all ending nodes
text(lines(i,4),lines(i,5), num2str(topology(i,2)),'color','white',...
    'HorizontalAlignment','center','FontSize',fs2,'Interpreter','Latex')
% Plots element numbers
text(avcoord(i,1),avcoord(i,2),num2str(i),'color','black',...
    'backgroundcolor','white','Margin',1,'HorizontalAlignment','center',...
    'FontSize',fs2,'Interpreter','Latex')




xticks(unique(coords(:,1)));
yticks(unique(coords(:,2)));

xticklabels(num2str(unique(coords(:,1))));
yticklabels(num2str(unique(coords(:,2))));
% if sum(abs(coords(:,1)))~=0
% xlim([0-max(coords(:,1)*.1) max(coords(:,1)*1.1)])
% end
% 
% if sum(abs(coords(:,2)))~=0
% ylim([0-max(abs(coords(:,2)*.1)) max(coords(:,2)*1.1)])
% end
% 
% if sum(abs(coords(:,3)))~=0
% zlim([0-max(coords(:,3)*.1) max(coords(:,3)*1.1)])
% end
end