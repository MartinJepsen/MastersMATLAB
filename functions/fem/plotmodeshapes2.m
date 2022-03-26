% function for plotting mode shapes. Only works with 2D beam elements
function plotmodeshapes2(input,meshdata,el_props,Phi,modenum,el_type,scale)
n_dof = meshdata.n_dof;

coords = input.coords;
topology = input.topology;
bc = input.bc;

L = el_props.L;

% Characteristic size of structure for scaling purposes
char_size = sqrt((max(coords(:,1))-min(coords(:,1)))^2+...
                 (max(coords(:,2))-min(coords(:,2)))^2+...
                 (max(coords(:,3))-min(coords(:,3)))^2);
             
disp(['MESSAGE: Plotting mode shape number ',num2str(modenum),'.'])
f = figure
f.Units = 'centimeters';
% f.Position(4) = 8;
fs = 12;
if modenum == 5; xlabel('$x$','Interpreter','latex','FontSize',fs); end
ylabel('$y$','Interpreter','latex','FontSize',fs)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Palatino Linotype','FontSize',10,'XColor','k')
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'FontName','Palatino Linotype','FontSize',10,'YColor','k')

xticks(unique(coords(:,1)));
yticks(unique(coords(:,2)));

xticklabels(unique(coords(:,1)));
yticklabels(unique(coords(:,2)));
hold on
axis equal
grid on


for i = 1:size(topology,1)
% makes x-y-z coordinates for the end and start node of the ith element as 
% [x1 y1 z1 x2 y2 z1]    
lines(i,:) = [coords(topology(i,1),:) coords(topology(i,2),:)];
end


% Creating eigenvector including constrained DOF
Phi_mode = zeros(n_dof,1);
ind = find(bc==0);

% For a 1D case, Phi only contains displacements in one direction. The
% following adds zeros for the other direction
if sum(abs(coords(:,2)))==0 && el_type == 1
    Phi = reshape([Phi(:,modenum)'; zeros(size(Phi(:,modenum)'))],[],1)
    Phi_mode(ind) = Phi;
elseif sum(abs(coords(:,1)))==0 && el_type == 1
    Phi = reshape([zeros(size(Phi(:,modenum)')); Phi(:,modenum)'],[],1);
    Phi_mode(ind) = Phi;
else
    Phi_mode(ind)=Phi(:,modenum);
end

% Excluding rotational DOF from mode shape
if el_type == 2
    Phi_mode(3:3:end)=[];
end
% Normalising mode shape by largest displacement
Phi_mode = Phi_mode/max(abs(Phi_mode));
% Adding mode shapes to nodal coordinates
modeshape = lines;

for i=1:size(topology,1)
    % DOF numbers of the start node of the i'th element
    ind1 = [(topology(i,1)*2-1):(topology(i,1)*2)];
    % DOF numbers of the end node of the i'th element
    ind2 = [(topology(i,2)*2-1):(topology(i,2)*2)];
    % Line start and end points in format [x1 y1 z1 x2 y2 z2]
    shape_coords(i,:)=[[Phi_mode(ind1)]' 0 [Phi_mode(ind2)]' 0];
end

modeshape = lines+shape_coords*(scale*char_size);


for i = 1:size(topology,1)
%     % Average coordinates for placing element numbers
%     avcoord(i,:) = [(lines(i,1)+lines(i,4))/2, (lines(i,2)+lines(i,5))/2,...
%                     (lines(i,3)+lines(i,6))/2];

    % Plots the lines and nodes
    plot3(lines(i,[1,4]),lines(i,[2,5]),lines(i,[3,6]),'.-k','MarkerSize',10,'linewidth',1.2)
end

%% Element and node numbering
for i = 1:size(topology,1)
    plot3(modeshape(i,[1,4]),modeshape(i,[2,5]),modeshape(i,[3,6]),'.-r','MarkerSize',10,'linewidth',1.2)
end

% Plots nodal numbers for all starting nodes
% text(lines(i,1),lines(i,2),lines(i,3), num2str(topology(i,1)),'color','white',...
%     'HorizontalAlignment','center','FontSize',fs*0.8,'Interpreter','Latex')
% % Plots nodal numbers for all ending nodes
% text(lines(i,4),lines(i,5),lines(i,6), num2str(topology(i,2)),'color','white',...
%     'HorizontalAlignment','center','FontSize',fs*0.8,'Interpreter','Latex')
% % Plots element numbers
% text(avcoord(i,1),avcoord(i,2),avcoord(i,3),num2str(i),'color','black'...
%     ,'backgroundcolor','white','FontSize',fs*0.8,'Interpreter','Latex')

if sum(abs(coords(:,1)))~=0
xlim([0-max(coords(:,1)*.05) max(coords(:,1)*1.05)])
end

% if sum(abs(coords(:,2)))~=0
% ylim([-max(abs(modeshape(:,2)*.4)) max(modeshape(:,2)*1.5)])
% end
ylim([-3, 12])

if sum(abs(coords(:,3)))~=0
zlim([0-max(coords(:,3)*.1) max(coords(:,3)*1.1)])
end


end