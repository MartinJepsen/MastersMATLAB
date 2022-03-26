% function for plotting mode shapes. Only works with 2D beam elements
function plotstaticdef(input,meshdata,el_props,el_type,scale,d)
n_dof = meshdata.n_dof;

coords = input.coords;
topology = input.topology;
bc = input.bc;

L = el_props.L;

disp(['MESSAGE: Plotting static deformation shape.'])
figure
title(['Deformed structure, deformations scaled by a factor ',num2str(scale),'.'])

for i = 1:size(topology,1)
% makes x-y-z coordinates for the end and start node of the ith element as 
% [x1 y1 z1 x2 y2 z1]    
lines(i,:) = [coords(topology(i,1),:) coords(topology(i,2),:)];
end

d_full = zeros(n_dof,1);
d_full(find(bc==0))=d_full(find(bc==0))+d;

% Excluding rotational DOF from mode shape
if el_type == 2
    d_full(3:3:end)=[];
end
% Normalising mode shape by largest displacement

% Adding mode shapes to nodal coordinates
def_shape = lines;



for i=1:size(topology,1)
    % DOF numbers of the start node of the i'th element
    ind1 = [(topology(i,1)*2-1):(topology(i,1)*2)];
    % DOF numbers of the end node of the i'th element
    ind2 = [(topology(i,2)*2-1):(topology(i,2)*2)];
    % Line start and end points in format [x1 y1 z1 x2 y2 z2]
    d_def(i,:)=[[d_full(ind1)]' 0 [d_full(ind2)]' 0];
end

def_shape = lines+d_def*scale;

hold on
axis equal
grid on
for i = 1:size(topology,1)
% Average coordinates for placing element numbers
avcoord(i,:) = [(lines(i,1)+lines(i,4))/2, (lines(i,2)+lines(i,5))/2,...
                (lines(i,3)+lines(i,6))/2];

% Plots the lines and nodes
plot3(lines(i,[1,4]),lines(i,[2,5]),lines(i,[3,6]),'.-k','MarkerSize',50)
plot3(def_shape(i,[1,4]),def_shape(i,[2,5]),def_shape(i,[3,6]),'.-r','MarkerSize',25)
%% Element and node numbering
% Plots nodal numbers for all starting nodes
text(lines(i,1),lines(i,2),lines(i,3), num2str(topology(i,1)),'color','white',...
    'HorizontalAlignment','center')
% Plots nodal numbers for all ending nodes
text(lines(i,4),lines(i,5),lines(i,6), num2str(topology(i,2)),'color','white',...
    'HorizontalAlignment','center')
% Plots element numbers
text(avcoord(i,1),avcoord(i,2),avcoord(i,3),num2str(i),'color','black','backgroundcolor','white')
xlabel('x')
ylabel('y')
zlabel('z')

if sum(abs(coords(:,1)))~=0
xlim([0-max(coords(:,1)*.1) max(coords(:,1)*1.1)])
end

if sum(abs(coords(:,2)))~=0
ylim([0-max(abs(coords(:,2)*.1)) max(coords(:,2)*1.1)])
end

if sum(abs(coords(:,3)))~=0
zlim([0-max(coords(:,3)*.1) max(coords(:,3)*1.1)])
end

end
end