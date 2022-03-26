function [eps] = barstrains(meshdata,input,d)
% EPS = barstrains(meshdata,input,d)
% Computes nodal strains for a structure of bar elements.
% meshdata: struct containing n_el, n_node, n_dof
% input: struct containing coords, topology and bc
% d: displacement field without constrained DOF


n_el = meshdata.n_el;
n_node = meshdata.n_node;
n_dof = meshdata.n_dof;
coords = input.coords;
topology = input.topology;
bc = input.bc;

% Includes constrained dof in d
d_ = zeros(n_dof,1);
d_(bc==0) = d;
d = d_;

for i=1:n_el
% element length and orientation
dx = coords(topology(i,1),1)-coords(topology(i,2),1);
dy = coords(topology(i,1),2)-coords(topology(i,2),2);
dz = coords(topology(i,1),3)-coords(topology(i,2),3);
rot =[atan(dy/dx) atan(dz/dx) atan(dz/dy)];

%%%% for bar element model %%%%
% if el_type==1

% Transformation matrix
T=[cos(rot(1)) sin(rot(1))  0           0;
  -sin(rot(1)) cos(rot(1))  0           0;
   0           0            cos(rot(1)) sin(rot(1));
   0           0           -sin(rot(1)) cos(rot(1))];
    
% Strain matrix in local system
B = [-1, 0, 1, 0];
    
% Determines the DOF in each node
for j=1:n_node
    dof_node(j,:)=[j*2-1 j*2];
end

% Determines which DOF belong to which elements, for indexing purposes
dof_element = [dof_node(topology(i,1),:), dof_node(topology(i,2),:)];

eps(i,:) = B*(T*d(dof_element));
end


eps = abs(eps)/max(abs(eps));
% %%%% for beam element model %%%
% if el_type==2
%     
%     % Determines the DOF in each node
%     for j=1:n_node
%         dof_node(j,:)=[j*3-2 j*3-1 j*3];
%     end
%     dof_element = [dof_node(topology(i,1),:), dof_node(topology(i,2),:)];
%     
%     T = [cos(rot(1)) sin(rot(1)) 0 0 0 0;
%         -sin(rot(1)) cos(rot(1)) 0 0 0 0;
%         0 0 1 0 0 0;
%         0 0 0 cos(rot(1)) sin(rot(1)) 0;
%         0 0 0 -sin(rot(1)) cos(rot(1)) 0;
%         0 0 0 0 0 1];
%     
%     B0 = [-1/L(i), 0, 0, 1/L(i), 0, 0;
%         0, (12*0)/L(i)^3 - 6/L(i)^2, (6*0)/L(i)^2 - 4/L(i), 0, 6/L(i)-(12*0)/L(i)^3, (6*0)/L(i)^2 - 2/L(i)];
%     BL = [-1/L(i), 0, 0, 1/L(i), 0, 0;
%         0, (12*L(i))/L(i)^3 - 6/L(i)^2, (6*L(i))/L(i)^2 - 4/L(i), 0, 6/L(i)^2-(12*L(i))/L(i)^3, (6*L(i))/L(i)^2 - 2/L(i)];
% 
%     eps(i,:) = BL*(T*d(dof_element));%-B0*(T*d(dof_element)));
%     
%     
%     
% end
%    
% end

%% Post processing of strains
% for i = 1:2 eps(:,i) = abs(eps(:,i)) / max(abs(eps(:,i))); end
% 
%     avg_ax = mean(abs(eps(:,1)));
%     avg_bend = mean(abs(eps(:,2)));
%     
%     threshold = 1/1000;
%     ax_dam = find(abs(eps(:,1))<= threshold*avg_ax);
%     bend_dam = find(abs(eps(:,2))<= threshold*avg_bend);
%     damaged_elements = NaN(n_el,2);
%     if exist('ax_dam')==1
%         damaged_elements(1:numel(ax_dam),1) = ax_dam;
%     end
%     if exist('bend_dam')==1
%         damaged_elements(1:numel(bend_dam),2) = bend_dam;
%     end
%         damaged_elements
    