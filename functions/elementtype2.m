function [input, el_type, meshdata] = elementtype(input)
% Defining parameters from input structure array
bc = input.bc;
topology = input.topology;
coords = input.coords;


freedof = numel(bc(find(bc==0)));
fixdof = numel(bc)-freedof;
n_el = size(topology,1);
n_node = size(coords,1);



el_type = 1;
if el_type ==1 && sum(abs(coords(:,1)))==0
    bc(1,:)=[];
elseif el_type ==1 && sum(abs(coords(:,2)))==0 && sum(abs(coords(:,3)))==0
    bc([2,3],:)=[];
elseif el_type ==1 && sum(abs(coords(:,3)))==0
    bc(3,:)=[];
end
    freedof = numel(bc(find(bc==0)));
    fixdof = numel(bc)-freedof;
    n_dof = numel(bc);
      if sum(abs(coords(:,3)))==0
        display(['MESSAGE: Creating a 2D'...
        ' bar element model with ',num2str(freedof),...
        ' free DOF, ',num2str(fixdof),' fixed DOF and ',num2str(n_el),' elements.'])
end
%         
%             
%     case 'Beam element (BE)'
%         if sum(abs(coords(:,3)))==0
%             n_dof = 3*n_node;
%         elseif sum(abs(coords(:,3)))~=0
%             n_dof = 6*n_node;
%         end
%         
%         display(['MESSAGE: Creating a ',num2str(size(coords,2)),'D '...
%             'Bernoulli-Euler beam element model with ',num2str(freedof),...
%             ' free DOF, ',num2str(fixdof),' fixed DOF and ',num2str(n_el),' elements.'])
%         el_type = 2;
        
% end



input = struct('coords',coords,'topology',topology,'bc',bc);
meshdata = struct('n_el',size(topology,1),'n_dof',n_dof,'n_node',n_node);
end
