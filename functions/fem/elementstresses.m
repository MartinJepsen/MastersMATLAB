% computes nodal sectional loads for all elements
function [R_el, d_el] = elementstresses(Kg,meshdata,input,el_type,el_props,d_DLV)
n_el = meshdata.n_el;
n_node = meshdata.n_node;
n_dof = meshdata.n_dof;
coords = input.coords;
topology = input.topology;
bc = input.bc;
A = el_props.A;
E = el_props.E;
I = el_props.I;
L = el_props.L;
rho = el_props.rho;
  
% unreduced displacement vector
d = zeros(n_dof,size(d_DLV,2));
d(bc==0,:) = d_DLV;

%%%%%%% iteration over the elements %%%%%%%%%
for i=1:n_el
% element length and orientation
dx = coords(topology(i,1),1)-coords(topology(i,2),1);
dy = coords(topology(i,1),2)-coords(topology(i,2),2);
dz = coords(topology(i,1),3)-coords(topology(i,2),3);
rot =[atan(dy/dx) atan(dz/dx) atan(dz/dy)];


%%%% for bar element model %%%%
if el_type==1
    
% Transformation matrix
T=[cos(rot(1)) sin(rot(1)) 0 0;
    -sin(rot(1)) cos(rot(1)) 0 0;
    0 0 cos(rot(1)) sin(rot(1));
     0 0 -sin(rot(1)) cos(rot(1))];
% Element stiffness matrix in local system
kel = A(i)*E(i)/L(i)*[1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0];    
    
% Determines the DOF in each node
    if el_type ==1 && sum(abs(coords(:,2)))==0 && sum(abs(coords(:,3)))==0
        for j=1:n_node
            dof_node(j,:)=j;
        end
        T([2,4],:)=[];
        T(:,[2,4])=[];
        kel([2,4],:)=[];
        kel(:,[2,4])=[];
    else
        for j=1:n_node
            dof_node(j,:)=[j*2-1 j*2];
        end
    end
% Determines which DOF belong to which elements, for indexing purposes
dof_element(i,:) = [dof_node(topology(i,1),:), dof_node(topology(i,2),:)];
   



% Element stiffness matrix in global system
% ke = T'*kel*T;

R_el(i,:) = kel*(T*d(dof_element(i,:)));
d_el(i,:) = T*d(dof_element(i,:));
end

%%%% for beam element model %%%
if el_type==2
% Determines the DOF in each node
    for j=1:n_node
        dof_node(j,:)=[j*3-2 j*3-1 j*3];
    end
dof_element(i,:) = [dof_node(topology(i,1),:), dof_node(topology(i,2),:)];
    
T = [cos(rot(1)) sin(rot(1)) 0 0 0 0;
    -sin(rot(1)) cos(rot(1)) 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 cos(rot(1)) sin(rot(1)) 0;
    0 0 0 -sin(rot(1)) cos(rot(1)) 0;
    0 0 0 0 0 1];

X1 = A(i)*E(i)/L(i);  % helping parameter
X2 = 12*E(i)*I(i)/(L(i)^3);
X3 = 6*E(i)*I(i)/(L(i)^2);
X4 = 4*E(i)*I(i)/L(i);
X5 = 2*E(i)*I(i)/L(i);

% element stiffness matrix in local coordinates
kel = [X1 0 0 -X1 0 0;
        0 X2 X3 0 -X2 X3;
        0 X3 X4 0 -X3 X5;
        -X1 0 0 X1 0 0;
        0 -X2 -X3 0 X2 -X3;
        0 X3 X5 0 -X3 X4];
% element stiffness matrix in global coordinates
ke = T'*kel*T;

R_el(i,:) = kel*(T*d(dof_element(i,:)));
d_el(i,:) = T*d(dof_element(i,:));
end
    
end

end

