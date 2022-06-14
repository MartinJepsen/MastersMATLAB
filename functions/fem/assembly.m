function [Kg, Mg] = assembly(input,meshdata,el_type,el_props)
% Defining input from struct arrays
n_node = meshdata.n_node;
n_dof = meshdata.n_dof;
n_el = meshdata.n_el;

% Assembly script assumes 2 dof per node for bar elements, regardsless of
% 1D case. This temporarily adds these "missing" dof. They are removed
% again at the end of the script.
if n_dof/n_node == 1
    n_dof = n_dof*2;
end

coords = input.coords;
topology = input.topology;
bc = input.bc;

A = el_props.A;
E = el_props.E;
I = el_props.I;
L = el_props.L;
rho = el_props.rho;

Kg = zeros(n_dof,n_dof);
Mg = zeros(n_dof,n_dof);
rot = zeros(n_el,3);


for i = 1:n_el
    
% compute element lengths 
dx = coords(topology(i,1),1)-coords(topology(i,2),1);
dy = coords(topology(i,1),2)-coords(topology(i,2),2);
dz = coords(topology(i,1),3)-coords(topology(i,2),3); 
% direction cosines
c = dx/L(i);
s = dy/L(i);

% rotation about [z y x] axes
% rot =[atan(dy/dx) atan(dz/dx) atan(dz/dy)];

%% Bar element
if el_type == 1
    % Determines the DOF in each node
    for j=1:n_node
        dof_node(j,:)=[j*2-1 j*2];
    end
    
% Determines which DOF belong to which elements, for indexing purposes
dof_element(i,:) = [dof_node(topology(i,1),:), dof_node(topology(i,2),:)];

% Transformation matrix
T=[c s 0 0;
  -s c 0 0;
   0 0 c s;
   0 0 -s c];

% Element stiffness matrix in local system
kel = A(i)*E(i)/L(i)*[1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0];
% Element stiffness matrix in global system
ke = T'*kel*T;

% Element mass matrix in local system
% mel = rho(i)*A(i)*L(i)/6*[2 0 1 0; 0 2 0 1; 1 0 2 0; 0 1 0 2];
% mel = rho(i)*A(i)*L(i)/6*[2 0 1 0; 0 0 0 0; 1 0 2 0; 0 0 0 0];

mel = rho(i)*A(i)*L(i)/2*eye(4);
% display('Using lumped bar element mass matrix')
% Element mass matrix in global system
me = T'*mel*T;

Kg(dof_element(i,:),dof_element(i,:)) =...
        Kg(dof_element(i,:),dof_element(i,:))+ke;
    
Mg(dof_element(i,:),dof_element(i,:)) =...
        Mg(dof_element(i,:),dof_element(i,:))+me;

%% Beam element
elseif el_type == 2
 % Determines the DOF in each node
 % NOTE: does not adapt to 3D case yet
    for j=1:n_node
        dof_node(j,:)=[j*3-2 j*3-1 j*3];
    end
dof_element(i,:) = [dof_node(topology(i,1),:), dof_node(topology(i,2),:)];
    
T = [c s 0 0 0 0;
    -s c 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 c s 0;
    0 0 0 -s c 0;
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

% Element mass matrix in local system
mel = (rho(i)*A(i)*L(i)/420)*[140 0 0 70 0 0;
                       0 156 22*L(i) 0 54 -13*L(i);
                       0 22*L(i) 4*L(i)^2 0 13*L(i) -3*L(i)^2;
                       70 0 0 140 0 0;
                       0 54 13*L(i) 0 156 -22*L(i);
                       0 -13*L(i) -3*L(i)^2 0 -22*L(i) 4*L(i)^2];

% Element mass matrix in global system
me = T'*mel*T;

Kg(dof_element(i,:),dof_element(i,:)) =...
        Kg(dof_element(i,:),dof_element(i,:))+ke;
    
Mg(dof_element(i,:),dof_element(i,:)) =...
        Mg(dof_element(i,:),dof_element(i,:))+me;
end

end
% Applies boundary conditions as specified by bc.
Kg(find(bc==1),:)=[];
Kg(:,find(bc==1))=[];
Mg(find(bc==1),:)=[];
Mg(:,find(bc==1))=[];

%% Removes rows/columns with no stiffness (relevant for bar elements)
% As of 17-09-2021, the element stiffness matrices are two-dimenional.
% For one-dimensional structures, this leads to zero stiffness in the other
% direction. This part of the script removes these columns in Kg.
for j = 1:size(Kg,1)            % iterates over the no. of columns in Kg
    if sum(abs(Kg(:,j)))==0     % if column j only contains zeros then...
       ind(j)=j;                % ...store column number j in ind
       ind(find(ind==0))=[];    % removes 0 indeces
    else
        continue
    end
end
exist_empty_rows = exist('ind','var');
if exist_empty_rows == 1
% disp(['ATTENTION: Zero row/column in system matrices detected.'...,
%             ' Removing rows/columns ',num2str(ind),'.'])
Kg(ind,:)=[];
Kg(:,ind)=[];
Mg(ind,:)=[];
Mg(:,ind)=[];
end

end