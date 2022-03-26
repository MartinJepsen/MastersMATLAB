function [el_props]=elementproperties(meshdata,input)
n_el = meshdata.n_el;
topology = input.topology;
coords = input.coords;

% Material parameters (BERNAL)
% d = 0.0906; % circular cross section diameter
% % A = zeros(n_el,1)+pi*d^2/4;
% A = zeros(n_el,1)+64.5e-4;
% E = zeros(n_el,1)+2e11;
% I = zeros(n_el,1)+pi*d^4/64;
% rho = zeros(n_el,1)+7850;

d = 0.2; % circular cross section diameter
A = zeros(n_el,1)+pi*(d/2)^2;
% A = zeros(n_el,1)+64.5e-4;
E = zeros(n_el,1)+2e11;
I = zeros(n_el,1)+pi*d^4/64;
rho = zeros(n_el,1)+7850;




for i=1:n_el
dx = coords(topology(i,1),1)-coords(topology(i,2),1);
dy = coords(topology(i,1),2)-coords(topology(i,2),2);
dz = coords(topology(i,1),3)-coords(topology(i,2),3); 
L(i) = norm([dx dy dz]);
end

el_props = struct('A',A,'E',E,'I',I,'rho',rho,'L',L);
end