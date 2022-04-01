
function [A, B, C, D] = n4sid_(yk, uk, ii)
% [A, B, C, D] = N4SID_(YK, UK, I)
% YK: system output
% UK: system input
% ii: system order



% Dimension check of captured measurements
[nu,Nu]=size(uk);
if Nu<nu
    uk=uk';
    [nu,Nu]=size(uk);
end
[ny,Ny]=size(yk);
if Ny < ny
    yk=yk';
    [ny,Ny]=size(yk);
end
if Nu~=Ny
    disp('Discrepancy in input and output sample lengths')
    return
end

% User-defined number of block rows
s=Nu;
m=size(uk,1);
l=size(yk,1);
j=s-2*ii+1; % Number of columns in block Hankel matrix

%% Block Hankel matrix
uk1=uk(:,ii+1:s)/sqrt(j);
yk1=yk(:,ii+1:s)/sqrt(j);
uk=uk/sqrt(j);
yk=yk/sqrt(j);
            
for k=1:ii 
    Up((k-1)*nu+1:k*nu,:)=uk(:,k:k+j-1);
    Yp((k-1)*ny+1:k*ny,:)=yk(:,k:k+j-1);
    Uf((k-1)*nu+1:k*nu,:)=uk1(:,k:k+j-1);
    Yf((k-1)*ny+1:k*ny,:)=yk1(:,k:k+j-1);
end

%% LQ Decomposition
U = [Up;Uf];
Y = [Yp;Yf];
L = triu(qr([U;Y]'))';
L = L(1:2*ii*(m+l),1:2*ii*(m+l));

%% Oblique projection
Lf = L((2*m+l)*ii+1:2*(m+l)*ii,:);              
Lp = [L(1:m*ii,:);L(2*m*ii+1:(2*m+l)*ii,:)];    
Lu=L(m*ii+1:2*m*ii,:);                         
Lfp = Lf*(eye(2*(m+l)*ii)-Lu'*pinv(Lu*Lu')*Lu);
Lpp = Lp*(eye(2*(m+l)*ii)-Lu'*pinv(Lu*Lu')*Lu);

Oi = Lfp*pinv(Lpp)*Lp;                          
 
%% SVD of Oblique projection
[U1,S,V1] = svd(Oi);

%% Plot singular values to determine the rank
ss=diag(S);
f = figure;
f.Units = 'centimeters';
f.Position = [2, 2, 15, 6];
fs = 12;
% bar((1:l*ii), ss/max(ss), 1,'FaceColor',[1, 1, 1]);

% plot((1:l*ii),ss/max(ss),'k-x')
% set(gca,'YScale','log')
% grid on
% xlabel('Model order','Interpreter','Latex','FontSize',fs);
% ylabel('Singular value (normalised)','Interpreter','Latex','FontSize',fs);
% zoom('on');
% YScale = 'log';
% % xticks([0:5:size(Oi,1)])
% ns=input('Provide system order via command window. Press Ctrl+C to cancel.');
% close all
ns = 8;

%% Retrieve the system matrices
U2=U1(:,1:ns);
S1=S(1:ns,1:ns);

gammai=U2*sqrt(S1);
gammaip1=gammai(1:l*(ii-1),:);

Lf=L((2*m+l)*ii+l+1:2*(m+l)*ii,:);              
Lp=[L(1:m*(ii+1),:);L(2*m*ii+1:(2*m+l)*ii+l,:)]; 
Lu=L(m*ii+m+1:2*m*ii,:);                        
Lfp=Lf*(eye(2*(m+l)*ii)-Lu'*pinv(Lu*Lu')*Lu);   
Lpp=Lp*(eye(2*(m+l)*ii)-Lu'*pinv(Lu*Lu')*Lu);    
 
Oip=Lfp*pinv(Lpp)*Lp;                            
 
Xi=pinv(gammai)*Oi;
Xip=pinv(gammaip1)*Oip;

Left=[Xip ; L((2*m+l)*ii+1:(2*m+l)*ii+l,:)];
Right=[Xi ; L(m*ii+1:m*(ii+1),:)]; 
sol=Left/Right;

A =sol(1:ns,1:ns);
B =sol(1:ns,ns+1:ns+m);
C =sol(ns+1:ns+l,1:ns);
D =sol(ns+1:ns+l,ns+1:ns+m);
disp('The estimated matrices have been returned')