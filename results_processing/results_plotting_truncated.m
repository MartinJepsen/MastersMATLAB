clear; close all

%% Set simulation variables
scheme = 1;
err = 0.00;
sensor = "dis";
poles = 1:2:13;
nsr = 0.05;
dam_ = 0.8;
pole_fac = 1.12;
modes = [14:-2:2];
for i = 1:numel(modes)
    mode = modes(i);
    poles = 1:2:(mode+1)/2;
    results = get_results(nsr, err, dam_, sensor, poles, pole_fac, scheme, mode);
    OL(:, i) = results.OL;
    CL(:, i) = results.CL;
    DEL(:, i) = results.CL - results.OL;
end

%% Process data
DEL = DEL';
OL = OL';
del_neg = find(DEL<0);
del_zero = find(DEL==0);
OL_2 = OL;
OL_2(del_neg) = OL_2(del_neg) + DEL(del_neg);
DEL_2 = abs(DEL);
DEL_2(del_zero) = DEL_2(del_zero) + 0.01;

%% Make figure
close all
fig = figure;
x = [1:numel(modes)];
y = 1:8;
z(:, :, 1) = OL_2;
z(:, :, 2) = DEL_2;
[y1,x1]=meshgrid(y,x);
ngroups = 2;
z1=zeros(size(z,1),size(z,2));    % initial 'zldata'
for i1=1:ngroups
    z2=z1;
    z1=z1+squeeze(z(:,:,i1));
    h(i1)=CREATESTACKEDMULTIBAR3d(x1, y1, z2, z1, i1.*ones(numel(z1(:)),1), 1, ngroups);
    hold on
end

a = gca;
a.DataAspectRatio = [1, 1, 35];
n_patches = prod(size(OL));

for ii = 1:floor(n_patches)
    a.Children(ii).FaceColor = [33,255,82]/255;
    a.Children(ii).FaceAlpha = 1;
    a.Children(ii).EdgeAlpha = 1;
end
for ii = (n_patches+1):2*n_patches
    a.Children(ii).FaceColor = [200,200,200]/255;
    a.Children(ii).EdgeColor = 'k';
    a.Children(ii).FaceAlpha = 0.3;
end

h_DDLV = a.Children(end);
handles = [h_DDLV];
labels ={'DDLV'};
h_neg = NaN;
h_zero = NaN;
h_plus = NaN;
idx_plus = find(DEL>1);
if ~isempty(idx_plus)
    idx_plus = idx_plus(1);
    h_plus = a.Children(idx_plus);
    handles = [handles, h_plus];
    labels = {labels{:}, 'CLDDLV (better)'};
end
if ~isempty(del_zero)
    for ii = 1:numel(del_zero)
        a.Children(n_patches+1-del_zero(ii)).FaceColor = 'b';
    end
    h_zero = a.Children(n_patches+1-del_zero(1));
    handles = [handles, h_zero];
    labels = {labels{:}, 'CLDDLV (same)'};
end
if ~isempty(del_neg)
    for ii = 1:numel(del_neg)
        a.Children(n_patches+1-del_neg(ii)).FaceColor = 'r';
    end
    h_neg = a.Children(n_patches+1-del_neg(1));
    handles = [handles, h_neg];
    labels = {labels{:}, 'CLDDLV (worse)'};
end

axis tight
xlabel('Model order','Interpreter','latex')
xticks(x);
xticklabels(string((modes)));
a.XLabel.Rotation = -19;
a.XLabel.VerticalAlignment = 'bottom';
a.XLabel.HorizontalAlignment = 'l';
a.XLabel.Position = [2,0,-30];

ylabel('Damage pattern')
yticks(y);
a.YLabel.Rotation = 19;
a.YLabel.VerticalAlignment = 'middle';

zlim([0, 100])
zticks([0:10:100])
a.ZAxis.MinorTick = 'on';
a.GridColor = 'k';
a.GridAlpha = 1;
a.XGrid = 'off';
a.YGrid = 'off';
a.ZGrid = 'on';
zlabel('POL (%)')

view(45,20)
box on


l = legend(handles,labels,'Orientation','vertical');
l.Position([1,2]) = [.07, .75];

% exportgraphics(fig, sprintf('D:/Programming/MastersLaTeX/figures/svals%d.pdf',scheme),'ContentType','image','Resolution',500)