clear; close all

%% Set simulation variables
err = 0.00;
sensor = "dis";
i = 1;
poles = 1:2:13;
pole_fac = 1.12;
nsr = 0.05;
scheme = 1;

% pole_fac = [1.01, 1.05, 1.12];
damages = [0.5, 0.6, 0.75, 0.8, 0.9, 0.95, 0.99, 1];
pole_fac = 1.12;
for dam_ = damages
    results = get_results(nsr, err, dam_, sensor, poles, pole_fac, scheme)
    OL(:, i) = results.OL;
    CL(:, i) = results.CL;
    DEL(:, i) = results.CL - results.OL;
    i = i + 1;
end

%%
DEL = DEL';
OL = OL';
del_neg = find(DEL<0);
del_zero = find(DEL==0);
OL_2 = OL;
OL_2(del_neg) = OL_2(del_neg) + DEL(del_neg);
DEL_2 = abs(DEL);
DEL_2(del_zero) = DEL_2(del_zero) + 0.01;
%%
close all
fig = figure;
x = 1:numel(damages);
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
n_patches = prod(size(OL));

for ii = 1:floor(n_patches)
    a.Children(ii).FaceColor = [33,255,82]/255;
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

a.DataAspectRatio = [1, 1, 35];
axis tight
xlabel('Damage (%)')
xticks(x);
xticklabels(string((1-damages)*100));
a.XLabel.Rotation = -16;
a.XLabel.VerticalAlignment = 'middle';

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
l.Position([1,2]) = [.04, .8];

exportgraphics(fig, sprintf("D:/Programming/MastersLaTeX/figures/damage_size%d.pdf", scheme),'ContentType','image','Resolution',500)