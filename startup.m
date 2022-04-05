wdir = "D:\OneDrive - Aalborg Universitet\Speciale\MastersMATLAB";

cd(wdir);

addpath(genpath(pwd));

%% Run get(groot, 'factoryObjectName') to get all changeable properties of ObjecName object.

%% set figure defaults

set(groot, "defaultFigureUnits", "centimeters")
set(groot, "defaultFigurePosition", [10, 10, 18, 8])
% set(groot, "defaultFigureFontSize", 12)

set(groot, "defaultAxesTickLabelInterpreter", "Latex")
set(groot, "defaultLegendInterpreter", "Latex")
set(groot, "defaultAxesFontName", "Palatino Linotype")
set(groot, "defaultAxesFontSize", 9)
set(groot, {"defaultAxesXColor", "defaultAxesYColor", "defaultAxeszColor"}, {'k', 'k', 'k'})

%% set default order of colors
set(groot,'defaultAxesColorOrder', [0, 0, 0; eye(3)])