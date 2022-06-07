clear
load('C:\Users\MAJP\MAJP_personal\Masters_thesis\MastersMATLAB\simulation\SYSID\model_error_002_dis\SetUp.mat')

Kg = ReferenceModels.Kg;
Cg = ReferenceModels.Cg;
Mg = ReferenceModels.Mg;
in_dof = GeneralParameters.in_dof;
out_dof = GeneralParameters.out_dof;

s = 1;

[~,v1,v3] = intersect(out_dof,in_dof);

unsorted_dof = [out_dof, setdiff(in_dof,out_dof)];
[~,order] = sort(unsorted_dof)

SS = StateSpaceModel();
SS.set_io(in_dof, out_dof);
SS.dt_from_FE(Kg, Cg, Mg, GeneralParameters.dt, "dis");
SS.get_modal_parameters();
SS.to_ct();

A_r = SS.A;
B_r = SS.B;
C_r = SS.C;

kappa = union(in_dof, out_dof);



%%

SS_f = StateSpaceModel();
SS_f.set_io(kappa, kappa);
SS_f.dt_from_FE(Kg, Cg, Mg, GeneralParameters.dt, "dis");
SS_f.get_modal_parameters();
B_f = SS_f.B;
C_f = SS_f.C;
% SS_f.to_ct();


[Bexp,Cexp]=ExpandBandC2(A_r,B_r,C_r,v1,v3,order);
Bexp = real(Bexp);
Cexp = real(Cexp);

SS_s = StateSpaceModel();
SS_s.set_io(kappa, kappa);
SS_s.dt_from_FE(Kg, Cg, Mg, GeneralParameters.dt, "dis");
SS_s.get_modal_parameters();
% SS_s.to_ct();
% SS_s.B = Bexp;
SS_s.C = Cexp;

t = GeneralParameters.t;
u = randn(numel(kappa),numel(t));


[u, y_s] = SS_s.time_response(u, t, 0, false);
[~, y_f] = SS_f.time_response(u, t, 0, false);

t = TRAC(y_s, y_f);
sum(round(t,3)==1,'all')

%%
close all
figure
hold on
plot(y_f(2,:),'r')
plot(y_s(7,:),'k--')