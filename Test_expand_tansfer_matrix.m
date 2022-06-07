clear
load('C:\Users\MAJP\MAJP_personal\Masters_thesis\MastersMATLAB\simulation\SYSID\model_error_002_dis\SetUp.mat')

Kg = ReferenceModels.Kg;
Cg = ReferenceModels.Cg;
Mg = ReferenceModels.Mg;
in_dof = GeneralParameters.in_dof;
out_dof = GeneralParameters.out_dof;

s = complex(1.2, 3.1);

[~,v1,v3] = intersect(out_dof,in_dof);
unsorted_dof = [out_dof, setdiff(in_dof,out_dof)];
[~,order] = sort(unsorted_dof);

kappa = union(in_dof, out_dof);

SS_f = StateSpaceModel();
SS_f.set_io(kappa, kappa);
SS_f.dt_from_FE(Kg, Cg, Mg, GeneralParameters.dt, "dis");
SS_f.get_modal_parameters();
SS_f.to_ct();
H_f = SS_f.transfer_matrix(s);

% [Bexp,Cexp]=ExpandBandC2(A_r,B_r,C_r,v1,v3,order);
% Bexp = real(Bexp);
% Cexp = real(Cexp);

SS_s = StateSpaceModel();
SS_s.set_io(in_dof, out_dof);
SS_s.dt_from_FE(Kg, Cg, Mg, GeneralParameters.dt, "dis");
SS_s.get_modal_parameters();
SS_s.to_ct();
SS_s.expand();
H_s = SS_s.transfer_matrix(s);
