
clear

if ~exist('damage', 'var')
    warning('No damage defined. Setting damaged configuration equal to the reference.')
    damage = [1, 1];
end
if ~exist('sensor', 'var')
    warning('No sensor type defined. Setting sensor to displacements.')
    sensor = "dis";
end
if ~exist('err', 'var')
    warning('No model error defined. Setting model error to 0.')
    err = 0;
end
if ~exist('nsr', 'var')
    warning('No nsr defined. Setting nsr to 0.')
    nsr = 0;
end
if ~exist('n_runs', 'var')
    warning('No number of runs defined. Setting number of runs to 4.')
    n_runs = 4;
end

dt = 0.0001;
GeneralParameters.in_dof = [1:12];
GeneralParameters.out_dof = [1:12];
GeneralParameters.r = numel(GeneralParameters.in_dof);
GeneralParameters.m = numel(GeneralParameters.out_dof);
GeneralParameters.dt = dt;
GeneralParameters.t = [0:dt:(20000 * dt - dt)]';          % time sequence
GeneralParameters.u = randn(numel(GeneralParameters.in_dof), numel(GeneralParameters.t));
GeneralParameters.blockrows = 60;
GeneralParameters.sensor = sensor;
GeneralParameters.nsr = nsr;
GeneralParameters.err = err;
GeneralParameters.n_runs = n_runs;

[ReferenceModels, GeneralParameters] = generate_reference_models(err, GeneralParameters);
[ReferenceModels, GeneralParameters] = generate_state_space_models(ReferenceModels, GeneralParameters);

function [ReferenceModels, GeneralParameters] = generate_state_space_models(ReferenceModels,...
                                                                            GeneralParameters);
    Kg = ReferenceModels.Kg;
    Cg = ReferenceModels.Cg;
    Mg = ReferenceModels.Mg;
    free_dof = GeneralParameters.free_dof;
    sensor = GeneralParameters.sensor;
    in_dof = GeneralParameters.in_dof;
    out_dof = GeneralParameters.out_dof;
    dt = GeneralParameters.dt;
    m = GeneralParameters.m;
    r = GeneralParameters.r;

    SS_exact = StateSpaceModel();
    SS_exact.set_io(1:free_dof, 1:free_dof, 2*free_dof);
    SS_exact.dt_from_FE(Kg, Cg, Mg, dt, sensor);
    SS_exact.get_modal_parameters();
    SS_exact.to_ct();
    Lambda = SS_exact.modal_parameters.Lambda;

%     % Exact, damaged
%     SS_exact_d = StateSpaceModel();
%     SS_exact_d.set_io(1:free_dof, 1:free_dof, 2*free_dof);
%     SS_exact_d.dt_from_FE(Kg_d, Cg_d, Mg, dt, sensor);
%     SS_exact_d.get_modal_parameters();
%     SS_exact_d.to_ct();
%     Lambda_d = SS_exact_d.modal_parameters.Lambda;

    cdis = zeros(m, free_dof);
    for ii=1:m
        cdis(ii, out_dof(ii)) = 1;
    end
    
    GeneralParameters.B2 = StateSpaceModel().set_io(in_dof, out_dof, 2 * free_dof);
    GeneralParameters.cdis = cdis;
    ReferenceModels.Lambda = Lambda;
%     DamagedModels.Lambda_d = Lambda_d;
end