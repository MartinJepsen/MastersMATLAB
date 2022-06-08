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
if ~exist('dt', 'var')
    warning('No dt defined. Setting dt to 1.')
    dt = 1;
end

GeneralParameters.in_dof = in_dof;
GeneralParameters.out_dof = out_dof;
GeneralParameters.r = numel(GeneralParameters.in_dof);
GeneralParameters.m = numel(GeneralParameters.out_dof);
GeneralParameters.dt = dt;
GeneralParameters.t = [0:dt:(n_samples * dt - dt)]';          % time sequence
GeneralParameters.blockrows = blockrows;
GeneralParameters.sensor = sensor;
GeneralParameters.nsr = nsr;
GeneralParameters.err = err;
GeneralParameters.n_runs = n_runs;

[ReferenceModels, GeneralParameters] = generate_reference_models(err, GeneralParameters);
[ReferenceModels, GeneralParameters] = generate_state_space_models(ReferenceModels, GeneralParameters);

t = GeneralParameters.t;
fs = GeneralParameters.dt^-1;
omega_hz = ReferenceModels.FE.modal_parameters.omega / (2*pi);
f0 = 0.8*omega_hz(1);
f1 = 1.2*omega_hz(truncated_mode);
t1 = t(end);

for in = 1:numel(in_dof)
    signal = zeros(1, numel(t));
    signal(:) = 50*chirp(t,f0,t1,f1,'linear',(in-1)*90);
    u(in, :) = signal;
    base_dir = sprintf("simulation/SYSID/model_error_%03d_%s", err*100, sensor);
end
GeneralParameters.u = u;

GeneralParameters.base_dir = base_dir;
if ~isfolder(base_dir)
    mkdir(base_dir)
    addpath(base_dir)
end

filepath = fullfile(base_dir, "SetUp.mat");
save(filepath, "GeneralParameters", "ReferenceModels")
copyfile(filepath, "gaindesign/01/SetUp.mat")
copyfile(filepath, "gaindesign/02/SetUp.mat")
copyfile(filepath, "gaindesign/03/SetUp.mat")

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
    ReferenceModels.SS_exact = SS_exact;
%     DamagedModels.Lambda_d = Lambda_d;
end