function [ReferenceModels, GeneralParameters] = generate_reference_models(err, GeneralParameters)
    %% Exact reference model
    n_el = 10;
    stiffness = [10:-1:1]'*1e3;
    mass = ones*3;
    damping_ratios = [0.01, 0.05];
    modes = [1, 6];
    FE = ChainSystem();
    FE.assembly(stiffness, mass)
    FE.rayleigh_damping(damping_ratios, modes)
    
    %% Error reference model
    stiffness_e = stiffness .* unifrnd(1-err, 1+err, [n_el, 1]);
    FE_e = ChainSystem();
    FE_e.assembly(stiffness_e, mass)
    FE_e.Cg = FE.modal_parameters.alpha * FE_e.Mg + FE.modal_parameters.beta * FE_e.Kg;
    FE_e.get_modal_parameters();
    
    ReferenceModels.FE = FE;
    ReferenceModels.FE_e = FE_e;
    ReferenceModels.Kg = FE.Kg;
    ReferenceModels.Cg = FE.Cg;
    ReferenceModels.Mg = FE.Mg;
    ReferenceModels.Kg_e = FE_e.Kg;
    ReferenceModels.Cg_e = FE_e.Cg;

    GeneralParameters.n_dof = FE.n_dof+1;
    GeneralParameters.free_dof = size(FE.Kg, 1);
    GeneralParameters.n_el = n_el;
    GeneralParameters.bc = 1;
    GeneralParameters.B_strain = FE.B;
    GeneralParameters.idx = setdiff(1:GeneralParameters.n_dof, GeneralParameters.bc);
end