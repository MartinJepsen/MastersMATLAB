function [ReferenceModels, GeneralParameters] = generate_reference_models(err, GeneralParameters)
    %% Exact reference model
    FE = FiniteElementModel();
    FE.from_xlsx('structures/paper_truss.xlsx');
    FE.assembly("bar", [1, 2, 9, 10]);
    FE.modal_damping(0.02);
    FE.strains_from_disp([]);
    n_el = size(FE.mesh.topology, 1);
    
    %% Error reference model
    FE_e = FiniteElementModel();
    FE_e.from_xlsx('structures/paper_truss.xlsx');
    FE_e.mesh.element_properties.E = FE.mesh.element_properties.E .* unifrnd(1-err, 1+err, [n_el, 1]);
    FE_e.assembly("bar", [1, 2, 9, 10]);
    FE_e.modal_damping(0.02);
    
    ReferenceModels.FE = FE;
    ReferenceModels.FE_e = FE_e;
    ReferenceModels.Kg = FE.Kg;
    ReferenceModels.Cg = FE.Cg;
    ReferenceModels.Mg = FE.Mg;
    ReferenceModels.Kg_e = FE_e.Kg;
    ReferenceModels.Cg_e = FE_e.Cg;

    GeneralParameters.n_dof = FE.n_dof;
    GeneralParameters.free_dof = size(FE.Kg, 1);
    GeneralParameters.n_el = n_el;
    GeneralParameters.bc = FE.mesh.bc_dof;
    GeneralParameters.B_strain = FE.B;
    GeneralParameters.idx = setdiff(1:FE.n_dof, FE.mesh.bc_dof);

end