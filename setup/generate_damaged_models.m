function [DamagedModels] = generate_damaged_models(FE, FE_e, damage)%% Exact reference model
    stiffness_d = FE.element_properties.k;
    stiffness_d(damage(:, 1)) = stiffness_d(damage(:, 1)) .* damage(:, 2);
    mass = FE.element_properties.m;

    FE_d = ChainSystem();
    FE_d.assembly(stiffness_d, mass);
    FE_d.Cg = FE.modal_parameters.alpha * FE_d.Mg + FE.modal_parameters.beta * FE_d.Kg;
    FE_d.get_modal_parameters();

    %% Error damaged model
    stiffness_de = FE_e.element_properties.k;
    stiffness_de(damage(:, 1)) = stiffness_de(damage(:, 1)) .* damage(:, 2);
    FE_de = ChainSystem();
    FE_de.assembly(stiffness_de, mass);
    FE_de.Cg = FE.modal_parameters.alpha * FE_de.Mg + FE.modal_parameters.beta * FE_de.Kg;
    FE_de.get_modal_parameters();
    
    DamagedModels.FE_d = FE_d;
    DamagedModels.FE_de = FE_de;
    DamagedModels.Kg_d = FE_d.Kg;
    DamagedModels.Cg_d = FE_d.Cg;
    DamagedModels.Kg_de = FE_de.Kg;
    DamagedModels.Cg_de = FE_de.Cg;

end