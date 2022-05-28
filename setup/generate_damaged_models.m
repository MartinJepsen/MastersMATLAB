
function [DamagedModels] = generate_damaged_models(FE, FE_e, damage)%% Exact reference model
    FE_d = FiniteElementModel();
    FE_d.from_xlsx('structures/paper_truss.xlsx');
    FE_d.mesh.element_properties.E(damage(:, 1)) = FE.mesh.element_properties.E(damage(:, 1)) * damage(:, 2);
    FE_d.assembly("bar", [1, 2, 9, 10]);
    FE_d.modal_damping(0.02);

    %% Error damaged model
    FE_de = FiniteElementModel();
    FE_de.from_xlsx('structures/paper_truss.xlsx');
    FE_de.mesh.element_properties.E(damage(:, 1)) = FE_e.mesh.element_properties.E(damage(:, 1)) * damage(:, 2);
    FE_de.assembly("bar", [1, 2, 9, 10]);
    FE_de.modal_damping(0.02);
    
    DamagedModels.FE_d = FE_d;
    DamagedModels.FE_de = FE_de;
    DamagedModels.Kg_d = FE_d.Kg;
    DamagedModels.Cg_d = FE_d.Cg;
    DamagedModels.Kg_de = FE_e.Kg;
    DamagedModels.Cg_de = FE_e.Cg;

end