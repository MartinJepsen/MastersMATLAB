%% Exact model
FE = FiniteElementModel();
FE.from_xlsx('structures/paper_truss.xlsx')
FE.assembly('bar', dam);
FE.apply_bc([1, 2, 9, 10]);
FE.modal_damping(0.02);
FE.strains_from_disp([])

n_dof = FE.n_dof;
free_dof = size(FE.Kg, 1);
n_el = size(FE.mesh.topology, 1);
B_strain = FE.B;
bc = FE.mesh.bc_dof;
Kg = FE.Kg;
Kg_d = FE.Kg_d;
Cg = FE.Cg;
Cg_d = FE.Cg_d;
Mg = FE.Mg;

%% Error model
FE_e = FiniteElementModel();
FE_e.from_xlsx('structures/paper_truss.xlsx')
FE_e.mesh.element_properties.E = FE.mesh.element_properties.E .* unifrnd(1-err, 1+err, [n_el, 1]);
FE_e.assembly('bar', dam);
FE_e.apply_bc([1, 2, 9, 10]);
FE_e.modal_damping(0.02);

Kg_e = FE_e.Kg;
Kg_de = FE_e.Kg_d;
Cg_e = FE_e.Cg;
Cg_de = FE_e.Cg_d;
Mg_e = FE_e.Mg;