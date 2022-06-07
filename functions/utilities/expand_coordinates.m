function GeneralParameters = expand_coordinates(GeneralParameters)

    O = union(GeneralParameters.in_dof, GeneralParameters.out_dof);
    r = numel(O);
    m = numel(O);
    in_dof = O;
    out_dof = O;
    free_dof = GeneralParameters.free_dof;
    B2 = zeros(free_dof, numel(in_dof));
    if isempty(in_dof) == 0
        for i = 1:numel(in_dof)
            B2(in_dof(i), i)=1;
        end
    end
    
    cdis = zeros(m, free_dof);
    for ii=1:m
        cdis(ii, out_dof(ii)) = 1;
    end
    
    GeneralParameters.in_dof = in_dof;
    GeneralParameters.out_dof = out_dof;
    GeneralParameters.cdis = cdis;
    GeneralParameters.B2 = B2;
    GeneralParameters.r = r;
    GeneralParameters.m = m;

end
