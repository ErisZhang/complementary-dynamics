function [W]=save_bbw_pnt(filename, V,F,C)

    [NF] = remesh_at_points(V,F,C);
    NV = [V;C];
    [b,bc] = boundary_conditions(NV,NF,C);
    W = biharmonic_bounded(NV,NF,b,bc,'OptType','quad');
    I=snap_points(V,NV);
    W=W(I,:);
    W = W./repmat(sum(W,2),1,size(W,2));
    writeDMAT(filename, W);
    
end