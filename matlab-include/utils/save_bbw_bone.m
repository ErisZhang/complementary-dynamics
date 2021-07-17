function [W]=save_bbw_bone(filename, V,F,C,P,BE)

    [NV,NF] = remesh_at_handles(V,F,C,P,BE,[]);
    [b,bc] = boundary_conditions(NV,NF,C,P,BE);
    W = biharmonic_bounded(NV,NF,b,bc,'OptType','quad');
    I=snap_points(V,NV);
    W=W(I,:);
    W = W./repmat(sum(W,2),1,size(W,2));
    writeDMAT(filename, W);
    
end