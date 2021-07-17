function Mxyz = massmatrix_xyz(V,F)
    M = massmatrix(V,F);
    Meles = diag(M);
    % 2D
    if size(V,2) == 2
        M2eles = zeros(2*size(M,1),1);
        M2eles(2*(1:size(M,1))) = Meles;
        M2eles(2*(1:size(M,1))-1) = Meles;
        Mxyz = sparse(1:2*size(M,1),1:2*size(M,1),M2eles,2*size(M,1),2*size(M,1));
    % 3D
    elseif size(V,2) == 3
        M3eles = zeros(3*size(M,1),1);
        M3eles(3*(1:size(M,1))) = Meles;
        M3eles(3*(1:size(M,1))-1) = Meles;
        M3eles(3*(1:size(M,1))-2) = Meles;
        Mxyz = sparse(1:3*size(M,1),1:3*size(M,1),M3eles,3*size(M,1),3*size(M,1));
    end
end