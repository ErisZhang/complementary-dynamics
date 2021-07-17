% phi default fudge factor to use
% Em indices of the points that are controlled by the orthogonality
% constraint
function [phi,Em] = default_D_matrix(V,T)

    dim = size(V,2);

    if dim == 2

        L = cotmatrix(V,T);
        Mr = massmatrix(V,T);
        Q = -L;
        l = Mr * ones(size(V,1),1);
        [phi_e,~] = min_quad_with_fixed(Q,l,unique(outline(T)),zeros(size(unique(outline(T)),1),1),[],[],[]);
        phi = zeros(size(V,1)*size(V,2));
        for i = 1:size(phi_e)
            phi(2*i,2*i) = phi_e(i);
            phi(2*i-1,2*i-1) = phi_e(i);
        end
        phi = sparse(phi);
        Em = find(phi_e ~= 0);

    elseif dim == 3

        F = boundary_faces(T);
        L = cotmatrix(V,T);
        Mr = massmatrix(V,T);
        Q = -L;
        l = Mr * ones(size(V,1),1);
        [phi_e,~] = min_quad_with_fixed(Q,l,unique(F),zeros(size(unique(F),1),1),[],[],[]);

        % phi = zeros(size(V,1)*size(V,2));
        % for i = 1:size(phi_e)
        %     phi(3*i,3*i) = phi_e(i);
        %     phi(3*i-1,3*i-1) = phi_e(i);
        %     phi(3*i-2,3*i-2) = phi_e(i);
        % end
        % phi = sparse(phi);
        % Em = find(phi_e ~= 0);

        phiI = [3*(1:size(phi_e,1))';3*(1:size(phi_e,1))'-1;3*(1:size(phi_e,1))'-2];
        phiJ = phiI;
        phiV = repmat(phi_e,3,1);
        phi = sparse(phiI,phiJ,phiV);
        Em = find(phi_e ~= 0);
        
    end

end