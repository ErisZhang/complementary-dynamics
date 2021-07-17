% xy xy xy

% S is sparse
function S = build_selection_matrix(Uc,Bd_idx)
    Bd_idx = reshape(Bd_idx,size(Bd_idx,1)*size(Bd_idx,2),1);
    I = zeros(size(Bd_idx));
    J = zeros(size(Bd_idx));
    V = zeros(size(Bd_idx));
    for i = 1:size(Bd_idx,1) 
        I(i) = i;
        J(i) = Bd_idx(i);
        V(i) = 1;
    end
    S = sparse(I,J,V,size(Bd_idx,1),size(Uc,1));
end