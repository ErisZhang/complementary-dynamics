% V vertex postions
% C cluster postions
% I indices of each vertex's closest cluster

function [I] = cluster_pnts_to_closest(V,C) 
    I = zeros(size(V,1),1);
    for i = 1:size(V,1)
        v = V(i,:);
        % euclidean distance
        D = sum((C-v).^2,2);
        idx = find(D == min(D));
        I(i) = idx(1);
    end
end