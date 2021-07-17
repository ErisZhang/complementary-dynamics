function [M] = lbs_matrix_xyz(V,W)
    
    n = size(V,1); % num of vertices
    m = size(W,2); % num of handles
    dim = size(V,2);
    
    rows = zeros(dim*(dim+1)*m*n,1);
    cols = zeros(dim*(dim+1)*m*n,1);
    vals = zeros(dim*(dim+1)*m*n,1);
    
    cur = 1;
    
    for x = 1:dim
        for j = 1:n
            for i = 1:m
                for c = 1:(dim+1)
                    value = W(j,i);
                    if c ~= (dim+1)
                        value = value*V(j,c);
                    end
%                     % xxx yyy zzz
%                     rows(cur) = (x-1)*n + j;
%                     cols(cur) = (x-1)*m + (c-1)*m*dim + i;
%                     vals(cur) = value;

                    % xyz xyz
                    rows(cur) = dim*(j-1)+x;
                    cols(cur) = (x-1)*m + (c-1)*m*dim + i;
                    vals(cur) = value;
                    
                    cur = cur + 1;
                end
            end
        end
    end
    
    M = sparse(rows,cols,vals,n*dim,m*dim*(dim+1));

end