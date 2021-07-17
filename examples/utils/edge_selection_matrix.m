% xxx yyy

% S is sparse
function S = edge_selection_matrix(V,E)

    dim = size(V,2);

    if dim == 2
        num_edges = size(E,1);
        num_vertices = size(V,1);

        % each edge introduces 4 entries to S
        I = zeros(2*2*num_edges,1);
        J = zeros(2*2*num_edges,1);
        Va = zeros(2*2*num_edges,1);

        for i = 1:num_edges

            I(4*(i-1)+1) = i;
            I(4*(i-1)+2) = i;
            I(4*(i-1)+3) = i+num_edges;
            I(4*(i-1)+4) = i+num_edges;

            J(4*(i-1)+1) = E(i,1);
            J(4*(i-1)+2) = E(i,2);
            J(4*(i-1)+3) = E(i,1)+num_vertices;
            J(4*(i-1)+4) = E(i,2)+num_vertices;

            Va(4*(i-1)+1) = 1;
            Va(4*(i-1)+2) = -1;
            Va(4*(i-1)+3) = 1;
            Va(4*(i-1)+4) = -1;

        end
        S = sparse(I,J,Va,2*num_edges,2*num_vertices);
        
    elseif dim == 3
        
        num_edges = size(E,1);
        num_vertices = size(V,1);

        % each edge introduces 6 entries to S
        I = zeros(3*2*num_edges,1);
        J = zeros(3*2*num_edges,1);
        Va = zeros(3*2*num_edges,1);

        for i = 1:num_edges

            I(6*(i-1)+1) = i;
            I(6*(i-1)+2) = i;
            I(6*(i-1)+3) = i+num_edges;
            I(6*(i-1)+4) = i+num_edges;
            I(6*(i-1)+5) = i+num_edges*2;
            I(6*(i-1)+6) = i+num_edges*2;

            J(6*(i-1)+1) = E(i,1);
            J(6*(i-1)+2) = E(i,2);
            J(6*(i-1)+3) = E(i,1)+num_vertices;
            J(6*(i-1)+4) = E(i,2)+num_vertices;
            J(6*(i-1)+5) = E(i,1)+num_vertices*2;
            J(6*(i-1)+6) = E(i,2)+num_vertices*2;

            Va(6*(i-1)+1) = 1;
            Va(6*(i-1)+2) = -1;
            Va(6*(i-1)+3) = 1;
            Va(6*(i-1)+4) = -1;
            Va(6*(i-1)+5) = 1;
            Va(6*(i-1)+6) = -1;

        end
        
        S = sparse(I,J,Va,3*num_edges,3*num_vertices);
        
        
    end
        
end