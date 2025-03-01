function isEdge = isEdgeCell(x, y, nRaw, nCol, jj, clus)

% description
% function of defining whether it is an edge cell

        % initialize
        isEdge = false;
        
        % find in rook neighborhood method
        directions = [-1 0; 1 0; 0 -1; 0 1];

        for d = 1:size(directions, 1)
            nxx = x + directions(d, 1);
            nyy = y + directions(d, 2);
            
            if nxx > 0 && nxx <= nCol && nyy > 0 && nyy <= nRaw
                if clus(nyy, nxx) ~= jj  % check whether it has different value = it is in different cluster
                    isEdge = true;
                    break;
                end
            end
        end
    end
