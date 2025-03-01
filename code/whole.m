function [whole, indirect] = whole(jj, clus, areas, nCol, nRaw, direct, cool_d)

% description
% function of calculating the whole extent and indirect cooling extent

    area_ha = areas(jj)/10000; % unit : ha
    [xGrid, yGrid] = meshgrid(1:nCol, 1:nRaw);  % x = col, y = row
    
    %% whole extent
    if area_ha >= 20
       
        % make empty matrix
        edge_cell = [];

        % calculate radius
        rad = floor(sqrt(area_ha/(2*pi)));  % eqaution from Tan and Li(2015)

        % get index of direct cooling     
        [Rows, Cols] = find(clus == jj);
                
        %% indirect cooling 범위 구하기
        for ii = 1:length(Rows)
            r = Rows(ii);
            c = Cols(ii);
            if isEdgeCell(c, r, nRaw, nCol, jj, clus) % verify if it is edge or not           
                circle = (xGrid - c).^2 + (yGrid - r).^2 <= rad^2; % draw a circle in the edge cell and save areas within the circle
                [rValues, cValues] = find(circle==1);
                for iii = 1:length(rValues)
                    linear_index = sub2ind([nRaw, nCol], rValues(iii), cValues(iii));
                    edge_cell(end+1,1) = linear_index;
                end
            end
        end

        % refine values
        indirect = unique(edge_cell); % remove non-unique values
        indirect = setdiff(indirect, cool_d);  % remove direct cooling area
        whole = [direct; indirect];

    else  % area_ha < 20 / 20 = indirect cooling threshold
        whole = direct;
        indirect = [];
    end

end

