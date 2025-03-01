function [areas] = extent(clus)

% description
% function of calculating the area of each cluster

    % get area (unit : ha)
    props = regionprops(clus, 'Area');
    
    % make empty matrix for saving
    numClusters = numel(props);
    areas = zeros(numClusters, 1);
    
    % save area value (unit : m^2)
    for i = 1:numClusters
        areas(i) = props(i).Area*10000;
    end
end
