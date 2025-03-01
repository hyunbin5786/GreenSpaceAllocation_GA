function [y] = cluster(q)

% description
% function for definig clusters
% cluster is define as green space (origianl green space + new green space)

    q(q == 3 | q == 8) = 10; 
    q(q ~= 10) = 0;
    y = bwlabel(q,8);
    
end

