function [yy] = fitness_vulclass(vulclass,nCol,nRaw,cooling_place)

% description
% fitness function of calculating the p65+ population in new green space

    vvulclass = reshape(vulclass,nCol,nRaw)';
    cooled_pop = vvulclass(cooling_place ~= 0);
    yy = sum(cooled_pop);

end


