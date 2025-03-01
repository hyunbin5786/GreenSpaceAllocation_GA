function [y] = fitness_allpop(allpop,nCol,nRaw,cooling_place)

% description
% fitness function of calculating the total population in new green space

    aallpop = reshape(allpop,nCol,nRaw)';
    cooled_pop = aallpop(cooling_place ~= 0);
    y = sum(cooled_pop);

end
