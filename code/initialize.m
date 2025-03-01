function [lulc_mod] = initialize(nCount, fixed, lulc, pActual)

% description
% function of making initial population

    lulc_mod = lulc;

    % find addable new green space
    actual_3 = find(lulc == 3);
    nactual_3 = numel(actual_3);
    add_3 = round(nCount * pActual) - nactual_3;

    % make new green space
    candidate = find(lulc ~= 3 & fixed == 0);
    d = randperm(numel(candidate), add_3);
    rand_idx = candidate(d);  % pick candidate of new green space

    % make new green space lulc into 8
    for i = 1:add_3
        idx = rand_idx(i);
        lulc_mod(idx) = 8; % 8 입력
    end


end
