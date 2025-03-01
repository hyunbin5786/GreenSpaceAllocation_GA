function q = Normal( p, maxValue, minValue, pp )

% description
% function of normallization of fitness value
% calculation is done based on the maximum and minimum value in the iteration

normalization=(p-minValue)/(maxValue-minValue);
q=normalization.^pp;

end

