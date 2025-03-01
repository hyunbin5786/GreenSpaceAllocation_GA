function f = tournament_selection(pop_new2, f, pool_size, tour, tot)

% description
% function of doing torurnament selection

can = zeros(pool_size,tour);  % initialize
candidate=randperm(tot, tot); % make random array

% save random array
for i=1:tot
    can(i) = candidate(i);
end


% collect information about the selected candidates

% initialize
c_TotalCost = zeros(pool_size,tour);
max_candidate = zeros(pool_size,1);

for j = 1 : pool_size
    j
    c_TotalCost(j,1) = pop_new2(can(j,1)).TotalCost;
    c_TotalCost(j,2) = pop_new2(can(j,2)).TotalCost;
    target = can(j,:);

    % do tournament selection (choose the max total value)
    if c_TotalCost(j,1) == c_TotalCost(j,2) % if the candidates have same value
        max_candidate(j) = target(1);
    else
        max_candidate(j) = target(c_TotalCost(j,:) == max(c_TotalCost(j,:))); 
    end
    
    f(j) = pop_new2(max_candidate(j));
end

end

