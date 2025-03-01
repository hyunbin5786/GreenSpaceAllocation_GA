function [yy, sum_dt] = fitness_temp(nClus, range, k, pop, nCol, nRaw, cool_d_ori, clus, original)

% description
% fitness function of calculating the cooling effect created by new green space

% make empty matrix for saving
sum_dtd = zeros([nRaw nCol]);
sum_dti = zeros([nRaw nCol]);
count_indir = zeros([nRaw nCol]);

% find new green space
l = pop(k).Position;
new = find(l==8);

% create cartesian grid
[xGrid, yGrid] = meshgrid(1:nCol, 1:nRaw);  % x = col, y = row

% find original cluster
clus_ori = zeros([nRaw nCol]);
clus_ori_idx = find(clus ~= 0 & original == 3);  % index of original cluster

% save original cluster index
for i = 1:length(clus_ori_idx)
    idx = clus_ori_idx(i);
    clus_ori(idx) = clus(idx);
end
   
% start calculating cooling effect
for j=1:nClus
    d = range(j).direct;
    i = range(j).indirect;

    %% 01 Direct Cooling
    area = length(d); % cluster area (unit : ha)

    % New Cooling Effect
    D_new = intersect(d, new);  % new cooling area (lulc = 8 & direct cooling)

    if ~isempty(D_new)  % if there is new cluster
        % new cooling effect
        nD_new = length(D_new);  % count new direct cooling
        dt_dnew = 1.6034+0.8560*log(area); % calculate new direct cooling effect
    
        % original cooling effect
        D_ori = setdiff(d, D_new);
        nD_ori = area-nD_new;  % nD_ori = original direct cooling area (unit : ha)
        dt_dori = 1.6034+0.8560*log(nD_ori); % cooling effect of original direct cooling area
        ddt_dori = dt_dnew-dt_dori; % additional cooling effect by new green space
        
    
        % save dt
        sum_dtd(D_new) = dt_dnew;
        sum_dtd(D_ori) = ddt_dori;
    

        %% 02 Indirect Cooling
        if nD_ori >= 20

            % make empty matrix for saving
            edge_cell = [];
    
            % calculate radius
            rad = floor(sqrt(nD_ori/(2*pi)));  % 내림, Tan and Li(2015)식 사용
    
            % find index of original direct cooling
            [Rows, Cols] = find(clus_ori == j);
                    
            % find the edge cell of the cluster
            for ii = 1:length(Rows)
                r = Rows(ii);
                c = Cols(ii);
                if isEdgeCell(c, r, nRaw, nCol, j, clus_ori) % verify whether it is edge or no              
                    circle = (xGrid - c).^2 + (yGrid - r).^2 <= rad^2;  % draw a circle in the edge cell and save areas within the circle
                    [rValues, cValues] = find(circle==1);
                    for iii = 1:length(rValues)
                        linear_index = sub2ind([nRaw, nCol], rValues(iii), cValues(iii));
                        edge_cell(end+1,1) = linear_index;
                    end
                end
            end

        % refining values
        I_ori = unique(edge_cell); % remove non-unique values
        I_ori = setdiff(I_ori, cool_d_ori);  % remove direct cooling area
    
        else  % nD < 20 / 20 = indirect cooling threshold
            I_ori = [];
        end
       
        %% 2-1. calculate new indirect cooling
        % find new neighbor
        I_new = setdiff(i, I_ori);  % new cooling area
        if ~isempty(I_new)
            % new cooling extent
            dt_inew = dt_dnew-(0.8284*log(area)+1.4030); % new indirect cooling dt           
            dt_iori = dt_dori-(0.8284*log(nD_ori)+1.4030);
            ddt_iori = dt_inew - dt_iori;  % dt from original indirect cooling

            % save dt
            sum_dti(I_new) = sum_dti(I_new) + dt_inew;
            sum_dti(I_ori) = sum_dti(I_ori) + ddt_iori;

            % count indirect cooling cell
            count_indir(I_new) = count_indir(I_new) + 1;
            count_indir(I_ori) = count_indir(I_ori) + 1;
           
        end
        
    end

end

%% 03 multi cooling
multi_cool = find(count_indir > 1);  % get index of multiple cooling

for i = 1:length(multi_cool)
    idx = multi_cool(i);
    sum_dti(idx) = sum_dti(idx)/count_indir(idx);  % arithmetic mean
end

%% 04 new cooling
sum_dt = sum_dtd + sum_dti;
yy = sum(sum_dt(:));

end