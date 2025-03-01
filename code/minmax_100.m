clc;
clear;
close all;

%% Problem Definition
% Defining the cell size
% nCol = number of Columns
% nRaw = number of Rows
nCol=80;
nRaw=115;
nVar=nCol*nRaw;
VarSize=[1 nVar];

%% Parameters
MaxIt=100000; % Maximum Iterations
nPop=50; % Population Size

% Maximum green space cover ratio (0~1)
pActual=0.486;

% Normalization Rate
pp=1;

%% load files (100m*100m)
lulc = readmatrix('lulc_100.txt');
lulc = reshape(lulc, 1, []);
nCount = nnz(lulc); % 0 = boundary

fixed = readmatrix('fixed_100.txt');
fixed = reshape(fixed, 1, []);
% 0 = can place new green space
% 1 = cannot place new green space

allpop = readmatrix('allpop_100.txt');
allpop = reshape(allpop, 1, []);

vulclass = readmatrix('vulclass_100.txt');
vulclass = reshape(vulclass, 1, []);

%% find addable new green space
actual_3 = find(lulc == 3);
nactual_3 = numel(actual_3);
add_3 = round(nCount * pActual) - nactual_3;

%% find min&max : allpop
% find candidate
conditionIndices = find(lulc(:) ~= 3 & fixed(:) == 0);
allpop_sort = allpop(conditionIndices)';
allpop_mat = [allpop_sort, conditionIndices];

% count fixed=0 cell
zzeros = allpop_mat(allpop_mat (:, 1) == 0, :);

% calculate max
allpop_des = sortrows(allpop_mat, 1, 'descend');  % sort in descending order
sorted_pop_max = allpop_des(1:add_3, :);  % choose top add_3 cells
sum_pop_max = sum(sorted_pop_max(:,1));  % sum

% calculate min
allpop_asc = sortrows(allpop_mat, 1, 'ascend');  % sort in ascending order
allpop_asc(allpop_asc(:, 1) == 0, :) = [];  % make 0 population cell into blank 
sorted_pop_min = allpop_asc(1:add_3, :);  % choose top add_3 cells
sum_pop_min = sum(sorted_pop_min(:,1));  % sum

%% find min&max : P65+
% find candidate
conditionIndices = find(lulc(:) ~= 3 & fixed(:) == 0);
vulclass_sort = vulclass(conditionIndices)';
vulclass_mat = [vulclass_sort, conditionIndices];

% count fixed=0 cell
zzeros = vulclass_mat(vulclass_mat (:, 1) == 0, :);

% calculate max
vulclass_des = sortrows(vulclass_mat, 1, 'descend');  % sort in descending order
sorted_vul_max = vulclass_des(1:add_3, :);  % choose top add_3 cells
sum_vul_max = sum(sorted_vul_max(:,1));  % sum

% calculate min
vulclass_asc = sortrows(vulclass_mat, 1, 'ascend');  % sort in ascending order
vulclass_asc(vulclass_asc(:, 1) == 0, :) = [];  % make 0 population cell into blank 
sorted_vul_min = vulclass_asc(1:add_3, :);  % choose top add_3 cells
sum_vul_min = sum(sorted_vul_min(:,1));  % sum


%% find max : cooling
% 1. allocate cells for max cooling
% 1-1. reshape data in ascii data format
lulc_ori = reshape(lulc, nRaw, []);
llulc = reshape(lulc, nRaw, []);
ffixed = reshape(fixed, nRaw, []);

llulc(llulc ~= 3 & ffixed == 0) = 10; % designate coolable poisition to 10
neighbour_3_count = zeros(nRaw, nCol);

lulc10 = find(llulc(:) == 10);


% 1-2. for lulc = 10, count lulc=3 by queen neighborhood method
for i = 1:nRaw
    for j = 1:nCol
        if llulc(i, j) == 10
            % calculate in queen neighborhood method
            for di = -1:1
                for dj = -1:1
                    if di == 0 && dj == 0  % skip self
                        continue;
                    end
                    ni = i + di;
                    nj = j + dj;
                    if ni >= 1 && ni <= nRaw && nj >= 1 && nj <= nCol
                        if llulc(ni, nj) == 3
                            neighbour_3_count(i, j) = neighbour_3_count(i, j) + 1;
                        end
                    end
                end
            end
        end
    end
end

% 1-3. input area value into cell which has nearby cluster
% allocate new value in position
lllulc = llulc;
lllulc(lllulc == 3) =  100 ; % allocate 100 in lulc=3 place 
lllulc(lllulc ~= 100) = 0; % designate other places as 0

% calculate cluster
clus = bwlabel(lllulc,8);
nClus = length(unique(clus))-1;

% calculate extent of neighbor
areas = extent(clus);  % unit : m^2

% put area value in clus matrix
for i=1:nClus
    aa = find(clus == i);
    clus(aa) = areas(i);
end

% calculate area of cluster in queen neighborhood method
for i = 1:nRaw
    for j = 1:nCol
        if neighbour_3_count(i, j) ~= 0  % for place where there is lulc=3 in queen neighborhood method
            % calculate in queen neighborhood method
            for di = -1:1
                for dj = -1:1
                    if di == 0 && dj == 0 % skip self
                        continue;
                    end
                    ni = i + di;
                    nj = j + dj;
                    if ni >= 1 && ni <= nRaw && nj >= 1 && nj <= nCol
                        if clus(ni, nj) ~= 0
                            clus(i, j) = clus(ni, nj);  % put the area of cluster value into clus matrix
                            break
                        end
                    end
                end
            end
        end
    end
end

% initial update : allocate 8(new green space) for maximum value
new_mat = clus.*neighbour_3_count; % make new matrix(considering area and nearby 3 count) with element-wise multiplication
max_neighbour_3 = max(new_mat(:));
idx_new = find(new_mat == max_neighbour_3);
llulc(idx_new) = 8;

% 1-4. for lulc = 10, count lulc=3 or 8 by queen neighborhood method
while true
    neighbour_3_8_count = zeros(nRaw, nCol);
    for i = 1:nRaw
        for j = 1:nCol
            if llulc(i, j) == 10
                % calculate in queen neighborhood method
                for di = -1:1
                    for dj = -1:1
                        if di == 0 && dj == 0  % skip self
                            continue;
                        end
                        ni = i + di;
                        nj = j + dj;
                        if ni >= 1 && ni <= nRaw && nj >= 1 && nj <= nCol
                            % add count for place with lulc = 3 or 8
                            if llulc(ni, nj) == 3 || llulc(ni, nj) == 8
                                neighbour_3_8_count(i, j) = neighbour_3_8_count(i, j) + 1;
                            end
                        end
                    end
                end
            end
        end
    end
    
    % 1-5. Allocation
    lllulc = llulc;
    lllulc(lllulc == 3 | lllulc == 8) =  100 ; % allcoate 100 for lulc = 3 or 8
    lllulc(lllulc ~= 100) = 0; % other places to 0

    % calculate cluster
    clus = bwlabel(lllulc,8);
    nClus = length(unique(clus))-1;
    
    % calculate extent of neighbor
    areas = extent(clus); % unit : m^2
    
    % put area value in clus matrix
    for i=1:nClus
        aa = find(clus == i);
        clus(aa) = areas(i);
    end

    % input area value in cells in queen neighborhood method
    for i = 1:nRaw
        for j = 1:nCol
            if neighbour_3_8_count(i, j) ~= 0
                % calculate in queen neighborhood method
                for di = -1:1
                    for dj = -1:1
                        if di == 0 && dj == 0  % skip self
                            continue;
                        end
                        ni = i + di;
                        nj = j + dj;
                        if ni >= 1 && ni <= nRaw && nj >= 1 && nj <= nCol
                            if clus(ni, nj) ~= 0
                                clus(i, j) = clus(ni, nj);
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    
    % update : allocate 8(new green space) for maximum value
    new_mat = clus.*neighbour_3_8_count; % make new matrix(considering area and nearby 3 | 8 count) with element-wise multiplication
    max_neighbour_38 = max(new_mat(:));
    idx_new = find(new_mat == max_neighbour_38);
    llulc(idx_new) = 8;
    
    % 1-6. final calculation
    % lulc = 8 If the number of cells is less than add_3, go to step 1-4
    % or end the loop if it is greater or equal

    num_8_cells = sum(llulc(:) == 8);
    % disp(num_8_cells)
    if num_8_cells >= add_3
        index_10 = find(llulc == 10);
        llulc(index_10) = lulc_ori(index_10);
        break;
    end
end


% 2. calculate cooling effect
% 2-1. make a cluster
% calculate cluster
clus = cluster(llulc);
nClus = length(unique(clus))-1;

% calculate extent of neighbor
areas = extent(clus);   % area 단위 : m^2
cool_d = find(clus~=0);  % direct cooling에 속하는 모든 영역
cool_d_ori = find(lulc_ori == 3);  % 원래 direct cooling에 속하는 모든 영역

% range matrix
empty.area=[];
empty.direct=[];
empty.indirect = [];
empty.whole = [];
range=repmat(empty,nClus,1);

 % calculate range of each cluster (unit : m)
for j=1:nClus
    range(j).area = areas(j);  % unit m^2
    range(j).direct = find(clus==j);
    [range(j).whole, range(j).indirect] = whole(j, clus, areas, nCol, nRaw, range(j).direct, cool_d);  % indirect cooling이 없는 지역에는 0넣기
end

% 2-2. calculate fitness
% create empty matrix
sum_dtd = zeros([nRaw nCol]);
sum_dti = zeros([nRaw nCol]);
count_indir = zeros([nRaw nCol]);

% find lulc = 8 area
l = llulc;
new = find(l==8);

% make cartesian grid
[xGrid, yGrid] = meshgrid(1:nCol, 1:nRaw);  % x = col, y = row

% find original Clus
clus_ori = zeros([nRaw nCol]);
clus_ori_idx = find(clus ~= 0 & lulc_ori == 3);  % original clus index
for i = 1:length(clus_ori_idx)
    idx = clus_ori_idx(i);
    clus_ori(idx) = clus(idx);
end
    
for j=1:nClus
    d = range(j).direct;
    i = range(j).indirect;

    % 2-2-1. Direct Cooling
    area = length(d); % cluster area (unit : ha)

    % New Cooling Effect
    D_new = intersect(d, new);  % new cooling area (lulc = 8 & direct cooling)

    if ~isempty(D_new)  % if there is new cluster
        % new cooling effect
        nD_new = length(D_new);  % new direct cooling area
        dt_dnew = 1.6034+0.8560*log(area); % new direct cooling effect
    
        % original cooling effect
        D_ori = setdiff(d, D_new);
        nD_ori = area-nD_new;   % nD_ori = original direct cooling area (unit : ha)
        dt_dori = 1.6034+0.8560*log(nD_ori);  % cooling effect of original direct cooling area
        ddt_dori = dt_dnew-dt_dori; % additional cooling effect by new green space
        
    
        % save dt
        sum_dtd(D_new) = dt_dnew;
        sum_dtd(D_ori) = ddt_dori;
    

        % 2-2-2. Indirect Cooling
        if nD_ori >= 20

            % make empty matrix for saving
            edge_cell = [];
    
            % calculate radius
            rad = floor(sqrt(nD_ori/(2*pi)));  % eqaution from Tan and Li(2015)
    
            % find index of original direct cooling
            [Rows, Cols] = find(clus_ori == j);
                    
            % find the edge cell of the cluster
            for ii = 1:length(Rows)
                r = Rows(ii);
                c = Cols(ii);
                if isEdgeCell(c, r, nRaw, nCol, j, clus_ori) % verify whether it is edge or not               
                    circle = (xGrid - c).^2 + (yGrid - r).^2 <= rad^2; % draw a circle in the edge cell and save areas within the circle
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
       
        % 2-2-2-1. calculate new indirect cooling
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

% 2-2-3 multi cooling
multi_cool = find(count_indir > 1);  % get index of multiple cooling

for i = 1:length(multi_cool)
    idx = multi_cool(i);
    sum_dti(idx) = sum_dti(idx)/count_indir(idx);  % arithmetic mean
end

% 2-2-4. new cooling calculation
sum_dt = sum_dtd + sum_dti;
sum_temp_max = sum(sum_dt(:));


%% find min : cooling
% assume that all new green space is seperated
nD_new = add_3;
dt_d = 1.6034+0.8560*log(1); % unit : ha
sum_temp_min = dt_d*nD_new;

%% Visulaization
% color map
lulc_raw = uint8([255 0 0
    238 233 7
    42 75 45
    42 75 45
    6 2 250
    89 206 202
    6 2 250]);

lulc_op = uint8([255 0 0
    238 233 7
    42 75 45
    42 75 45
    6 2 250
    89 206 202
    6 2 250
    57 153 38]);

% visualize lulc
figure(1);
lulc_= reshape(lulc, nRaw, []);
lulc_(lulc_==0)=7; %change background to waterbody
imagesc(lulc_);
colormap(lulc_raw);
title('lulc');

figure(2);
lulc_= llulc;
lulc_(lulc_==0)=7; %change background to waterbody
imagesc(lulc_);
colormap(lulc_op);
title('lulc modified');