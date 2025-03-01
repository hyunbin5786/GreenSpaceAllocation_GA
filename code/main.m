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

% Number of Objectives
nObj=3;


%% Parameters
MaxIt=100000; % Maximum Iterations
nPop=50; % Population Size

% Maximum green space cover ratio (0~1)
pActual=0.486;

% Normalization Rate
pp=1;

%%Weight of each objective
Weight_obj1=1/4; % objective1 : Maximization of Cooling
Weight_obj2=1/4; % objective2 : Maximization of Cooled POP
Weight_obj3=1/2; % objective3 : Maximization of Cooled P65+

% Global Max and Min (calculated at minmax_100.m)
max_obj1 = 1.975208976711999e+03;  % cooling 
min_obj1 = 3.367140000000000e+02;

max_obj2 = 119675;  % allpop
min_obj2 = 8273;

max_obj3 = 29693;  % P65+
min_obj3 = 4848;


%% load files (100m*100m)
lulc = readmatrix('lulc_100.txt');
lulc = reshape(lulc, 1, []);
nCount = nnz(lulc); % 0 = boundary
original= reshape(lulc, nRaw, []);

fixed = readmatrix('fixed_100.txt');
fixed = reshape(fixed, 1, []);
% 0 = can place new green space
% 1 = cannot place new green space

allpop = readmatrix('allpop_100.txt');
allpop = reshape(allpop, 1, []);

vulclass = readmatrix('vulclass_100.txt');
vulclass = reshape(vulclass, 1, []);

%% colormap
% 1 = urban
% 2 = agriculture
% 3 = green space
% 6 = bare land
% 7 = water
% 8 = new green space

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


% initiallization for repmat
empty_individual.Position=[];
empty_individual.Cost1=[];
empty_individual.Cost2=[];
empty_individual.Cost3=[];
empty_individual.NormCost=[];
empty_individual.TotalCost=[];

pop=repmat(empty_individual,nPop,1);

%% initial children
for k=1:nPop
    k
    % make 50 children pop
    [pop(k).Position] = initialize(nCount, fixed, lulc, pActual);

    q = pop(k).Position;
    q = reshape(q,nRaw, []);

    % calculate cluster
    clus = cluster(q);
    nClus = length(unique(clus))-1;  % number of cluster

    % calculate extent of neighbor
    areas = extent(clus);   % area / unit : m^2
    cool_d = find(clus~=0);  % direct cooling extent
    cool_d_ori = find(original == 3);  % original direct cooling

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
        [range(j).whole, range(j).indirect] = whole(j, clus, areas, nCol, nRaw, range(j).direct, cool_d);
    end

    % calculate fitness function
    [pop(k).Cost1, cooling_place] = fitness_temp(nClus, range, k, pop, nCol, nRaw, cool_d_ori, clus, original); % cooling
    pop(k).Cost2 = fitness_allpop(allpop, nCol, nRaw, cooling_place); % allpop
    pop(k).Cost3 = fitness_vulclass(vulclass, nCol, nRaw, cooling_place);  % P65+
end


% Cost
cost_obj1=[pop.Cost1];  % cooling
cost_obj2=[pop.Cost2];  % allpop
cost_obj3=[pop.Cost3];  % P65+


% Normalization
for i=1:nPop
    pop(i).NormCost=[Normal(cost_obj1(i), max_obj1, min_obj1, pp), Normal(cost_obj2(i), max_obj2, min_obj2, pp), Normal(cost_obj3(i), max_obj3, min_obj3, pp)]';
end

% Total Cost
for i=1:nPop
    pop(i).TotalCost=Weight_obj1*pop(i).NormCost(1,1)+Weight_obj2*pop(i).NormCost(2,1)+Weight_obj3*pop(i).NormCost(3,1);
end

%% start making children : crossover and mutation
empty_individuall.Position=[];
empty_individuall.Cost1=[];
empty_individuall.Cost2=[];
empty_individuall.Cost3=[];
empty_individuall.new_Norm=[];
empty_individuall.new_TotCost=[];

% matrix containing the value of the best children
Fitness_Info=repmat(empty_individuall,MaxIt,1);

consecutive_count = 0; % initialize consecutive count

for it=1:MaxIt
    tic
    it

    %% Crossover
     
    nCrossover = nPop*0.6; % 60% data is for crossover
    em.Position = [];
    popc = repmat(em, nCrossover, 1);  % matrix for crossover

    for s = 1:nCrossover/2  % we are going to input data in pairs
        random_num = datasample(1:nPop, 2, 'Replace', false); % pick 2 numbers with sampling without replacement
        p = random_num(1); % crossover candidate 1
        q = random_num(2); % crossover candidate 2
        c1 = pop(p).Position; 
        c2 = pop(q).Position;

        % create candidate positon for crossover
        c1_nh = datasample(find(c1 == 8 & fixed == 0 & lulc ~=3), round(numel(find(c1 == 8& fixed == 0))/5), 'Replace', false);
        c2_nh = datasample(find(c2 == 8 & fixed == 0 & lulc ~=3), round(numel(find(c2 == 8& fixed == 0))/5), 'Replace', false);
    
        c1(c1_nh(:)) = original(c1_nh(:));
        c2(c2_nh(:)) = original(c2_nh(:));
        c1(c2_nh(:)) = 8; 
        c2(c1_nh(:)) = 8; 
    
        % if green space exceeds nCount * pActual, delete it randomly
        limit_c1 = numel(find(c1 == 3 & c1 == 8)); 
        limit_c2 = numel(find(c2 == 3 & c2 == 8)); 
    
        if limit_c1 < nCount * pActual  % the case which doesn't exceed the limit
            popc(s+(s-1)).Position = c1;
        else
            c1_t = find(c1 == 8); % find new green space
            num_to_remove = max(0, limit_c1 - round(nCount * pActual)); % number of positions to delete
            if num_to_remove > 0
                c1_lo = datasample(c1_t, num_to_remove, 'Replace', false); % pick position to delete
                c1_v = c1;
                c1_v(c1_lo) = original(c1_lo); % make it to original landcover
                popc(s+(s-1)).Position = c1_v; % input data in first pair
            else
                popc(s+(s-1)).Position = c1;
            end
        end
    
        if limit_c2 < nCount * pActual
            popc(s+s).Position = c2;
        else
            c2_t = find(c2 == 8);
            num_to_remove = max(0, limit_c2 - round(nCount * pActual));
            if num_to_remove > 0
                c2_lo = datasample(c2_t, num_to_remove, 'Replace', false);
                c2_v = c2;
                c2_v(c2_lo) = original(c2_lo);
                popc(s+s).Position = c2_v;
            else
                popc(s+s).Position = c2;
            end
        end
    end

%% Mutation
nMutation=nPop*0.4; % 40% data is for mutation
popm=repmat(em,nMutation,1);
original_m= reshape(lulc, nRaw, []);  % matrix for crossover

     for k = 1:nMutation
         t=datasample(1:nPop,1); % pick one number from 1~nPop
         m1=pop(t).Position; % pick mutation candidate
         num_m = round(numel(find(m1 == 8 & fixed == 0 & lulc ~=3))*0.01); % calculate (mutable position number) * 0.01 and round it
         nMu_remove=datasample(find(m1 == 8 & fixed == 0 & lulc ~=3),num_m); % pick remove candidate
         nMu_new=datasample(find(m1 ~= 8 & fixed == 0 & lulc ~=3),num_m);  % pick new green space candiate
         m1(nMu_remove)=original_m(nMu_remove); m1(nMu_new)=8; % make remove candiate to original landcover and make new green space candidate to landcover 8
         popm(k).Position=m1; % put mutated landcover into mutation matrix
     end

    %% calculate the cost of new population
    pop_new=repmat(empty_individual,nPop,1);
    pop_change=[popc;popm]; % matrix with crossover and mutation candidate

    
    for k=1:nPop
        %k
        pop_new(k).Position=pop_change(k).Position;
        x = pop_new(k).Position;
        x = reshape(x,nCol,nRaw)';
    
        % calculate cluster
        clus = cluster(x);
        nClus = length(unique(clus))-1;
    
        % calculate extent of neighbor
        areas = extent(clus);    % area / unit : m^2
        cool_d = find(clus~=0);  % direct cooling extent
        cool_d_ori = find(clus ~= 0 & original == 3);  % original direct cooling

        % create range matrix
        range_new=repmat(empty,nClus,1);    
    
        % calculate range of each cluster (unit : m)
        for j=1:nClus
            range_new(j).area = areas(j); % unit m^2
            range_new(j).direct = find(clus==j);
            [range_new(j).whole, range_new(j).indirect] = whole(j, clus, areas, nCol, nRaw, range_new(j).direct, cool_d);  % indirect cooling이 없는 지역에는 0넣기
        end
    
        % calculate fitness function
        [pop_new(k).Cost1, cooling_place] = fitness_temp(nClus, range_new, k, pop_new, nCol, nRaw, cool_d_ori, clus, original);  % cooling
        pop_new(k).Cost2 = fitness_allpop(allpop, nCol, nRaw, cooling_place);  % allpop
        pop_new(k).Cost3 = fitness_vulclass(vulclass, nCol, nRaw, cooling_place);  % P65+
    end
    

pop_new2=[pop;pop_new];
cost_obj1=[pop_new2.Cost1];  % cooling
cost_obj2=[pop_new2.Cost2];  % allpop
cost_obj3=[pop_new2.Cost3];  % p65+

    
% Normalization
for i=1:nPop*2
    pop_new2(i).NormCost=[Normal(cost_obj1(i), max_obj1, min_obj1, pp), Normal(cost_obj2(i), max_obj2, min_obj2, pp), Normal(cost_obj3(i), max_obj3, min_obj3, pp)]';
end

% Total Cost
for i=1:nPop*2
    pop_new2(i).TotalCost=Weight_obj1*pop_new2(i).NormCost(1,1)+Weight_obj2*pop_new2(i).NormCost(2,1)+Weight_obj3*pop_new2(i).NormCost(3,1);
end



    %% Tournament selection
    % pool - size of the mating pool.
    % tour - Tournament size.
    pool_size = nPop;
    tour = 2;
    f=repmat(empty_individual, pool_size,1); 
    f=tournament_selection(pop_new2, f, pool_size, tour, nPop*2); % choose pop with higher total cost
    pop=f; % pop after tournament selection
    
    tt_Costs=[pop.TotalCost];
    
    % find the best position
    [~,b]=sort(tt_Costs,'descend');
    a=b(1); % save the best position as a

    % put best position data in Fitness_Info matrix
    Fitness_Info(it).Position = pop(a).Position;
    Fitness_Info(it).Cost1 = pop(a).Cost1;
    Fitness_Info(it).Cost2 = pop(a).Cost2;
    Fitness_Info(it).Cost3 = pop(a).Cost3;

    
    %% Plotting
    if it >= 2
        cost_obj1_fit=[Fitness_Info.Cost1];  % cooling
        cost_obj2_fit=[Fitness_Info.Cost2];  % allpop
        cost_obj3_fit=[Fitness_Info.Cost3];  % P65+
        
        
        % Normalization
        for i=1:it
            Fitness_Info(i).new_Norm=[Normal(cost_obj1_fit(i), max_obj1, min_obj1, pp), Normal(cost_obj2_fit(i), max_obj2, min_obj2, pp), Normal(cost_obj3_fit(i), max_obj3, min_obj3, pp)]';
        end
        
        % Total Cost
        for i=1:it
            Fitness_Info(i).new_TotCost=Weight_obj1*Fitness_Info(i).new_Norm(1,1)+Weight_obj2*Fitness_Info(i).new_Norm(2,1)+Weight_obj3*Fitness_Info(i).new_Norm(3,1);
        end
        
        % Draw figure
        figure(1);
        PlotCosts([Fitness_Info.new_TotCost],1:it);
        drawnow;

        cost_obj1_disp = round(Fitness_Info(it).Cost1,2);  % cooling
        cost_obj2_disp = round(Fitness_Info(it).Cost2,2);  % allpop
        cost_obj3_disp = round(Fitness_Info(it).Cost3,2);  % P65+

        


        %% Early stopping (for it>=100)
        current_total_cost = Fitness_Info(it).new_TotCost;
        pre = it - 1;
        previous_total_cost = Fitness_Info(pre).new_TotCost;
        
        % print the result
        if current_total_cost == previous_total_cost
            consecutive_count = consecutive_count + 1;
            disp(['횟수 : ' , num2str(it), '회, ', num2str(consecutive_count), '번 / 값 : ', num2str(cost_obj1_disp), ',', num2str(cost_obj2_disp), ',',num2str(cost_obj3_disp)]);
        else
            consecutive_count = 1;
        end
        
        
        % stop the loop if consecutive_count > it*0.2
        threshold = it * 0.2;
        if it > 1000
            if consecutive_count >= threshold
                disp(['루프를 중단합니다.']);
                break;
            end
            if it >= 5000  
                if mod(it,5000) == 0              
                    save('4-3_1117.mat', 'Fitness_Info', "-v7.3");  % myStruct를 myStruct.mat 파일로 저장
                    save('4-3_1117_pop.mat', 'pop', "-v7.3");
                    % load('1-1_0926.mat');  % myStruct.mat 파일 불러오기
                end
            end
        end
        
    end




toc
end

%% Result
% print the optimized plan
figure(7);
optimum=Fitness_Info(it).Position;
optimum(optimum==0)=7; % change background to waterbody
optimum_ = reshape(optimum, nRaw, []);
imagesc(optimum_);
colormap(lulc_op);
title('the result of optimization');

% print the original LULC
figure(8);
lulc_= reshape(lulc, nRaw, []);
lulc_(lulc_==0)=7; % change background to waterbody
imagesc(lulc_);
colormap(lulc_raw);
title('current LULC');

%% submit the result
% submit optimized plan
data = reshape(Fitness_Info(it).Position, nRaw, []);  % Convert data to 115x90
fid = fopen('4-3_1117.txt', 'w');  % Open the file in writing mode
for row = 1:size(data, 1)
    fprintf(fid, '%i ', data(row, :));  % write by rows
    fprintf(fid, '\n');  % add line break at the end of the row
end
fclose(fid);  % close the file


% submit Fitness_Info matrix
save('4-3_1117.mat', 'Fitness_Info', "-v7.3");
% load('1-1_0926.mat');  % loading code

