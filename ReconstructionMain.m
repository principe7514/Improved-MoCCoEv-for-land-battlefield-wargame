clc
clear
close all

warning('off','MATLAB:mesh');
%% % 定义一全局初值
%武器对车辆的攻击等级
Vrank=[10,10,10,10,10,9,9,9,9,8,8,8,7,7,6,5,4,3,2,0,0];
%攻击等级高度差修正
Amod=[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
      -2,-2,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
      -2,-2,-2,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
      -3,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
      -3,-3,-3,-2,-2,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
      -4,-3,-3,-3,-2,-2,-2,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0;
      -4,-4,-4,-3,-3,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0;
      -5,-4,-4,-4,-3,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0;
      -5,-5,-5,-4,-4,-3,-3,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1, 0];
%车辆战斗结果
VRes=[1/36, 6/36, 9/36,12/36,15/36,18/36,24/36,30/36,39/36,42/36;
      3/36, 9/36,15/36,18/36,24/36,30/36,36/36,42/36,51/36,58/36;
      6/36,12/36,18/36,24/36,36/36,42/36,51/36,58/36,66/36,72/36;
      6/36,12/36,18/36,30/36,39/36,51/36,58/36,66/36,78/36,91/36;
      6/36,12/36,30/36,39/36,58/36,66/36,72/36,91/36,101/36,115/36];
%车辆战损修正
Vmod=[-12/36,  6/36,14/36,33/36,42/36;
      -17/36,  3/36, 7/36,24/36,36/36;
      -28/36, -2/36, 3/36,15/36,31/36;
      -39/36, -7/36, 0/36,10/36,26/36;
      -51/36,-13/36,-3/36, 6/36,21/36;
      -64/36,-22/36,-7/36, 2/36,15/36;
      -64/36,-22/36,-7/36, 2/36,15/36;
      -64/36,-22/36,-7/36, 2/36,15/36];
%地形对应修正量
Smod=[0, 1, 2, 2];
%% 初始化随机初始值
%随机生成一个地形
maprow=20;
mapcol=20;
%写死一个20*20的地形
Map=[4,3,3,3,2,2,2,1,1,1,1,0,0,0,0,1,1,1,2,2;
     3,3,3,3,2,2,2,1,1,1,0,0,0,0,1,1,1,2,2,3;
     3,3,3,2,2,2,1,1,1,0,0,0,0,0,0,1,1,1,2,2;
     3,3,2,2,2,1,1,1,0,0,0,0,0,0,0,0,1,1,1,2;
     3,2,2,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1;
     2,2,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1;
     2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;
     1,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0;
     2,1,1,1,1,0,0,0,1,2,3,1,0,0,0,0,0,0,1,1;
     2,2,1,1,0,0,0,0,1,2,2,1,0,0,0,0,0,1,1,1;
     3,3,2,2,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,2;
     3,2,2,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,2,3,
     3,2,1,1,0,0,0,0,0,1,1,0,0,0,0,0,1,2,3,4;
     4,3,2,1,1,0,0,1,1,1,2,1,0,0,0,0,0,1,2,3;
     5,4,3,2,1,1,1,2,2,2,2,2,1,0,0,0,0,0,1,2;
     5,4,3,2,1,1,0,1,2,3,3,2,1,0,0,0,0,0,0,1;
     4,3,3,2,1,0,0,0,1,2,3,2,2,1,1,0,0,0,0,0;
     4,3,2,1,0,0,0,0,1,2,3,2,1,1,0,0,0,0,0,1;
     3,2,1,0,0,0,0,1,2,3,4,3,2,1,1,1,0,0,1,1;
     2,1,0,0,0,0,1,1,2,3,3,2,1,1,0,0,0,1,1,2;];
% maphprob=[0.4,0.28,0.15,0.10,0.5,0.2];
% Map=zeros(maprow,mapcol);
% Map(1,1)=0;
% for i=1:maprow
%     for j=1:mapcol
%         if i==1 && j==1
%             continue
%         elseif j==1
%             prev=Map(i-1,mapcol);
%         elseif i==1
%             prev=Map(i,j-1);
%         else
%             prev=ceil((Map(i,j-1)+Map(i-1,j))/2);
%         end
%         possiblevalue=max(0,prev-1):min(5,prev+1);
%         possibleprob=maphprob(possiblevalue+1);
%         possibleprob=possibleprob/sum(possibleprob);
%         Map(i,j)=randsample(possiblevalue,1,true,possibleprob);
%     end
% end
%绘制地形图
% [X,Y]=meshgrid(1:size(Map,2),1:size(Map,1));
% surf(X,Y,Map);
% zlim([0,20]);

%随机生成一个地貌
% shapenum=[1,2,3,4];
% probabilities=[0.4,0.1,0.3,0.2];
% Shape=randsample(shapenum,mapcol*maprow,true,probabilities);
% Shape=reshape(Shape,[maprow mapcol]);
Shape=[3,3,3,3,3,3,1,1,4,1,1,1,1,1,1,1,1,1,1,1;
       3,3,3,3,3,3,3,1,4,1,1,1,1,1,1,1,1,3,3,3;
       3,3,3,3,3,3,1,1,4,1,1,1,1,1,1,1,1,3,3,3;
       3,3,3,3,3,1,1,4,1,3,3,3,1,1,1,1,1,1,3,3;
       3,3,3,3,1,1,4,1,1,3,3,3,3,1,1,1,1,1,1,3;
       3,3,3,1,1,1,4,3,3,3,3,3,1,1,1,1,3,3,3,3;
       3,3,1,1,1,1,1,4,3,3,3,1,1,1,1,1,3,3,3,1;
       3,1,1,1,1,1,4,1,1,1,1,1,1,1,1,1,3,3,1,1;
       3,3,1,1,1,1,4,1,1,2,2,1,1,1,1,1,3,1,1,3;
       3,3,3,1,1,1,4,1,1,2,2,1,1,1,1,1,1,1,1,3;
       3,3,3,3,1,1,1,4,1,1,1,1,1,1,1,1,1,1,1,3;
       3,3,3,3,3,1,1,1,4,4,4,4,4,4,4,4,1,1,1,1;
       3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1,4,1,1,1;
       3,3,1,1,1,1,1,1,1,3,3,3,1,1,1,1,4,2,2,1;
       3,1,1,1,1,1,1,1,3,3,3,3,1,1,1,1,4,2,1,1;
       1,1,1,1,1,1,1,3,3,1,1,3,1,1,1,4,1,1,1,1;
       1,1,1,1,1,1,1,3,3,1,1,3,1,1,1,4,1,1,1,1;
       1,1,1,1,1,1,1,1,3,1,1,1,1,1,1,1,4,1,1,1;
       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,4,1,1;
       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,4,1,1,1;];



% imagesc(Shape);
% colormap([1 1 0;1 1 1;0 1 0;0 0 1;]);
% colorbar;
% colormatrix=zeros(20,20,3);
% for i=1:20
%     for j=1:20
%         colormatrix(i,j,:)=colormap(shape(i,j),:);
%     end
% end
% surf(X,Y,Map,colormatrix);
% zlim([0,20]);


%生成插值函数句柄以减少代码运行时间
findh=griddedInterpolant(Map);
finds=griddedInterpolant(Shape,'nearest');
findr=griddedInterpolant(0:20,Vrank);
findam=griddedInterpolant(Amod);
findvr=griddedInterpolant(VRes);
findvm=griddedInterpolant(Vmod);

%生成初始的agent状态
ag1num=4;
ag2num=4;
stat1=zeros(4,ag1num);%决策变量每个包含位置x，位置y，剩余班数量，护甲类型
stat2=zeros(4,ag2num);
% stat1(1,:)=randi([1 mapcol/2],1,ag1num);
% stat1(2,:)=randi([1 mapcol],1,ag1num);
% stat1(3,:)=3;
% stat1(4,:)=1;
% stat2(1,:)=randi([mapcol/2 mapcol],1,ag2num);
% stat2(2,:)=randi([1 mapcol],1,ag2num);
% stat2(3,:)=4;
% stat2(4,:)=1;
stat1=[3,4,5,6;
       9,5,3,7;
       3,3,3,3;
       1,1,1,1;];
stat2=[10,12,13,12;
       10,12,15,11;
       3,3,3,3;
       1,1,1,1;];
save('stat1.mat','stat1');
save('stat2.mat','stat2');

%绘制初始化后的示意图
% scatter(stat1(1,:),stat1(2,:),'blue');
% hold on;
% scatter(stat2(1,:),stat2(2,:),'red');
% xlim([0,100]);
% ylim([0 100]);

%% 问题定义

% VarSize1 = [4 ag1num];   % Size of Decision Variables Matrix
% VarSize2 = [4 ag2num];

%搜索范围限制
VarMin = -3;          % Lower Bound of Variables
VarMax = 3;          % Upper Bound of Variables

CostFunction=@(x,y,stat1,stat2) COMOP0(x,y,stat1,stat2,Map,Shape,Vrank,Amod,VRes,Vmod,Smod,findh,finds,findr,findam,findvr,findvm);

%% 参数定义

tao1=3;   %两个种群在大进化代间的小进化次数
tao2=3;   

MaxIt = 100;      % Maximum Number of Iterations

nPop = 50;        % Population Size

pCrossover = 0.7;                         % Crossover Percentage
nCrossover = 2*round(pCrossover*nPop/2);  % Number of Parnets (Offsprings)

pMutation = 0.3;                          % Mutation Percentage
nMutation = round(pMutation*nPop);        % Number of Mutants

mu = 0.02;                    % Mutation Rate

sigma = 0.1*(VarMax-VarMin);  % Mutation Step Size

%% 初始化

empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Rank = [];
empty_individual.DominationSet = [];
empty_individual.DominatedCount = [];
empty_individual.CrowdingDistance = [];

pop1 = repmat(empty_individual, nPop, 1);
pop2 = repmat(empty_individual, nPop, 1);

for i = 1:nPop

    %生成决策变量
    Agent1=zeros(4,ag1num);
    Agent2=zeros(4,ag2num);
    %第一个Agent进行随机搜索取位置
    % Agent1(1:2,:)=randi([VarMin VarMax],2,ag1num);
     Agent1(1:2,:)=unifrnd(VarMin,VarMax,2,ag1num);
    % rp1((rp1+stat1(1:2,:))<1)=1;%剔除非法值
    % rp1((rp1+stat1(1:2,:))>mapcol)=-stat1(1:2,:)+mapcol;%剔除非法值
    Agent1(3,:)=randi([1 ag2num],1,ag1num);
    Agent1(4,:)=ones(1,ag1num);
    %进行第二个agent的初始化
    % Agent2(1:2,:)=randi([VarMin VarMax],2,ag2num);
    Agent2(1:2,:)=unifrnd(VarMin,VarMax,2,ag2num);
    % Agent2(Agent2<1)=1;%剔除非法值
    % Agent2(Agent2>mapcol)=mapcol;%剔除非法值
    Agent2(3,:)=randi([1 ag1num],1,ag2num);
    Agent2(4,:)=ones(1,ag2num);
        
    pop1(i).Position = Agent1;

    pop2(i).Position = Agent2;
    
    % pop1(i).Cost = CostFunction(pop1(i).Position,pop2(i).Position,stat1,stat2);
    % 
    % pop2(i).Cost = CostFunction(pop2(i).Position,pop1(i).Position,stat2,stat1);
    
end
%在初始化时计算初始种群的函数值
for i=1:nPop
    fitness1s=zeros(nPop,2);
    fitness2s=zeros(nPop,2);

    fitness1=0;
    fitness2=0;
    for j=1:nPop
        % fitness1s(j,:)=CostFunction(pop1(i).Position,pop2(j).Position,stat1,stat2);
        % fitness2s(j,:)=CostFunction(pop2(i).Position,pop1(j).Position,stat2,stat1);
        fitness1=fitness1+CostFunction(pop1(i).Position,pop2(j).Position,stat1,stat2);
        fitness2=fitness2+CostFunction(pop2(i).Position,pop1(j).Position,stat2,stat1);
    end
    % atk1=fitness1s(:,1);
    % atk2=fitness2s(:,1);
    % [~,mindex1]=min(atk1);
    % [~,mindex2]=min(atk2);
    % pop1(i).Cost=fitness1s(mindex1,:);
    % pop2(i).Cost=fitness2s(mindex2,:);
    % pop1(i).Cost=
end

% Non-Dominated Sorting
[pop1, F] = NonDominatedSorting(pop1);

[pop2, G] = NonDominatedSorting(pop2);

% Calculate Crowding Distance
pop1 = CalcCrowdingDistance(pop1, F);

pop2=CalcCrowdingDistance(pop2,G);

% Sort Population
[pop1, F] = SortPopulation(pop1);

[pop2, G] = SortPopulation(pop2);

%% NSGA-II Main Loop

t1=0;
t2=0;
fileID=fopen('output.txt','w');
clf
%计算初始种群
for i=1:nPop
    fitness2s=zeros(nPop,2);
    fitness2=0;
    for j=1:nPop
        fitness2=fitness2+CostFunction(pop2(i).Position,pop1(j).Position,stat2,stat1);
    end
    pop2(i).Cost=fitness2/nPop;
end
for i=1:nPop
    fitness1s=zeros(nPop,2);
    fitness1=0;
    for j=1:nPop
        fitness1 = fitness1+CostFunction(pop1(i).Position,pop2(j).Position,stat1,stat2);
    end
    pop1(i).Cost=fitness1/nPop;
end
for it = 1:MaxIt
    %计算对手种群改变以后的父代函数值
    for i=1:nPop
         fitness1s=zeros(nPop,2);
         fitness1=0;
         for j=1:nPop
             fitness1 = fitness1+CostFunction(pop1(i).Position,pop2(j).Position,stat1,stat2);
         end
         pop1(i).Cost=fitness1/nPop;
    end

    for t=t1+1:t1+tao1
        % Crossover
        popc = repmat(empty_individual, nCrossover/2, 2);
        for k = 1:nCrossover/2

            i1 = randi([1 nPop]);
            p1 = pop1(i1);

            i2 = randi([1 nPop]);
            p2 = pop1(i2);

            [popc(k, 1).Position, popc(k, 2).Position] = Crossoveri(p1.Position, p2.Position);

            fitness1=0;
            fitness2=0;

            fitness1s=zeros(nPop,2);
            fitness2s=zeros(nPop,2);
            for y=1:nPop
                fitness1 = fitness1 + CostFunction(popc(k, 1).Position,pop2(y).Position,stat1,stat2);
                fitness2 = fitness2 + CostFunction(popc(k, 2).Position,pop2(y).Position,stat1,stat2);
                % fitness1s(y,:)=CostFunction(popc(k, 1).Position,pop2(y).Position,stat1,stat2);
                % fitness2s(y,:)=CostFunction(popc(k, 2).Position,pop2(y).Position,stat1,stat2);
            end
            % atk1=fitness1s(:,1);
            % atk2=fitness2s(:,1);
            % [~,mindex1]=min(atk1);
            % [~,mindex2]=min(atk2);
            % popc(k,1).Cost=fitness1s(mindex1,:);
            % popc(k,2).Cost=fitness2s(mindex2,:);
            popc(k, 1).Cost = fitness1/nPop;
            popc(k, 2).Cost = fitness2/nPop;

        end
        popc = popc(:);

        % Mutation
        popm = repmat(empty_individual, nMutation, 1);
        for k = 1:nMutation

            i = randi([1 nPop]);
            p = pop1(i);

            range1=[VarMin VarMax;VarMin VarMax; 1 ag2num;1 1];

            m1 = Mutatei(p.Position, mu, sigma);
            % %处理越界的情况
            p1=m1(1:2,:);
            p1(p1<1)=1;
            p1(p1>mapcol)=mapcol;
            m1(1:2,:)=p1;
            tg1=m1(3,:);
            tg1(tg1<1)=1;
            tg1(tg1>ag2num)=ag2num;
            m1(3,:)=tg1;
            m1(4,:)=1;

            popm(k).Position=m1;

            fitness=0;
            fitness1s=zeros(nPop,2);
            for y=1:nPop
                fitness = fitness + CostFunction(popm(k).Position,pop2(y).Position,stat1,stat2);
                % fitness1s(y,:)=CostFunction(popm(k).Position,pop2(y).Position,stat1,stat2);
            end
            popm(k).Cost = fitness/nPop;
            % atk1=fitness1s(:,1);
            % [~,mindex1]=min(atk1);
            % popm(k).Cost=fitness1s(mindex1,:);
        end

        % Merge
        pop1 = [pop1
             popc
             popm]; %#ok

        % Non-Dominated Sorting
        [pop1, F] = NonDominatedSorting(pop1);

        % Calculate Crowding Distance
        pop1 = CalcCrowdingDistance(pop1, F);

        % Sort Population
        pop1 = SortPopulation(pop1);

        % Truncate
        pop1 = pop1(1:nPop);

        % Non-Dominated Sorting
        [pop1, F] = NonDominatedSorting(pop1);

        % Calculate Crowding Distance
        pop1 = CalcCrowdingDistance(pop1, F);

        % Sort Population
        [pop1, F] = SortPopulation(pop1);

        % Store F1
        F1 = pop1(F{1});

        % Show Iteration Information
        disp(['Iteration ' num2str(it) ':t ' num2str(t) ': Number of F1 Members = ' num2str(numel(F1))]);
        fprintf(fileID,'iteration');
        fprintf(fileID,num2str(it));
        fprintf(fileID,'\t it');
        fprintf(fileID,num2str(t));
        fprintf(fileID,'\t F');
        fprintf(fileID,'\t');
        fprintf(fileID,'%.5f ',F1.CrowdingDistance);
        fprintf(fileID,'\n');
        % Plot F1 Costs
        % figure(1);
        % PlotCosts(F1);
        % pause(0.01);
        % figure(2);
        % figure(2);
        % Plot1Agent(F1);
    end

    t1=t1+tao1;

    %计算对手种群改变以后的父代函数值
    for i=1:nPop
         fitness2s=zeros(nPop,2);
         fitness2=0;
         for j=1:nPop
             fitness2=fitness2+CostFunction(pop2(i).Position,pop1(j).Position,stat2,stat1);
         end
         pop2(i).Cost=fitness2/nPop;
     end

    for t=t2+1:t2+tao2
        % Crossover
        popc = repmat(empty_individual, nCrossover/2, 2);
        for k = 1:nCrossover/2

            i1 = randi([1 nPop]);
            p1 = pop2(i1);

            i2 = randi([1 nPop]);
            p2 = pop2(i2);

            [popc(k, 1).Position, popc(k, 2).Position] = Crossoveri(p1.Position, p2.Position);

            fitness1=0;
            fitness2=0;
            fitness1s=zeros(nPop,2);
            fitness2s=zeros(nPop,2);
            for x=1:nPop
                 fitness1=fitness1+CostFunction(popc(k, 1).Position,pop1(x).Position,stat2,stat1);
                 fitness2=fitness2+CostFunction(popc(k, 2).Position,pop1(x).Position,stat2,stat1);
                %fitness1s(x,:)=CostFunction(popc(k, 1).Position,pop1(x).Position,stat2,stat1);
                %fitness2s(x,:)=CostFunction(popc(k, 2).Position,pop1(x).Position,stat2,stat1);
            end

            popc(k, 1).Cost = fitness1/nPop;
             popc(k, 2).Cost = fitness2/nPop;

            % atk1=fitness1s(:,1);
            % atk2=fitness2s(:,1);
            % [~,mindex1]=min(atk1);
            % [~,mindex2]=min(atk2);
            % popc(k,1).Cost=fitness1s(mindex1,:);
            % popc(k,2).Cost=fitness2s(mindex2,:);

        end
        popc = popc(:);

        % Mutation
        popm = repmat(empty_individual, nMutation, 1);
        for k = 1:nMutation

            i = randi([1 nPop]);
            p = pop2(i);

            range2=[VarMin VarMax;VarMin VarMax; 1 ag1num;1 1];
            m2 = Mutatei(p.Position, mu, sigma);
            % %处理越界的情况
            p2=m2(1:2,:);
            p2(p2<1)=1;
            p2(p2>mapcol)=mapcol;
            m2(1:2,:)=p2;
            tg2=m2(3,:);
            tg2(tg2<1)=1;
            tg2(tg2>ag2num)=ag2num;
            m2(3,:)=tg2;
            m2(4,:)=1;

            popm(k).Position=m2;

            fitness=0;
            fitness1s=zeros(nPop,2);
            for x=1:nPop
                 fitness = fitness+CostFunction(popm(k).Position,pop1(x).Position,stat2,stat1);
                % fitness1s(y,:)=CostFunction(popm(k).Position,pop1(x).Position,stat2,stat1);
            end
            popm(k).Cost = fitness/nPop;
            % atk1=fitness1s(:,1);
            % [~,mindex1]=min(atk1);
            % popm(k).Cost=fitness1s(mindex1,:);
        end

        % Merge
        pop2 = [pop2
             popc
             popm]; %#ok

        % Non-Dominated Sorting
        [pop2, G] = NonDominatedSorting(pop2);

        % Calculate Crowding Distance
        pop2 = CalcCrowdingDistance(pop2, G);

        % Sort Population
        pop2 = SortPopulation(pop2);

        % Truncate
        pop2 = pop2(1:nPop);

        % Non-Dominated Sorting
        [pop2, G] = NonDominatedSorting(pop2);

        % Calculate Crowding Distance
        pop2 = CalcCrowdingDistance(pop2, G);

        % Sort Population
        [pop2, G] = SortPopulation(pop2);

        % Store F1
        G1 = pop2(G{1});

        % Show Iteration Information
        disp(['Iteration ' num2str(it) ':t ' num2str(t) ': Number of G1 Members = ' num2str(numel(G1))]);
        fprintf(fileID,'iteration');
        fprintf(fileID,num2str(it));
        fprintf(fileID,'\t it');
        fprintf(fileID,num2str(t));
        fprintf(fileID,'\t G');
        fprintf(fileID,'\t');
        fprintf(fileID,'%.5f ',G1.CrowdingDistance);
        fprintf(fileID,'\n');
        % Plot F1 Costs
        % figure(1);
        % PlotCosts(G1);
        % pause(0.01);
        % figure(3);
        % Plot1Agent(G1);
    end
    t2=t2+tao2;
    %绘制两个agent的图像
    figure(1);
    pause(0.1);
    Plot2AgentR(F1,G1);
    % Costs1 = [F1.Cost];
    % Costs1(1,:)=-Costs1(1,:);
    % Costs1(2,:)=Costs1(2,:)+5;
    % plot(Costs1(1, :), Costs1(2, :), 'r*', 'MarkerSize', 8);
    % xlabel('1^{st} Objective');
    % ylabel('2^{nd} Objective');
    % title('Non-dominated Solutions (F_{1})');
    % clf
    % scatter(F1(1).Position(1,:)+stat1(1,:),F1(1).Position(2,:)+stat1(2,:),'red');
    % hold on;
    % scatter(G1(1).Position(1,:)+stat2(1,:),G1(1).Position(2,:)+stat2(2,:),'blue');
    % xlim([1 20]);
    % ylim([1 20]);
    % pause(0.5)
end
fclose(fileID);
% OG=G1;
%% Results
% save('newF1.mat','F1');
% save('newG1.mat','G1');
% save('OG.mat','OG');
NF=F1;
NG=G1;
save('NF.mat','NF');
save('NG.mat','NG');


