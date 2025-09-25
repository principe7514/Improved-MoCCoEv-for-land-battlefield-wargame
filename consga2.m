clc;
clear;
close all;

%% Problem Definition

CostFunction = @(x,y) COMOP1(x,y);      % Cost Function

nVar1 = 3;             % Number of Decision Variables
nVar2 = 4;

VarSize1 = [1 nVar1];   % Size of Decision Variables Matrix
VarSize2 = [1 nVar2];

VarMin = -5;          % Lower Bound of Variables
VarMax = 5;          % Upper Bound of Variables

% Number of Objective Functions
nObj1 = numel(CostFunction(unifrnd(VarMin, VarMax, VarSize1),unifrnd(VarMin, VarMax, VarSize2)));
nObj2 = numel(CostFunction(unifrnd(VarMin, VarMax, VarSize1),unifrnd(VarMin, VarMax, VarSize2)));

%% NSGA-II Parameters

tao1=5;   %两个种群在大进化代间的小进化次数
tao2=5;   

MaxIt = 50;      % Maximum Number of Iterations

nPop = 50;        % Population Size

pCrossover = 0.7;                         % Crossover Percentage
nCrossover = 2*round(pCrossover*nPop/2);  % Number of Parnets (Offsprings)

pMutation = 0.4;                          % Mutation Percentage
nMutation = round(pMutation*nPop);        % Number of Mutants

mu = 0.02;                    % Mutation Rate

sigma = 0.1*(VarMax-VarMin);  % Mutation Step Size

%% Initialization

empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Rank = [];
empty_individual.DominationSet = [];
empty_individual.DominatedCount = [];
empty_individual.CrowdingDistance = [];

pop1 = repmat(empty_individual, nPop, 1);
pop2 = repmat(empty_individual, nPop, 1);

for i = 1:nPop
    
    pop1(i).Position = unifrnd(VarMin, VarMax, VarSize1);

    pop2(i).Position = unifrnd(VarMin,VarMax,VarSize2);
    
    pop1(i).Cost = CostFunction(pop1(i).Position,pop2(i).Position);

    pop2(i).Cost = CostFunction(pop1(i).Position,pop2(i).Position);
    
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

for it = 1:MaxIt

    for t=t1+1:t1+tao1
        % Crossover
        popc = repmat(empty_individual, nCrossover/2, 2);
        for k = 1:nCrossover/2
            
            i1 = randi([1 nPop]);
            p1 = pop1(i1);
            
            i2 = randi([1 nPop]);
            p2 = pop1(i2);
            
            [popc(k, 1).Position, popc(k, 2).Position] = Crossover(p1.Position, p2.Position);
            
            fitness1=0;
            fitness2=0;

            for y=1:nPop
                fitness1 = fitness1 + CostFunction(popc(k, 1).Position,pop2(y).Position);
                fitness2 = fitness2 + CostFunction(popc(k, 2).Position,pop2(y).Position);
            end

            popc(k, 1).Cost = fitness1/nPop;
            popc(k, 2).Cost = fitness2/nPop;
            
        end
        popc = popc(:);
        
        % Mutation
        popm = repmat(empty_individual, nMutation, 1);
        for k = 1:nMutation
            
            i = randi([1 nPop]);
            p = pop1(i);
            
            popm(k).Position = Mutate(p.Position, mu, sigma);

            fitness=0;

            for y=1:nPop
                fitness = fitness + CostFunction(popm(k).Position,pop2(y).Position);
            end
            
            popm(k).Cost = fitness/nPop;
            
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
        disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);
        
        % Plot F1 Costs
        % figure(1);
        % PlotCosts(F1);
        % pause(0.01);
    end

    t1=t1+tao1;

    for t=t2+1:t2+tao2
        % Crossover
        popc = repmat(empty_individual, nCrossover/2, 2);
        for k = 1:nCrossover/2
            
            i1 = randi([1 nPop]);
            p1 = pop2(i1);
            
            i2 = randi([1 nPop]);
            p2 = pop2(i2);
            
            [popc(k, 1).Position, popc(k, 2).Position] = Crossover(p1.Position, p2.Position);

            fitness1=0;
            fitness2=0;

            for x=1:nPop
                fitness1=fitness1+CostFunction(pop1(x).Position,popc(k, 1).Position);
                fitness2=fitness2+CostFunction(pop1(x).Position,popc(k, 2).Position);
            end
            
            popc(k, 1).Cost = fitness1/nPop;
            popc(k, 2).Cost = fitness2/nPop;
            
        end
        popc = popc(:);
        
        % Mutation
        popm = repmat(empty_individual, nMutation, 1);
        for k = 1:nMutation
            
            i = randi([1 nPop]);
            p = pop2(i);
            
            popm(k).Position = Mutate(p.Position, mu, sigma);

            fitness=0;

            for x=1:nPop
                fitness = fitness+CostFunction(pop1(x).Position,popm(k).Position);
            end
            
            popm(k).Cost = fitness/nPop;
            
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
        disp(['Iteration ' num2str(it) ': Number of G1 Members = ' num2str(numel(G1))]);
        
        % Plot F1 Costs
        % figure(1);
        % PlotCosts(G1);
        % pause(0.01);
    end
    t2=t2+tao2;
    
    %绘制两个agent的图像
    figure(1);
    Plot2Agent(F1,G1);
    pause(0.5)
end

%% Results