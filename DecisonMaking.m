load('newF1.mat');
load('newG1.mat');
load("stat1.mat");
load("stat2.mat");

F1p=zeros(50,4,4);
% F1pn=zeros(50,4,4);
for i=1:50
    F1p(i,:,:)=round(F1(i).Position);
    % F1pn(i,:,:)=F1(i).Position;
end
costmapg=zeros(50,2);
for i=1:50
    costmapg(i,:)=G1(i).Cost;
end
costmapf=zeros(50,2);
for i=1:50
    costmapf(i,:)=F1(i).Cost;
end
% F1p=round(F1.Position);
G1p=zeros(50,4,4);

% G1pn=zeros(50,4,4);
for i=1:50
    G1p(i,:,:)=round(G1(i).Position);
    % G1pn(i,:,:)=G1(i).Position;
end
% G1p=round(G1.Position);
% palateofont=zeros(size(F1pn,1),2);
% for j=1:size(G1pn,1)
%     fitness=0;
%     for i=1:size(F1pn,1)
%         fitness=fitness+CostFunction(reshape(F1pn(i,:,:),4,4),reshape(G1pn(j,:,:),4,4));
%     end
%     fitness=fitness/size(F1pn,1);
%     palateofont(j,:)=fitness;
% end


%检查离散化后是否有重复的决策变量，如果有的话去掉
uniqueF1p=[];
for i=1:size(F1p,1)
    currentmatrix=F1p(i,:,:);
    dup=false;
    for j=1:size(uniqueF1p,1)
        if isequal(currentmatrix,uniqueF1p(j,:,:))
            dup=true;
            break;
        end
    end
    if ~dup
        uniqueF1p=cat(1,uniqueF1p,currentmatrix);
    end
end
uniqueG1p=[];
for i=1:size(G1p,1)
    currentmatrix=G1p(i,:,:);
    dup=false;
    for j=1:size(uniqueG1p,1)
        if isequal(currentmatrix,uniqueG1p(j,:,:))
            dup=true;
            break;
        end
    end
    if ~dup
        uniqueG1p=cat(1,uniqueG1p,currentmatrix);
    end
end
size1=size(uniqueF1p,1);
size2=size(uniqueG1p,1);
cost1=zeros(size1,2);
cost2=zeros(size2,2);
% %离散化后的目标函数比较
% for i=1:size1
%     fitness1s=zeros(nPop,2);
%     fitness1=0;
%     for j=1:size2
%         fitness1 = fitness1+CostFunction(reshape(uniqueF1p(i,:,:),4,4),reshape(uniqueG1p(j,:,:),4,4),stat1,stat2);
%     end
%     cost1(i,:)=fitness1/nPop;
% end
% for i=1:size2
%     fitness2s=zeros(nPop,2);
%     fitness2=0;
%     for j=1:size1
%         fitness2=fitness2+CostFunction(reshape(uniqueG1p(i,:,:),4,4),reshape(uniqueF1p(j,:,:),4,4),stat2,stat1);
%     end
%     cost2(i,:)=fitness2/nPop;
% end
% figure(1)
% palateofont(:,1)=-costmapg(:,1);
% palateofont(:,2)=costmapg(:,2);
% palateofont1(:,1)=-costmapf(:,1);
% palateofont1(:,2)=costmapf(:,2);
% plot(palateofont1(:, 1), palateofont1(:, 2), 'r*',palateofont(:, 1), palateofont(:, 2), 'b*',cost1(:,1),cost1(:,2),'g*', 'MarkerSize', 8);
% xlabel('1^{st} Objective');
% ylabel('2^{nd} Objective');
% title('Corresponding Inital Strategy');
% legend('Pareto front A','Pareto front B','Discreat strategy set A','Location','northeast');

%选取第一个值
alldiv=zeros(size(uniqueF1p,1),1);
for i=1:size(uniqueF1p,1)
    divmatrix=zeros(size(uniqueG1p,1),2);
    for j=1:size(uniqueG1p,1)
        divmatrix(j,:)=CostFunction(reshape(uniqueF1p(i,:,:),4,4),reshape(uniqueG1p(j,:,:),4,4),stat1,stat2);
    end
    colmeans=mean(divmatrix);
    colstds=std(divmatrix);
    normaldiv=(divmatrix-repmat(colmeans,size(divmatrix,1),1))./repmat(colstds,size(divmatrix,1),1);
    avenormstd=std(normaldiv(:));
    alldiv(i,1)=avenormstd;
end
[~,minindex]=min(alldiv);
costmap=zeros(size(uniqueF1p,1),2);
for i=1:size(uniqueF1p,1)
    fitness=0;
    for j=1:size(uniqueG1p,1)
        fitness=fitness+CostFunction(reshape(uniqueF1p(i,:,:),4,4),reshape(uniqueG1p(j,:,:),4,4),stat1,stat2);
    end
    fitness=fitness/size(uniqueG1p,1);
    costmap(i,:)=fitness;
end
%绘图
costmap(:,1)=-costmap(:,1);
costmap(:,2)=costmap(:,2);
plot(costmap(:, 1), costmap(:, 2), 'r*',costmap(minindex, 1), costmap(minindex, 2), 'k*', 'MarkerSize', 8);
xlabel('1^{st} Objective');
ylabel('2^{nd} Objective');
title('Initial Strategy Selection');
legend('discretize Non-dominate Set','Lowest Normalized Standard Deviation','Location','southeast');
xlim([0 5]);
ylim([0 5]);
% grid on;
% hold on;
% plot(costmap(minindex, 1), costmap(minindex, 2), 'k*', 'MarkerSize', 8);
% grid on;
hold off;

des1set=[];
des1set(end+1)=minindex;
des2set=[];
havecir=false;
tfindex=minindex;
while havecir~=true
    divmatrix2=zeros(size(uniqueG1p,1),2);
    for j=1:size(uniqueG1p,1)
        % sp1=CostFunction(reshape(uniqueF1p(minindex,:,:),4,4),reshape(uniqueG1p(j,:,:),4,4));
        sp2=CostFunction(reshape(uniqueG1p(j,:,:),4,4),reshape(uniqueF1p(tfindex,:,:),4,4),stat2,stat1);
        divmatrix2(j,:)=sp2;
    end
    tfindex=TradeOffValue(divmatrix2);
    [isinarray,index]=ismember(tfindex,des2set);
    if isinarray
        des2set(end+1)=tfindex;
        break;
    else
        des2set(end+1)=tfindex;
    end
    divmatrix=zeros(size(uniqueF1p,1),2);
    for i=1:size(uniqueF1p,1)
        sp1=CostFunction(reshape(uniqueF1p(i,:,:),4,4),reshape(uniqueG1p(tfindex,:,:),4,4),stat1,stat2);
        divmatrix(i,:)=sp1;
    end
    tfindex=TradeOffValue(divmatrix);
    [isinarray,index]=ismember(tfindex,des1set);
    if isinarray
        des1set(end+1)=tfindex;
        break;
    else
        des1set(end+1)=tfindex;
    end
end

va1b6=CostFunction(reshape(uniqueF1p(1,:,:),4,4),reshape(uniqueG1p(6,:,:),4,4),stat1,stat2);
va1b1=CostFunction(reshape(uniqueF1p(1,:,:),4,4),reshape(uniqueG1p(1,:,:),4,4),stat1,stat2);
va10b1=CostFunction(reshape(uniqueF1p(10,:,:),4,4),reshape(uniqueG1p(1,:,:),4,4),stat1,stat2);
va10b6=CostFunction(reshape(uniqueF1p(10,:,:),4,4),reshape(uniqueG1p(6,:,:),4,4),stat1,stat2);

Apoints=zeros(4,2);
Apoints(1,:)=va1b6;Apoints(2,:)=va1b1;Apoints(3,:)=va10b1;Apoints(4,:)=va10b6;


figure(1)
divmatrix(:,1)=-divmatrix(:,1);
divmatrix(:,2)=divmatrix(:,2);
divmatrix2(:,1)=-divmatrix2(:,1);
divmatrix2(:,2)=divmatrix2(:,2);
palateofont(:,1)=-costmapg(:,1);
palateofont(:,2)=costmapg(:,2);
palateofont1(:,1)=-costmapf(:,1);
palateofont1(:,2)=costmapf(:,2);
plot(palateofont1(:, 1), palateofont1(:, 2), 'r*',palateofont(:, 1), palateofont(:, 2), 'b*',Apoints(:,1),Apoints(:,2),'g*', 'MarkerSize', 8);
xlabel('1^{st} Objective');
ylabel('2^{nd} Objective');
title('Corresponding Inital Strategy');
legend('Pareto front A','Pareto front B','Conter strategy A','Location','northeast');



% % 计算种群中每个个体针对初始策略的回应（计算目标函数）
% divmatrix=zeros(size(uniqueG1p,1),2);
% for j=1:size(uniqueG1p,1)
% 
%     % sp1=CostFunction(reshape(uniqueF1p(minindex,:,:),4,4),reshape(uniqueG1p(j,:,:),4,4));
%     sp2=CostFunction(reshape(uniqueG1p(j,:,:),4,4),reshape(uniqueF1p(minindex,:,:),4,4),stat2,stat1);
%     divmatrix(j,:)=sp2;
% end
% %计算最高的权衡值
% avgrset=zeros(size(divmatrix,1),1);
% for i=1:size(divmatrix,1)
%     avgr=0;
%     avgn=0;
%     % s=zeros(size(divmatrix,1),2);
%     for j=1:size(divmatrix,1)
%         if j~=i
%             s=divmatrix(i,:)-divmatrix(j,:);
%             loss=0;
%             gain=0;
%             if s~=0
%                 for k=1:size(s,2)
%                     if s(k)<0
%                         loss=loss+abs(s(k));
%                     else
%                         gain=gain+abs(s(k));
%                     end
%                 end
%                 rij=gain/loss;
%             else
%                 rij=0;
%             end
%             if rij~=0 && ~isinf(rij)
%                 avgr=avgr+rij;
%                 avgn=avgn+1;
%             end
%         end
%     end
%     avgr=avgr/avgn;
%     avgrset(i)=avgr;
% end
% [~,mindex]=min(avgrset);
% divmatrix(:,1)=-divmatrix(:,1);
% divmatrix(:,2)=divmatrix(:,2);
% 
% 
% % 绘图
% palateofont(:,1)=-costmapg(:,1);
% palateofont(:,2)=costmapg(:,2);
% plot(divmatrix(:, 1), divmatrix(:, 2), 'r*',palateofont(:, 1), palateofont(:, 2), 'b*',divmatrix(mindex,1),divmatrix(mindex,2),'k*', 'MarkerSize', 8);
% xlabel('1^{st} Objective');
% ylabel('2^{nd} Objective');
% title('Corresponding Inital Strategy');
% legend('strategy set','Pareto front','highest trade-off','Location','southeast');
% xlim([0 5]);
% ylim([0 5]);
% hold on;
% grid on;
% plot([3.2],[0.3],'r*', 'MarkerSize', 8);
% hold on;
% grid on;
% plot([3.5], [0.25], 'g*', 'MarkerSize', 8);