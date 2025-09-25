load('newF1.mat');
load('newG1.mat');
load("OF.mat");
load("OG.mat");

    % Plot2AgentR(F1,G1);
Costs1 = [F1.Cost];
Costs1(1,:)=-Costs1(1,:);
Costs1(2,:)=Costs1(2,:);
Costs2 = [G1.Cost];
Costs2(1,:)=-Costs2(1,:);
Costs2(2,:)=Costs2(2,:);
Costs3 = [OF.Cost];
Costs3(1,:)=-Costs3(1,:);
Costs3(2,:)=Costs3(2,:);
Costs4 = [OG.Cost];
Costs4(1,:)=-Costs4(1,:);
Costs4(2,:)=Costs4(2,:);
plot(Costs1(1, :), Costs1(2, :), 'r*',Costs2(1, :), Costs2(2, :), 'b*' ...
    ,Costs3(1, :), Costs3(2, :), 'm*',Costs4(1, :), Costs4(2, :), 'c*', 'MarkerSize', 8);
xlabel('1^{st} Objective');
ylabel('2^{nd} Objective');
title('Non-dominated Solutions (F_{1})');
legend('PF A with competitive co-evluation','PF B with competitive co-evluation','PF A with independently optimal','PF B with independently optimal','Location','northeast');
grid on;
hold off;
xlim([0,5]);
ylim([0,5]);