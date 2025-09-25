
function Plot2AgentR(pop1,pop2)
    num=numel(pop1);
    Costs1=zeros(num,2);
    for i=1:num
        Costs1(i,:)=pop1(i).Cost;
    end
    num2=numel(pop2);
    Costs2=zeros(num2,2);
    for i=1:num2
        Costs2(i,:)=pop2(i).Cost;
    end
    Costs1 = [pop1.Cost];
    Costs1(1,:)=-Costs1(1,:);
    Costs1(2,:)=Costs1(2,:);
    Costs2 = [pop2.Cost];
    Costs2(1,:)=-Costs2(1,:);
    Costs2(2,:)=Costs2(2,:);

    plot(Costs1(1, :), Costs1(2, :), 'r*',Costs2(1, :), Costs2(2, :), 'b*', 'MarkerSize', 8);
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    title('Non-dominated Solutions (F_{1})');
    legend('Pareto front A','Pareto front B','Location','northwest');
    grid on;
    hold off;
    xlim([0,5]);
    ylim([0,5]);
end