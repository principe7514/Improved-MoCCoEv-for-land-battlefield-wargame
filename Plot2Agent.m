
function Plot2Agent(pop1,pop2)

    Costs1 = [pop1.Cost];
    Costs1(1,:)=-Costs1(1,:);
    Costs1(2,:)=Costs1(2,:)+5;
    Costs2 = [pop2.Cost];
    Costs2(1,:)=-Costs2(1,:);
    Costs2(2,:)=Costs2(2,:)+5;

    plot(Costs1(1, :), Costs1(2, :), 'r*',Costs2(1, :), Costs2(2, :), 'b*', 'MarkerSize', 8);
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    title('Non-dominated Solutions (F_{1})');
    legend('Pareto front A','Pareto front B','Location','northeast');
    grid on;
    hold off;
end
