function Plot1Agent(pop1)
    num=numel(pop1);
    Costs1=zeros(num,2);
    for i=1:num
        Costs1(i,:)=pop1(i).Cost;
    end
    Costs1(:,1)=-Costs1(:,1);
    Costs1(:,2)=Costs1(:,2);
    plot(Costs1(:, 1), Costs1(:, 2), 'r*', 'MarkerSize', 8);
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    title('Non-dominated Solutions (F_{1})');
end

