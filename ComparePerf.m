load('newF1.mat');
load('newG1.mat');
load("NF.mat");
load('NG.mat');

mato=zeros(50);
for i=1:50
    for j=1:50
        mato(i,j)=CostFunction(F1(i).Position,G1(j).Position,stat1,stat2);
    end
end
