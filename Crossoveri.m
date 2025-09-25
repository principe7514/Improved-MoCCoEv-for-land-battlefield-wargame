


function [y1, y2] = Crossoveri(x1, x2)

    alpha = randi([0 1],size(x1));
    beta = rand(size(x1));
    alpha(1:2,:)=beta(1:2,:);
    y1 = alpha.*x1+(1-alpha).*x2;
    y2 = alpha.*x2+(1-alpha).*x1;
    
end

