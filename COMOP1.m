
function z=COMOP1(x,y)
        
    a = 0.8;
    
    b = 3;

    n = numel(x);
    
    z1 = 1-exp(-sum((x-1/sqrt(n)).^2))+sum(-10*exp(-0.2*sqrt(y(1:end-1).^2+y(2:end).^2)));
    
    z2 = 1-exp(-sum((x+1/sqrt(n)).^2))+sum(abs(y).^a+5*(sin(y)).^b);
    
    z = [z1 z2]';

end
