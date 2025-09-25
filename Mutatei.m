

function y = Mutatei(x, mu, sigma)

    nVar = numel(x);
    
    nMu = ceil(mu*nVar);

    j = randsample(nVar, nMu);
    
    % col=size(x,2);

    if numel(sigma)>1
        sigma = sigma(j);
    end
    [sc,sr]=ind2sub([size(x,2),2],j);
    % sr=fix(j/size(x,2))+1;
    % sc=j-(sr-1)*size(x,2);
    y = x;
    if sr<3
        r=sigma.*randn(size(j));
        % r=round(sigma.*randn(size(j)));
        y(sr,sc) = x(sr,sc)+r ;
    else 
        if sr==3
            y(sr,sc)=randi([1 4]);
        end
    end
    
    % [jx jy]=ind2sub(size(x),j);
    % ad=range(jx,:);
    % for i=1:size(j)
    %     y(jx,jy)=randi([ad(i,1),ad(i,2)],1);
    % end


end
