function mindex = TradeOffValue(divmatrix)
%TRADEOFFVALUE 此处显示有关此函数的摘要
%   此处显示详细说明
avgrset=zeros(size(divmatrix,1),1);
    for i=1:size(divmatrix,1)
        avgr=0;
        avgn=0;
        % s=zeros(size(divmatrix,1),2);
        for j=1:size(divmatrix,1)
            if j~=i
                s=divmatrix(i,:)-divmatrix(j,:);
                loss=0;
                gain=0;
                if s~=0
                    for k=1:size(s,2)
                        if s(k)<0
                            loss=loss+abs(s(k));
                        else
                            gain=gain+abs(s(k));
                        end
                    end
                    rij=gain/loss;
                else
                    rij=0;
                end
                if rij~=0 && ~isinf(rij)
                    avgr=avgr+rij;
                    avgn=avgn+1;
                end
            end
        end
        avgr=avgr/avgn;
        avgrset(i)=avgr;
    end
    [~,mindex]=min(avgrset);
end

