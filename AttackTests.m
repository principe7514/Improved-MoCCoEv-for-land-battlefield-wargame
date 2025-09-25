function z=AttackTests(d,h,s,w,an,tg,td,ts,Vrank,Amod,VRes,Vmod,findr,findam,findvr,findvm)
%AttackTests d距离，h高差，s状态，w采用武器，an攻击者班数，tg目标地形，td目标护甲类型，ts目标状态
    %处理超过20格的射击
    overdis=all(d<=20,1);
    d=overdis.*d;
    d(d<1)=1;
    % rank=overdis.*(Vrank(sub2ind(size(Vrank),w,d+1)));%计算攻击等级
    
    % rank=interp1(0:1:20,Vrank,d)
    rank=findr(d);%计算攻击等级（插值）
    % modrank=rank+overdis.*(Amod(sub2ind(size(Amod),h+1,d)));%附加高差修正
    % modrank=rank+interp2(1:1:20,0:1:8,Amod,d,h)%附加高差修正（插值）
    % modrank=rank+findam(h,d);%附加高差修正（插值）
    modrank=rank;%无高程测试
    %处理攻击等级为0或更低
    zeroat=~all(modrank<=1,1);
    modrank(modrank<=1)=1;
    %进行主要程序
    % at=VRes(sub2ind(size(VRes),an,modrank));%计算战斗结果
    % at=interp2(1:1:10,1:1:5,VRes,modrank,an)%计算战斗结果（插值）
    at=findvr(an,modrank);%计算战斗结果（插值）
    % a=s+tg+ts*2-1;
    % if any(a(:)==0)
    %         sp=0;
    % end
    % vm=Vmod(sub2ind(size(Vmod),s+tg+ts*2+1,td));%计算修正结果
    % vm=interp2(1:1:5,0:1:7,Vmod,td,s+tg+ts*2)%计算修正结果（插值）
    as=s+tg+ts.*2-1;
    vm=findvm(as',td')';
    z=zeroat.*(at+vm);%得到最终的值
    z(z<0)=0;%0值处理
end

