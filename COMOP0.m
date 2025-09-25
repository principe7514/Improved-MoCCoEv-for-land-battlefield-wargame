function z = COMOP0(x,y,stat1,stat2,Map,Shape,Vrank,Amod,VRes,Vmod,Smod,findh,finds,findr,findam,findvr,findvm)
%COMOP 此处显示有关此函数的摘要
%   此处显示详细说明
    px=x(1:2,:)+stat1(1:2,:);
    py=y(1:2,:)+stat2(1:2,:);
    px(px<1)=1;%剔除非法值
    px(px>size(Map,1))=size(Map,1);%剔除非法值
    py(py<1)=1;%剔除非法值
    py(py>size(Map,1))=size(Map,1);%剔除非法值
    d1=sqrt(sum((px-py(:,x(3,:))).^2,1));%计算我方单位和目标的距离
    
    h1=findh(py(2,x(3,:))',py(1,x(3,:))')'-findh(px(2,:)',px(1,:)')';
    % h11=interp2(1:1:size(Map,1),1:1:size(Map,2),Map,px(1,:),px(2,:))
    % h12=interp2(1:1:size(Map,1),1:1:size(Map,2),Map,py(1,x(3,:)),py(2,x(3,:)))
    % h1=h12-h11;
    % h1=Map(py(1,x(3,:))+(py(2,x(3,:))-1)*size(Map,1))-Map(px(1,:)+(px(2,:)-1)*size(Map,1));%计算我方单位与目标的高差

    h1(h1<0)=0;
    as=~all(abs(px-stat1(1:2,:))<=0.5,1);%计算我方该回合是否行动（我方状态）
    ts=~all(abs(py-stat2(1:2,:))<=0.5,1);%计算敌方该回合是否行动（敌方状态）
    td1=stat2(4,x(3,:));%读取目标护甲类型
    % tg1=Smod(Shape(py(1,x(3,:))+(py(2,x(3,:))-1)*size(Shape,1)));%获取目标地形
    tg1=finds(py(2,x(3,:))',py(1,x(3,:))')';
    % tg1=interp2(1:1:size(Map,1),1:1:size(Map,2),Map,py(1,x(3,:)),py(2,x(3,:)),'nearest')
    z1=sum(AttackTests(d1,h1,as,x(4,:),stat1(3,:),tg1,td1,ts(x(3,:)),Vrank,Amod,VRes,Vmod,findr,findam,findvr,findvm));


    d2=sqrt(sum((py(:,:)-px(:,y(3,:))).^2,1));%敌方计算到我的距离
    h2=findh(px(2,y(3,:))',px(1,y(3,:))')'-findh(py(2,:)',py(1,:)')';
    % h21=interp2(1:1:size(Map,1),1:1:size(Map,2),Map,py(1,:),py(2,:));
    % h22=interp2(1:1:size(Map,1),1:1:size(Map,2),Map,px(1,y(3,:)),px(2,y(3,:)));
    % h2=h22-h21;
    % h2=Map(px(1,y(3,:))+(px(2,y(3,:))-1)*size(Map,1))-Map(py(1,:)+(py(2,:)-1)*size(Map,1));%敌方计算高差

    h2(h2<0)=0;
    td2=stat1(4,y(3,:));%敌方读取护甲类型
    % tg2=Smod(Shape(px(1,y(3,:))+(px(2,y(3,:))-1)*size(Shape,1)));%敌方获取目标地形
    % tg2=interp2(1:1:size(Map,1),1:1:size(Map,2),Map,px(1,y(3,:)),px(2,y(3,:)),'nearest');
    tg2=finds(px(2,y(3,:))',px(1,y(3,:))')';
    z2=sum(AttackTests(d2,h2,ts,y(4,:),stat2(3,:),tg2,td2,as(y(3,:)),Vrank,Amod,VRes,Vmod,findr,findam,findvr,findvm));
    z=[-z1,z2]';
end

