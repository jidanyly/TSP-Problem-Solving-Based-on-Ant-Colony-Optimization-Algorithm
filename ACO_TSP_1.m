%%%%%%%%%%%%%%% 蚁群算法解决 TSP 问题I %%%%%%%%%%%%%%%
%% 1.清空环境变量
clear;
close all;
clc;
%% 2.导入数据
citys = xlsread('Problem_I');
%% 3.计算城市间相互距离
n=size(citys,1);
D=zeros(n,n);
for i=1:n
    for j=1:n
        if i~=j
            D(i,j)=sqrt(sum((citys(i,:)-citys(j,:)).^2));
        else
            D(i,j)=eps;
        end
    end
end
%% 4.初始化参数
m=40;                         %蚂蚁数量
alpha=1;                      %信息素重要程度因子
beta=5;                       %启发函数重要程度因子
rho=0.2;                      %信息素挥发因子
Q=10;                         %常系数
Eta=1./D;                     %启发函数
Tau=ones(n,n);                %信息素矩阵
Table=zeros(m,n);             %路径记录表
iter=1;                       %迭代次数初值
itermax=200;                  %最大迭代次数
Route_best=zeros(itermax,n);  %各代最佳路径
Length_best=zeros(itermax,1); %各代最佳路径的长度
Length_ave=zeros(itermax,1);  %各代路径的平均长度
%% 5.迭代寻找最佳路径
tic    
while iter<=itermax
    % 5.1随机产生各个蚂蚁的起点城市
start=zeros(m,1);
    for i=1:m
        temp=randperm(n);
        start(i,:)=temp(1);
    end
    Table(:,1)=start;
    citys_index=1:n;
    % 5.2产生解空间（路径表）
    %逐个蚂蚁路径选择
    for i=1:m
        %逐个城市路径选择
        for j=2:n
            tabu=Table(i,:); %已访问城市
            allow_index=~ismember(citys_index,tabu);
            allow=citys_index(allow_index);
            P=allow;
            % 5.3计算城市间转移概率
            for k=1:length(allow) 
                P(k)=Tau(Table(i,j-1),allow(k))^alpha * Eta(Table(i,j-1),allow(k))^beta;
            end
            P=P/sum(P);
            % 5.4轮盘赌法选择下一个访问城市
            Pc=cumsum(P);
            target_index=find(Pc>=rand);
            target=allow(target_index(1));
            Table(i,j)=target;
        end
    end
    % 5.5计算各个蚂蚁的路径距离
    Length=zeros(m,1);
    for i=1:m
        for j=1:n-1
            Length(i,1)=Length(i,1)+D(Table(i,j),Table(i,j+1));
        end
        Length(i,1)=Length(i,1)+D(Table(i,end),Table(i,1));
    end
    % 5.6计算最短路径及平均距离
    if iter==1
        [min_Length,min_index]=min(Length);
        Length_best(iter)=min_Length;
        Route_best(iter,:)=Table(min_index,:);
        Length_ave(iter)=mean(Length);
    else
        [min_Length,min_index]=min(Length);
        Length_best(iter)=min(Length_best(iter-1),min_Length);
        Length_ave(iter)=mean(Length);
        if Length_best(iter-1)>min_Length
            Route_best(iter,:)=Table(min_index,:);
        else
            Route_best(iter,:)=Route_best(iter-1,:);
        end
    end
    % 5.7更新信息素
    Delta_Tau=zeros(n,n);
    %逐个蚂蚁计算
    for i=1:m
        %逐个城市计算
        for j=1:n-1
            Delta_Tau(Table(i,j),Table(i,j+1))=Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
        end
        Delta_Tau(Table(i,end),Table(i,1))=Delta_Tau(Table(i,end),Table(i,1)) + Q/Length(i);
    end
    Tau=(1-rho)*Tau+Delta_Tau;
    % 5.8迭代次数加1，清空路径表
    iter=iter+1;
    Table=zeros(m,n);
    [iter min_Length]
end
toc
%% 6.结果显示
best_route=Route_best(end,:);
best_length=Length_best(end,:);
disp(['最短距离: ' num2str(best_length)]);
disp(['最短路径顺序: ' num2str(best_route)]);
%% 7.绘图
figure(1)
plot([citys(best_route,1);citys(best_route(1),1)],[citys(best_route,2);citys(best_route(1),2)],'g-','linewidth',1.5)
hold on;
plot([citys(best_route,1);citys(best_route(1),1)],[citys(best_route,2);citys(best_route(1),2)],'bo','MarkerFaceColor','b')
hold on;
plot(citys(best_route(1),1),citys(best_route(1),2),'ro','MarkerFaceColor','r');
for i=1:size(citys,1)
    text(citys(i,1),citys(i,2),[' ' num2str(i)],'color','k');
end
text(citys(best_route(1),1),citys(best_route(1),2),'    起点','color','r');
text(citys(best_route(end),1),citys(best_route(end),2),'    终点','color','r');
xlabel('城市位置横坐标')
ylabel('城市位置纵坐标')
title(['蚁群算法TSP问题优化路径(最短距离：' num2str(best_length) ')'])
figure(2)
plot(1:itermax,Length_ave,'b-',1:itermax,Length_best,'r-','linewidth',1)
legend('平均距离','最短距离')
xlabel('迭代次数')
ylabel('目标函数值 & 路程距离')
title('各代最短距离与平均距离对比')
