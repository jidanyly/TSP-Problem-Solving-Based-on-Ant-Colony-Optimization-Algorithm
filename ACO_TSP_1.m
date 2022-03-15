%%%%%%%%%%%%%%% ��Ⱥ�㷨��� TSP ����I %%%%%%%%%%%%%%%
%% 1.��ջ�������
clear;
close all;
clc;
%% 2.��������
citys = xlsread('Problem_I');
%% 3.������м��໥����
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
%% 4.��ʼ������
m=40;                         %��������
alpha=1;                      %��Ϣ����Ҫ�̶�����
beta=5;                       %����������Ҫ�̶�����
rho=0.2;                      %��Ϣ�ػӷ�����
Q=10;                         %��ϵ��
Eta=1./D;                     %��������
Tau=ones(n,n);                %��Ϣ�ؾ���
Table=zeros(m,n);             %·����¼��
iter=1;                       %����������ֵ
itermax=200;                  %����������
Route_best=zeros(itermax,n);  %�������·��
Length_best=zeros(itermax,1); %�������·���ĳ���
Length_ave=zeros(itermax,1);  %����·����ƽ������
%% 5.����Ѱ�����·��
tic    
while iter<=itermax
    % 5.1��������������ϵ�������
start=zeros(m,1);
    for i=1:m
        temp=randperm(n);
        start(i,:)=temp(1);
    end
    Table(:,1)=start;
    citys_index=1:n;
    % 5.2������ռ䣨·����
    %�������·��ѡ��
    for i=1:m
        %�������·��ѡ��
        for j=2:n
            tabu=Table(i,:); %�ѷ��ʳ���
            allow_index=~ismember(citys_index,tabu);
            allow=citys_index(allow_index);
            P=allow;
            % 5.3������м�ת�Ƹ���
            for k=1:length(allow) 
                P(k)=Tau(Table(i,j-1),allow(k))^alpha * Eta(Table(i,j-1),allow(k))^beta;
            end
            P=P/sum(P);
            % 5.4���̶ķ�ѡ����һ�����ʳ���
            Pc=cumsum(P);
            target_index=find(Pc>=rand);
            target=allow(target_index(1));
            Table(i,j)=target;
        end
    end
    % 5.5����������ϵ�·������
    Length=zeros(m,1);
    for i=1:m
        for j=1:n-1
            Length(i,1)=Length(i,1)+D(Table(i,j),Table(i,j+1));
        end
        Length(i,1)=Length(i,1)+D(Table(i,end),Table(i,1));
    end
    % 5.6�������·����ƽ������
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
    % 5.7������Ϣ��
    Delta_Tau=zeros(n,n);
    %������ϼ���
    for i=1:m
        %������м���
        for j=1:n-1
            Delta_Tau(Table(i,j),Table(i,j+1))=Delta_Tau(Table(i,j),Table(i,j+1)) + Q/Length(i);
        end
        Delta_Tau(Table(i,end),Table(i,1))=Delta_Tau(Table(i,end),Table(i,1)) + Q/Length(i);
    end
    Tau=(1-rho)*Tau+Delta_Tau;
    % 5.8����������1�����·����
    iter=iter+1;
    Table=zeros(m,n);
    [iter min_Length]
end
toc
%% 6.�����ʾ
best_route=Route_best(end,:);
best_length=Length_best(end,:);
disp(['��̾���: ' num2str(best_length)]);
disp(['���·��˳��: ' num2str(best_route)]);
%% 7.��ͼ
figure(1)
plot([citys(best_route,1);citys(best_route(1),1)],[citys(best_route,2);citys(best_route(1),2)],'g-','linewidth',1.5)
hold on;
plot([citys(best_route,1);citys(best_route(1),1)],[citys(best_route,2);citys(best_route(1),2)],'bo','MarkerFaceColor','b')
hold on;
plot(citys(best_route(1),1),citys(best_route(1),2),'ro','MarkerFaceColor','r');
for i=1:size(citys,1)
    text(citys(i,1),citys(i,2),[' ' num2str(i)],'color','k');
end
text(citys(best_route(1),1),citys(best_route(1),2),'    ���','color','r');
text(citys(best_route(end),1),citys(best_route(end),2),'    �յ�','color','r');
xlabel('����λ�ú�����')
ylabel('����λ��������')
title(['��Ⱥ�㷨TSP�����Ż�·��(��̾��룺' num2str(best_length) ')'])
figure(2)
plot(1:itermax,Length_ave,'b-',1:itermax,Length_best,'r-','linewidth',1)
legend('ƽ������','��̾���')
xlabel('��������')
ylabel('Ŀ�꺯��ֵ & ·�̾���')
title('������̾�����ƽ������Ա�')
