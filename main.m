clear;
clc;
close all;
% 
% % 启动机器人工具箱
% startup_rvc

global duration;

%DH参数
global alpha;
global d;
global ar;
global tool_robot;

%参数赋值
% duration = 6*pi;%螺旋线
duration = 6*pi;%花瓣线
%duration = 3;%path圆形线


%-------------------------------------------------------------UR5机械臂参数%
% ar = [0, -0.4251, -0.3921, 0, 0, 0];
% d = [0.0746, 0, 0, 0.11, 0.0947, 0.075];
d = [0.0746, 0, 0, 0.11, 0, 0];
% alpha = [pi/2, 0, 0, pi/2, -pi/2, 0];
% % alpha = [0,-pi/2,-pi/2, 0, pi/2,0];

ar = [0, -0.4251, -0.3921, 0, 0, 0];
% d = [0.0746, 0, 0, 0.11, 0.0947, 0.075];
alpha = [pi/2, 0, 0, pi/2, -pi/2, 0];


% 建立UR5机械臂模型
L(1) = Link('d', d(1),  'a', ar(1), 'alpha', alpha(1),  'standard');
L(2) = Link('d', d(2),  'a', ar(2), 'alpha', alpha(2),  'standard');
L(3) = Link('d', d(3),  'a', ar(3), 'alpha', alpha(3),  'standard');
L(4) = Link('d', d(4),  'a', ar(4), 'alpha', alpha(4),  'standard');
L(5) = Link('d', d(5),  'a', ar(5), 'alpha', alpha(5),  'standard');
L(6) = Link('d', d(6),  'a', ar(6), 'alpha', alpha(6),  'standard');
tool_robot = SerialLink(L , 'name', 'UR5');

%--------------------------------------------------------------求解znn模型%
options=odeset('RelTol',10^(-5),'AbsTol',10^(-5)); %'Refine',6,
tspan = [0 duration];

%初始误差为0
%x0 = [0,-pi/6,-pi/6,-pi/6,pi/2,0]; %path 圆形线
% x0 = [0,-2*pi/3,-pi/3,0,pi/2,0];%path1螺旋线
% x0 = [0,-pi/6,-pi/6,-pi/6,pi/2,0];%花瓣
x0=[0,0,0,0,0,0];
[t,x] = ode45(@ZNNmodel,tspan,x0,options);%矩阵返回值ode45
%记录抽样前数据
t_original = t;
x_original = x;

%-------------------------------------对t和x抽样%，适用于矩阵返回向量值ode45%
jj = 0;
epsilon = 0;
interval = 0.025;
for ii=1:length(t)
    if(t(ii,1)>=epsilon) 
        jj=jj+1;
        tn(jj,1)=t(ii,1);
        xn(jj,:)=x(ii,:);
        epsilon=tn(jj,1)+interval;  %抽样间隔
    elseif(ii==length(t)) %确保最后一个数值被抽取
        jj=jj+1;
        tn(jj,1)=t(ii,1);
        xn(jj,:)=x(ii,:);
    end
end
clear t x;
t=tn;
x=xn(:,1:6);
%---------------------------------------------------------------------画图%
% % 
% % 根据给定关节变量绘制空间中机器人模型图
% figure(1);
% tool_robot.plot(x,'trail','p');
% View_Set();
% hold on;
% 绘制末端轨迹对比



p_end_a = zeros(length(t),3);
p_end_d = zeros(length(t),3);
for i=1:length(t)
    para = x(i,:);
    para1 = t(i);
    %获取每一点的末端坐标
    p_a = fk_num(para);
    p_d = path1(para1);
    p_end_a(i,:)=p_a';
    p_end_d(i,:)=p_d';
end

figure(2);
hold on
plot3(p_end_d(:,1),p_end_d(:,2),p_end_d(:,3),'-','Color',[0 1 1],'Linewidth', 1);
plot3(p_end_a(:,1),p_end_a(:,2),p_end_a(:,3),':p','Color',[1.00 0.60 0.07],...
      'Linewidth', 1.7 ,'MarkerIndices',1:10:length(p_end_a(:,2)));
set(legend('Actual path',...
           'Desired path',...
           'Location','SouthEast'),'Interpreter','latex','FontSize',10,'fontname','Times New Roman')
xlabel('X(m)','rotation',0,'Interpreter','latex','fontname','Times New Roman');
ylabel('Y(m)','rotation',0,'Interpreter','latex','fontname','Times New Roman');
% zlabel('Z(m)','rotation',0,'Interpreter','latex','fontname','Times New Roman');
zlabel('Z(m)','Interpreter','latex','fontname','Times New Roman');
view(3); % 设置为三维视角
%绘制关节角变化曲线
figure (3);
hold on
box on
plot(t,x(:,1));
plot(t,x(:,2));
plot(t,x(:,3));
plot(t,x(:,4));
plot(t,x(:,5));
plot(t,x(:,6));
set(legend({'$\theta_1(t)$',...
           '$\theta_2(t)$',...
           '$\theta_3(t)$',...
           '$\theta_4(t)$',...
           '$\theta_5(t)$',...
           '$\theta_6(t)$'},...
           'Location','SouthEast'),'Interpreter','latex','FontSize',12,...
           'fontname','Times New Roman')
xlabel('$t(s)$','Interpreter','latex','fontname','Times New Roman');
ylabel('$\theta(t)(rad)$','rotation',0,'Interpreter','latex','fontname','Times New Roman');
       
       
%计算斜率估计关节角速度
dtheta = zeros(length(t),6);
x_vari = zeros(length(t),6);
for i1 = 1 : length(t)-2
    x_vari(i1,:) = x(i1+1,:) - x(i1,:);
    dtheta(i1,:) = x_vari(i1,:)/(t(i1+1) - t(i1));
end
dtheta(i1+1,:) = (x(i1+2,:) - x(i1+1,:))/(t(i1+2)- t(i1+1));
dtheta(i1+2,:) = dtheta(i1+1,:);

figure (4);
hold on
box on
plot(t,dtheta(:,1));
plot(t,dtheta(:,2));
plot(t,dtheta(:,3));
plot(t,dtheta(:,4));
plot(t,dtheta(:,5));
plot(t,dtheta(:,6));
set(legend({'$\dot{\theta}_1(t)$',...
           '$\dot{\theta}_2(t)$',...
           '$\dot{\theta}_3(t)$',...
           '$\dot{\theta}_4(t)$',...
           '$\dot{\theta}_5(t)$',...
           '$\dot{\theta}_6(t)$'},...
           'Location','SouthEast'),'Interpreter','latex',...
           'FontSize',12,'fontname','Times New Roman')
xlabel('$t(s)$','Interpreter','latex','fontname','Times New Roman');
ylabel('$\dot{\theta}(t)(rad/s)$','rotation',0,'Interpreter','latex','fontname','Times New Roman');

%误差范数
err = zeros(3,length(t));
nerr = zeros(length(t));
for j=1:length(t)
    T=t(j);
    err(:,j)=errfunction(T,x(j,:));
    nerr(j)=norm(err(:,j));
end
figure(5);
box on
plot(t,nerr);
set(legend('$\vert\vert \mathit e(t) \vert\vert_2$',...
           'Location','SouthEast'),'Interpreter','latex','FontSize',12,'fontname','Times New Roman')
xlabel('$t(s)$','Interpreter','latex','fontname','Times New Roman');
ylabel('$\vert\vert \mathit e(t) \vert\vert_2$','rotation',0,'Interpreter','latex','fontname','Times New Roman');

%末端轨迹在x,y,z轴上的误差
figure(6); 
box on
plot(t,err);
%ylim([-1e-3 1e-3]);
set(legend('$e_{x}(t)$',...
           '$e_{y}(t)$',...
           '$e_{z}(t)$',...
           'Location','SouthEast'),'Interpreter','latex','FontSize',15,'fontname','Times New Roman')
xlabel('$t(s)$','Interpreter','latex','fontname','Times New Roman');
ylabel('$e(t)(m)$','rotation',0,'Interpreter','latex','fontname','Times New Roman');



