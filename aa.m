clear;
clc;
close all;

global duration;
%DH参数
global alpha;
global ar;
global d;
global tool_robot;

%参数赋值
duration = 6*pi; % 花瓣线/nengpao

%-------------------------------------------------------------UR5机械臂参数%
ar = [0, -0.4251, -0.3921, 0, 0, 0];
d = [0.0746, 0, 0, 0.11, 0, 0];
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
options=odeset('RelTol',10^(-5),'AbsTol',10^(-5)); 
tspan = [0 duration];

% === 修改点 1：使用 fsolve 自动求解完美的初始关节角消除初始误差 ===
disp('正在计算匹配路径起点的初始关节角...');
p_start = path1(0); % 获取路径的物理起点
q_guess = [0, -pi/3, -pi/3, -pi/6, pi/2, 0]; % 提供一个接近的初始猜测值避免奇异点
options_fsolve = optimset('Display','off','TolFun',1e-8,'TolX',1e-8);
x0 = fsolve(@(q) fk_num(q) - p_start, q_guess, options_fsolve);
% x0=[0,0,0,0,0,0]; % 注释掉原硬编码

[t,x] = ode45(@ZNNmodel,tspan,x0,options);%矩阵返回值ode45

%记录抽样前数据
t_original = t;
x_original = x;

%-------------------------------------对t和x抽样%
jj = 0;
epsilon = 0;
interval = 0.025;
for ii=1:length(t)
    if(t(ii,1)>=epsilon) 
        jj=jj+1;
        tn(jj,1)=t(ii,1);
        xn(jj,:)=x(ii,:);
        epsilon=tn(jj,1)+interval;  
    elseif(ii==length(t)) 
        jj=jj+1;
        tn(jj,1)=t(ii,1);
        xn(jj,:)=x(ii,:);
    end
end
clear t x;
t=tn;
x=xn(:,1:6);

%---------------------------------------------------------------------画图%
p_end_a = zeros(length(t),3);
p_end_d = zeros(length(t),3);
for i=1:length(t)
    para = x(i,:);
    para1 = t(i);
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
set(legend('Actual path','Desired path','Location','SouthEast'),...
    'Interpreter','latex','FontSize',10,'fontname','Times New Roman')
xlabel('X(m)','rotation',0,'Interpreter','latex','fontname','Times New Roman');
ylabel('Y(m)','rotation',0,'Interpreter','latex','fontname','Times New Roman');
zlabel('Z(m)','Interpreter','latex','fontname','Times New Roman');
view(3); 

figure (3);
hold on
box on
plot(t,x(:,1)); plot(t,x(:,2)); plot(t,x(:,3));
plot(t,x(:,4)); plot(t,x(:,5)); plot(t,x(:,6));
set(legend({'$\theta_1(t)$','$\theta_2(t)$','$\theta_3(t)$',...
           '$\theta_4(t)$','$\theta_5(t)$','$\theta_6(t)$'},...
           'Location','SouthEast'),'Interpreter','latex','FontSize',12,'fontname','Times New Roman')
xlabel('$t(s)$','Interpreter','latex','fontname','Times New Roman');
ylabel('$\theta(t)(rad)$','rotation',0,'Interpreter','latex','fontname','Times New Roman');
       
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
plot(t,dtheta(:,1)); plot(t,dtheta(:,2)); plot(t,dtheta(:,3));
plot(t,dtheta(:,4)); plot(t,dtheta(:,5)); plot(t,dtheta(:,6));
set(legend({'$\dot{\theta}_1(t)$','$\dot{\theta}_2(t)$','$\dot{\theta}_3(t)$',...
           '$\dot{\theta}_4(t)$','$\dot{\theta}_5(t)$','$\dot{\theta}_6(t)$'},...
           'Location','SouthEast'),'Interpreter','latex','FontSize',12,'fontname','Times New Roman')
xlabel('$t(s)$','Interpreter','latex','fontname','Times New Roman');
ylabel('$\dot{\theta}(t)(rad/s)$','rotation',0,'Interpreter','latex','fontname','Times New Roman');

err = zeros(3,length(t));
nerr = zeros(length(t),1);
for j=1:length(t)
    T=t(j);
    err(:,j)=errfunction(T,x(j,:));
    nerr(j)=norm(err(:,j));
end

figure(5);
box on
plot(t,nerr);
set(legend('$\vert\vert \mathit e(t) \vert\vert_2$','Location','SouthEast'),...
    'Interpreter','latex','FontSize',12,'fontname','Times New Roman')
xlabel('$t(s)$','Interpreter','latex','fontname','Times New Roman');
ylabel('$\vert\vert \mathit e(t) \vert\vert_2$','rotation',0,'Interpreter','latex','fontname','Times New Roman');

figure(6); 
box on
plot(t,err);
set(legend('$e_{x}(t)$','$e_{y}(t)$','$e_{z}(t)$','Location','SouthEast'),...
    'Interpreter','latex','FontSize',15,'fontname','Times New Roman')
xlabel('$t(s)$','Interpreter','latex','fontname','Times New Roman');
ylabel('$e(t)(m)$','rotation',0,'Interpreter','latex','fontname','Times New Roman');


function p = path1(t)
% nengpao 曲线 (作为当前激活的路径)
p=[-0.8172;-0.15+0.15*cos((pi/3)*t)*(2-2^(sin((5*pi/3)*t)));0.15*sin((pi/3)*t)*(2-2^(sin((5*pi/3)*t)))+0.0201];
end

function dp = dpath1(t)
% nengpao 曲线导数
dp=[0;-(pi*sin((pi*t)/3)*(2^sin((5*pi*t)/3) - 2))/20 - (2^sin((5*pi*t)/3)*pi*log(2)*cos((pi*t)/3)*cos((5*pi*t)/3))/4;(pi*cos((pi*t)/3)*(2^sin((5*pi*t)/3) - 2))/20 - (2^sin((5*pi*t)/3)*pi*log(2)*cos((5*pi*t)/3)*sin((pi*t)/3))/4];
end

function J = jaco_num(theta)
theta1 = theta(1);theta2 = theta(2);theta3 = theta(3);
theta4 = theta(4);theta5 = theta(5);theta6 = theta(6);
J = ...
[ (11*cos(theta1))/100 + (3*cos(theta1)*cos(theta5))/40 + (4251*cos(theta2)*sin(theta1))/10000 - (947*cos(theta4)*(cos(theta2)*sin(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2)))/10000 + (947*sin(theta4)*(sin(theta1)*sin(theta2)*sin(theta3) - cos(theta2)*cos(theta3)*sin(theta1)))/10000 - (3*sin(theta5)*(cos(theta4)*(sin(theta1)*sin(theta2)*sin(theta3) - cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta2)*sin(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2))))/40 - (3921*sin(theta1)*sin(theta2)*sin(theta3))/10000 + (3921*cos(theta2)*cos(theta3)*sin(theta1))/10000, (4251*cos(theta1)*sin(theta2))/10000 - (947*cos(theta4)*(cos(theta1)*sin(theta2)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))/10000 - (947*sin(theta4)*(cos(theta1)*cos(theta2)*sin(theta3) + cos(theta1)*cos(theta3)*sin(theta2)))/10000 + (3*sin(theta5)*(cos(theta4)*(cos(theta1)*cos(theta2)*sin(theta3) + cos(theta1)*cos(theta3)*sin(theta2)) - sin(theta4)*(cos(theta1)*sin(theta2)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3))))/40 + (3921*cos(theta1)*cos(theta2)*sin(theta3))/10000 + (3921*cos(theta1)*cos(theta3)*sin(theta2))/10000, (3*sin(theta5)*(cos(theta4)*(cos(theta1)*cos(theta2)*sin(theta3) + cos(theta1)*cos(theta3)*sin(theta2)) - sin(theta4)*(cos(theta1)*sin(theta2)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3))))/40 - (947*sin(theta4)*(cos(theta1)*cos(theta2)*sin(theta3) + cos(theta1)*cos(theta3)*sin(theta2)))/10000 - (947*cos(theta4)*(cos(theta1)*sin(theta2)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))/10000 + (3921*cos(theta1)*cos(theta2)*sin(theta3))/10000 + (3921*cos(theta1)*cos(theta3)*sin(theta2))/10000, (3*sin(theta5)*(cos(theta4)*(cos(theta1)*cos(theta2)*sin(theta3) + cos(theta1)*cos(theta3)*sin(theta2)) - sin(theta4)*(cos(theta1)*sin(theta2)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3))))/40 - (947*sin(theta4)*(cos(theta1)*cos(theta2)*sin(theta3) + cos(theta1)*cos(theta3)*sin(theta2)))/10000 - (947*cos(theta4)*(cos(theta1)*sin(theta2)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))/10000, (3*cos(theta5)*(cos(theta4)*(cos(theta1)*sin(theta2)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta1)*cos(theta2)*sin(theta3) + cos(theta1)*cos(theta3)*sin(theta2))))/40 - (3*sin(theta1)*sin(theta5))/40, 0;...
  (11*sin(theta1))/100 - (4251*cos(theta1)*cos(theta2))/10000 + (3*cos(theta5)*sin(theta1))/40 + (947*cos(theta4)*(cos(theta1)*cos(theta2)*sin(theta3) + cos(theta1)*cos(theta3)*sin(theta2)))/10000 - (947*sin(theta4)*(cos(theta1)*sin(theta2)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))/10000 + (3*sin(theta5)*(cos(theta4)*(cos(theta1)*sin(theta2)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta1)*cos(theta2)*sin(theta3) + cos(theta1)*cos(theta3)*sin(theta2))))/40 + (3921*cos(theta1)*sin(theta2)*sin(theta3))/10000 - (3921*cos(theta1)*cos(theta2)*cos(theta3))/10000, (4251*sin(theta1)*sin(theta2))/10000 - (947*cos(theta4)*(sin(theta1)*sin(theta2)*sin(theta3) - cos(theta2)*cos(theta3)*sin(theta1)))/10000 - (947*sin(theta4)*(cos(theta2)*sin(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2)))/10000 + (3*sin(theta5)*(cos(theta4)*(cos(theta2)*sin(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2)) - sin(theta4)*(sin(theta1)*sin(theta2)*sin(theta3) - cos(theta2)*cos(theta3)*sin(theta1))))/40 + (3921*cos(theta2)*sin(theta1)*sin(theta3))/10000 + (3921*cos(theta3)*sin(theta1)*sin(theta2))/10000, (3*sin(theta5)*(cos(theta4)*(cos(theta2)*sin(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2)) - sin(theta4)*(sin(theta1)*sin(theta2)*sin(theta3) - cos(theta2)*cos(theta3)*sin(theta1))))/40 - (947*sin(theta4)*(cos(theta2)*sin(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2)))/10000 - (947*cos(theta4)*(sin(theta1)*sin(theta2)*sin(theta3) - cos(theta2)*cos(theta3)*sin(theta1)))/10000 + (3921*cos(theta2)*sin(theta1)*sin(theta3))/10000 + (3921*cos(theta3)*sin(theta1)*sin(theta2))/10000, (3*sin(theta5)*(cos(theta4)*(cos(theta2)*sin(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2)) - sin(theta4)*(sin(theta1)*sin(theta2)*sin(theta3) - cos(theta2)*cos(theta3)*sin(theta1))))/40 - (947*sin(theta4)*(cos(theta2)*sin(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2)))/10000 - (947*cos(theta4)*(sin(theta1)*sin(theta2)*sin(theta3) - cos(theta2)*cos(theta3)*sin(theta1)))/10000, (3*cos(theta1)*sin(theta5))/40 + (3*cos(theta5)*(cos(theta4)*(sin(theta1)*sin(theta2)*sin(theta3) - cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta2)*sin(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2))))/40, 0;...
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                     (3921*sin(theta2)*sin(theta3))/10000 - (3921*cos(theta2)*cos(theta3))/10000 - (4251*cos(theta2))/10000 - (3*sin(theta5)*(cos(theta4)*(cos(theta2)*cos(theta3) - sin(theta2)*sin(theta3)) - sin(theta4)*(cos(theta2)*sin(theta3) + cos(theta3)*sin(theta2))))/40 + (947*cos(theta4)*(cos(theta2)*sin(theta3) + cos(theta3)*sin(theta2)))/10000 + (947*sin(theta4)*(cos(theta2)*cos(theta3) - sin(theta2)*sin(theta3)))/10000,                                                                                                                         (3921*sin(theta2)*sin(theta3))/10000 - (3921*cos(theta2)*cos(theta3))/10000 - (3*sin(theta5)*(cos(theta4)*(cos(theta2)*cos(theta3) - sin(theta2)*sin(theta3)) - sin(theta4)*(cos(theta2)*sin(theta3) + cos(theta3)*sin(theta2))))/40 + (947*cos(theta4)*(cos(theta2)*sin(theta3) + cos(theta3)*sin(theta2)))/10000 + (947*sin(theta4)*(cos(theta2)*cos(theta3) - sin(theta2)*sin(theta3)))/10000,                                                                                                 (947*cos(theta4)*(cos(theta2)*sin(theta3) + cos(theta3)*sin(theta2)))/10000 - (3*sin(theta5)*(cos(theta4)*(cos(theta2)*cos(theta3) - sin(theta2)*sin(theta3)) - sin(theta4)*(cos(theta2)*sin(theta3) + cos(theta3)*sin(theta2))))/40 + (947*sin(theta4)*(cos(theta2)*cos(theta3) - sin(theta2)*sin(theta3)))/10000,                                                                                 -(3*cos(theta5)*(cos(theta4)*(cos(theta2)*sin(theta3) + cos(theta3)*sin(theta2)) + sin(theta4)*(cos(theta2)*cos(theta3) - sin(theta2)*sin(theta3))))/40, 0];
end

function position = fk_num(theta)
theta1 = theta(1);theta2 = theta(2);theta3 = theta(3);
theta4 = theta(4);theta5 = theta(5);theta6 = theta(6);
p1 = (11*sin(theta1))/100 - (4251*cos(theta1)*cos(theta2))/10000 + (3*cos(theta5)*sin(theta1))/40 + (947*cos(theta4)*(cos(theta1)*cos(theta2)*sin(theta3) + cos(theta1)*cos(theta3)*sin(theta2)))/10000 - (947*sin(theta4)*(cos(theta1)*sin(theta2)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)))/10000 + (3*sin(theta5)*(cos(theta4)*(cos(theta1)*sin(theta2)*sin(theta3) - cos(theta1)*cos(theta2)*cos(theta3)) + sin(theta4)*(cos(theta1)*cos(theta2)*sin(theta3) + cos(theta1)*cos(theta3)*sin(theta2))))/40 + (3921*cos(theta1)*sin(theta2)*sin(theta3))/10000 - (3921*cos(theta1)*cos(theta2)*cos(theta3))/10000;
p2 = (947*cos(theta4)*(cos(theta2)*sin(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2)))/10000 - (3*cos(theta1)*cos(theta5))/40 - (4251*cos(theta2)*sin(theta1))/10000 - (11*cos(theta1))/100 - (947*sin(theta4)*(sin(theta1)*sin(theta2)*sin(theta3) - cos(theta2)*cos(theta3)*sin(theta1)))/10000 + (3*sin(theta5)*(cos(theta4)*(sin(theta1)*sin(theta2)*sin(theta3) - cos(theta2)*cos(theta3)*sin(theta1)) + sin(theta4)*(cos(theta2)*sin(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2))))/40 + (3921*sin(theta1)*sin(theta2)*sin(theta3))/10000 - (3921*cos(theta2)*cos(theta3)*sin(theta1))/10000;
p3 = (947*sin(theta4)*(cos(theta2)*sin(theta3) + cos(theta3)*sin(theta2)))/10000 - (3921*cos(theta2)*sin(theta3))/10000 - (3921*cos(theta3)*sin(theta2))/10000 - (3*sin(theta5)*(cos(theta4)*(cos(theta2)*sin(theta3) + cos(theta3)*sin(theta2)) + sin(theta4)*(cos(theta2)*cos(theta3) - sin(theta2)*sin(theta3))))/40 - (947*cos(theta4)*(cos(theta2)*cos(theta3) - sin(theta2)*sin(theta3)))/10000 - (4251*sin(theta2))/10000 + 373/5000;
position = [p1;p2;p3];
end

function err = errfunction(t,theta)
p = fk_num(theta);
err = p - path1(t);
end

function output = ZNNmodel(t,x)
z=x(1:6);
J = jaco_num(z);
dp = dpath1(t);
err = errfunction(t,z);

ksi=0.01;
absErr = abs(err);
absErr2 = absErr;
signErr2 = 1000*atan(absErr2/ksi);
lambda = signErr2(1:end);
lambda = diag(lambda);

alpha = 0.5;
gamma1 = lambda*tanh(alpha*abs(err));

% === 修改点 2：修复 ZNN 演化方程的符号，正确叠加期望速度 ===
% 原代码导致系统背离轨迹: dth = -pinv(J)*(gamma1.*AFMnew_3(err) + dp);
dth = pinv(J) * (dp - gamma1 .* AFMnew_3(err)); 

output = dth;
end

function y = AFMnew_3(L_ij)
    y= arrayfun(@(x) exp(abs(x).^0.5).*((log(abs(x) + 1)).^0.5).*sign(x)./0.5 + x+sign(x),L_ij);
end