function output = ZNNmodel(t,x)

z=x(1:6);%关节变量
%y=x(7:9);%积分项

%雅克比矩阵
J = jaco_num(z);

%预设末端轨迹 螺旋
dp = dpath1(t);

%误差r3
err = errfunction(t,z);%计算误差
ksi=0.01;
% 向量化形式变参vec
absErr = abs(err);
absErr2 = absErr;
signErr2 = 1000*atan(absErr2/ksi);
% signErr2 = ksi*cosh(absErr2/ksi);
lambda = signErr2(1:end);
lambda = diag(lambda);

%CMR变参
% lambda = 1;
% alpha = 0.1;
% p = 0.5;
% eta = 2;

%LJW变参
% lambda = 1000;%1000
alpha = 0.5;%1

gamma1 = lambda*tanh(alpha*abs(err));

dth = -pinv(J)*(gamma1.*AFMnew_3(err) + dp);%+1干扰项

%输出
output = dth;

t

end