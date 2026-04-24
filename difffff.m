% syms t;
% 
% % 定义函数 p
% p1 = 0.6796;
% p2 = 0.1 * sin(t) * (exp(cos(t)) - sin(t/12)^5 - 2*cos(4*t)) - 0.11;
% p3 = 0.1 * cos(t) * (exp(cos(t)) - sin(t/12)^5 - 2*cos(4*t)) + 0.4174;
% 
% % 对各个分量求导
% dp1_dt = diff(p1, t);  % p1的导数
% dp2_dt = diff(p2, t);  % p2的导数
% dp3_dt = diff(p3, t);  % p3的导数
% 
% % 将导数存入dp向量
% dp = [dp1_dt; dp2_dt; dp3_dt];
% 
% % 显示结果
% disp('dp/dt = ');
% disp(dp);
