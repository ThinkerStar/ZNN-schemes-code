% % 定义符号变量 t
% syms t;
% 
% % 定义 p 向量
% % p1 = 0;
% % p2 = 0.15 * cos((pi/3) * t) * (2 - 2^(sin((5*pi/3) * t)));
% % p3 = 0.15 * sin((pi/3) * t) * (2 - 2^(sin((5*pi/3) * t)));
% p1=0.6796;
% p2=0.1*sin(t)*(exp(cos(t))-sin(t/12)^5-2*cos(4*t))-0.11;
% p3=0.1*cos(t)*(exp(cos(t))-sin(t/12)^5-2*cos(4*t))+0.4174;
% % 创建 t 的数值范围从 0 到 6pi
% t_vals = linspace(0, 6*pi, 1000);  % 生成 1000 个点的 t 值
% 
% % 计算 p1, p2, p3 的值
% p1_vals = double(subs(p1, t, t_vals));
% p2_vals = double(subs(p2, t, t_vals));
% p3_vals = double(subs(p3, t, t_vals));
% 
% % 绘制 3D 曲线
% figure;
% plot3(p2_vals, p3_vals, p1_vals, 'LineWidth', 2);
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('Plot of p(t) from 0 to 6\pi');
% grid on;
