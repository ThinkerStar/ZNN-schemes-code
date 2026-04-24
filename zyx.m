
delta = input('请输入 delta 的值: '); % 用户输入 delta
t1 = linspace(0, 1, 1000); % 在 [0, 1] 区间生成 100 个 t1 值

f_t1 =(1./(t1+1)) .* ((25*t1 + 71 + 39*delta) / 51) .* ((t1 + 11 - 27*delta) / 17) ...
             - (1 ./ (2*(t1+1))) .* ((t1 + 11 - 27*delta) / 17).^2 ...
             + (1 ./ (t1+1)) .* ((5*t1 + 4 + 18*delta) / 17).^2 ...
             + 1 - (1 ./ t1);


% 绘制函数图像
figure;
plot(t1, f_t1, 'LineWidth', 2);
xlabel('t1');
ylabel('f(t1)');
title(['f(t1) 的曲线，\delta = ' num2str(delta)]);
hold on
plot(t1+2, f_t1, 'LineWidth', 2);
plot(t1+1, f_t1, 'LineWidth', 2);
hold off
figure;
plot(t1+2, f_t1, 'LineWidth', 2);
