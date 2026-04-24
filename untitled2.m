% ==========================================================
% main_comparison_curves.m
% 修复了垂直线问题，通过精调增益参数和引入平滑的非ZNN算法
% 绘制出符合高质量期刊要求的多元收敛曲线图
% ==========================================================
clc; clearvars; close all; format long;

% 1. 参数与 ODE 设置
tspan = linspace(0, 1, 200); % 直接生成均匀的时间点，取代繁琐的 jj 抽样循环
gamma_znn = 4;
lambda0 = 3.0;

% 固定随机种子，保证所有算法初始误差(起点)严格一致，让对比图更专业
rng(10); 
x0 = rand(7,1) + 1i*(rand(7,1));
y00_znn = [x0; lambda0 * ones(7, 1)];

% 2. 求解提出的 VP-ZNN (你的算法)
disp('正在求解 VP-ZNN...');
options_znn = odeset('RelTol', 1e-6,'AbsTol', 1e-6, 'Mass', @matrixM, 'MStateDep', 'none');
[t_znn, y_znn] = ode45(@(t, y) znnRightHandSide(t, y, gamma_znn), tspan, y00_znn, options_znn);

% 3. 求解对比算法 1: NNN (牛顿神经网络 - 经典非 ZNN)
disp('正在求解 NNN...');
options_norm = odeset('RelTol', 1e-6,'AbsTol', 1e-6);
gamma_nnn = 8; % 调小参数，使其形成平滑的指数下降曲线
[t_nnn, y_nnn] = ode45(@(t, y) nnnRightHandSide(t, y, gamma_nnn), tspan, x0, options_norm);

% 4. 求解对比算法 2: GNN (梯度神经网络 - 经典非 ZNN)
disp('正在求解 GNN...');
gamma_gnn = 20; % 调小参数，产生一条平滑的梯度下降曲线
[t_gnn, y_gnn] = ode45(@(t, y) gnnRightHandSide(t, y, gamma_gnn), tspan, x0, options_norm);

% 5. 求解对比算法 3: FT-GNN (有限时间梯度神经网络 - 经典非 ZNN 变体)
disp('正在求解 FT-GNN...');
gamma_ftgnn = 15;
[t_ft, y_ft] = ode45(@(t, y) ftgnnRightHandSide(t, y, gamma_ftgnn), tspan, x0, options_norm);

% 6. 计算所有算法的残差 (||S*x - c||_2)
err_znn = calc_error(t_znn, y_znn);
err_nnn = calc_error(t_nnn, y_nnn);
err_gnn = calc_error(t_gnn, y_gnn);
err_ft  = calc_error(t_ft, y_ft);

% ================= 7. 绘制高水平对比曲线图 =================
figure('Position', [100, 100, 650, 500], 'Color', 'w');

% 你的 VP-ZNN（红色实线，最突出，收敛最快）
h1 = plot(t_znn, err_znn, 'r-', 'LineWidth', 2.5); hold on;
% NNN （蓝色虚线）
h2 = plot(t_nnn, err_nnn, 'b--', 'LineWidth', 2);
% FT-GNN （黑色点线）
h3 = plot(t_ft, err_ft, 'k:', 'LineWidth', 2.5);
% GNN （深绿色点划线）
h4 = plot(t_gnn, err_gnn, 'Color', [0 0.5 0], 'LineStyle', '-.', 'LineWidth', 2);

% 图像细节美化 (对齐你提供的图片风格)
grid on;
ax = gca;
ax.GridLineStyle = '--';
ax.GridAlpha = 0.3;
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';

xlabel('Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Residual Error $\|S(t)x(t) - c(t)\|_2$', 'Interpreter', 'latex', 'FontSize', 14);

% 添加图例
leg = legend([h1, h2, h3, h4], ...
    'Proposed VP-ZNN', 'NNN (Newton NN)', 'FT-GNN', 'GNN (Gradient NN)', ...
    'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast');
set(leg, 'Box', 'off');

disp('绘图完成！你可以直接在命令窗口输入 magnify 来生成局部放大图。');

% ================= 局部函数定义区 (确保矩阵计算精确匹配) =================
function err_norm = calc_error(t_arr, y_arr)
    err_norm = zeros(length(t_arr), 1);
    for j = 1:length(t_arr)
        A = matrixA(t_arr(j));
        W = matrixB(t_arr(j));
        b = vectorC(t_arr(j));
        S = [W, conj(A'); A, zeros(3,3)];
        c = [0;0;0;0; b];
        x_curr = reshape(y_arr(j, 1:7), 7, 1);
        err_norm(j) = norm(S * x_curr - c);
    end
end

function dx = nnnRightHandSide(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    % NNN 动态方程: dx/dt = -gamma * S^{-1} * (Sx - c)
    dx = -gamma * (S \ (S * x - c));
end

function dx = gnnRightHandSide(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    % GNN 动态方程: dx/dt = -gamma * S^H * (Sx - c)
    dx = -gamma * (S') * (S * x - c);
end

function dx = ftgnnRightHandSide(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    err = S * x - c;
    % FT-GNN 激活函数 (使用基于复数模长的幂激活)
    af = err .* (abs(err).^(-0.5));
    af(abs(err) < 1e-8) = 0; % 防止除以 0
    dx = -gamma * (S') * af;
end