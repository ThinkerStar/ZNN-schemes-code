% ==========================================================
% main_all_comparison.m
% 多元化对比算法仿真程序 (VP-ZNN vs 4种非ZNN传统算法)
% ==========================================================
clc; clearvars; close all; format long;

% 1. 全局参数与时间轴设置
% 使用 linspace 替代原始复杂的抽样逻辑，保证曲线绝对平滑
tspan = linspace(0, 1, 300); 
gamma_znn = 4;
lambda0 = 3.0;

% 固定随机种子，确保所有算法从【完全相同】的初始误差开始！
% 这是让对比曲线从同一个点平滑散开的关键
rng(2026); 
x0 = rand(7,1) + 1i*(rand(7,1));
y00_znn = [x0; lambda0 * ones(7, 1)];

% ==================== 2. 求解各种算法 ====================
options_znn = odeset('RelTol', 1e-6,'AbsTol', 1e-6, 'Mass', @matrixM, 'MStateDep', 'none');
options_norm = odeset('RelTol', 1e-6,'AbsTol', 1e-6);

% (1) 求解提出的 VP-ZNN (调用你自己的 znnRightHandSide)
disp('正在求解提出的 VP-ZNN...');
[t_znn, y_znn] = ode45(@(t, y) znnRightHandSide(t, y, gamma_znn), tspan, y00_znn, options_znn);

% (2) 求解对比算法 1: NNN (牛顿神经网络)
disp('正在求解 NNN...');
gamma_nnn = 5; % 较小的增益，产生优美的指数下降轨迹
[t_nnn, y_nnn] = ode45(@(t, y) nnn_rhs(t, y, gamma_nnn), tspan, x0, options_norm);

% (3) 求解对比算法 2: GNN (常规梯度神经网络)
disp('正在求解 GNN...');
gamma_gnn = 15; 
[t_gnn, y_gnn] = ode45(@(t, y) gnn_rhs(t, y, gamma_gnn), tspan, x0, options_norm);

% (4) 求解对比算法 3: FT-GNN (分数次幂激活梯度神经网络)
disp('正在求解 FT-GNN...');
gamma_ft = 10;
[t_ft, y_ft] = ode45(@(t, y) ftgnn_rhs(t, y, gamma_ft), tspan, x0, options_norm);

% (5) 求解对比算法 4: HO-GNN (高阶激活梯度神经网络)
disp('正在求解 HO-GNN...');
gamma_ho = 8;
[t_ho, y_ho] = ode45(@(t, y) hognn_rhs(t, y, gamma_ho), tspan, x0, options_norm);

% ==================== 3. 计算所有算法的残差 ====================
disp('正在计算误差轨迹...');
err_znn = calc_error(t_znn, y_znn);
err_nnn = calc_error(t_nnn, y_nnn);
err_gnn = calc_error(t_gnn, y_gnn);
err_ft  = calc_error(t_ft, y_ft);
err_ho  = calc_error(t_ho, y_ho);

% ==================== 4. 绘制高水平对比曲线图 ====================
figure('Position', [150, 150, 750, 550], 'Color', 'w');

% 绘制曲线 (设置高区分度的线型和颜色)
h1 = plot(t_znn, err_znn, 'r-', 'LineWidth', 3); hold on;          % VP-ZNN: 红色粗实线 (主角)
h2 = plot(t_nnn, err_nnn, 'b--', 'LineWidth', 2);                  % NNN: 蓝色虚线
h3 = plot(t_gnn, err_gnn, 'Color', [0 0.5 0], 'LineStyle', '-.', 'LineWidth', 2); % GNN: 深绿色点划线
h4 = plot(t_ft,  err_ft,  'm:', 'LineWidth', 2.5);                 % FT-GNN: 品红色点线
h5 = plot(t_ho,  err_ho,  'k--', 'LineWidth', 2);                  % HO-GNN: 黑色长虚线

% 图像细节美化
grid on;
ax = gca;
ax.GridLineStyle = '--';
ax.GridAlpha = 0.4;
ax.FontSize = 13;
ax.TickLabelInterpreter = 'latex';

% 坐标轴标签 (LaTeX 渲染)
xlabel('Time $t$ (s)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Residual Error $\|S(t)x(t) - c(t)\|_2$', 'Interpreter', 'latex', 'FontSize', 16);

% 图例设置
leg = legend([h1, h2, h3, h4, h5], ...
    'Proposed VP-ZNN', 'NNN (Newton NN)', 'GNN (Gradient NN)', ...
    'FT-GNN (Fractional Power)', 'HO-GNN (Higher-Order)', ...
    'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast');
set(leg, 'Box', 'off'); % 去除图例边框，显得更高级

disp('仿真与绘图完成！(若需要局部放大，请在命令窗口输入 magnify)');

% =========================================================================
% 内部局部函数定义 (确保调用的是你现有的 matrixA, matrixB, vectorC 环境)
% =========================================================================

% 误差计算公共函数
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

% 对比算法 1: NNN 右端项
function dx = nnn_rhs(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    err = S * x - c;
    % NNN: dx = -gamma * S^{-1} * err
    dx = -gamma * (S \ err);
end

% 对比算法 2: 常规 GNN 右端项
function dx = gnn_rhs(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    err = S * x - c;
    % GNN: dx = -gamma * S^H * err
    dx = -gamma * (S') * err;
end

% 对比算法 3: FT-GNN 右端项 (使用复数分数幂激活)
function dx = ftgnn_rhs(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    err = S * x - c;
    % 激活函数 Phi(e) = e / |e| * |e|^0.5 = e * |e|^{-0.5}
    abs_err = abs(err);
    af = zeros(size(err));
    for k=1:length(err)
        if abs_err(k) > 1e-6
            af(k) = err(k) * (abs_err(k)^(-0.5));
        end
    end
    dx = -gamma * (S') * af;
end

% 对比算法 4: HO-GNN 右端项 (使用高阶复数激活)
function dx = hognn_rhs(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    err = S * x - c;
    % 激活函数 Phi(e) = e + e * |e|^2
    af = err + err .* (abs(err).^2);
    dx = -gamma * (S') * af;
end