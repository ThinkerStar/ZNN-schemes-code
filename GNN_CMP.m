% ==========================================================
% main_comprehensive_comparison.m
% 多维对比算法仿真程序 (VP-ZNN vs 4种非ZNN传统算法)
% 包含：线性残差、对数残差、状态轨迹跟踪、状态跟踪误差
% ==========================================================
clc; clearvars; close all; format long;

% ==================== 1. 全局参数设置 ====================
tspan = linspace(0, 1, 300); % 统一时间轴，保证各算法画图点完全对齐
gamma_znn = 4;
lambda0 = 3.0;

% 固定随机种子，保证公平对比 (统一初始起跑线)
rng(2026); 
x0 = rand(7,1) + 1i*(rand(7,1));
y00_znn = [x0; lambda0 * ones(7, 1)];

% 增益参数微调，保证各传统算法曲线具有代表性且平滑
gamma_nnn = 5; 
gamma_gnn = 15; 
gamma_ft  = 10;
gamma_ho  = 8;

options_znn = odeset('RelTol', 1e-6,'AbsTol', 1e-6, 'Mass', @matrixM, 'MStateDep', 'none');
options_norm = odeset('RelTol', 1e-6,'AbsTol', 1e-6);

% ==================== 2. 求解所有算法 ====================
disp('正在求解提出的 VP-ZNN...');
[t_znn, y_znn] = ode45(@(t, y) znnRightHandSide(t, y, gamma_znn), tspan, y00_znn, options_znn);

disp('正在求解非 ZNN 对比算法...');
[~, y_nnn] = ode45(@(t, y) nnn_rhs(t, y, gamma_nnn), tspan, x0, options_norm);
[~, y_gnn] = ode45(@(t, y) gnn_rhs(t, y, gamma_gnn), tspan, x0, options_norm);
[~, y_ft]  = ode45(@(t, y) ftgnn_rhs(t, y, gamma_ft), tspan, x0, options_norm);
[~, y_ho]  = ode45(@(t, y) hognn_rhs(t, y, gamma_ho), tspan, x0, options_norm);

% ==================== 3. 计算多维评价指标 ====================
disp('正在计算多维评价指标...');
len = length(tspan);
err_znn = zeros(len, 1); err_nnn = zeros(len, 1); err_gnn = zeros(len, 1);
err_ft  = zeros(len, 1); err_ho  = zeros(len, 1);

track_znn = zeros(len, 1); track_nnn = zeros(len, 1); track_gnn = zeros(len, 1);
track_ft  = zeros(len, 1); track_ho  = zeros(len, 1);

% 提取所有算法的状态 x(t) (前 7 个维度)
x_znn = y_znn(:, 1:7); x_nnn = y_nnn(:, 1:7); x_gnn = y_gnn(:, 1:7);
x_ft  = y_ft(:, 1:7);  x_ho  = y_ho(:, 1:7);

x_exact = zeros(len, 7); % 理论精确解

for j = 1:len
    t = tspan(j);
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)];
    c = [0;0;0;0; b];
    
    % 1. 计算理论绝对精确解 x*(t) = S^{-1}(t)c(t)
    x_star = pinv(S) * c; 
    x_exact(j, :) = x_star.';
    
    % 2. 计算残差 ||S*x - c||_2
    err_znn(j) = norm(S * x_znn(j, :).' - c);
    err_nnn(j) = norm(S * x_nnn(j, :).' - c);
    err_gnn(j) = norm(S * x_gnn(j, :).' - c);
    err_ft(j)  = norm(S * x_ft(j, :).' - c);
    err_ho(j)  = norm(S * x_ho(j, :).' - c);
    
    % 3. 计算状态跟踪误差 ||x(t) - x*(t)||_2
    track_znn(j) = norm(x_znn(j, :).' - x_star);
    track_nnn(j) = norm(x_nnn(j, :).' - x_star);
    track_gnn(j) = norm(x_gnn(j, :).' - x_star);
    track_ft(j)  = norm(x_ft(j, :).' - x_star);
    track_ho(j)  = norm(x_ho(j, :).' - x_star);
end

% ==================== 4. 绘制四维高水平对比图 ====================
disp('正在生成图表...');

% 颜色与线型配置
c_znn = 'r'; ls_znn = '-';  lw_znn = 2.5; % 你的算法
c_nnn = 'b'; ls_nnn = '--'; lw_nnn = 2;
c_gnn = [0 0.5 0]; ls_gnn = '-.'; lw_gnn = 2;
c_ft  = 'm'; ls_ft  = ':';  lw_ft  = 2.5;
c_ho  = 'k'; ls_ho  = '--'; lw_ho  = 2;

% --- 图 1: 常规残差误差 (Linear Residual) ---
figure('Position', [100, 550, 500, 400], 'Color', 'w');
plot(tspan, err_znn, c_znn, 'LineStyle', ls_znn, 'LineWidth', lw_znn); hold on;
plot(tspan, err_nnn, c_nnn, 'LineStyle', ls_nnn, 'LineWidth', lw_nnn);
plot(tspan, err_gnn, 'Color', c_gnn, 'LineStyle', ls_gnn, 'LineWidth', lw_gnn);
plot(tspan, err_ft,  c_ft,  'LineStyle', ls_ft,  'LineWidth', lw_ft);
plot(tspan, err_ho,  c_ho,  'LineStyle', ls_ho,  'LineWidth', lw_ho);
format_plot('Time $t$ (s)', 'Residual Error $\|S(t)x(t) - c(t)\|_2$', 'Linear Residual Comparison');

% --- 图 2: 对数残差误差 (Logarithmic Residual) - 凸显你的精度 ---
figure('Position', [620, 550, 500, 400], 'Color', 'w');
semilogy(tspan, err_znn, c_znn, 'LineStyle', ls_znn, 'LineWidth', lw_znn); hold on;
semilogy(tspan, err_nnn, c_nnn, 'LineStyle', ls_nnn, 'LineWidth', lw_nnn);
semilogy(tspan, err_gnn, 'Color', c_gnn, 'LineStyle', ls_gnn, 'LineWidth', lw_gnn);
semilogy(tspan, err_ft,  c_ft,  'LineStyle', ls_ft,  'LineWidth', lw_ft);
semilogy(tspan, err_ho,  c_ho,  'LineStyle', ls_ho,  'LineWidth', lw_ho);
format_plot('Time $t$ (s)', 'Log Residual $\|S(t)x(t) - c(t)\|_2$', 'Logarithmic Residual Comparison');

% --- 图 3: 状态轨迹跟踪 (以第一维度的实部 Re(x1) 为例) ---
figure('Position', [100, 70, 500, 400], 'Color', 'w');
% 画出理论参考线 (灰色粗底线)
plot(tspan, real(x_exact(:, 1)), 'Color', [0.7 0.7 0.7], 'LineWidth', 4); hold on;
plot(tspan, real(x_znn(:, 1)), c_znn, 'LineStyle', ls_znn, 'LineWidth', lw_znn);
plot(tspan, real(x_nnn(:, 1)), c_nnn, 'LineStyle', ls_nnn, 'LineWidth', lw_nnn);
plot(tspan, real(x_gnn(:, 1)), 'Color', c_gnn, 'LineStyle', ls_gnn, 'LineWidth', lw_gnn);
format_plot('Time $t$ (s)', 'State Trajectory Re$(x_1(t))$', 'State Trajectory Tracking');
% 特殊图例
legend('Exact Solution $x^*(t)$', 'Proposed EFTCVPZNN', 'NNN', 'GNN', 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');

% --- 图 4: 状态跟踪误差 (State Tracking Error ||x - x*||) ---
figure('Position', [620, 70, 500, 400], 'Color', 'w');
plot(tspan, track_znn, c_znn, 'LineStyle', ls_znn, 'LineWidth', lw_znn); hold on;
plot(tspan, track_nnn, c_nnn, 'LineStyle', ls_nnn, 'LineWidth', lw_nnn);
plot(tspan, track_gnn, 'Color', c_gnn, 'LineStyle', ls_gnn, 'LineWidth', lw_gnn);
plot(tspan, track_ft,  c_ft,  'LineStyle', ls_ft,  'LineWidth', lw_ft);
plot(tspan, track_ho,  c_ho,  'LineStyle', ls_ho,  'LineWidth', lw_ho);
format_plot('Time $t$ (s)', 'State Error $\|x(t) - x^*(t)\|_2$', 'State Tracking Error');

disp('所有绘图已完成！');

% ==================== 局部辅助与公式定义函数 ====================

% 图像统一样式格式化函数
function format_plot(x_lbl, y_lbl, fig_title)
    grid on; ax = gca; ax.GridLineStyle = '--'; ax.GridAlpha = 0.4;
    ax.FontSize = 13; ax.TickLabelInterpreter = 'latex';
    xlabel(x_lbl, 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(y_lbl, 'Interpreter', 'latex', 'FontSize', 15);
    title(fig_title, 'Interpreter', 'latex', 'FontSize', 16);
    % 除了图3外，其它调用默认图例
    if ~contains(fig_title, 'Trajectory')
        leg = legend('Proposed EFTCVPZNN', 'NNN (Newton NN)', 'GNN (Gradient NN)', ...
            'FT-GNN (Fractional Power)', 'HO-GNN (Higher-Order)', ...
            'Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast');
        set(leg, 'Box', 'off');
    end
end

% 对比算法 1: NNN 右端项
function dx = nnn_rhs(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    err = S * x - c;
    dx = -gamma * (S \ err);
end

% 对比算法 2: 常规 GNN 右端项
function dx = gnn_rhs(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    err = S * x - c;
    dx = -gamma * (S') * err;
end

% 对比算法 3: FT-GNN 右端项
function dx = ftgnn_rhs(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    err = S * x - c;
    abs_err = abs(err);
    af = zeros(size(err));
    for k=1:length(err)
        if abs_err(k) > 1e-6
            af(k) = err(k) * (abs_err(k)^(-0.5));
        end
    end
    dx = -gamma * (S') * af;
end

% 对比算法 4: HO-GNN 右端项
function dx = hognn_rhs(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    err = S * x - c;
    af = err + err .* (abs(err).^2);
    dx = -gamma * (S') * af;
end