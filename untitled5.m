% ==========================================================
% main_advanced_comparison.m
% 顶刊级别多维对比仿真 (Proposed VP-ZNN vs 先进 ZNN 变体 vs GNN)
% 包含算法: VP-ZNN(你的), FP-ZNN, SBP-ZNN, VZN, GNN
% ==========================================================
clc; clearvars; close all; format long;

% ==================== 1. 全局参数设置 ====================
tspan = linspace(0, 1, 300); 
gamma_znn = 4;
lambda0 = 3.0;

% 固定随机种子保证统一起跑线
rng(2026); 
x0 = rand(7,1) + 1i*(rand(7,1));
y00_znn = [x0; lambda0 * ones(7, 1)];

% 传统算法的常数增益
gamma_base = 10; 

% ODE 设置
% 你的 VP-ZNN 依然调用你自己编写的包含 mass matrix 的右端项
options_znn = odeset('RelTol', 1e-6,'AbsTol', 1e-6, 'Mass', @matrixM, 'MStateDep', 'none');
% 新的对比 ZNN 算法被我优化为了标准显式 ODE，速度极快
options_norm = odeset('RelTol', 1e-6,'AbsTol', 1e-6);

% ==================== 2. 求解所有算法 ====================
disp('正在求解提出的 VP-ZNN (Proposed)...');
[t_znn, y_znn] = ode45(@(t, y) znnRightHandSide(t, y, gamma_znn), tspan, y00_znn, options_znn);

disp('正在求解 FP-ZNN (传统固定参数 ZNN)...');
[~, y_fp] = ode45(@(t, y) fp_znn_rhs(t, y, gamma_base), tspan, x0, options_norm);

disp('正在求解 SBP-ZNN (先进符号双幂有限时间 ZNN)...');
[~, y_sbp] = ode45(@(t, y) sbp_znn_rhs(t, y, gamma_base), tspan, x0, options_norm);

disp('正在求解 VZN (指数变参数 ZNN)...');
[~, y_vzn] = ode45(@(t, y) vzn_rhs(t, y, gamma_base), tspan, x0, options_norm);

disp('正在求解 GNN (传统梯度神经网络)...');
[~, y_gnn] = ode45(@(t, y) gnn_rhs(t, y, 20), tspan, x0, options_norm);

% ==================== 3. 计算多维评价指标 ====================
disp('正在计算评价指标...');
len = length(tspan);
err_znn = zeros(len, 1); err_fp = zeros(len, 1); err_sbp = zeros(len, 1);
err_vzn = zeros(len, 1); err_gnn = zeros(len, 1);

track_znn = zeros(len, 1); track_fp = zeros(len, 1); track_sbp = zeros(len, 1);
track_vzn = zeros(len, 1); track_gnn = zeros(len, 1);

x_znn = y_znn(:, 1:7); x_exact = zeros(len, 7);

for j = 1:len
    t = tspan(j);
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)];
    c = [0;0;0;0; b];
    
    % 计算理论绝对精确解
    x_star = pinv(S) * c; 
    x_exact(j, :) = x_star.';
    
    % 残差
    err_znn(j) = norm(S * x_znn(j, :).' - c);
    err_fp(j)  = norm(S * y_fp(j, :).' - c);
    err_sbp(j) = norm(S * y_sbp(j, :).' - c);
    err_vzn(j) = norm(S * y_vzn(j, :).' - c);
    err_gnn(j) = norm(S * y_gnn(j, :).' - c);
    
    % 跟踪误差
    track_znn(j) = norm(x_znn(j, :).' - x_star);
    track_fp(j)  = norm(y_fp(j, :).' - x_star);
    track_sbp(j) = norm(y_sbp(j, :).' - x_star);
    track_vzn(j) = norm(y_vzn(j, :).' - x_star);
    track_gnn(j) = norm(y_gnn(j, :).' - x_star);
end

% ==================== 4. 绘制四维高水平对比图 ====================
disp('正在生成图表...');

c_znn = 'r'; ls_znn = '-';  lw_znn = 3;   % 你的 VP-ZNN: 红色粗实线
c_fp  = 'b'; ls_fp  = '--'; lw_fp  = 2;   % FP-ZNN: 蓝色虚线
c_sbp = 'm'; ls_sbp = '-.'; lw_sbp = 2.5; % SBP-ZNN: 品红点划线
c_vzn = 'k'; ls_vzn = ':';  lw_vzn = 2.5; % VZN: 黑色点线
c_gnn = [0 0.5 0]; ls_gnn = '--'; lw_gnn = 1.5; % GNN: 深绿细虚线

% 图 1: 常规残差
figure('Position', [100, 550, 500, 400], 'Color', 'w');
plot(tspan, err_znn, c_znn, 'LineStyle', ls_znn, 'LineWidth', lw_znn); hold on;
plot(tspan, err_sbp, c_sbp, 'LineStyle', ls_sbp, 'LineWidth', lw_sbp);
plot(tspan, err_vzn, c_vzn, 'LineStyle', ls_vzn, 'LineWidth', lw_vzn);
plot(tspan, err_fp,  c_fp,  'LineStyle', ls_fp,  'LineWidth', lw_fp);
plot(tspan, err_gnn, 'Color', c_gnn, 'LineStyle', ls_gnn, 'LineWidth', lw_gnn);
format_plot('Time $t$ (s)', 'Residual Error $\|S(t)x(t) - c(t)\|_2$', 'Linear Residual Comparison');

% 图 2: 对数残差 (决定性图表，看谁先撞底)
figure('Position', [620, 550, 500, 400], 'Color', 'w');
semilogy(tspan, err_znn, c_znn, 'LineStyle', ls_znn, 'LineWidth', lw_znn); hold on;
semilogy(tspan, err_sbp, c_sbp, 'LineStyle', ls_sbp, 'LineWidth', lw_sbp);
semilogy(tspan, err_vzn, c_vzn, 'LineStyle', ls_vzn, 'LineWidth', lw_vzn);
semilogy(tspan, err_fp,  c_fp,  'LineStyle', ls_fp,  'LineWidth', lw_fp);
semilogy(tspan, err_gnn, 'Color', c_gnn, 'LineStyle', ls_gnn, 'LineWidth', lw_gnn);
format_plot('Time $t$ (s)', 'Log Residual $\|S(t)x(t) - c(t)\|_2$', 'Logarithmic Residual Comparison');

% 图 3: 状态轨迹跟踪 (Re(x1))
figure('Position', [100, 70, 500, 400], 'Color', 'w');
plot(tspan, real(x_exact(:, 1)), 'Color', [0.7 0.7 0.7], 'LineWidth', 4); hold on;
plot(tspan, real(x_znn(:, 1)), c_znn, 'LineStyle', ls_znn, 'LineWidth', lw_znn);
plot(tspan, real(y_sbp(:, 1)), c_sbp, 'LineStyle', ls_sbp, 'LineWidth', lw_sbp);
plot(tspan, real(y_fp(:, 1)),  c_fp,  'LineStyle', ls_fp,  'LineWidth', lw_fp);
format_plot('Time $t$ (s)', 'State Trajectory Re$(x_1(t))$', 'State Trajectory Tracking');
legend('Exact Solution $x^*(t)$', 'Proposed VP-ZNN', 'SBP-ZNN', 'FP-ZNN', 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');

% 图 4: 状态跟踪误差
figure('Position', [620, 70, 500, 400], 'Color', 'w');
plot(tspan, track_znn, c_znn, 'LineStyle', ls_znn, 'LineWidth', lw_znn); hold on;
plot(tspan, track_sbp, c_sbp, 'LineStyle', ls_sbp, 'LineWidth', lw_sbp);
plot(tspan, track_vzn, c_vzn, 'LineStyle', ls_vzn, 'LineWidth', lw_vzn);
plot(tspan, track_fp,  c_fp,  'LineStyle', ls_fp,  'LineWidth', lw_fp);
plot(tspan, track_gnn, 'Color', c_gnn, 'LineStyle', ls_gnn, 'LineWidth', lw_gnn);
format_plot('Time $t$ (s)', 'State Error $\|x(t) - x^*(t)\|_2$', 'State Tracking Error');

disp('所有绘图已完成！');

% ==================== 局部辅助与 RHS 核心计算区 ====================
function format_plot(x_lbl, y_lbl, fig_title)
    grid on; ax = gca; ax.GridLineStyle = '--'; ax.GridAlpha = 0.4;
    ax.FontSize = 13; ax.TickLabelInterpreter = 'latex';
    xlabel(x_lbl, 'Interpreter', 'latex', 'FontSize', 15); ylabel(y_lbl, 'Interpreter', 'latex', 'FontSize', 15);
    title(fig_title, 'Interpreter', 'latex', 'FontSize', 16);
    if ~contains(fig_title, 'Trajectory')
        leg = legend('Proposed VP-ZNN', 'SBP-ZNN (Finite-Time)', 'VZN (Exponential Varying)', ...
            'FP-ZNN (Fixed Parameter)', 'GNN (Gradient NN)', 'Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast');
        set(leg, 'Box', 'off');
    end
end

% ------ 高效数值求导模块 (替代耗时的 syms) ------
function [dS, dc] = get_derivatives(t)
    dt = 1e-5;
    A_p = matrixA(t+dt); W_p = matrixB(t+dt); b_p = vectorC(t+dt);
    S_p = [W_p, conj(A_p'); A_p, zeros(3,3)]; c_p = [0;0;0;0; b_p];
    A_m = matrixA(t-dt); W_m = matrixB(t-dt); b_m = vectorC(t-dt);
    S_m = [W_m, conj(A_m'); A_m, zeros(3,3)]; c_m = [0;0;0;0; b_m];
    dS = (S_p - S_m) / (2*dt);
    dc = (c_p - c_m) / (2*dt);
end

% 对比算法 1: FP-ZNN (固定参数线性 ZNN)
function dx = fp_znn_rhs(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    [dS, dc] = get_derivatives(t);
    err = S * x - c;
    % 核心动态方程: S*dx = -dS*x + dc - gamma * err
    dx = S \ (-dS * x + dc - gamma * err);
end

% 对比算法 2: SBP-ZNN (有限时间符号双幂 ZNN)
function dx = sbp_znn_rhs(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    [dS, dc] = get_derivatives(t);
    err = S * x - c;
    
    % SBP 激活函数 (复数域推广: e * |e|^(p-1) + e * |e|^(q-1))
    abs_e = abs(err);
    af = zeros(size(err));
    for k = 1:length(err)
        if abs_e(k) > 1e-8
            af(k) = err(k)*(abs_e(k)^(-0.5)) + err(k)*(abs_e(k)^(0.5));
        end
    end
    dx = S \ (-dS * x + dc - gamma * af);
end

% 对比算法 3: VZN (纯指数变参数 ZNN)
function dx = vzn_rhs(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    [dS, dc] = get_derivatives(t);
    err = S * x - c;
    % 时变增益参数
    gamma_t = gamma + exp(t);
    dx = S \ (-dS * x + dc - gamma_t * err);
end

% 经典 GNN 算法 (保留作反面对比)
function dx = gnn_rhs(t, x, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    err = S * x - c;
    dx = -gamma * (S') * err;
end