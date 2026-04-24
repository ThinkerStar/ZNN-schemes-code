% ==========================================================
% main_architecture_comparison.m
% 跨架构多维对比仿真 (VP-ZNN vs PDNN vs PI-RNN vs HNN)
% 彻底脱离了 GNN 和 NNN 的架构限制
% ==========================================================
clc; clearvars; close all; format long;

% ==================== 1. 全局参数与设置 ====================
tspan = linspace(0, 1, 300); 
gamma_znn = 4;
lambda0 = 3.0;

% 固定随机种子，确保严格公平对比
rng(2026); 
x0 = rand(7,1) + 1i*(rand(7,1));

% 初始化设置
y00_znn = [x0; lambda0 * ones(7, 1)]; % VP-ZNN 初始状态 (14维)
y00_pi  = [x0; zeros(7, 1)];          % PI-RNN 初始状态 (前7维x, 后7维积分u)
y00_pd  = [x0; zeros(7, 1)];          % PDNN 初始状态 (前7维x, 后7维对偶v)

% 求解器设置
options_znn  = odeset('RelTol', 1e-6,'AbsTol', 1e-6, 'Mass', @matrixM, 'MStateDep', 'none');
options_norm = odeset('RelTol', 1e-6,'AbsTol', 1e-6);

% ==================== 2. 求解各类网络架构 ====================
disp('1. 求解提出的 VP-ZNN (基于时变参数前馈预测架构)...');
[t_znn, y_znn] = ode45(@(t, y) znnRightHandSide(t, y, gamma_znn), tspan, y00_znn, options_znn);

disp('2. 求解 PDNN (原始-对偶鞍点优化架构)...');
alpha_pd = 5; beta_pd = 50; 
[~, y_pd] = ode45(@(t, y) pdnn_rhs(t, y, alpha_pd, beta_pd), tspan, y00_pd, options_norm);

disp('3. 求解 PI-RNN (比例积分记忆控制架构)...');
alpha_pi = 20; beta_pi = 50;
[~, y_pi] = ode45(@(t, y) pirnn_rhs(t, y, alpha_pi, beta_pi), tspan, y00_pi, options_norm);

disp('4. 求解 HNN (霍普菲尔德漏积分仿生架构)...');
eta_hnn = 2; gamma_hnn = 30;
[~, y_hnn] = ode45(@(t, y) hnn_rhs(t, y, eta_hnn, gamma_hnn), tspan, x0, options_norm);

% ==================== 3. 提取状态与计算误差 ====================
disp('正在计算多维评价指标...');
len = length(tspan);

% 提取真正的决策变量 x(t) (均为前7维)
x_znn = y_znn(:, 1:7);
x_pd  = y_pd(:, 1:7);
x_pi  = y_pi(:, 1:7);
x_hnn = y_hnn(:, 1:7);

x_exact = zeros(len, 7);
err_znn = zeros(len, 1); err_pd = zeros(len, 1); err_pi = zeros(len, 1); err_hnn = zeros(len, 1);
trk_znn = zeros(len, 1); trk_pd = zeros(len, 1); trk_pi = zeros(len, 1); trk_hnn = zeros(len, 1);

for j = 1:len
    t = tspan(j);
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)];
    c = [0;0;0;0; b];
    
    % 绝对理论解析解
    x_star = pinv(S) * c; 
    x_exact(j, :) = x_star.';
    
    % 计算残差 ||Sx - c||
    err_znn(j) = norm(S * x_znn(j, :).' - c);
    err_pd(j)  = norm(S * x_pd(j, :).' - c);
    err_pi(j)  = norm(S * x_pi(j, :).' - c);
    err_hnn(j) = norm(S * x_hnn(j, :).' - c);
    
    % 计算跟踪误差 ||x - x*||
    trk_znn(j) = norm(x_znn(j, :).' - x_star);
    trk_pd(j)  = norm(x_pd(j, :).' - x_star);
    trk_pi(j)  = norm(x_pi(j, :).' - x_star);
    trk_hnn(j) = norm(x_hnn(j, :).' - x_star);
end

% ==================== 4. 生成学术对比图表 ====================
disp('正在生成图表...');

c_znn = 'r'; ls_znn = '-';  lw_znn = 3;   % 你的算法 (红色)
c_pi  = 'b'; ls_pi  = '--'; lw_pi  = 2;   % PI-RNN (蓝色)
c_pd  = 'm'; ls_pd  = '-.'; lw_pd  = 2;   % PDNN (品红)
c_hnn = 'k'; ls_hnn = ':';  lw_hnn = 2.5; % HNN (黑色)

% --- 图 1: 常规线性残差 ---
figure('Position', [100, 550, 500, 400], 'Color', 'w');
plot(tspan, err_znn, c_znn, 'LineStyle', ls_znn, 'LineWidth', lw_znn); hold on;
plot(tspan, err_pi,  c_pi,  'LineStyle', ls_pi,  'LineWidth', lw_pi);
plot(tspan, err_pd,  c_pd,  'LineStyle', ls_pd,  'LineWidth', lw_pd);
plot(tspan, err_hnn, c_hnn, 'LineStyle', ls_hnn, 'LineWidth', lw_hnn);
format_plot('Time $t$ (s)', 'Residual Error $\|S(t)x(t) - c(t)\|_2$', 'Linear Residual Comparison');

% --- 图 2: 对数残差 (高精度区分) ---
figure('Position', [620, 550, 500, 400], 'Color', 'w');
semilogy(tspan, err_znn, c_znn, 'LineStyle', ls_znn, 'LineWidth', lw_znn); hold on;
semilogy(tspan, err_pi,  c_pi,  'LineStyle', ls_pi,  'LineWidth', lw_pi);
semilogy(tspan, err_pd,  c_pd,  'LineStyle', ls_pd,  'LineWidth', lw_pd);
semilogy(tspan, err_hnn, c_hnn, 'LineStyle', ls_hnn, 'LineWidth', lw_hnn);
format_plot('Time $t$ (s)', 'Log Residual $\|S(t)x(t) - c(t)\|_2$', 'Logarithmic Residual Comparison');

% --- 图 3: 状态轨迹跟踪 (Re(x1)) ---
figure('Position', [100, 70, 500, 400], 'Color', 'w');
plot(tspan, real(x_exact(:, 1)), 'Color', [0.7 0.7 0.7], 'LineWidth', 4); hold on;
plot(tspan, real(x_znn(:, 1)), c_znn, 'LineStyle', ls_znn, 'LineWidth', lw_znn);
plot(tspan, real(x_pi(:, 1)),  c_pi,  'LineStyle', ls_pi,  'LineWidth', lw_pi);
plot(tspan, real(x_pd(:, 1)),  c_pd,  'LineStyle', ls_pd,  'LineWidth', lw_pd);
format_plot('Time $t$ (s)', 'State Trajectory Re$(x_1(t))$', 'State Trajectory Tracking');
legend('Exact Solution $x^*(t)$', 'Proposed EFTCVPZNN', 'PI-RNN', 'PDNN', 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');

% --- 图 4: 状态空间误差 ---
figure('Position', [620, 70, 500, 400], 'Color', 'w');
plot(tspan, trk_znn, c_znn, 'LineStyle', ls_znn, 'LineWidth', lw_znn); hold on;
plot(tspan, trk_pi,  c_pi,  'LineStyle', ls_pi,  'LineWidth', lw_pi);
plot(tspan, trk_pd,  c_pd,  'LineStyle', ls_pd,  'LineWidth', lw_pd);
plot(tspan, trk_hnn, c_hnn, 'LineStyle', ls_hnn, 'LineWidth', lw_hnn);
format_plot('Time $t$ (s)', 'State Error $\|x(t) - x^*(t)\|_2$', 'State Tracking Error');

disp('所有绘图已完成！');

% ==================== 局部辅助与 RHS 动态方程区 ====================
function format_plot(x_lbl, y_lbl, fig_title)
    grid on; ax = gca; ax.GridLineStyle = '--'; ax.GridAlpha = 0.4;
    ax.FontSize = 13; ax.TickLabelInterpreter = 'latex';
    xlabel(x_lbl, 'Interpreter', 'latex', 'FontSize', 15); ylabel(y_lbl, 'Interpreter', 'latex', 'FontSize', 15);
    title(fig_title, 'Interpreter', 'latex', 'FontSize', 16);
    if ~contains(fig_title, 'Trajectory')
        leg = legend('Proposed EFTCVPZNN ', 'PI-RNN ', ...
            'PDNN ', 'HNN ', 'Interpreter', 'latex', 'FontSize', 11, 'Location', 'northeast');
        set(leg, 'Box', 'off');
    end
end

% 1. PI-RNN 右端项 (包含积分状态 u)
function dy = pirnn_rhs(t, y, alpha, beta)
    x = y(1:7); u = y(8:14); % u 为积分项
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    err = S * x - c;
    grad = S' * err;
    dx = -alpha * grad - beta * u;
    du = grad; % du/dt = grad => u 是 grad 的积分
    dy = [dx; du];
end

% 2. PDNN 右端项 (包含对偶状态 v)
function dy = pdnn_rhs(t, y, alpha, beta)
    x = y(1:7); v = y(8:14); % v 为对偶变量
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    err = S * x - c;
    dx = -alpha * x - S' * v; % 原始变量沿拉格朗日负梯度
    dv = beta * err;          % 对偶变量沿拉格朗日正梯度更新
    dy = [dx; dv];
end

% 3. HNN (Hopfield/Leaky NN) 右端项
function dx = hnn_rhs(t, x, eta, gamma)
    A = matrixA(t); W = matrixB(t); b = vectorC(t);
    S = [W, conj(A'); A, zeros(3,3)]; c = [0;0;0;0; b];
    err = S * x - c;
    dx = -eta * x - gamma * (S') * err; % 包含 -eta*x 的衰减(漏电)项
end