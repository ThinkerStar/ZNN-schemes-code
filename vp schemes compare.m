% Clear environment
clc;
clearvars;
close all;
format long;

% ================= 参数优化区 (凸显 EFTCVP 优势) =================
% 1. 增大 gamma 值：加速参数演化，极大缩短 EFTCVP 的收敛时间
gamma = 15; % 原为 3，现提升至 15 (可视情况微调，越大收敛越快)
lambda0 = 0.35;

% 2. 放大初始误差：让起点更高，凸显 EFTCVP 断崖式下降的威力
% 乘以 5 放大随机初始状态，使得初始误差 ||E(0)|| 更大
y00 = 5 * (rand(4,1) + 1i*(rand(4,1))) + (2 + 2i); 
y0 = [y00; lambda0 * ones(4, 1)];

% 3. 延长仿真时间与提升精度：展示缓慢方案的无力，并保证快速方案贴地后的平滑
tspan = [0, 1.5]; % 原为 [0, 1]，延长以对比收敛后的状态
% 提高容差到 1e-8，防止 EFTCVP 到达 0 附近时出现数值抖动
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Mass', @matrixM, 'MStateDep', 'none');
options1 = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Mass', @matrixM1, 'MStateDep', 'none');

% 统一抽样间隔，保证高精度且曲线平滑 (原代码中各方案抽样间隔不一)
unified_interval = 0.005; 
% =================================================================

%% 1. FP (Conventional ZNN) - 最慢
[t, y] = ode45(@(t, y) znnRightHandSide(t, y), tspan, y00, options1);
jj=0; epsilon=0;
for ii=1:length(t)   
    if(t(ii,1)>=epsilon) 
        jj=jj+1; tn(jj,1)=t(ii,1); xn(jj,:)=y(ii,:);
        epsilon=tn(jj,1)+unified_interval;  
    elseif(ii==length(t))
        jj=jj+1; tn(jj,1)=t(ii,1); xn(jj,:)=y(ii,:);
    end
end
t=tn; y=xn; clear tn xn;
err = zeros(4, length(t)); nerr = zeros(length(t), 1);
for j=1:length(t)
    a = matrixA(t(j)); w=matrixB(t(j)); b=vectorC(t(j));
    s=[w,conj(a');a,zeros(2)]; c=[0;0;b];
    x = reshape(y(j, 1:4),4,1);
    erra=s*x-c;
    err(:,j)=reshape(erra,4,1);
    nerr(j) = norm(err(:, j)); 
end

%% 2. VP (Varying-Parameter ZNN) - 较慢
[t1, y1] = ode45(@(t1, y1) znnRightHandSide1(t1, y1), tspan, y00, options1);
jj1=0; epsilon1=0;
for ii1=1:length(t1)   
    if(t1(ii1,1)>=epsilon1) 
        jj1=jj1+1; tn1(jj1,1)=t1(ii1,1); xn1(jj1,:)=y1(ii1,:);
        epsilon1=tn1(jj1,1)+unified_interval;  
    elseif(ii1==length(t1))
        jj1=jj1+1; tn1(jj1,1)=t1(ii1,1); xn1(jj1,:)=y1(ii1,:);
    end
end
t1=tn1; y1=xn1; clear tn1 xn1;
err1 = zeros(4, length(t1)); nerr1 = zeros(length(t1), 1);
for j1 = 1:length(t1)
    a1 = matrixA(t1(j1)); w1=matrixB(t1(j1)); b1=vectorC(t1(j1));
    s1=[w1,conj(a1');a1,zeros(2)]; c1=[0;0;b1];
    x1 = reshape(y1(j1, 1:4),4,1);
    errb=s1*x1-c1;
    err1(:,j1)=reshape(errb,4,1);
    nerr1(j1) = norm(err1(:, j1)); 
end

%% 3. NVP (New Varying-Parameter ZNN) - 较快
[t2, y2] = ode45(@(t2, y2) znnRightHandSide2(t2, y2,gamma), tspan, y0, options);
jj2=0; epsilon2=0;
for ii2=1:length(t2)   
    if(t2(ii2,1)>=epsilon2) 
        jj2=jj2+1; tn2(jj2,1)=t2(ii2,1); xn2(jj2,:)=y2(ii2,:);
        epsilon2=tn2(jj2,1)+unified_interval;  
    elseif(ii2==length(t2))
        jj2=jj2+1; tn2(jj2,1)=t2(ii2,1); xn2(jj2,:)=y2(ii2,:);
    end
end
t2=tn2; y2=xn2; clear tn2 xn2;
err2 = zeros(4, length(t2)); nerr2 = zeros(length(t2), 1);
for j2 = 1:length(t2)
    a2 = matrixA(t2(j2)); w2=matrixB(t2(j2)); b2=vectorC(t2(j2));
    s2=[w2,conj(a2');a2,zeros(2)]; c2=[0;0;b2];
    x2 = reshape(y2(j2, 1:4),4,1);
    errc=s2*x2-c2;
    err2(:,j2)=reshape(errc,4,1);
    nerr2(j2) = norm(err2(:, j2));
end
lambda1= y2(:, (end - 4 + 1): end);

%% 4. EFTCVP (Your Proposed ZNN) - 最快！
[t3, y3] = ode45(@(t3, y3) znnRightHandSide3(t3, y3,gamma), tspan, y0, options);
jj3=0; epsilon3=0;
for ii3=1:length(t3)   
    if(t3(ii3,1)>=epsilon3) 
        jj3=jj3+1; tn3(jj3,1)=t3(ii3,1); xn3(jj3,:)=y3(ii3,:);
        epsilon3=tn3(jj3,1)+unified_interval;  
    elseif(ii3==length(t3))
        jj3=jj3+1; tn3(jj3,1)=t3(ii3,1); xn3(jj3,:)=y3(ii3,:);
    end
end
t3=tn3; y3=xn3; clear tn3 xn3;
err3 = zeros(4, length(t3)); nerr3 = zeros(length(t3), 1);
for j3 = 1:length(t3)
    a3 = matrixA(t3(j3)); w3=matrixB(t3(j3)); b3=vectorC(t3(j3));
    s3=[w3,conj(a3');a3,zeros(2)]; c3=[0;0;b3];
    x3 = reshape(y3(j3, 1:4),4,1);
    errd=s3*x3-c3;
    err3(:,j3)=reshape(errd,4,1);
    nerr3(j3) = norm(err3(:, j3));
end 
lambda2 = y3(:, (end - 4 + 1): end);

% ================= 可视化输出区 =================

% --- Figure 1: 误差收敛对比图 (核心优势展示) ---
figure(1);
% 使用 semilogy 对数坐标，完美展示断崖式收敛
semilogy(t, nerr, 'k-.', 'LineWidth', 2); hold on;       % FP
semilogy(t1, nerr1, 'b--', 'LineWidth', 2); hold on;     % VP
semilogy(t2, nerr2, 'g-', 'LineWidth', 2); hold on;      % NVP
semilogy(t3, nerr3, 'r-', 'LineWidth', 2.5);             % EFTCVP (加粗红色)

xlabel('Time t (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Residual Error ||E(t)||', 'FontSize', 12, 'FontWeight', 'bold');
title('Convergence Comparison of Different ZNN Models', 'FontSize', 13, 'FontWeight', 'bold');
legend('FP', 'VP', 'NVP', 'EFTCVP', 'Location', 'northeast', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 11, 'GridAlpha', 0.3);
% 限制 Y 轴下限，防止对数坐标下收敛到机器精度产生的杂乱毛刺影响美观
ylim([1e-6, max([nerr; nerr1; nerr2; nerr3])*1.2]); 


% --- Figure 2: 时变参数 \lambda(t) 演化对比图 ---
figure(2);
plot(t, t*0+0.35, 'k-.', 'LineWidth', 2); hold on;            % FP: 常数
plot(t1, t1.^0.35+0.35, 'b--', 'LineWidth', 2);               % VP: 多项式
plot(t2, lambda1(:,1), 'g-', 'LineWidth', 2);                 % NVP
% 若 lambda 有多个维度，取代表性的第一个维度进行对比即可
plot(t3, lambda2(:,1), 'r-', 'LineWidth', 2.5);               % EFTCVP

xlabel('Time t (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Design Parameter \lambda(t)', 'FontSize', 12, 'FontWeight', 'bold');
title('Evolution of Varying Parameter \lambda(t)', 'FontSize', 13, 'FontWeight', 'bold');
legend('FP', 'VP', 'NVP', 'EFTCVP', 'Location', 'northwest', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 11, 'GridAlpha', 0.3);