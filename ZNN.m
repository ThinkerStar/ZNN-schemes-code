
% Clear environment
clc;
clearvars;
close all;
format long;

% Values of gamma and lambda(0) 
gamma =3; % gamma = 0: Conventional ZNN, gamma > 0: Varying-parameter ZNN
lambda0 = 0.35;

% ODE
% setting
tspan =[0,1];                                                                                                                           
options = odeset('RelTol', 1e-6,'AbsTol', 1e-6, 'Mass', @matrixM, 'MStateDep', 'none');
options1 = odeset('RelTol', 1e-6,'AbsTol', 1e-6, 'Mass', @matrixM1, 'MStateDep', 'none');
% Initialization
y00=rand(4,1)+1i*(rand(4,1));
y0 = [y00; lambda0 * ones(4, 1)];

% Solving using ode45 solver
[t, y] = ode45(@(t, y) znnRightHandSide(t, y), tspan, y00, options1);
jj=0;
epsilon=0;
%interval=0.00000005;
interval=0.0025;  %primal

   sndsizey=size(y);
for ii=1:length(t)   
    if(t(ii,1)>=epsilon) 
        jj=jj+1;
        tn(jj,1)=t(ii,1);
        xn(jj,:)=y(ii,:);
        epsilon=tn(jj,1)+interval; 
    elseif(ii==length(t))
        jj=jj+1;
        tn(jj,1)=t(ii,1);
        xn(jj,:)=y(ii,:);
    end
end
clear t y;
t=tn;
y=xn;

err = zeros(4, length(t));
nerr = zeros(length(t), 1);
 for j=1:length(t)
     a = matrixA(t(j));
     w=matrixB(t(j));
     b=vectorC(t(j));
     s=[w,conj(a');a,zeros(2)];
     c=[0;0;b];
      x = y(j, 1:4);
      x=reshape(x,4,1);
     erra=s*x-c;
    err(:,j)=reshape(erra,4,1);
    nerr(j) = norm(err(:, j)); 
 end

[t1, y1] = ode45(@(t1, y1) znnRightHandSide1(t1, y1), tspan, y00, options1);

 jj1=0;
epsilon1=0;
%interval=0.00000005;
interval1=0.0125;  %primal

   sndsizey1=size(y1);
for ii1=1:length(t1)   
    if(t1(ii1,1)>=epsilon1)
        jj1=jj1+1;
        tn1(jj1,1)=t1(ii1,1);
        xn1(jj1,:)=y1(ii1,:);
        epsilon1=tn1(jj1,1)+interval1;  %抽样间隔
    elseif(ii1==length(t1))
        jj1=jj1+1;
        tn1(jj1,1)=t1(ii1,1);
        xn1(jj1,:)=y1(ii1,:);
    end
end
clear t1 y1;
t1=tn1;
y1=xn1;
%Data visualization
err1 = zeros(4, length(t1));
nerr1 = zeros(length(t1), 1);
for j1 = 1:length(t1)
    a1 = matrixA(t1(j1));
     w1=matrixB(t1(j1));
     b1=vectorC(t1(j1));
     s1=[w1,conj(a1');a1,zeros(2)];
     c1=[0;0;b1];
      x1 = y1(j1, 1:4);
      x1=reshape(x1,4,1);
     errb=s1*x1-c1;
    err1(:,j1)=reshape(errb,4,1);
    nerr1(j1) = norm(err1(:, j1)); 
end
% 
[t2, y2] = ode45(@(t2, y2) znnRightHandSide2(t2, y2,gamma), tspan,y0, options);

 jj2=0;
epsilon2=0;
%interval=0.00000005;
interval2=0.0125;  %primal

   sndsizey2=size(y2);
for ii2=1:length(t2)   
    if(t2(ii2,1)>=epsilon2)
        jj2=jj2+1;
        tn2(jj2,1)=t2(ii2,1);
        xn2(jj2,:)=y2(ii2,:);
        epsilon2=tn2(jj2,1)+interval2;  
    elseif(ii2==length(t2))
        jj2=jj2+1;
        tn2(jj2,1)=t2(ii2,1);
        xn2(jj2,:)=y2(ii2,:);
    end
end
clear t2 y2;
t2=tn2;
y2=xn2;
%Data visualization
err2 = zeros(4, length(t2));
nerr2 = zeros(length(t2), 1);
for j2 = 1:length(t2)
     a2 = matrixA(t2(j2));
     w2=matrixB(t2(j2));
     b2=vectorC(t2(j2));
     s2=[w2,conj(a2');a2,zeros(2)];
     c2=[0;0;b2];
      x2 = y2(j2, 1:4);
      x2=reshape(x2,4,1);
     errc=s2*x2-c2;
    err2(:,j2)=reshape(errc,4,1);
    nerr2(j2) = norm(err2(:, j2));
end
lambda1= y2(:, (end - 4 + 1): end);
% 
[t3, y3] = ode45(@(t3, y3) znnRightHandSide3(t3, y3,gamma), tspan, y0, options);

 jj3=0;
epsilon3=0;
%interval=0.00000005;
interval3=0.0125;  %primal

   sndsizey3=size(y3);
for ii3=1:length(t3)   
    if(t3(ii3,1)>=epsilon3) 
        jj3=jj3+1;
        tn3(jj3,1)=t3(ii3,1);
        xn3(jj3,:)=y3(ii3,:);
        epsilon3=tn3(jj3,1)+interval3;  %抽样间隔
    elseif(ii3==length(t3))
        jj3=jj3+1;
        tn3(jj3,1)=t3(ii3,1);
        xn3(jj3,:)=y3(ii3,:);
    end
end
clear t3 y3;
t3=tn3;
y3=xn3;
%Data visualization
err3 = zeros(4, length(t3));
nerr3 = zeros(length(t3), 1);
for j3 = 1:length(t3)
     a3 = matrixA(t3(j3));
     w3=matrixB(t3(j3));
     b3=vectorC(t3(j3));
     s3=[w3,conj(a3');a3,zeros(2)];
     c3=[0;0;b3];
      x3 = y3(j3, 1:4);
      x3=reshape(x3,4,1);
     errd=s3*x3-c3;
    err3(:,j3)=reshape(errd,4,1);
    nerr3(j3) = norm(err3(:, j3));
end 

lambda2 = y3(:, (end - 4 + 1): end);



figure (1);
plot(t, t*0+0.35,'c','LineWidth',2);hold on;
plot(t1,t1.^0.35+0.35,'b','LineWidth',2);hold on;
plot(t2, lambda1(:,1),'y','LineWidth',2);hold on;
plot(t2,lambda1(:,2),'r','LineWidth',2);hold on;
plot(t3, lambda2(:,1),'g','LineWidth',2);hold on;
plot(t3,lambda2(:,2),'m','LineWidth',2);