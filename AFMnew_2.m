function y = AFMnew_2(x,t)
%AFMSIGNBIPOWER 此处显示有关此函数的摘要
%   此处显示详细说明
% 
%     p=0.3;
%     r=2;
%     y = zeros(length(x), 1);
%     i=find(x>0);
%     y(i)=r/(p).*exp((abs(x(i))).^(p)).*abs(x(i)).^(1-p);
%     j=find(x<=0);
%     y(j)=-r/(p).*exp((abs(x(j))).^(p)).*abs(x(j)).^(1-p);

y= arrayfun(@(u) exp(abs(u).^0.5).*((log(abs(u) + 1)).^0.5).*sign(u)./0.5 + u+sign(u),x);
% p=0.5;
% y = zeros(length(x), 1);
% i=find(x>0);
% y(i)=1/3*abs(x(i)).^p+1/3*abs(x(i)).^(1/p) + 1/3*x(i).^3+1;
% j=find(x<=0);
% y(j)=-1/3*abs(x(j)).^p-1/3*abs(x(j)).^(1/p) + 1/3*x(j).^3-1;
end

