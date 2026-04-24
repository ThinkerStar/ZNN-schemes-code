function y = AFMnew_3(L_ij)
%AFMSIGNBIPOWER 此处显示有关此函数的摘要
%   此处显示详细说明

%     p=0.5;
%     y = zeros(length(x), 1);
%     i=find(x>0);
%     y(i)=1/3*abs(x(i)).^p+1/3*abs(x(i)).^(1/p) + 1/3*x(i).^3+1/3;
%     j=find(x<=0);
%     y(j)=-1/3*abs(x(j)).^p-1/3*abs(x(j)).^(1/p) + 1/3*x(j).^3-1/3;
    y= arrayfun(@(x) exp(abs(x).^0.5).*((log(abs(x) + 1)).^0.5).*sign(x)./0.5 + x+sign(x),L_ij);
end

