function y = AFMnew(x,t)
%AFMSIGNBIPOWER 此处显示有关此函数的摘要
%   此处显示详细说明

    p=0.1;
    r=1;
    y = zeros(length(x), 1);
    i=find(x>0);
    y(i)=r/(p).*exp((abs(x(i))).^(p)).*abs(x(i)).^(1-p);
    j=find(x<=0);
    y(j)=-r/(p).*exp((abs(x(j))).^(p)).*abs(x(j)).^(1-p);
end

