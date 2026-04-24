function dp = dpath1(t)

% dp = [-0.015; 0.12*cos(t);-0.2*sin(t)];%ÂÝÐýÏß
% 
% dp = [0; - (sin(t)*((5*cos(t/12)*sin(t/12)^4)/12 - 8*sin(4*t) + exp(cos(t))*sin(t)))/10 - (cos(t)*(2*cos(4*t) - exp(cos(t)) + sin(t/12)^5))/10;(sin(t)*(2*cos(4*t) - exp(cos(t)) + sin(t/12)^5))/10 - (cos(t)*((5*cos(t/12)*sin(t/12)^4)/12 - 8*sin(4*t) + exp(cos(t))*sin(t)))/10];%ÂÝÐýÏß
%ºûµû
% dp=[0;(2^sin((5*pi)/(3*t))*pi*log(2)*cos(pi/(3*t))*cos((5*pi)/(3*t)))/(4*t^2)...
%     - (pi*sin(pi/(3*t))*(2^sin((5*pi)/(3*t)) - 2))/(20*t^2);...
%     (pi*cos(pi/(3*t))*(2^sin((5*pi)/(3*t)) - 2))/(20*t^2) +...
%     (2^sin((5*pi)/(3*t))*pi*log(2)*cos((5*pi)/(3*t))*sin(pi/(3*t)))/(4*t^2)];
% »¨°ê
%nengpao
% dp=[0;(pi*sin((pi*t)/3)*(2^sin((5*pi*t)/3) - 2))/20 - (2^sin((5*pi*t)/3)*pi*log(2)*cos((pi*t)/3)*cos((5*pi*t)/3))/4;- (pi*cos((pi*t)/3)*(2^sin((5*pi*t)/3) - 2))/20 - (2^sin((5*pi*t)/3)*pi*log(2)*cos((5*pi*t)/3)*sin((pi*t)/3))/4];


dp=[0;-(pi*sin((pi*t)/3)*(2^sin((5*pi*t)/3) - 2))/20 - (2^sin((5*pi*t)/3)*pi*log(2)*cos((pi*t)/3)*cos((5*pi*t)/3))/4;(pi*cos((pi*t)/3)*(2^sin((5*pi*t)/3) - 2))/20 - (2^sin((5*pi*t)/3)*pi*log(2)*cos((5*pi*t)/3)*sin((pi*t)/3))/4];
