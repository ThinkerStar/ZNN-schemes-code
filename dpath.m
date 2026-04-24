function dp = dpath(t)

global a;
global b;
cof = 2*pi/3;

dp = [a*cof*sin(cof*t); b*cof*cos(cof*t);   0];