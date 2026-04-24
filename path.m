function p = path(t)

a = 0.3;
b = 0.3;
px = -0.3589;
py = -0.11;
pz = 0.7017;
cof = 2*pi/3;

p = [px-a*cos(cof*t);py+b*sin(cof*t);pz];

