function output = matrixB(t)

output =[cos(t)+1i*sin(t),cos(3*t),         2-1i*sin(2*t),1-1i*cos(3*t);
         cos(3*t),         1i*sin(t),        1+1i*cos(t),   sin(2*t)-1i*cos(3*t);
          2-1i*sin(2*t),    1+1i*cos(t),      1-1i*cos(t),  sin(t);
        1-1i*cos(3*t),     sin(2*t)-1i*cos(t),  sin(t),     2+1i*sin(t)];

end
