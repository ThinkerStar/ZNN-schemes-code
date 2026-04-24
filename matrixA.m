% Matrix A should be updated accordingly for other problems

function output = matrixA(t)
output=[1+1i*cos(3*t),2+1i*sin(2*t),sin(2*t)-1i*cos(t),2-1i*cos(3*t);
    1+1i*sin(t),cos(3*t)-1i*sin(t),3*cos(2*t),sin(2*t);
    sin(t)+1i*cos(t),6i*cos(t),cos(3*t)-1i*sin(t),cos(t)+1i];
end