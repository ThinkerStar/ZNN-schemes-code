function output = znnRightHandSide1(t, y)

disp(t);

% Update x and Lambda

A1 = matrixA(t);
W1 = matrixB(t);
b1=vectorC(t);
syms u ;
dotA = diff(matrixA(u));
u =t;
dA = eval(dotA);
syms v ;
dotW = diff(matrixB(v));
v=t;
dW = eval(dotW);
syms w ;
dotb = diff(vectorC(w));
w=t;
db = eval(dotb);
S1=[W1,conj(A1');A1,zeros(2)];
c=[0;0;b1];
dotS=[dW,conj(dA');dA,zeros(2)];

dotc=[0;0;db];

err = S1*y-c;

vecD= dotc-(t^0.35+0.35)*(AFMnew_2(abs(err)).*exp(1i*angle(err)));
% Right handside of the ODE
output = -dotS*y + vecD;

end
