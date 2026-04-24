function output = znnRightHandSide(t, y, gamma)

disp(t);

% Update x and Lambda
x = y(1:7);
lambda = y((end - 7 + 1): end);
Lambda = diag(lambda);
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
S1=[W1,conj(A1');A1,zeros(3,3)];
c=[0;0;0;0;b1];
dotS=[dW,conj(dA');dA,zeros(3,3)];
dotM=[dotS,zeros(7);zeros(7),zeros(7)];
dotc=[0;0;0;0;db];

err = S1*x-c;
absErr1 = abs(real(err));
absErr2=abs(imag(err));
err1=err;
% tanh(absErr2/ksi) is a continuous function to mock up sign(absErr2)
% errTol denotes the error tolerance
errTol = 1e-6;

for k = 1:length(err)
    if (absErr1(k) <= errTol && absErr2(k)<=errTol)
        err1(k) = 0;
    end
end
ksi = 0.001;
signErr2=2/pi*atan(abs(err1)/ksi);
vecD= [dotc-Lambda*(AFMnew(abs(err)).*exp(1i*(angle(err))));gamma*exp(t)*signErr2];
% Right handside of the ODE
output = -dotM*y + vecD;

end
