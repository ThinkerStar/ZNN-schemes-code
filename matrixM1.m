% Matrix M should be updated accordingly for other problems

function output = matrixM1(t)

A = matrixA(t);
W=matrixB(t);
S=[W,conj(A'); A, zeros(2)];

output = S;

end