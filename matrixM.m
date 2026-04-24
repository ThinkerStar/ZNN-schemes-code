% Matrix M should be updated accordingly for other problems

function output = matrixM(t)

A = matrixA(t);
W=matrixB(t);
S=[W,conj(A'); A, zeros(3,3)];

output = [S,zeros(7, 7);zeros(7,7), eye(7,7)];

end