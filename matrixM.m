% Matrix M should be updated accordingly for other problems

function output = matrixM(t)

A = matrixA(t);
W=matrixB(t);
S=[W,conj(A'); A, zeros(2)];

output = [S,zeros(4, 4);zeros(4, 4), eye(4,4)];

end