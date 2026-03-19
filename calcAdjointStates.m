function [b_adj, s_adj, eta_adj] = calcAdjointStates(Uadj, omega, timeMesh, u0Laplace, u0tt, u0Sqtt)

adjointState = calcSolution(timeMesh, Uadj, omega);

s_adj = adjointState.*u0Laplace;
s_adj = sum(s_adj .* (timeMesh(2) - timeMesh(1)),1).';

b_adj = adjointState.*u0tt;
b_adj = (-1).*sum(b_adj .* (timeMesh(2) - timeMesh(1)),1).';

eta_adj = adjointState.*u0Sqtt;
eta_adj = sum(eta_adj .* (timeMesh(2) - timeMesh(1)),1).';


end

