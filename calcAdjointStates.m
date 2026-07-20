function [b_adj, s_adj, eta_adj] = calcAdjointStates(Uadj, omega, timeMesh, u0Laplace, u0tt, u0Sqtt)

adjointState = calcSolution(timeMesh, Uadj, omega);
avgT = 1/(timeMesh(1,end) - timeMesh(1,1));
s_adj = adjointState.*u0Laplace;
s_adj = avgT.*trapz(timeMesh.', s_adj,1).';

b_adj = adjointState.*u0tt;
b_adj = -avgT.*trapz(timeMesh.', b_adj,1).';

eta_adj = adjointState.*u0Sqtt;
eta_adj = avgT.*trapz(timeMesh.', eta_adj,1).';

end

