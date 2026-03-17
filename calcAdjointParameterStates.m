function [db,ds,deta] = calcAdjointParameterStates(db_a, ds_a, deta_a, timeMesh, omega)

db_t = calcSolution(timeMesh, db_a, omega);
db = sum(db_t .* timeMesh(2) - timeMesh(1),1).'; % TODO: use trapezoid rule

ds_t = calcSolution(timeMesh, ds_a, omega);
ds = sum(ds_t .* timeMesh(2) - timeMesh(1),1).'; % TODO: use trapezoid rule

deta_t = calcSolution(timeMesh, deta_a, omega);
deta = sum(deta_t .* timeMesh(2) - timeMesh(1),1).'; % TODO: use trapezoid rule

end

