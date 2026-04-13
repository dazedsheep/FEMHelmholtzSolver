function [x, k, stopres] = landweber(elements, A, b, x0, tol, maxit)
% This function assumes A* = A
x = x0;
% estimate step size by power iteration
x = scalarMulParameters(1/sqrt(calcInnerProductParameters(x,x, elements)), x);
maxIt = 3;
for k = 1:maxIt
    Ax = A(x);
    x = A(x);
    x = scalarMulParameters(1/sqrt(calcInnerProductParameters(x,x, elements)), x);
end
omega = 0.9 * 1/calcInnerProductParameters(Ax,Ax,elements);

x = x0;
stopres = zeros(maxit,1);
for k = 1:maxit
    
    Az = A(x);
    
    res = minusParameters(Az, b);
    
    Ar =  A(res);
   
    x = minusParameters(x, scalarMulParameters(omega,Ar));

    stopres(k) = sqrt(calcInnerProductParameters(res, res, elements));

    if stopres(k) < tol
        return
    end
end
end