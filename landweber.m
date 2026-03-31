function [x, k, stopres] = landweber(elements, A, b, x0, omega, tol, maxit)
% This function assumes A* = A
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