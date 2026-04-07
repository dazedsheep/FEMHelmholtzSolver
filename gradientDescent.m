function [x, iter, stopres] = gradientDescent(elements, A, b, x0, tol, maxit)
x = x0;
res = minusParameters(b,A(x)); %Az = Axn = y

stopres = zeros(maxit,1);
betak = zeros(maxit,1);
d = betak;
for iter = 1:maxit
    Ar = A(res);
    rr = calcInnerProductParameters(res,res, elements);
    d(iter) = rr / (calcInnerProductParameters(res, Ar ,elements));
    x = addParameters(x, scalarMulParameters(d(iter), res));
    stopres(iter) = sqrt(rr);

    if stopres(iter) < tol
        break;
    end

    res = minusParameters(res, scalarMulParameters(d(iter), res));
end
end
