function [x, iter, stopres] = conjugateGradient(elements, A, b, x0, tol, maxit)
x = x0;
res = minusParameters(b,A(x)); %Az = Axn = y
pk = res;
stopres = zeros(maxit,1);
betak = zeros(maxit,1);
d = betak;
for iter = 1:maxit
    Apk = A(pk);
    rr = calcInnerProductParameters(res,res, elements);
    d(iter) = rr / (calcInnerProductParameters(pk, Apk ,elements));
    x = addParameters(x, scalarMulParameters(d(iter), pk));
    % the following line may propagate numerical errors
    resNew = minusParameters(res, scalarMulParameters(d(iter), Apk));
    %resNew = minusParameters(b, A(x));
    rrN = calcInnerProductParameters(resNew, resNew, elements);
    stopres(iter) = sqrt(rrN);

    if stopres(iter) < tol
        break;
    end

    betak(iter) = rrN/rr;

    pk = addParameters(resNew, scalarMulParameters(betak(iter), pk));
    res = resNew;
end
end

