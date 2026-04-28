function [x, iter, stopres] = conjugateGradient(elements, A, b, x0, tol, maxit, innerProduct, add, minus, scalarMul, plot_handles, newtonIter)
x = x0;
res = minus(b,A(x)); %Ax = Axn = y
pk = res;
stopres = zeros(maxit,1);
betak = zeros(maxit,1);
d = betak;
for iter = 1:maxit
    Apk = A(pk);
    rr = innerProduct(res,res);
    d(iter) = rr / (innerProduct(pk, Apk));
    x = add(x, scalarMul(d(iter), pk));
    % the following line may propagate numerical errors
    resNew = minus(res, scalarMul(d(iter), Apk));
    %resNew = minusParameters(b, A(x));
    rrN = innerProduct(resNew, resNew);
    stopres(iter) = sqrt(rrN);

    if stopres(iter) < tol
        break;
    end
    
    updateVideoFigure(x.refState_1.s0, x.refState_2.b0, projectParameters(x).refState_3.eta0, plot_handles, newtonIter, iter); % update figure in every iteration

    betak(iter) = rrN/rr;

    pk = add(resNew, scalarMul(betak(iter), pk));
    res = resNew;
end
end

