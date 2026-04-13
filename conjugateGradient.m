function [x, iter, stopres] = conjugateGradient(elements, A, b, x0, tol, maxit, innerProduct, add, minus, scalarMul, precond)
x = x0;
res = minus(b,A(x)); 
zk = precond(res);
pk = zk;
stopres = zeros(maxit,1);
betak = zeros(maxit,1);
d = betak;
for iter = 1:maxit
    Apk = A(pk);
    rr = innerProduct(res,zk); % (r_{k}, z_{k})
    d(iter) = rr / (innerProduct(pk, Apk)); % \alpha_{k+1}  = (r_{k}, z_{k}) / (p_k, A(p_k))
    x = add(x, scalarMul(d(iter), pk)); % x_{k+1} = x_k + \alpha_{k+1}*p_k 
    % the following line may propagate numerical errors
    resNew = minus(res, scalarMul(d(iter), Apk)); % r_{k+1} = r_k - \alpha_{k+1} * A(p_k)
    %resNew = minusParameters(b, A(x));
    zk = precond(resNew);
    rrN = innerProduct(resNew, zk);
    stopres(iter) = sqrt(rrN);

    if stopres(iter) < tol
        break;
    end
    betak(iter) = innerProduct(resNew, zk)/rr; % \beta_{k+1} = (r_{k+1}, z_{k+1}) / (r_k, z_k)
    pk = add(zk, scalarMul(betak(iter), pk)); % p_{k+1} = z_{k+1} + \beta_{k+1} p_{k}
    res = resNew;
end
end
