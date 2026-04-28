function [dx] = GNStep(elements, xk, res, ref, A, timeMeshh, observation, alpha, nCG)

rhs = adjointLinearisedForwardOperatorAllAtOnce(elements, xk, res, A, timeMeshh, observation);

rhs.ds   = -rhs.ds;
rhs.db   = -rhs.db;
rhs.deta = -rhs.deta;
rhs.du1  = -rhs.du1;
rhs.du2  = -rhs.du2;
rhs.du3  = -rhs.du3;

dx = zeroLike(rhs);
r  = rhs;
p  = r;

for k=1:nCG

    Hp = applyGN(xk, p, ref, A, timeMeshh, observation, alpha);

    rr = innerProd(r, r);
    alphaCG = rr / innerProd(elements, timeMeshh, p, Hp);

    dx = addStruct(dx, scale(p, alphaCG));
    r  = addStruct(r, scale(Hp, -alphaCG));

    if sqrt(innerProd(elements, timeMeshh, r, r)) < 1e-6
        break;
    end

    beta = innerProd(elements, timeMeshh, r, r)/rr;
    p = addStruct(r, scale(p, beta));
end

end


