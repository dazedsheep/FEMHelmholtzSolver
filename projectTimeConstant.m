function [x] = projectTimeConstant(timeMeshh, x)
    T = timeMeshh*(size(x.du,1));
    x.deta = sum( x.deta .* timeMeshh, 1)./T;
    x.db = sum( x.db .* timeMeshh, 1)./T;
    x.ds = sum( x.ds .* timeMeshh, 1)./T;
    
end