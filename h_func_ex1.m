function [val] = h_func_ex1(x, y)

k = 20;
a = k*cos(pi/2);
b = k*sin(pi/2);
beta = 1;


origu = @(x,y,a,b) (cos(a*x+b*y)+1i*sin(a*x + b*y));
gradu = @(x, y, a, b) [-a*sin(a*x + b*y) + 1i*a*cos(a*x + b*y), -b*sin(a*x + b*y) + 1i*b*cos(a*x + b*y)];

%deal with the corner nodes separately 
% if (x==0 && y==0)
%     val = gradu(x,y,a,b)*[-1;0] + 1i*beta*k*origu(x,y,a,b);
%     return
% end
% 
% if (x==1 && y==0)
%     val = gradu(x,y,a,b)*[1;0] + 1i*beta*k*origu(x,y,a,b);
%     return
% end
% 
% if (x==0 && y==1)
%     val = gradu(x,y,a,b)*[-1;0]  + 1i*beta*k*origu(x,y,a,b);
%     return
% end
% 
% if (x==1 && y==1)
%     val = gradu(x,y,a,b)*[1;0] + 1i*beta*k*origu(x,y,a,b);
%     return
% end

if (y==0 && x >= 0)
    val = gradu(x,y,a,b)*[0;-1] + 1i*beta*k*origu(x,y,a,b);
    return
end

if (y==1 && x >= 0)
    val = gradu(x,y,a,b)*[0;1] + 1i*beta*k*origu(x,y,a,b);
    return
end

if (x==1 && y >= 0)
    val = gradu(x,y,a,b)*[1;0] + 1i*beta*k*origu(x,y,a,b);
    return
end

if (x==0 && y >= 0)
    val = gradu(x,y,a,b)*[-1;0] + 1i*beta*k*origu(x,y,a,b);
    return
end

val = 0;


end

