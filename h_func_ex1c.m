function [val] = h_func_ex1c(x, y, beta, origu, gradu, omega, gamma)

% if(x==0 && y==0)
%     val = (gradu(x,y).') * [-1;-1]./sqrt(2) + (1i*beta*omega + gamma)*origu(x,y);
%     return
% end
% 
% if(x==0 && y==1)
%     val = (gradu(x,y).') * [-1;1]./sqrt(2) + (1i*beta*omega + gamma)*origu(x,y);
%     return
% end
% 
% if(x==1 && y==1)
%     val = (gradu(x,y).') * [1;1]./sqrt(2) + (1i*beta*omega + gamma)*origu(x,y);
%     return
% end
% 
% 
% if(x==1 && y==0)
%     val = (gradu(x,y).') * [1;-1]./sqrt(2) + (1i*beta*omega + gamma)*origu(x,y);
%     return
% end

if (y==0 && x >= 0)
    val = gradu(x,y).'*[0;-1] + (1i*beta*omega + gamma)*origu(x,y);
    return
end

if (y==1 && x >= 0)
    val = gradu(x,y).'*[0;1] + (1i*beta*omega + gamma)*origu(x,y);
    return
end

if (x==1 && y >= 0)
    val = gradu(x,y).'*[1;0] + (1i*beta*omega + gamma)*origu(x,y);
    return
end

if (x==0 && y >= 0)
    val = gradu(x,y).'*[-1;0] + (1i*beta*omega + gamma)*origu(x,y);
    return
end


val = 0;

end

