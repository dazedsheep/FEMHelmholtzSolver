function [df] = numGradient(f,x)
h = 0.001;
for i=1:size(x,2)
    v = zeros(1,size(x,2));
    v(i) = 1;
    df(i) = 1/(2*h) *(f(x-h*v) + f(x+h*v));
   
end

end

