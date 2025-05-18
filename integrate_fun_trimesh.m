function [area,integral] = integrate_fun_trimesh(p,t,f)

%--------------------------------------------------------------------
%
% Integrate a function f evaluated at the nodes p of 2D linear
% triangulation t, using linear interpolation of the function. 
%
% The triangulation can be arbitrary, comprising disconnected or 
% multiply-connnected domains.
%
% Inputs:  p is a 2xN matrix of [x;y] nodal coordinates.
%          t is a 3xE node connectivity matrix for the nodes in p.
%          f is a N-vector of values of the function you wish to 
%            integrate over the triangulation, given at the nodes p.
%
% Outputs: area is the area of the triangulation.
%          integral is the integral of f over the triangulation.
%
% Author: Tom Montenegro-Johnson
% http://www.damtp.cam.ac.uk/user/tdj23
% Release date: Nov 2015
%
%--------------------------------------------------------------------

%Corner point indices
it1 = t(1,:);
it2 = t(2,:);
it3 = t(3,:);

%Affine transform to canonical element
detadx = p(2,it1) - p(2,it2);
detady = p(1,it2) - p(1,it1);
dpsidx = p(2,it3) - p(2,it1);
dpsidy = p(1,it1) - p(1,it3);

%Jacobian of transform
jac = dpsidx.*detady-dpsidy.*detadx;
jac = 0.5*jac;

%Gauss points on canonical triangle
[psi,eta,w] = gausspoints;
weights = w(:)*ones(1,length(jac));

%Shape functions on linear triangle
[phi,~,~] = shapes(psi,eta);

%Corner points of function
f1 = f(it1);
f2 = f(it2);
f3 = f(it3);
f_out = phi(1,:)'*f1 + phi(2,:)'*f2 + phi(3,:)'*f3;

%Do the integration
area = sum(jac);
elm_integrals = sum(f_out.*weights,1).*jac;
integral = sum(elm_integrals);

end

function [psi,eta,w] = gausspoints

    %Degree of precision 2, sufficient for linear interpolation
    psi=[2/3 1/6 1/6];
	eta=[1/6 1/6 2/3];
	w=[1/3 1/3 1/3];

end

function [phi,dphidx,dphidy] = shapes(x,y)

    %Linear shape functions on canonical triangle
    phi    = [1-x-y  ; x     ; y];
    dphidx = [-1+0*x ; 1+0*x ; 0*x];
    dphidy = [-1+0*x ; 0*x   ; 1+0*x];

end