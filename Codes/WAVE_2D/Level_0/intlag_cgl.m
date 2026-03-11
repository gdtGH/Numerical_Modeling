function [a]=intlag_cgl(x_cgl, x_new);
% INTLAG_CGL  Computes matrix "a" to evaluate 1D Lagrange interpolant at CGL 
%
%    (Chebyshev Gauss Lobatto)  nodes in [-1,1] at another mesh x_new in [-1,1]
%             (formula (2.4.30)  pag. 88 CHQZ2)
%
%  [a]=intlag_cgl(x_cgl, x_new) 
%
% Input: x_cgl = array of np_cgl CGL nodes in [-1,1] (ordered from left to
%                   right)
%        x_new = array of np_new another set in [-1,1] (ordered from left to
%                   right)
%
% Output: a = matrix of size (np_new,np_cgl) 
%
%
%       If u_cgl is the array with evaluations  of a function u at nodes
%          x_cgl, it holds:  u_new = a * u_cgl, i.e. u_new=u(x_new)
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

x_cgl=-x_cgl;
np_cgl=length(x_cgl);
np_new=length(x_new);
n=np_cgl-1;
a=zeros(np_new,np_cgl);

[pn,pn1]=jacobi_eval(x_new,n,-0.5,-0.5);
pn1(:,1)=pn1(:,1)*2^(2*n)*(prod(1:n))^2/(prod(1:2*n));
c=[2;ones(n-1,1);2];

for i=1:np_new
for l=1:np_cgl
if abs(x_cgl(l)-x_new(i))>1.e-14
a(i,np_cgl+1-l)=(-1)^(l)*(1-x_new(i)^2)*pn1(i,1)/(c(l)*n*n*(x_new(i)-x_cgl(l)));
else
a(i,np_cgl+1-l)=1;   
end
end
end

return
