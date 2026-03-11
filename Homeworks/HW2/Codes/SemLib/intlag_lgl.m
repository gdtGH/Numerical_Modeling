function [a]=intlag_lgl(x_lgl, x_new);
% INTLAG_LGL  Computes matrix "a" to evaluate 1D Lagrange interpolant at LGL
% (Legendre Gauss Lobatto)  nodes in [-1,1] at another mesh x_new in [-1,1]
%             (formula (1.2.55)  pag. 17 CHQZ2)
%
%  [a]=intlag_lgl(x_lgl, x_new)
%
% Input: x_lgl = array of np_lgl LGL nodes in [-1,1]
%        x_new = array of np_new another set in [-1,1]
%
% Output: a = matrix of size (np_new,np_lgl) 
%
%       If u_lgl is the array with evaluations  of a function u at nodes
%          x_lgl, it holds:  u_new = a * u_lgl, i.e. u_new=u(x_new)
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

np_lgl=length(x_lgl);
np_new=length(x_new);
n_lgl=np_lgl-1;
a=zeros(np_new,np_lgl);
nnp1=1/(n_lgl*np_lgl);

pn=pnleg(x_lgl,n_lgl);
pn1=pnleg1(x_new,n_lgl);
for i=1:np_new
    for j=1:np_lgl
if abs(x_lgl(j)-x_new(i))>1.e-14
a(i,j)=nnp1*(1-x_new(i)^2)*pn1(i)/((x_lgl(j)-x_new(i))*pn(j));
else
 a(i,j)=1;   
end
    end
end

return
