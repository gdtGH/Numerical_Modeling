function [a]=intlag_lg(x_lg, w_lg, x_new);
% INTLAG_LG  Computes matrix "a" to evaluate 1D Lagrange interpolant at LG 
%
%  (Legendre Gauss) nodes in [-1,1] at another mesh x_new in [-1,1]
%
%  [a]=intlag_lg(x_lg, w_lg, x_new)
%
% Input: x_lg = array of np_lg LG nodes in [-1,1]
%        w_lg = array of np_lg LG weights in [-1,1]
%        x_new = array of np_new another set in [-1,1]
%
% Output: a = matrix of size (np_new,np_lg) 
%
%       If u_lg is the array with evaluations  of a function u at nodes
%          x_lg, it holds:  u_new = a * u_lg, i.e. u_new=u(x_new)
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

np_lg=length(x_lg);
np_new=length(x_new);
n_lg=np_lg-1;
a=zeros(np_new,np_lg);

pn=pnleg_all([x_lg;x_new],n_lg);
pn1=pn(np_lg+1:np_lg+np_new,:);
pn=pn(1:np_lg,:);
nor=0.5+(0:n_lg);
for i=1:np_new
for j=1:np_lg
a(i,j)=sum(nor.*pn(j,:).*pn1(i,:))*w_lg(j);
end
end

return
