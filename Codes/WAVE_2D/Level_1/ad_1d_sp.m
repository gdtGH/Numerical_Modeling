function [A]=ad_1d_sp(nu,beta,wx,dx,jacx)
% AD_1D_SP Computes 1D local spectral matrix associated to ad operator -nu u'' + beta u'
%
%    advection-diffusion operator (one spectral element)  -nu u'' + beta u'
%
%    [A]=ad_1d_sp(nu,beta,wx,dx,jacx) produces the matrix 
%        A_{ij}=(nu phi_j', phi'_i)_N +beta (phi_j' ,phi_i)_N
%
% Input: nu   = viscosity (constant>0)
%        beta = coefficient of first order term (constant)
%        wx = npdx LGL weigths in [-1,1],
%           (produced by calling [x,w]=xwlgl(npdx))
%           (npdx = number of nodes, =n+1, if n=polynomial degree)
%        dx = first derivative LGL matrix (produced by calling d=derlgl(x,npdx))
%        jacx =  jacobian of the map F:[-1,1]---->[xa_ie,b_ie]
%
% Output: A = matrix (npdx,npdx) defined above
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

n=length(wx);
coef1=nu/jacx;
A=(coef1*dx'+beta*speye(n))*spdiags(wx,0,n,n)*dx;

