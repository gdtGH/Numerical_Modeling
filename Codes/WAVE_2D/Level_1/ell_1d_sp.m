function [A]=ell_1d_sp(nu,beta,gam,wx,dx,jacx)
% ELL_1D_SP   Computes 1D local SEM matrix associated to -nu* u''+beta u' +gam u
%
%     [A]=ell_1d_sp(nu,beta,gam,wx,dx,jacx) computes stiffness matrix A:
%         A_{ij} =nu*(phi_j', phi'_i)_N+beta(phi_j',phi_i)_N
%
% Input: 
%        nu   = viscosity (constant>0)
%        beta  = coefficient of first order term (constant)
%        gam  = coefficient of zero-order term (constant>=0)
%        wx= npdx LGL weigths in [-1,1],  
%            (produced by calling [x,wx]=xwlgl(npdx))
%            (npdx = number of nodes, =n+1, if n=polynomial degree)
%        dx= first derivative LGL matrix (produced by calling dx=derlgl(x,npdx))
%        jacx =  jacobian of the map F:[-1,1]---->[a,b]
%
% Output: A = matrix (npdx,npdx) 
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

  n=length(wx);
  coef1=nu/jacx; coef2=gam*jacx;
  A=(coef1*dx'+beta*speye(n))*spdiags(wx,0,n,n)*dx+coef2*spdiags(wx,0,n,n);

return
