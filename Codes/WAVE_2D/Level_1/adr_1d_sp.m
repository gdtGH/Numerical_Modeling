function [A]=adr_1d_sp(wx,dx,jacx,nu,b,gam)
% ADR_1D_SP Computes 1D local spectral matrix associated to adr operator -(nu u' + b(x) u)'+gam u
%
%          advection-diffusion-raction operator (one spectral element)
%
%  [A]=adr_1d_sp(wx,dx,jacx,nu,b,gam) produces the matrix 
%
% Input:
%      A_{ij}=(nu phi_j' + b phi_j, phi'_i)_N +(gam phi_j,phi_i)_N
%
%      wx = npdx LGL weigths in [-1,1],
%          (produced by calling [x,w]=xwlgl(npdx))
%          (npdx = number of nodes, =n+1, if n=polynomial degree)
%      dx = first derivative LGL matrix (produced by calling d=derlgl(x,npdx))
%      jacx =  jacobian of the map F:[-1,1]---->[xa_ie,b_ie]
%      nu   = viscosity (constant>0)
%      b = column vector with evaluation of b(x) at LGL node of the spectral
%          element
%      gam  = coefficient of zeroth order term (constant>0)
%
% Output: A = matrix (npdx,npdx) defined above
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

  coef1=nu/jacx; coef2=gam*jacx;
  A=dx'*(diag(wx)*(coef1*dx+diag(b)))+coef2*diag(wx);

