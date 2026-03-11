function [A]=stiff_1d_sp(w,d,jac)
% STIFF_1D_SP    Computes local stiffness matrix for 1D problem
%
%     [A]=stiff_1d_sp(w,d,jac) computes stiffness matrix A:
%         A_{ij} = ( (phi_j)',(phi_i)')_N
%         by Galerkin-Numerical Integration
%
% Input: 
%        w = npdx LGL weigths in [-1,1],  
%            (produced by calling [x,w]=xwlgl(npdx))
%            (npdx = number of nodes, =n+1, if n=polynomial degree)
%        d = first derivative LGL matrix (produced by calling d=derlgl(x,npdx))
%        jac =  jacobian of the map F:[-1,1]---->[a,b]
%
% Output: A = matrix (npdx,npdx) with the stiffness contribution
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

A=d'*diag(w)*d;
A=A/jac;
return
