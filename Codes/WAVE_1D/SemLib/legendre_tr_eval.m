function [u_int]=legendre_tr_eval(x,u,x_int);
% LEGENDRE_TR_EVAL     Evaluates Discrete Legendre Transform
%
%            formula (4.4.17)  pag. 118, Quarteroni-Valli, Springer 1997
%            with coefficients computed with formula (4.4.18)  
%                   pag. 118, Quarteroni-Valli, Springer 1997
%
%  [u_int]=legendre_tr_eval(x,u,x_int) returns the evaluation at nodes x_int
%       of the expansion of (I_N u)(x) with respect to Legendre polynomials
%
%       (I_n u)(x)=\Sum_{k=0}^n u_k L_k(x)             (*)
%
%       where L_k(x) are the Legendre polynomials
%
% Input: x = array of np Legendre Gauss Lobatto (LGL) nodes on [-1,1]
%        u = array of np values of u at LGL nodes x : u_j=u(x_j)
%            size(u)=[np,nc], where nc is the number of transforms which
%            should be computed.
%        If u is a 2-indices array, then every column of u is transformed
%        following formula (*)
%        x_int = array of a second set of nodes in [-1,1]
%
% Output: u_int = array of the values of  the
%              expansion of I_Nu(x) with respect to Legendre polynomials.
%
% Reference: A. Quarteroni, A. Valli:
%            "Numerical Approximation of Partial Differential Equations"
%            Springer Verlag / 1997 (2nd Ed)
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


n1=length(x_int);
np=length(x); n=np-1;

uk=legendre_tr_coef(x,u);
pn=pnleg_all(x_int,n);
u_int=pn*uk;

return
