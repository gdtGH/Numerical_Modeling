function [uk]=legendre_tr_coef(x,u);
% LEGENDRE_TR_COEF   Computes Discrete Legendre Transform  coefficients
%
%            formula (4.4.18)  pag. 118, Quarteroni-Valli, Springer 1997
%
%  [uk]=legendre_tr(x,u) returns the discrete coefficients uk of 
%       the expansion of (I_N u)(x) with respect to Legendre polynomials:
%       (I_n u)(x)=\Sum_{k=0}^n u_k L_k(x)             (*)
%
%       where L_k(x) are the Legendre polynomials
%
% Input: x = array of np LGL (Legendre Gauss Lobatto) nodes in [-1,1]
%        u = array of np values of u at LGL nodes x : u_j=u(x_j)
%            size(u)=[np,nc], where nc is the number of transforms which
%            should be computed.
%  If u is a 2-indices array, then every column of u is transformed 
%  following formula (*)
%
% Output: uk = array containing discrete coefficients of the
%              expansion of I_Nu(x) with respect to Legendre polynomials.
%              size(uk)=[np,nc]=size(u)
%
% Reference: A. Quarteroni, A. Valli:
%            "Numerical Approximation of Partial Differential Equations"
%            Springer Verlag / 1997 (2nd Ed) 
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$
%
[np,nc]=size(u); n=np-1;
uk=zeros(np,nc);
nnp1=1/(n*(n+1));

p=pnleg_all(x,n);
ln2=u./((p(:,np)).^2*ones(1,nc));
coef=(1:2:2*np)'*ones(1,nc);
uk=coef.*(p'*ln2);
uk=uk*nnp1;
s=1./p(:,np);
uk(np,:)=(u'*s)'/np;

return
