function [a]=legendre_tr_matrix(x);
% LEGENDRE_TR_MATRIX  Computes matrix to produce Discrete Legendre Transform 
%
%            formula (4.4.18)  pag. 118, Quarteroni-Valli, Springer 1997
%
%  [a]=legendre_tr_matrix(x) returns the matrix which maps the vector
%         of physical uknowns of u at LGL nodes, to the vector of frequency
%         unknowns with respect to Legendre expansion.
%
%       (I_n u)(x)=\Sum_{k=0}^n u^*_k L_k(x)             (*)
%
%       where L_k(x) are the Legendre polynomials
%
% Input: x = array of np (Legendre Gauss Lobatto) LGL nodes on [-1,1]
%
% Output: a = matrix (np,np)  associated to the map u_j ---> u^*_k
%
% Reference: A. Quarteroni, A. Valli:
%            "Numerical Approximation of Partial Differential Equations"
%            Springer Verlag / 1997 (2nd Ed)
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

np=length(x);
n=np-1;
a=zeros(np);
nnp1=1/(n*np);

pn=pnleg_all(x,n);
pn2=(pn(:,np)).^2;
for k=0:n-1
a(k+1,:)=(2*k+1)*nnp1*(pn(:,k+1))'./pn2';
end
a(np,:)=1./(np.*(pn(:,np))');

return
