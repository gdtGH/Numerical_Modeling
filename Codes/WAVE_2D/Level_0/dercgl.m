function [d] = dercgl(x,np) 
%DERCGL      Spectral (Chebyshev Gauss Lobatto) derivative matrix 
%
%    [D]=DERCGL(X,NP) returns the spectral derivative
%    matrix D at the NP CGL nodes X (in [-1,1]). NP-1 is the polynomial
%    degree used. Chebyshev Gauss Lobatto (CGL) grid ORDERED FROM LEFT TO RIGHT
%
%
% Input: x = array of CGL nodes on [-1,1] (computed by xwcgl)
%        np = number of CGL nodes (=n+1, n=polynomial interpolation degree)
%
% Output: d = spectral derivative matrix: formula  (2.4.31), pag. 89 CHQZ2
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$
  

n=np-1; 
x=x(:);

if n==0 d=0;  return; end; 

c=[2;ones(n-1,1);2];
for j=2:n
d(j,j)=-x(j)/(2*(1-x(j)^2));
end
for j=1:np
xj=x(j);
for l=1:j-1
d(j,l)=c(j)/c(l)*(-1)^(l+j)/(xj-x(l));
end
for l=j+1:np
d(j,l)=c(j)/c(l)*(-1)^(l+j)/(xj-x(l));
end
end
d(np,np)=(2*n^2+1)/6;
d(1,1)=-d(np,np);
d(1,np)=(-1)^(np+1)/(x(1)-x(np));
d(np,1)=-d(1,np);
return
