function [d] = der2cgl(x,np) 
% DER2CGL    Spectral (Chebyshev Gauss Lobatto) second derivative matrix 
%
%    [d]=der2cgl(x,np) returns the spectral second derivative
%    matrix d at the np CGL nodes X (in [-1,1]). np-1 is the polynomial
%    degree used. Chebyshev Gauss Lobatto (CGL) grid ORDERED FROM LEFT TO RIGHT
%
%
% Input: x = array of CGL nodes on [-1,1] (computed by xwcgl)
%        np = number of CGL nodes (=n+1, n=polynomial interpolation degree)
%
% Output: d = spectral derivative matrix: formula  (2.4.32), pag. 89 CHQZ2
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

  

n=np-1; 
x=x(:);

if n==0 d=0;  return; end; 

d=zeros(np,np);
dueterzi=2/3;
c=[2;ones(n-1,1);2];
duenp1=2*n*n+1;
d(1,1)=(n^4-1)/15;
for l =2:np
d(1,l)=dueterzi*(-1)^(l-1)/c(l)*(duenp1*(1+x(l))-6)/(1+x(l))^2;
end
for j = 2:n
for l =1:np
  if j ~= l          
  d(j,l) = (-1)^(j+l-2)/c(l)*(x(j)*(x(j)+x(l))-2)/((1-x(j)^2)*(x(j)-x(l))^2);
  else
  d(j,j) = -((n*n-1)*(1-x(j)^2)+3)/(3*(1-x(j)^2)^2);
  end       
end
end
for l=1:n
d(np,l)=dueterzi*(-1)^(l-1+n)/c(l)*(duenp1*(1-x(l))-6)/(1-x(l))^2;
end
d(np,np)=d(1,1);
 
return
