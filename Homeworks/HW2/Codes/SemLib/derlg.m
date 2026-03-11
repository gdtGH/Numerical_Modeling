function [d] = derlg(x,np) 
%DERLG      Spectral (Legendre Gauss) derivative matrix 
%
%    [d]=derlg(x,np) returns the spectral Legendre Gauss derivative
%    matrix d at the np LG nodes X (in (-1,1)). np-1 is the polynomial
%    degree used. Legendre Gauss (LG) grid
%
%
% Input: x = array of LG nodes on (-1,1) (computed by xwlg)
%        np = number of LG nodes (=n+1, n=polynomial interpolation degree)
%
% Output: d = spectral derivative matrix
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

 
n=np-1; 
d=zeros(np);
 [p,pd] = jacobi_eval(x,np,0,0);
for i=1:np
for j=1:np
if(i~=j)
    d(i,j)=pd(i,1)/(pd(j,1)*(x(i)-x(j)));
elseif(i==j)
d(i,j)=x(i)/(1-x(i)*x(i));
end
end
end
 
return
