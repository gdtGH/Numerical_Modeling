function [d] = der2lgl(x,np) 
% DER2LGL      Spectral (Legendre Gauss Lobatto) second derivative matrix 
%
%    [d]=der2lgl(x,np) returns the spectral derivative
%    matrix d at the np LGL nodex X (in [-1,1]). np-1 is the polynomial
%    degree used. Legendre Gauss Lobatto (LGL) grid
%
%
% Input: x = array of LGL nodes on [-1,1] (computed by xwlgl)
%        np = number of LGL nodes (=n+1, n=polynomial interpolation degree)
%
% Output: d = spectral second derivative matrix: formula (2.3.29), pag. 80 CHQZ2
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


 

n=np-1; 
d=zeros(np,np);
for l =1:np
   lnxl = pnleg(x(l),n);    
if l~=1
d(1,l)=(-1)^n/lnxl*(n*np*(1+x(l))-4)*0.5/(1+x(l))^2;
end
   for j = 2:n
      if j ~= l          
         lnxj = pnleg(x(j),n);          
         d(j,l) = -2*lnxj/(lnxl*(x(j)-x(l))^2);          
      else          
         d(j,j) = pnleg2(x(j),n)/(3*lnxl);
      end       
   end 
if l~= np
d(np,l)=1/lnxl*(n*np*(1-x(l))-4)*0.5/(1-x(l))^2;
end
end 
d(1,1)=n*(n+1)*(n^2+n-2)/24;
d(np,np)=d(1,1);
 
return
