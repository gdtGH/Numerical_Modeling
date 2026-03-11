function [p] = pnleg_all(x,n) 
%PNLEG_ALL      Evaluates Legendre polynomials, from degree 0 to n
%
%    [p]=pnleg_all(x,n) returns the evaluation of the Legendre
%    polynomials of degree 0,1,...,N at x, using the three
%    term relation (2.3.3), pag. 75 CHQZ2.
%
% Input: x = scalar or array
%        n = degree of Legendre polynomial to be evalueted
%
% Output: p = 1-index (size(p)=size(x)) or 
%             2-indexes (size(p)=[length(x),n+1])
%             array containing the 
%             evaluation of L_0(x), L_1(x), .... L_n(x)
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


nn=size(x);
if nn==1
if n==0
p=1;
elseif n==1
    p=[1,x];
else
p1=1;  p2=x;  p3=p2; p=[p1,p2];
for k =1:n-1    
   p3=((2*k+1)*x*p2-k*p1)/(k+1); 
   p1=p2; 
   p2=p3; 
   p=[p,p3]; 
end 
end

else
if n==0
p=ones(nn);
elseif n==1
    p=[ones(nn),x];
else
p1=ones(nn);  p2=x;  p3=p2; p=[p1,p2];
for k =1:n-1    
   p3=((2*k+1)*x.*p2-k*p1)/(k+1); 
   p1=p2; 
   p2=p3; 
   p=[p,p3]; 
end 
end
 
end
return 
