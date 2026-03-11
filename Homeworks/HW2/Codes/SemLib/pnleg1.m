function [p1,p] = pnleg1 (x, n) 
% PNLEG1    Evaluates the first derivative of Legendre polynomial
%           of degree n 
%
%    [p1,p]=pnleg1(x,n)  evaluates (L_n)'(x), L_n(x)
%    at the node(s) x, using the three term relation  (2.3.3), pag. 75 CHQZ2.
%
% Input: x = scalar or column array
%        n = degree of Legendre polynomial
%
% Output: p1 = scalar or column array (same dimension as x)
%             with the evaluation of (L_n)'(x)
%         p  = scalar or column array (same dimension as x)
%             with the evaluation of L_n(x)
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


nn=size(x);

if nn==1
    % x is a scalar
    p=0; p1=0;
if n==0
p=1;p1=0;
elseif n==1
p=x;p1=1;
else
p1=0;p2=0;p11=0;p21=0;p3=0;p31=0;
p1=1;  p2=x; 
p11=0; p21=1;
for k =1:n-1
   duekp1=2*k+1;kp1=1/(k+1);
   p3=(duekp1*x*p2-k*p1)*kp1; 
   p31=(duekp1*(x*p21+p2)-k*p11)*kp1; 
   p11=p21; p21=p31; 
   p1=p2; p2=p3; 
end 
p1=p31; p=p3;
end

else
x=x(:);
    % x is an array
p=zeros(nn); p1=p;
if n==0
p=ones(nn);p1=zeros(nn);
elseif n==1
p=x;p1=ones(nn);
else
p1=p;p2=p;p11=p;p21=p;p3=p;p31=p;
p1=1;  p2=x; 
p11=0; p21=1;
for k =1:n-1
   duekp1=2*k+1;kp1=1/(k+1);
   p3=(duekp1*x.*p2-k*p1)*kp1; 
   p31=(duekp1*(x.*p21+p2)-k*p11)*kp1; 
   p11=p21; p21=p31; 
   p1=p2; p2=p3; 
end 
p1=p31; p=p3;
end
end
return 
