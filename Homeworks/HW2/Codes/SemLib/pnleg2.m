function [p2,p1,p] = pnleg2 (x, n) 
% PNLEG2    Evaluates the second derivative of Legendre polynomial
%           of degree n. 
%
%    [p2,p1,p]=pnleg2(x,n)  evaluates (L_n)''(x), (L_n)'(x), L_n(x)
%    at the node(s) x, using the three term relation  (2.3.3), pag. 75 CHQZ2.
%
% Input: x = scalar or column array
%        n = degree of Legendre polynomial
%
% Output: p2 = scalar or column array (same dimension as x)
%             with the evaluation of (L_n)''(x)
%         p1 = scalar or column array (same dimension as x)
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
p=0;p1=p;p2=p;
if n<=1
p=0;p1=p;p2=p;
elseif n==2
p=1;p1=p;p2=p;
else
p1=1;  p2=x;  
p11=0; p21=1; 
p12=0; p22=0;
for k =1:n-1    
   p3=((2*k+1)*x*p2-k*p1)/(k+1); 
   p31=((2*k+1)*(x*p21+p2)-k*p11)/(k+1); 
   p32=((2*k+1)*(x*p22+p21*2)-k*p12)/(k+1); 
   p1=p2; p2=p3; 
   p11=p21; p21=p31; 
   p12=p22; p22=p32; 
end 
p2=p32; p1=p31;p=p3;
end

% x is an array
else
x=x(:);
p=zeros(nn,1);p1=p;p2=p;
if n<=1
p=zeros(nn,1);p1=p;p2=p;
elseif n==2
p=ones(nn,1);p1=p;p2=p;
else
p1=ones(nn,1);  p2=x;  
p11=zeros(nn,1); p21=ones(nn,1); 
p12=zeros(nn,1); p22=zeros(nn,1);
for k =1:n-1    
   p3=((2*k+1)*x.*p2-k*p1)/(k+1); 
   p31=((2*k+1)*(x.*p21+p2)-k*p11)/(k+1); 
   p32=((2*k+1)*(x.*p22+p21*2)-k*p12)/(k+1); 
   p1=p2; p2=p3; 
   p11=p21; p21=p31; 
   p12=p22; p22=p32; 
end 
p2=p32; p1=p31;p=p3;
end
end
return 
