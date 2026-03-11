function [x,w] = xwcgl(np,a,b)
%XWCGL  Computes nodes and weights of the Chebyshev-Gauss-Lobatto quadrature formula.
%
%    [x,w]=xwcgl(np) returns the np weigths and nodes 
%    of the corresponding Chebyshev Gauss-Lobatto quadrature 
%    formula in the reference interval [-1,1].
%
%    [x,w]=xwcgl(np,a,b) returns the NP weigths and the nodes 
%    of the corresponding Chebyshev Gauss-Lobatto quadrature 
%    formula in the  interval [a,b].
%
% Input: np = number of nodes
%        a, b = extrema of the interval
%
% Output: x(np,1) = CGL nodes  (CHQZ2, (2.4.14), pag. 86)
%         w(np,1) = CGL weigths (CHQZ2, (2.4.14), pag. 86)
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


if np <= 1 
  es='The number of the quadrature nodes should be greater than 1';
  disp(es); x = []; w = [];
  return 
end 

n=np-1;
ww=pi/n;
w=ww*ones(np,1);
x(1:np,1)=-cos((0:n)'*ww);
w(1)=w(1)*.5;
w(np)=w(1);

if nargin ==3
   bma=(b-a)*.5;
   bpa=(b+a)*.5;
   x=bma*x+bpa;
   w=w*bma;
end

return
