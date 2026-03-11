function [x,w] = xwlg(np,a,b)
%XWLG  Computes nodes and weights of the Legendre-Gauss  quadrature formula.
%
%    [x,w]=xwlg(np) returns the np weigths and nodes 
%    of the corresponding Legendre Gauss quadrature 
%    formula in the reference interval (-1,1).
%
%    [x,w]=xwlg(np,a,b) returns the np weigths and the nodes 
%    of the corresponding Legendre Gauss quadrature 
%    formula in the  interval (a,b).
%
% Input: np = number of nodes
%        a, b = extrema of the interval
%
% Output: x(np,1) = LG nodes  (CHQZ2, (2.3.10), pag. 76)
%         w(np,1) = LG weigths (CHQZ2, (2.3.10), pag. 76)
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


n=np-1;
if np<=1
  x=0;w=2;
  return
end
x=jacobi_roots(np,0,0);
w=2./(pnleg1(x,np).^2.*(1-x.^2));



%
% map on (a,b)
%
if nargin == 3
  bma=(b-a)*.5;
  bpa=(b+a)*.5;
  x=bma*x+bpa;
  w=w*bma;
end
return
