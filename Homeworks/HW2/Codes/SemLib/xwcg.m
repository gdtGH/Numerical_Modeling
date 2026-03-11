function [x,w] = xwcg(np,a,b)
%XWCG     Computes  nodes and weights of the Chebyshev-Gauss quadrature formula
%
%    [x,w]=xwcg(np) returns the np weigths and nodes 
%    of the corresponding Chebyshev Gauss quadrature 
%    formula in the reference interval (-1,1).
%
%    [x,w]=xwcg(np,a,b) REturns the np weigths and the nodes 
%    of the corresponding Chebyshev Gauss quadrature 
%    formula in the  interval (a,b).
%
% Input: np = number of nodes
%        a, b = extrema of the interval
%
% Output: x(np,1) = CG nodes  (CHQZ2, (2.4.12), pag. 85)
%         w(np,1) = CG weigths (CHQZ2, (2.4.12), pag. 85)
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


n=np-1;
if np<=1
  x=0;w=pi*(b-a)*.5;
  return
end
den=pi/(2*np);
w=pi/np*ones(np,1);
x(1:np,1)=-cos((2*(0:n)'+1)*den);

%
% mappatura su (a,b)
%
if nargin == 3
  bma=(b-a)*.5;
  bpa=(b+a)*.5;
  x=bma*x+bpa;
  w=w*bma;
end
return
