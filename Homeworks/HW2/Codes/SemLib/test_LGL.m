%  TEST Script for testing functions of this directory
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio and Ilario Mazzieri
%   $Date: 2007/04/01$

close all; clear; clc;

n=10; % polynomial degree
np=n+1; % number of nodes


% Legendre Gauss nodes
[x,w] = xwlg(np);
wwlg = sum(w);
figure(1);
plot(x,zeros(np,1),'kx','Markersize',8); hold on;

% Legendre Gauss Lobatto nodes
[x,w] = xwlgl(np);
wwlgl = sum(w);
plot(x,zeros(np,1),'ro','Markersize',8)
legend('LG nodes','LGL nodes');

%fprintf('LGL: weights sum is %13.6e \n', ww)
%fprintf('LG: weights sum is %13.6e \n', ww)

% Basis Functions in [-1,1]
figure(2);
plot(x,zeros(np,1),'ro','Markersize',8); hold on;

x_int = linspace(-1,1);
for i = 1 : np
    [p1,p] = pnleg1(x_int,n);
    [pi]   = pnleg(x(i), n); 
    y_int = - 1/(n*(n+1)) * ((1-x_int.^2).*p1') ./ ((x_int-x(i))*pi);
    if(i == 1) y_int(1) = 1; end
    if(i == np) y_int(end) = 1; end
    plot(x_int,y_int); hold on;
end


% % first derivative LGL matrix
% u=sin(x);
% [d]=derlgl(x,np);
% u1=d*u;
% err=norm(u1-cos(x),inf);
% fprintf('LGL: Error on the 1st order pseudospectral derivative %13.6e \n', err)
% 
% % %second derivative LGL matrix
% u=sin(x);
% [d2]=der2lgl(x,np);
% u2=d2*u;
% err=norm(u2+sin(x),inf);
% fprintf('LGL: Error on the 2nd order pseudospectral derivative %13.6e \n', err)


