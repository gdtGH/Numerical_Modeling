%  TEST Script for testing functions of this directory
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


n=5; % polynomial degree
np=n+1; % number of nodes

% Chebyshev Gauss Lobatto nodes
[x,w]=xwcgl(np);
ww=sum(w);
fig=figure(...
'Name','CGL nodes',...
'Visible', 'on');
plot(x,zeros(np,1),'kx','Markersize',8)
fprintf('CGL: weights sum is %13.6e \n', ww)

% Lagrange interpolation at CGL nodes
u=sin(x);
x_int=(linspace(-1,1,50))';
a=intlag_cgl(x,x_int);
u_int=a*u;
u1=sin(x_int); 
fig=figure(...
'Name','Lagrange interpolation at CGL nodes',...
'Visible', 'on');
plot(x,u,'ro',x_int,u_int,'b',x_int,u1,'k')
legend('physical unknowns','interpolant','original function')



% first derivative CGL matrix
u=sin(x);
[d]=dercgl(x,np);
u1=d*u;
err=norm(u1-cos(x),inf);
fprintf('CGL: Error on the 1st order pseudospectral derivative %13.6e \n', err)

% second derivative CGL matrix
u=sin(x);
[d2]=der2cgl(x,np);
u2=d2*u;
err=norm(u2+sin(x),inf);
fprintf('CGL: Error on the 2nd order pseudospectral derivative %13.6e \n', err)

% Chebyshev Gauss nodes
[x,w]=xwcg(np);
ww=sum(w);
fig=figure(...
'Name','CG nodes',...
'Visible', 'on');
plot(x,zeros(np,1),'kx','Markersize',8)
fprintf('CG: weights sum is %13.6e \n', ww)

pause
% Legendre Gauss Lobatto nodes
[x,w]=xwlgl(np);
ww=sum(w);
fig=figure(...
'Name','LGL nodes',...
'Visible', 'on');
plot(x,zeros(np,1),'kx','Markersize',8)
fprintf('LGL: weights sum is %13.6e \n', ww)

% Lagrange interpolation at LGL nodes
u=sin(x);
x_int=(linspace(-1,1,50))';
a=intlag_lgl(x,x_int);
u_int=a*u;
u1=sin(x_int);
fig=figure(...
'Name','Lagrange interpolation at LGL nodes',...
'Visible', 'on');
plot(x,u,'ro',x_int,u_int,'b',x_int,u1,'k')
legend('physical unknowns','interpolant','original function')

% Lagrange interpolation at LGL nodes (2nd way)
u=sin(x);
x_int=(linspace(-1,1,50))';
[u_int]=legendre_tr_eval(x,u,x_int);
u1=sin(x_int);
fig=figure(...
'Name','Lagrange interpolation at LGL nodesi (2nd)',...
'Visible', 'on');
plot(x,u,'ro',x_int,u_int,'b',x_int,u1,'k')
legend('physical unknowns','interpolant','original function')


% first derivative LGL matrix
u=sin(x);
[d]=derlgl(x,np);
u1=d*u;
err=norm(u1-cos(x),inf);
fprintf('LGL: Error on the 1st order pseudospectral derivative %13.6e \n', err)

% second derivative LGL matrix
u=sin(x);
[d2]=der2lgl(x,np);
u2=d2*u;
err=norm(u2+sin(x),inf);
fprintf('LGL: Error on the 2nd order pseudospectral derivative %13.6e \n', err)

% Legendre Gauss nodes
[x,w]=xwlg(np);
ww=sum(w);
fig=figure(...
'Name','LG nodes',...
'Visible', 'on');
plot(x,zeros(np,1),'kx','Markersize',8)
fprintf('LG: weights sum is %13.6e \n', ww)

% first derivative LG matrix
u=sin(x);
[d]=derlg(x,np);
u1=d*u;
err=norm(u1-cos(x),inf);
fprintf('LG: Error on the 1st order pseudospectral derivative %13.6e \n', err)


