function [uex,uexx,uexy,ff,g,h,gam]=setfun_lap_2d
%  SETFUN_LAP_2D  Sets functions and coefficients for lap_2d
%
%    [uex,uexx,uexy,ff,gam]=setfun_lap_2d
%
% Output: uex =@(x,y)[....] function handle to the expression of exact
%               solution
%         uexx =@(x,y)[....] function handle to the expression of the first 
%               x-derivative of the exact solution
%         uexy =@(x,y)[....] function handle to the expression of the first 
%               y-derivative of the exact solution
%         ff =@(x,y)[....] function handle to the expression of function 
%              at right hand side
%         g =@(x,y)[....] function handle to the expression of Dirichlet 
%              boundary data
%         h =@(x,y)[....] function handle to the expression of Neumann 
%              boundary data. It is a vector of 4 functions, each for any
%              side. 
%         gam = coefficient of zero-order term (constant >=0)
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


syms  x y
gam=0;
%Uex=sin(2*pi*x)*cos(2*pi*y); % exact solution
Uex=sin(pi*x)*cos(pi*y);
Uexx=diff(Uex,x); % first x-derivative of exact solution
Uexy=diff(Uex,y); % first y-derivative of exact solution
Uexx2=diff(Uexx,x); % second x-derivative of exact solution
Uexy2=diff(Uexy,y); % second y-derivative of exact solution
Ff=(gam*(Uex)-(Uexx2+Uexy2)); % right hand side
H1=-diff(Uex,y);  % Neumann data on Side 1
H2=diff(Uex,x);  % Neumann data on Side 2
H3=diff(Uex,y);  % Neumann data on Side 3
H4=-diff(Uex,x);  % Neumann data on Side 4

uex=vectorize(char(Uex));
uexx=vectorize(char(Uexx));
uexy=vectorize(char(Uexy));
ff=vectorize(char(Ff));
h1=vectorize(char(H1));
h2=vectorize(char(H2));
h3=vectorize(char(H3));
h4=vectorize(char(H4));

uex=@(x,y)[eval(uex)];
uexx=@(x,y)[eval(uexx)];
uexy=@(x,y)[eval(uexy)];
ff=@(x,y)[eval(ff)];
g=@(x,y)[eval(uex)];
h=@(x,y)[eval(h1);eval(h2);eval(h3);eval(h4)];

