function [uex,uexx,ff,nu,beta,gam]=setfun_ell_1d
%  SETFUN_ELL_1D  Sets functions and coefficients for ellprecofem_1d and ell_1d
%
%    [uex,uexx,ff,nu,beta,gam]=setfun_ell_1d
%
% Output: uex =@(x)[....] function handle to the expression of exact
%               solution
%         uexx =@(x)[....] function handle to the expression of the first 
%               derivative of the exact solution
%         ff =@(x)[....] function handle to the expression of function 
%              at right hand side
%         nu = viscosity coefficient (constant>0)
%         beta = coefficient of first-order term (constant)
%         gam = coefficient of zero-order term (constant >=0)
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

syms  x
nu=1;gam=0; beta=0;
Uex=cos((1+x)*pi*3).*sin((0.5+x)*pi/5)+sin(pi/10); % exact solution
Uexx=diff(Uex,x); % first derivative of exact solution
Ff=(gam*(Uex)-diff(nu*(Uexx))+beta*Uexx); % right hand side

% Ff=1+0*x;
% Uex=0+0*x;
uex=vectorize(char(Uex));
uexx=vectorize(char(Uexx));
ff=vectorize(char(Ff));

uex=@(x)[eval(uex)];
uexx=@(x)[eval(uexx)];
ff=@(x)[eval(ff)];

