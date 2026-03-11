function [uex,uexx,ff,nu,beta,gam]=setfun_adr_1d
% SETFUN_ADR_1D Sets functions and coefficients for adr_1d.m
% 
% [uex,uexx,ff,nu,beta,gam]=setfun_adr_1d
%
% Output: uex =@(x)[....] function handle to the expression of exact
%               solution
%         uexx =@(x)[....] function handle to the expression of the first
%               derivative of the exact solution
%         ff =@(x)[....] function handle to the expression of function
%              at right hand side
%         nu = viscosity coefficient (constant>0)
%         beta =@(x)[....]  function handle to the expression of 
%                transport coefficient function
%         gam = coefficient of zero-order term (constant >=0)
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


syms  x
gam=1; nu=1;
Beta=cos((1+x)*pi*.25);
Uex=cos((1+x)*pi*3).*sin((0.5+x)*pi/5)+sin(pi/10); % exact solution
Uexx=diff(Uex,x); % first derivative of exact solution
Ff=gam*(Uex)-diff(nu*(Uexx)+Beta*(Uex)); % right hand side

beta=strrep(char(Beta),'*','.*');
uex=strrep(char(Uex),'*','.*');
uexx=strrep(char(Uexx),'*','.*');
ff=strrep(char(Ff),'*','.*');

beta=strrep(char(beta),'/','./');
uex=strrep(char(uex),'/','./');
uexx=strrep(char(uexx),'/','./');
ff=strrep(char(ff),'/','./');

beta=strrep(char(beta),'^','.^');
uex=strrep(char(uex),'^','.^');
uexx=strrep(char(uexx),'^','.^');
ff=strrep(char(ff),'^','.^');


beta=@(x)[eval(beta)];
uex=@(x)[eval(uex)];
uexx=@(x)[eval(uexx)];
ff=@(x)[eval(ff)];


