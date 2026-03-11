% JACOBI_POL  Script for plotting some Jacobi polynomials for n=4
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

nx=100;
x=linspace(-1,1,nx);
x=x';
n=4;

% Jacobi alpha=beta=0 (Legendre)
[pj] = jacobi_eval(x,n,0,0) ;
pl=pj(:,1);

% Jacobi alpha=beta=-0.5 (Chebyshev)
[pj] = jacobi_eval(x,n,-0.5,-0.5) ;
pc=pj(:,1)*2^(2*n)*(prod(1:n))^2/prod(1:2*n);

% Gegenbauer nu=-3/4;
nu=-3/4;
[pj] = jacobi_eval(x,n,nu-0.5,nu-0.5) ;
pg1=pj(:,1)*gamma(nu+0.5)*gamma(2*nu+n)/(gamma(nu+n+0.5)*gamma(2*nu));


% Gegenbauer nu=0;
nu=0;
[pj] = jacobi_eval(x,n,nu-0.5,nu-0.5) ;
pg2=pj(:,1)*gamma(nu+0.5)*gamma(2*nu+n)/(gamma(nu+n+0.5)*prod(1:2*nu));



figure(1);
clf
set(gca,'Fontname','Times','Fontsize',16);
plot(x,pl,'k',x,pc,'k--','Linewidth',1);
hold on
plot(x,pg1,'k-.','Linewidth',1);
plot(x,pg2,'k:','Linewidth',2);
xlabel('x','Fontname','Times','Fontsize',16);
l=legend('Legendre','Chebyshev','Gegenbauer, \nu =-3/4',...
    'Gegenbauer, \nu =0');
set(gca,'PlotBoxAspectRatio',[3 2 1],...
    'Ylimmode','manual','Ylim',[-1.2,1.2],...
   'Xtickmode','manual','Xtick',[-1,-0.5,0,0.5,1],...
    'Xgrid','on','Xminorgrid','off',...
    'Ygrid','on','Yminorgrid','off',...
    'Fontname','Times','Fontsize',16);
set(l,'Position',[0.364,0.15,0.309,0.279]);
