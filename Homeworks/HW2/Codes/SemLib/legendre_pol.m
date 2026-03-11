% LEGENDRE_POL  plots Legendre polynomials for n=0,...,4
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

nx=100;
x=linspace(-1,1,nx);
x=x';

% Legendre
n=4;
pn=pnleg_all(x,n);

figure(1);
clf
set(gca,'Fontname','Times','Fontsize',16);
plot(x,pn(:,1),'k--','Linewidth',1);
hold on
plot(x,pn(:,2),'k',x,pn(:,3),'k--','Linewidth',1);
plot(x,pn(:,4),'k-.','Linewidth',1);
plot(x,pn(:,5),'k:','Linewidth',2);
xlabel('x','Fontname','Times','Fontsize',16);
set(gca,'PlotBoxAspectRatio',[3 2 1],...
    'Ylimmode','manual','Ylim',[-1.2,1.2],...
   'Xtickmode','manual','Xtick',[-1,-0.5,0,0.5,1],...
    'Xgrid','on','Xminorgrid','off',...
    'Ygrid','on','Yminorgrid','off',...
    'Fontname','Times','Fontsize',16);
