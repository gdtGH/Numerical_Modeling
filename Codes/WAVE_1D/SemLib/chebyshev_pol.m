% CHEBYSHEV_POL  plots Chebyshev polynomials for n=0,...,4
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

nx=100;
x=linspace(-1,1,nx);
x=x';

pc=zeros(nx,5);
% Chebyshev
for n=1:4
[pj] = jacobi_eval(x,n,-0.5,-0.5) ;
pc(:,n)=pj(:,1)*2^(2*n)*(prod(1:n))^2/prod(1:2*n);
end

figure(1);
clf
set(gca,'Fontname','Times','Fontsize',16);
plot(x,pc(:,1),'k',x,pc(:,2),'k--','Linewidth',1);
hold on
plot(x,pc(:,3),'k-.','Linewidth',1);
plot(x,pc(:,4),'k:','Linewidth',2);
plot(x,ones(nx,1),'k--','Linewidth',1);
xlabel('x','Fontname','Times','Fontsize',16);
set(gca,'PlotBoxAspectRatio',[3 2 1],...
    'Ylimmode','manual','Ylim',[-1.2,1.2],...
   'Xtickmode','manual','Xtick',[-1,-0.5,0,0.5,1],...
    'Xgrid','on','Xminorgrid','off',...
    'Ygrid','on','Yminorgrid','off',...
    'Fontname','Times','Fontsize',16);
