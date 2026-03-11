function [err_l2]=normal2_1d(fdq, nq, errtype, u, uex, ...
 x, wx, xx, jacx, xy,nov)
% NORMAL2_1D   Computes L2-norm in 1D
%
%  [err_h1]=normah1_1d(fdq, nq, errtype, u, uex, ...
%            x, wx, xx, jacx, xy,nov)
%
% Input : fdq = 0  uses Legendre Gauss quadrature formulas with nq nodes
%               in each element (exactness degree = 2*nq+1)
%             = 1   uses Legendre Gauss Lobatto quadrature formulas on the 
%               nodes in which the discrete function is knonwn.
%         nq = nodes (in each elemetn) for GL quadrature formulas. Not used if
%              fdq=1
%         errtype = 0 for absolute error ||u-u_ex||
%                   1 for relative error ||u-u_ex||/||u_ex||
%         u = numerical solution, column array of length noe
%         uex  = exact solution (uex=@(x)[uex(x)], with .*, .^, ./)
%         x = column array  with LGL nodes in [-1,1]
%         wx= column array  with LGL weigths in [-1,1]
%         xx = 2-indexes array of size (2,ne): xx(1:2,ie)=[xa_ie;xb_ie]
%         jacx = column array of length = ne, containing 
%                jacobians of the maps F_ie:[-1,1]---->[xa_ie,xb_ie]
%         xy = mesh, column array of length noe
%         nov = local to global map. 2-indexes array, size(nov)=[nov,ne]
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


npdx=length(x); [ldnov,ne]=size(nov); noe=nov(npdx,ne);
num=0; den=0;
err_l2=0;

if fdq==0

% Legendre Gauss nodes and weigths

[xg,wxg] = xwlg(nq,-1,1);

% Evaluation of Lagrange basis polynomials at quadrature nodes xg
%
[phix]= intlag_lgl (x,xg);

% Loop on spectral elements
for ie=1:ne
u_loc=u(nov(1:npdx,ie));
[norma_loc1,norma_loc2]=normal2_ex_loc(errtype,npdx,uex,...
 u_loc,x,jacx(ie),xx(1:2,ie),xg,wxg,nq,phix);
num=num+norma_loc1;
den=den+norma_loc2;
end


elseif fdq==1
for ie=1:ne
u_loc=u(nov(1:npdx,ie));
xy_loc=xy(nov(1:npdx,ie));
[norma_loc1,norma_loc2]=normal2_loc(errtype,npdx,uex,...
 u_loc,xy_loc, wx,jacx(ie));
num=num+norma_loc1;
den=den+norma_loc2;
end
end

if errtype==0
    err_l2=sqrt(num);
elseif errtype==1
    if abs(den)>1.d-14; err_l2=sqrt(num/den); end
end

return

function [norma_loc1,norma_loc2]=normal2_ex_loc(errtype,npdx,uex,...
 u, x,jacx,xx,xg,wxg,nq,phix);

% High degree Legendre Gaussian  formulas to compute L2-norm error

% mapping quadrature nodes on element ie

xxg=xg*jacx+(xx(2)+xx(1))*.5;

% evaluation of exact solution at quadrature nodes.

U=uex(xxg);


% evaluate numerical solution and its derivative at quadrature nodes.

u_i=phix*u;
%  figure(1)
%  mesh(xxg,yyg,reshape(U,npx,npy)); 
%  figure(2) 
%  mesh(xxg,yyg,uloc_i)

% compute the sum

norma_loc1=sum((U-u_i).^2.*wxg)*jacx;

if errtype==0
    norma_loc2=0;
else
    norma_loc2=sum(U.^2.*wxg)*jacx;
end

return

function [norma_loc1,norma_loc2]=normal2_loc(errtype,npdx,uex,uloc,...
 xy_loc, wx,jacx); 

% LGL quadrature formulas on npdx nodes to compute L2-norm error

% evaluation of exact solution at quadrature nodes.

U=uex(xy_loc);


% compute the sum

norma_loc1=sum((U-uloc).^2.*wx)*jacx;
if errtype==0
    norma_loc2=0;
else
    norma_loc2=sum(U.^2.*wx)*jacx;
end

return

