function [err_inf,err_h1,err_l2]=errors_1d(nx,ne,xa,xb,un,uex,uexx,param);
% ERRORS_1D Computes errors  for 1D b.v.p.
%
% Input: 
%      nx = polynomial degree in each element (the same in each element)
%             to be set only if p=4, otherwise, nx=p;
%      ne = number of elements (equally spaced)
%      xa, xb = extrema of computational domain Omega=(xa,xb)
%      un = numerical solution produced, e.g., by lap_1d, ellprecofem_1d,..
%      uex  = exact solution (uex=@(x)[uex(x)], with .*, .^, ./)
%      uexx = first derivative of exact solution (uexx=@(x)[uexx(x)], 
%             with .*, .^, ./)
%      param(1) = 1: compute errors (L^inf-norm, L2-norm, H1-norm)
%                      on the exact solution
%                 2: no  errors are computed
%      param(2) = 0: LG quadrature formulas with high precision degree are
%                      used to compute norms (exact norms)
%                 1: LGL quadrature formulas with npdx,npdy nodes are
%                      used to compute norms (discrete norms)
%                   (used only if param(1) == 1)
%      param(3) = number of nodes for high degree quadrature formula,
%                   (used only if param(2) == 0 & param(1) == 1)
%      param(4) = 0: absolute errors are computed
%                 1: relative errors are computed
%                   (used only if param(1) == 1)
%
% Output: 
%         err_inf = ||u_ex-u_N|| with respect to discrete maximum-norm
%         err_h1 = ||u_ex-u_N|| with respect to discrete H1-norm
%         err_l2 = ||u_ex-u_N|| with respect to discrete L2-norm
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

npdx=nx+1;

[x,wx]=xwlgl(npdx); dx=derlgl(x,npdx);

% nov
ldnov=npdx; nov=zeros(ldnov,ne);
[nov]=cosnov_1d(npdx,ne,nov);
noe=nov(npdx,ne);

% Uniform decomposition of  Omega in ne elements
[xx,jacx,xy,ww]=mesh_1d(xa,xb,ne,npdx,nov,x,wx);

if(param(1)==1)
% Evaluate exact solution to compute errors
nq=param(3); fdq=param(2); err_type=param(4);

u=uex(xy);

% 

err=u-un;

% ||u_ex-u_N|| with respect to discrete maximum-norm

if err_type==0
err_inf=norm(err,inf);
else
err_inf=norm(err,inf)/norm(u,inf);
end    


% ||u_ex-u_N|| with respect to H1 norm
% fdq=0 LG quadrature formula with nq nodes in each element
% fdq=1 LGL quadrature formula with npdx nodes in each element
[err_h1]=normah1_1d(fdq, nq, err_type, un, uex, uexx,...
 x, wx, dx, xx, jacx, xy,nov);

% ||u_ex-u_N|| with respect to L2 norm

[err_l2]=normal2_1d(fdq, nq, err_type, un, uex, ...
 x, wx, xx, jacx, xy,nov);
end

return
