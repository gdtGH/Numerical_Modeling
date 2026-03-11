function [A,M]=matrices_1d(xa,xb,cb,ne,nx);
% MATRICES_1D   Assembles SEM stiffness and mass matrices for 1D b.v.p.
%
%     M_{ij}=(phi_j,phi_i)_N
%     A_{ij}=(grad phi_j, grad phi_i)_N
%    
% SEM  Numerical Integration with LGL quadrature formulas.
%
%  [A,M]=matrices_1d(xa,xb,cb,ne,nx);
%
% Input: xa, xb = extrema of computational domain Omega=(xa,xb)
%      cb = string containing boundary conditions, for example
%           cb='dn' imposes Dirichlet in xa and Neumann in xb 
%      ne = number of elements (equally spaced)
%      nx = polynomial degree in each element (the same in each element)
%             to be set only if p=4, otherwise, nx=p;
%
% Output: A = stiffness matrix
%         M = Mass matrix
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


%
npdx=nx+1; 

[x,wx]=xwlgl(npdx); dx=derlgl(x,npdx);

% nov
ldnov=npdx; nov=zeros(ldnov,ne);
[nov]=cosnov1d(npdx,ne,nov);
noe=nov(npdx,ne);

% Uniform decomposition of  Omega in ne elements

[xx,jacx,xy,ww]=mesh1d(xa,xb,ne,npdx,nov,x,wx);

% ww contains the diagonal of the mass matrix

% Stiffness Matrix assembling
A=stiff_1d_se(npdx,ne,nov,wx,dx,jacx); 
M=spdiags(ww,0,noe,noe);

% Setting boundary conditions on the matrix

if cb(1)=='d' & cb(2)=='d'; 
lint=(2:noe-1); 
elseif cb(1)=='n' & cb(2)=='d';
lint=(1:noe-1);
elseif cb(1)=='d' & cb(2)=='n';
lint=(2:noe);
elseif cb(1)=='n' & cb(2)=='n';
lint=(1:noe);
end

A=A(lint,lint); M=M(lint,lint);

return
