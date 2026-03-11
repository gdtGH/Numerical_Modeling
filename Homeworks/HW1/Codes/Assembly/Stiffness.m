function [K_loc] = Stiffness(Grad, w_1D, nln, BJ, mu_loc)
%% [K_loc] = Stiffness(Grad,w_1D,nln,BJ,mu_loc)
%==========================================================================
% Build the local stiffness matrix for the term mu*u'*v'
%==========================================================================
%    INPUT:
%          Grad        : (array real) evaluation of the derivative on q.p.
%          w_1D        : (array real) quadrature weights
%          nln         : (integer) number of local unknowns
%          BJ          : (array real) Jacobian of the map
%          mu_loc      : (array real, optional) mu evaluated at q.p.
%                        Default = 1 (homogeneous medium)
%    OUTPUT:
%          K_loc       : (array real) Local stiffness matrix

% Default: mu = 1 if not provided (backward compatible)
if nargin < 5 || isempty(mu_loc)
    mu_loc = ones(length(w_1D), 1);
end

K_loc = zeros(nln, nln);

for i = 1:nln
    for j = 1:nln
        for k = 1:length(w_1D)
            Binv = 1./BJ;
            Jdet = BJ;
            K_loc(i,j) = K_loc(i,j) + mu_loc(k) .* (Jdet.*w_1D(k)) .* ...
                          ( (Grad(k,:,i) * Binv) * (Grad(k,:,j) * Binv)' );
        end
    end
end