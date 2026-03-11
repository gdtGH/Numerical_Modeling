function [M_loc] = Mass(dphiq, w_1D, nln, BJ, ro_loc)
%% [M_loc] = Mass(dphiq,w_1D,nln,BJ,ro_loc)
%==========================================================================
% Build the local mass matrix for the term ro(x)*u*v
%==========================================================================
%    INPUT:
%          dphiq       : (array) basis functions at quadrature nodes
%          w_1D        : (array) quadrature weights
%          nln         : (integer) number of local unknowns
%          BJ          : Jacobian of the map
%          ro_loc      : (array, optional) rho evaluated at q.p.
%                        Default = 1 (homogeneous medium)
%    OUTPUT:
%          M_loc       : Local mass matrix

if nargin < 5 || isempty(ro_loc)
    ro_loc = ones(length(w_1D), 1);
end

M_loc = zeros(nln, nln);

for i = 1:nln
    for j = 1:nln
        for k = 1:length(w_1D)
            Jdet = BJ;
            M_loc(i,j) = M_loc(i,j) + ro_loc(k) .* (Jdet.*w_1D(k)) .* ...
                          dphiq(1,k,i) .* dphiq(1,k,j);
        end
    end
end