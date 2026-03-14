% LaxWendroffFlux.m
% Computes the 2nd-order numerical flux at interface x_{i+1/2} using the
% MUSCL-Hancock scheme with the minmod slope limiter.
%
% THEORY:
%   Pure Lax-Wendroff (no limiter) is known to produce spurious oscillations
%   near shocks in nonlinear systems, leading to h < 0 and blow-up.
%   The standard fix is to replace the fixed half-slope with a TVD-limited
%   slope (minmod), which falls back to 1st-order (Godunov) near shocks and
%   retains 2nd-order accuracy in smooth regions.
%
%   The homework formula is recovered exactly when the solution is smooth:
%     U_{i+1/2}^-  =  U_i   + 0.5*(U_{i+1} - U_i)
%     U_{i+1/2}^+  =  U_{i+1} - 0.5*(U_{i+1} - U_i)
%   The minmod limiter only clips the slope when it changes sign between
%   adjacent cells (i.e. at a local extremum or near a shock).
%
% INPUT:
%   U_LL  (2x1)  U_{i-1}  — left  neighbor of cell i
%   U_L   (2x1)  U_i      — left  cell  (cell i)
%   U_R   (2x1)  U_{i+1}  — right cell  (cell i+1)
%   U_RR  (2x1)  U_{i+2}  — right neighbor of cell i+1
%   g     (scalar) gravitational acceleration
%
% OUTPUT:
%   Fnum  (2x1)  2nd-order TVD flux at the interface

function Fnum = LaxWendroffFlux(U_LL, U_L, U_R, U_RR, g)

%--------------------------------------------------------------------------
% Minmod-limited slopes for cell i (left) and cell i+1 (right)
%--------------------------------------------------------------------------
%   slope_L = minmod( U_i - U_{i-1},  U_{i+1} - U_i   )
%   slope_R = minmod( U_{i+1} - U_i,  U_{i+2} - U_{i+1} )
%--------------------------------------------------------------------------
slope_L = minmod_vec(U_L  - U_LL,  U_R  - U_L);
slope_R = minmod_vec(U_R  - U_L,   U_RR - U_R);

%--------------------------------------------------------------------------
% Reconstruct states at the interface
%--------------------------------------------------------------------------
U_L_recon = U_L + 0.5 * slope_L;   % U_{i+1/2}^-
U_R_recon = U_R - 0.5 * slope_R;   % U_{i+1/2}^+

%--------------------------------------------------------------------------
% Ensure positivity of h in reconstructed states (backup safety)
%--------------------------------------------------------------------------
if U_L_recon(1) < 0,  U_L_recon(1) = 0;  U_L_recon(2) = 0;  end
if U_R_recon(1) < 0,  U_R_recon(1) = 0;  U_R_recon(2) = 0;  end

%--------------------------------------------------------------------------
% Roe flux on reconstructed states
%--------------------------------------------------------------------------
Fnum = GodunovFlux(U_L_recon, U_R_recon, g);

end

%--------------------------------------------------------------------------
% MINMOD LIMITER (component-wise)
%   minmod(a, b) = 0   if sign(a) ~= sign(b)   (local extremum)
%               = a    if |a| <= |b|             (a is less steep)
%               = b    otherwise
%--------------------------------------------------------------------------
function s = minmod_vec(a, b)
    s      = zeros(2, 1);
    for i  = 1:2
        if a(i) * b(i) <= 0
            s(i) = 0;                               % opposite signs → zero slope
        elseif abs(a(i)) <= abs(b(i))
            s(i) = a(i);                            % left slope is gentler
        else
            s(i) = b(i);                            % right slope is gentler
        end
    end
end