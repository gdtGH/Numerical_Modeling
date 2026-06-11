% LaxWendroffFlux.m
% Second-order numerical flux at interface x_{i+1/2} using a PURE
% Lax-Wendroff (unlimited) linear reconstruction — exactly as required by
% the HW3 assignment.  NO slope limiter is applied.
%
% THEORY:
%   The interface states are reconstructed with the fixed (unlimited) slope
%   (U_{i+1} - U_i):
%       U_{i+1/2}^-  =  U_i     + 0.5*(U_{i+1} - U_i)
%       U_{i+1/2}^+  =  U_{i+1}  - 0.5*(U_{i+1} - U_i)
%   (both collapse to the cell-interface average 0.5*(U_i + U_{i+1})), and the
%   flux is the exact Godunov flux evaluated on those reconstructed states.
%   With Forward Euler in time this is the classic second-order centred
%   scheme: it is dispersive and is KNOWN to produce spurious oscillations
%   and h<0 (blow-up) near shocks, because no TVD/limiter mechanism clips the
%   unphysical over/undershoots.  Demonstrating that failure is the point of
%   the assignment — see the blow-up guard in MainFVM.m.
%
%   (The previous version of this file used a MUSCL reconstruction with the
%    minmod limiter; that TVD machinery has been removed on purpose.)
%
% SIGNATURE NOTE:
%   The 4-state signature (U_LL, U_L, U_R, U_RR) is kept UNCHANGED so the call
%   in MainFVM.m needs no edit.  The pure Lax-Wendroff stencil only uses the
%   two interface-adjacent cells U_L = U_i and U_R = U_{i+1}; U_LL and U_RR
%   are intentionally ignored (no wider stencil / no limiter).
%
% INPUT:
%   U_LL  (2x1)  U_{i-1}  — UNUSED (kept for signature compatibility)
%   U_L   (2x1)  U_i      — left  cell  (cell i)
%   U_R   (2x1)  U_{i+1}  — right cell  (cell i+1)
%   U_RR  (2x1)  U_{i+2}  — UNUSED (kept for signature compatibility)
%   g     (scalar) gravitational acceleration
%
% OUTPUT:
%   Fnum  (2x1)  pure (unlimited) Lax-Wendroff flux at the interface

function Fnum = LaxWendroffFlux(U_LL, U_L, U_R, U_RR, g) %#ok<INUSL>

%--------------------------------------------------------------------------
% Pure Lax-Wendroff reconstruction (fixed, unlimited slope = U_R - U_L)
%--------------------------------------------------------------------------
dU        = U_R - U_L;
U_L_recon = U_L + 0.5 * dU;     % = 0.5*(U_L + U_R)
U_R_recon = U_R - 0.5 * dU;     % = 0.5*(U_L + U_R)

%--------------------------------------------------------------------------
% MINIMAL POSITIVITY SAFEGUARD (NOT a slope limiter / NOT a TVD fix):
% if a reconstructed depth becomes non-positive, fall back to the
% (un-reconstructed) cell state for that interface side.  This only prevents
% sqrt of a negative number / division by zero inside the flux — it does NOT
% suppress the oscillations the pure scheme is meant to exhibit.
%--------------------------------------------------------------------------
if U_L_recon(1) <= 0,  U_L_recon = U_L;  end
if U_R_recon(1) <= 0,  U_R_recon = U_R;  end

%--------------------------------------------------------------------------
% Exact Godunov flux on the reconstructed states
%--------------------------------------------------------------------------
Fnum = GodunovFlux(U_L_recon, U_R_recon, g);

end
