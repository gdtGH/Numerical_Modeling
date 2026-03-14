% LaxWendroffFlux.m
% Computes the 2nd-order numerical flux at a single interface x_{i+1/2}
% using piecewise linear reconstruction (Lax-Wendroff approach) followed
% by the Roe approximate Riemann solver.
%
% Reconstruction formulas (as given in the HW3 assignment):
%   U_{i+1/2}^-  =  U_i   + 0.5 * (U_{i+1} - U_i)
%   U_{i+1/2}^+  =  U_{i+1} - 0.5 * (U_{i+1} - U_i)
%
% The reconstructed states are then passed to GodunovFlux (Roe solver).
%
% INPUT:
%   U_L  (2x1)  left  cell average  U(:,i)    = [h_i;   q_i  ]
%   U_R  (2x1)  right cell average  U(:,i+1)  = [h_i+1; q_i+1]
%   g    (scalar) gravitational acceleration
%
% OUTPUT:
%   Fnum (2x1)  2nd-order numerical flux at the interface

function Fnum = LaxWendroffFlux(U_L, U_R, g)

%--------------------------------------------------------------------------
% Piecewise linear reconstruction at the interface
%--------------------------------------------------------------------------
dU = U_R - U_L;

U_L_recon = U_L + 0.5 * dU;   % U_{i+1/2}^-  (left  side of interface)
U_R_recon = U_R - 0.5 * dU;   % U_{i+1/2}^+  (right side of interface)

% Positivity fix: the linear reconstruction can produce h < 0 near shocks.
% Clamp h to zero and zero out q consistently to avoid NaN in the Roe solver.
if U_L_recon(1) < 0
    U_L_recon(1) = 0;
    U_L_recon(2) = 0;
end
if U_R_recon(1) < 0
    U_R_recon(1) = 0;
    U_R_recon(2) = 0;
end

%--------------------------------------------------------------------------
% Compute numerical flux using the Roe solver on reconstructed states
%--------------------------------------------------------------------------
Fnum = GodunovFlux(U_L_recon, U_R_recon, g);

end