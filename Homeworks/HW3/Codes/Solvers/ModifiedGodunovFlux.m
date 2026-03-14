% ModifiedGodunovFlux.m
% Computes the numerical flux at interface x_{i+1/2} for the 1D SWE
% with a MODIFIED PRESSURE term (logarithmic equation of state).
%
% Modified flux:
%   F(U) = [ q ;  q^2/h + g*h*log(h/h_ref) ]
%
% Modified wave celerity (linearisation of the modified pressure):
%   c_mod = sqrt( g * max(log(h/h_ref) + 1, eps) )
%
%   NOTE: c_mod is only real for h > h_ref/e (~0.368 when h_ref=1).
%   For smaller h the argument is clamped to eps to keep the scheme
%   well-defined (the modified model is not physically meaningful in
%   that regime anyway).
%
% Riemann solver: Local Lax-Friedrichs (Rusanov)
%   F_num = 0.5*(F_L + F_R) - 0.5*s_max*(U_R - U_L)
%   s_max = max( |u_L|+c_mod_L,  |u_R|+c_mod_R )
%
% INPUT:
%   U_L   (2x1)  left  state [h_L; q_L]
%   U_R   (2x1)  right state [h_R; q_R]
%   g     (scalar) gravitational acceleration
%   h_ref (scalar) reference water height
%
% OUTPUT:
%   Fnum  (2x1)  numerical flux at the interface

function Fnum = ModifiedGodunovFlux(U_L, U_R, g, h_ref)

%--------------------------------------------------------------------------
% Extract states — guard h > 0
%--------------------------------------------------------------------------
h_L = max(U_L(1), eps);
q_L = U_L(2);
h_R = max(U_R(1), eps);
q_R = U_R(2);

u_L = q_L / h_L;
u_R = q_R / h_R;

%--------------------------------------------------------------------------
% Modified wave celerities
%   c_mod = sqrt( g*(log(h/h_ref) + 1) )
%   Clamp argument to eps to avoid sqrt of negative number.
%--------------------------------------------------------------------------
c_mod_L = sqrt(g * max(log(h_L / h_ref) + 1, eps));
c_mod_R = sqrt(g * max(log(h_R / h_ref) + 1, eps));

%--------------------------------------------------------------------------
% Modified physical fluxes
%   F = [ q ;  q^2/h + g*h*log(h/h_ref) ]
%--------------------------------------------------------------------------
F_L = [ q_L ;
        q_L^2 / h_L + g * h_L * log(h_L / h_ref) ];

F_R = [ q_R ;
        q_R^2 / h_R + g * h_R * log(h_R / h_ref) ];

%--------------------------------------------------------------------------
% Rusanov (Local Lax-Friedrichs) flux
%   s_max = max wave speed across the interface
%--------------------------------------------------------------------------
s_max = max(abs(u_L) + c_mod_L,  abs(u_R) + c_mod_R);

Fnum = 0.5 * (F_L + F_R) - 0.5 * s_max * (U_R - U_L);

end
