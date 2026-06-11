% GodunovFlux.m
% Numerical flux at a single interface x_{i+1/2} for the 1D nonlinear
% Shallow Water Equations, computed with the EXACT Riemann solver
% (Godunov's first-order scheme — Toro, "Shock-Capturing Methods for
% Free-Surface Shallow Flows"; course notes LES12).
%
% This is a TRUE Godunov flux: the exact self-similar solution of the
% local Riemann problem is sampled along the cell interface (xi = x/t = 0)
% and the physical flux is evaluated on that sampled state.  (The previous
% version of this file was a Roe APPROXIMATE solver — it has been replaced.)
%
% SWE in conservative form:
%   U_t + F(U)_x = 0,   U = [h; q],   q = h*u
%   F(U) = [ q ;  q^2/h + 0.5*g*h^2 ]
%   Two genuinely-nonlinear fields with eigenvalues u -+ c,  c = sqrt(g*h).
%
% EXACT SOLVER OUTLINE
%   1. Toro wave functions (depth change produced by the k-wave):
%        f_k(h,h_k) = 2*( sqrt(g*h) - sqrt(g*h_k) )           if h <= h_k  (rarefaction)
%                   = (h - h_k)*sqrt( 0.5*g*(1/h + 1/h_k) )   if h >  h_k  (shock)
%   2. Star depth h* solves   f_L(h*,h_L) + f_R(h*,h_R) + (u_R - u_L) = 0
%      (the function is monotone increasing in h => unique positive root,
%       found with fzero on a robust bracket).
%   3. Star velocity          u* = 0.5*(u_L+u_R) + 0.5*( f_R(h*,h_R) - f_L(h*,h_L) ).
%   4. Sample the self-similar solution at xi = 0.  The sign of u* selects
%      the left (xi<u*) or right (xi>u*) wave family; inside that family the
%      state at xi=0 is the data state, the star state, or a transonic
%      rarefaction value, depending on the wave speeds.
%
% DRY-STATE GUARD: if both depths are below a dry tolerance the interface is
% treated as a wall (zero flux).  A vacuum / single dry side (h<=tol) is
% handled with the exact dry-bed rarefaction sampling.  None of these arise
% in the HW3 wet-bed tests; they are only numerical safeguards.
%
% INPUT:
%   U_L  (2x1)  left  state  [h_L; q_L]
%   U_R  (2x1)  right state  [h_R; q_R]
%   g    (scalar) gravitational acceleration
%
% OUTPUT:
%   Fnum (2x1)  exact Godunov flux at the interface

function Fnum = GodunovFlux(U_L, U_R, g)

tol = 1e-6;     % dry-state threshold

%--------------------------------------------------------------------------
% Extract states (do NOT clamp h here — the dry logic below needs raw h)
%--------------------------------------------------------------------------
h_L = U_L(1);   q_L = U_L(2);
h_R = U_R(1);   q_R = U_R(2);

%--------------------------------------------------------------------------
% Consistency short-circuit: identical states => exact physical flux F(U_L).
% Guarantees GodunovFlux(U,U) == F(U) exactly (no numerical viscosity added
% where there is no wave), and avoids an unnecessary root solve.
%--------------------------------------------------------------------------
if h_L == h_R && q_L == q_R
    Fnum = phys_flux(h_L, q_L, g);
    return
end

%--------------------------------------------------------------------------
% Both sides dry -> reflecting/zero flux
%--------------------------------------------------------------------------
if h_L <= tol && h_R <= tol
    Fnum = [0; 0];
    return
end

u_L = q_L / max(h_L, tol);
u_R = q_R / max(h_R, tol);
c_L = sqrt(g * max(h_L, 0));
c_R = sqrt(g * max(h_R, 0));

%--------------------------------------------------------------------------
% DRY / VACUUM cases (exact dry-bed sampling at xi = 0)
%   - left side dry            -> single right rarefaction
%   - right side dry           -> single left  rarefaction
%   - dry generated in middle  -> two rarefactions with vacuum, when
%                                 (u_R - u_L) >= 2*(c_L + c_R)
%--------------------------------------------------------------------------
dry_left   = h_L <= tol;
dry_right  = h_R <= tol;
dry_middle = (u_R - u_L) >= 2 * (c_L + c_R);

if dry_left || dry_right || dry_middle
    [h0, u0] = sample_dry(h_L, u_L, c_L, h_R, u_R, c_R, g, dry_left, dry_right);
    Fnum = phys_flux(h0, h0 * u0, g);
    return
end

%--------------------------------------------------------------------------
% WET BED: solve for the star depth h* with fzero on a robust bracket.
% phi(h) is strictly increasing, so a sign change brackets the unique root.
%--------------------------------------------------------------------------
phi = @(h) wavefun(h, h_L, g) + wavefun(h, h_R, g) + (u_R - u_L);

h_lo = tol;                          % phi(h_lo) < 0 for a wet bed
h_hi = max(h_L, h_R);                % grow until phi(h_hi) > 0
while phi(h_hi) < 0
    h_hi = 2 * h_hi;
end

h_star = fzero(phi, [h_lo, h_hi]);
u_star = 0.5 * (u_L + u_R) + 0.5 * (wavefun(h_star, h_R, g) - wavefun(h_star, h_L, g));

%--------------------------------------------------------------------------
% Sample the self-similar solution at xi = x/t = 0.
% The contact (here the star streamline) moves at u_star: its sign decides
% whether xi=0 falls in the left (1-) or right (2-) wave family.
%--------------------------------------------------------------------------
if u_star >= 0
    %=================== LEFT FAMILY ==================================
    if h_star >= h_L
        % --- left SHOCK ---
        S_L = u_L - sqrt(0.5 * g * h_star * (h_star + h_L)) / sqrt(h_L);
        if 0 <= S_L                 % xi=0 left of the shock -> left data
            h0 = h_L;   u0 = u_L;
        else                        % xi=0 in the star region
            h0 = h_star; u0 = u_star;
        end
    else
        % --- left RAREFACTION ---
        S_HL = u_L     - c_L;                 % head  (fastest, left-most)
        S_TL = u_star  - sqrt(g * h_star);    % tail
        if 0 <= S_HL                % left of the fan -> left data
            h0 = h_L;   u0 = u_L;
        elseif 0 >= S_TL            % right of the fan -> star
            h0 = h_star; u0 = u_star;
        else                        % TRANSONIC fan: u - sqrt(g h) = 0 at xi=0
            % Riemann invariant  u_L + 2 c_L = u + 2 c  with  u = c  here
            c0 = (u_L + 2 * c_L) / 3;
            u0 = c0;
            h0 = c0^2 / g;
        end
    end
else
    %=================== RIGHT FAMILY =================================
    if h_star >= h_R
        % --- right SHOCK ---
        S_R = u_R + sqrt(0.5 * g * h_star * (h_star + h_R)) / sqrt(h_R);
        if 0 >= S_R                 % xi=0 right of the shock -> right data
            h0 = h_R;   u0 = u_R;
        else                        % star region
            h0 = h_star; u0 = u_star;
        end
    else
        % --- right RAREFACTION ---
        S_HR = u_R     + c_R;                 % head  (right-most)
        S_TR = u_star  + sqrt(g * h_star);    % tail
        if 0 >= S_HR                % right of the fan -> right data
            h0 = h_R;   u0 = u_R;
        elseif 0 <= S_TR            % left of the fan -> star
            h0 = h_star; u0 = u_star;
        else                        % TRANSONIC fan: u + sqrt(g h) = 0 at xi=0
            % Riemann invariant  u_R - 2 c_R = u - 2 c  with  u = -c  here
            c0 = (-u_R + 2 * c_R) / 3;
            u0 = -c0;
            h0 = c0^2 / g;
        end
    end
end

%--------------------------------------------------------------------------
% Godunov flux = physical flux of the sampled interface state
%--------------------------------------------------------------------------
Fnum = phys_flux(h0, h0 * u0, g);

end

%==========================================================================
% LOCAL FUNCTIONS
%==========================================================================

%--- physical flux F(U) = [q; q^2/h + 0.5 g h^2] -------------------------
function F = phys_flux(h, q, g)
    if h <= 0
        F = [0; 0];
    else
        F = [q; q^2 / h + 0.5 * g * h^2];
    end
end

%--- Toro wave function f_k(h, h_k) --------------------------------------
function f = wavefun(h, h_k, g)
    if h <= h_k
        f = 2 * (sqrt(g * h) - sqrt(g * h_k));            % rarefaction
    else
        f = (h - h_k) * sqrt(0.5 * g * (1 / h + 1 / h_k)); % shock
    end
end

%--- exact dry-bed sampling at xi = 0 ------------------------------------
%  Returns the (h,u) state on the interface for the three dry-bed
%  configurations.  Uses the rarefaction-fan formulas:
%    left  fan: u = (u_L + 2 c_L + 2*xi)/3 ,  c = ( u_L + 2 c_L - xi)/3
%    right fan: u = (u_R - 2 c_R + 2*xi)/3 ,  c = (-u_R + 2 c_R + xi)/3
%  evaluated at xi = 0.
function [h0, u0] = sample_dry(h_L, u_L, c_L, h_R, u_R, c_R, g, dry_left, dry_right)

    if dry_left && ~dry_right
        % ---- left dry: single right rarefaction into the dry bed ----
        S_star = u_R - 2 * c_R;     % wet/dry front
        if 0 >= u_R + c_R           % right of head -> right data
            h0 = h_R; u0 = u_R;
        elseif 0 <= S_star          % in the dry region
            h0 = 0;   u0 = 0;
        else                        % inside the right fan
            c0 = (-u_R + 2 * c_R) / 3;
            u0 = -c0;
            h0 = c0^2 / g;
        end

    elseif dry_right && ~dry_left
        % ---- right dry: single left rarefaction into the dry bed ----
        S_star = u_L + 2 * c_L;     % wet/dry front
        if 0 <= u_L - c_L           % left of head -> left data
            h0 = h_L; u0 = u_L;
        elseif 0 >= S_star          % in the dry region
            h0 = 0;   u0 = 0;
        else                        % inside the left fan
            c0 = (u_L + 2 * c_L) / 3;
            u0 = c0;
            h0 = c0^2 / g;
        end

    else
        % ---- vacuum generated between two wet states (two rarefactions) ----
        S_starL = u_L + 2 * c_L;
        S_starR = u_R - 2 * c_R;
        if 0 <= u_L - c_L           % left data
            h0 = h_L; u0 = u_L;
        elseif 0 < S_starL          % inside left fan
            c0 = (u_L + 2 * c_L) / 3;
            u0 = c0;  h0 = c0^2 / g;
        elseif 0 < S_starR          % vacuum
            h0 = 0;   u0 = 0;
        elseif 0 < u_R + c_R        % inside right fan
            c0 = (-u_R + 2 * c_R) / 3;
            u0 = -c0; h0 = c0^2 / g;
        else                        % right data
            h0 = h_R; u0 = u_R;
        end
    end
end
