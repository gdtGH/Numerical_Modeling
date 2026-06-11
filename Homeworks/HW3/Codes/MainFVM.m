% MainFVM.m
% Main FVM solver for the 1D nonlinear Shallow Water Equations.
%
% Solves:  U_t + F(U)_x = 0
%   Standard flux:  F(U) = [q;  q^2/h + 0.5*g*h^2]
%   Modified flux:  F(U) = [q;  q^2/h + g*h*log(h/h_ref)]
%
% Boundary conditions: REFLECTING WALL enforced via ghost cells.
% Time integration:    explicit Forward Euler, adaptive CFL-based dt.
%
% INPUT:
%   Data    (struct)  test case parameters — see DataTest.m
%   scheme  (char)    'godunov' | 'laxwendroff' | 'modified_godunov'
%   h_ref   (scalar, optional)  reference depth for modified flux (default: [])
%                               Required when scheme = 'modified_godunov'.
%
% OUTPUT:
%   Solutions (struct):  .x  .h  .q  .u  .t

function Solutions = MainFVM(Data, scheme, h_ref)

% Default: standard flux (no modified pressure)
if nargin < 3
    h_ref = [];
end

fprintf('============================================================\n');
if isempty(h_ref)
    fprintf('Test: %s  |  Scheme: %s  |  M = %d  |  CFL = %.2f\n', ...
            Data.name, scheme, Data.M, Data.CFL);
else
    fprintf('Test: %s  |  Scheme: %s  |  h_ref = %.3f  |  M = %d  |  CFL = %.2f\n', ...
            Data.name, scheme, h_ref, Data.M, Data.CFL);
end

%==========================================================================
% 1. GRID
%==========================================================================
x0 = Data.domain(1);
xL = Data.domain(2);
M  = Data.M;
g  = Data.g;

dx = (xL - x0) / M;
x  = (x0 + dx/2 : dx : xL - dx/2)';   % cell centres, M x 1

%==========================================================================
% 2. INITIAL CONDITIONS
%==========================================================================
U = zeros(2, M);
U(1,:) = Data.h0(x)';
U(2,:) = Data.q0(x)';

%==========================================================================
% 3. TIME LOOP
%==========================================================================
h_dry    = 1e-6;
t        = 0;
iter     = 0;
max_iter = 50000;     % standard: ~1500 steps; modified flux near h->0: up to ~20000
diverged = false;     % set true if a NaN/Inf or h<0 blow-up is detected
t_blowup = NaN;       % time at which the (rejected) blow-up step was attempted

while t < Data.T * (1 - 1e-12)

    iter = iter + 1;
    if iter > max_iter
        warning('MainFVM: max iterations (%d) at t = %.6f. Stopping.', max_iter, t);
        break
    end

    %----------------------------------------------------------------------
    % Snapshot of the last VALID state, used to roll back on a blow-up
    % (see the blow-up guard after the update below).
    %----------------------------------------------------------------------
    U_prev = U;
    t_prev = t;

    %----------------------------------------------------------------------
    % 3a. Dry-cell treatment
    %----------------------------------------------------------------------
    dry      = U(1,:) < h_dry;
    U(1,dry) = 0;
    U(2,dry) = 0;

    %----------------------------------------------------------------------
    % 3b. Adaptive dt — CFL condition
    %     Standard scheme:  c = sqrt(g*h)
    %     Modified scheme:  c_mod = sqrt(g*(log(h/h_ref)+1))  [clamped]
    %----------------------------------------------------------------------
    h_safe = max(U(1,:), eps);
    u_abs  = abs(U(2,:) ./ h_safe);

    if strcmp(scheme, 'modified_godunov')
        % Modified wave celerity for CFL
        c_wave = sqrt(g .* max(log(h_safe / h_ref) + 1, eps));
    else
        % Standard wave celerity
        c_wave = sqrt(g .* h_safe);
    end

    c_max = max(u_abs + c_wave);

    if ~isfinite(c_max) || c_max < eps
        warning('MainFVM: c_max = %g at t = %.6f (iter %d). Stopping.', ...
                c_max, t, iter);
        break
    end

    dt = Data.CFL * dx / c_max;
    dt = min(dt, Data.T - t);

    if dt < 1e-14
        t = Data.T;
        break
    end

    %----------------------------------------------------------------------
    % 3c. Extended array — TWO ghost layers each side (M+4 columns)
    %     Layout: [ g2 | g1 | cell_1 ... cell_M | g1 | g2 ]
    %     Reflecting wall: h mirrored, q negated
    %----------------------------------------------------------------------
    U_ext = zeros(2, M + 4);
    U_ext(:, 3:M+2) = U;
    U_ext(1,2) =  U(1,1);   U_ext(2,2) = -U(2,1);
    U_ext(1,1) =  U(1,1);   U_ext(2,1) = -U(2,1);
    U_ext(1,M+3) =  U(1,M); U_ext(2,M+3) = -U(2,M);
    U_ext(1,M+4) =  U(1,M); U_ext(2,M+4) = -U(2,M);

    %----------------------------------------------------------------------
    % 3d. Numerical fluxes at all M+1 interfaces
    %----------------------------------------------------------------------
    F_all = zeros(2, M+1);

    for k = 1 : M+1
        ULL = U_ext(:, k);
        UL  = U_ext(:, k+1);
        UR  = U_ext(:, k+2);
        URR = U_ext(:, k+3);

        switch scheme
            case 'godunov'
                F_all(:,k) = GodunovFlux(UL, UR, g);
            case 'laxwendroff'
                F_all(:,k) = LaxWendroffFlux(ULL, UL, UR, URR, g);
            case 'modified_godunov'
                F_all(:,k) = ModifiedGodunovFlux(UL, UR, g, h_ref);
            otherwise
                error('MainFVM: unknown scheme "%s".', scheme);
        end
    end

    %----------------------------------------------------------------------
    % 3e. Conservative FVM update
    %----------------------------------------------------------------------
    U = U - (dt / dx) * (F_all(:, 2:M+1) - F_all(:, 1:M));

    %----------------------------------------------------------------------
    % 3e-bis. BLOW-UP GUARD
    %     The pure Lax-Wendroff scheme (no limiter) is expected to develop
    %     unbounded oscillations / negative depths near shocks.  Rather than
    %     spinning to max_iter on NaN/Inf garbage, detect the blow-up, roll
    %     back to the last valid state, and stop.  This DOCUMENTS the failure
    %     (it does not mask it) so the report can plot the last stable time.
    %----------------------------------------------------------------------
    if any(~isfinite(U(:))) || min(U(1,:)) < -1e-6
        diverged = true;
        t_blowup = t + dt;          % time of the rejected step
        U = U_prev;                 % restore last valid state
        t = t_prev;                 % and its time
        if strcmp(scheme, 'laxwendroff')
            warning(['Lax-Wendroff puro: oscillazioni non limitate / ' ...
                     'blow-up vicino allo shock (atteso senza limitatore). ' ...
                     'Blow-up a t = %.6f; ultimo stato valido a t = %.6f.'], ...
                     t_blowup, t_prev);
        else
            warning(['MainFVM: blow-up (NaN/Inf o h<0) con scheme "%s" ' ...
                     'a t = %.6f; ultimo stato valido a t = %.6f.'], ...
                     scheme, t_blowup, t_prev);
        end
        break
    end

    %----------------------------------------------------------------------
    % 3f. Advance time
    %----------------------------------------------------------------------
    t = t + dt;

end

fprintf('Reached t = %.8f  (target T = %.6f)  |  %d steps\n', t, Data.T, iter);
fprintf('============================================================\n');

%==========================================================================
% 4. OUTPUT
%==========================================================================
Solutions.x        = x;
Solutions.h        = U(1,:)';
Solutions.q        = U(2,:)';
Solutions.u        = U(2,:)' ./ max(U(1,:)', eps);
Solutions.t        = t;           % last VALID time reached
Solutions.diverged = diverged;    % true if the run stopped on a blow-up
Solutions.t_blowup = t_blowup;    % time of the rejected blow-up step (NaN if none)

end