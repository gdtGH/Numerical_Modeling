% MainFVM.m
% Main FVM solver for the 1D nonlinear Shallow Water Equations.
%
% Solves:  U_t + F(U)_x = 0
%   U = [h; q],   F(U) = [q; q^2/h + 0.5*g*h^2]
%
% on the domain [x0, xL] with REFLECTING WALL boundary conditions
% (q = 0 at both walls) enforced via ghost cells.
%
% Time integration: explicit Forward Euler with adaptive CFL-based dt.
% Spatial discretization: FVM on a uniform grid of M cells.
%
% INPUT:
%   Data    (struct)  test case parameters — see DataTest.m
%   scheme  (char)    'godunov' or 'laxwendroff'
%
% OUTPUT:
%   Solutions (struct) with fields:
%     .x   (M x 1) cell-center coordinates
%     .h   (M x 1) water height  at time T
%     .q   (M x 1) discharge     at time T
%     .u   (M x 1) velocity      at time T  (= q/h)
%     .t   (scalar) actual final time reached

function Solutions = MainFVM(Data, scheme)

fprintf('============================================================\n');
fprintf('Test: %s  |  Scheme: %s  |  M = %d  |  CFL = %.2f\n', ...
        Data.name, scheme, Data.M, Data.CFL);

%==========================================================================
% 1. GRID CONSTRUCTION
%==========================================================================
x0 = Data.domain(1);
xL = Data.domain(2);
M  = Data.M;
g  = Data.g;

dx = (xL - x0) / M;
x  = (x0 + dx/2 : dx : xL - dx/2)';   % cell centers, M x 1

% Minimum depth threshold for wet/dry treatment
h_dry = 1e-6;

%==========================================================================
% 2. INITIAL CONDITIONS  (state matrix U is 2 x M)
%==========================================================================
U = zeros(2, M);
U(1,:) = Data.h0(x)';   % row 1: h
U(2,:) = Data.q0(x)';   % row 2: q

%==========================================================================
% 3. TIME LOOP — explicit, adaptive dt
%==========================================================================
t        = 0;
iter     = 0;
max_iter = 2e6;   % hard safety cap

while t < Data.T

    iter = iter + 1;
    if iter > max_iter
        warning('MainFVM: max iterations (%d) reached. Stopping at t = %.6f.', ...
                max_iter, t);
        break
    end

    %----------------------------------------------------------------------
    % WET/DRY TREATMENT — applied first, every step.
    %
    % Root cause of infinite loop: after the FVM update, h can go slightly
    % negative near shocks (Lax-Wendroff) or dry fronts. Clamping h to 0
    % but leaving q untouched means u = q/h_clamped ~ q/eps -> enormous
    % velocity -> c_max ~ 1e16 -> dt ~ 1e-19 -> t never reaches T.
    %
    % Fix: zero out BOTH h and q on any cell where h < h_dry.
    %----------------------------------------------------------------------
    dry       = U(1,:) < h_dry;
    U(1, dry) = 0;
    U(2, dry) = 0;

    %----------------------------------------------------------------------
    % 3a. Adaptive dt from CFL condition
    %----------------------------------------------------------------------
    h_cells = max(U(1,:), eps);
    u_cells = U(2,:) ./ h_cells;
    c_cells = sqrt(g .* h_cells);
    c_max   = max(abs(u_cells) + c_cells);

    % Guard: if c_max is not finite, something exploded upstream
    if ~isfinite(c_max) || c_max < eps
        warning('MainFVM: c_max = %g at t = %.6f (iter %d). Stopping early.', ...
                c_max, t, iter);
        break
    end

    dt = Data.CFL * dx / c_max;
    dt = min(dt, Data.T - t);   % clip last step to land exactly at T

    %----------------------------------------------------------------------
    % 3b. Ghost cells for reflecting wall BC  (2 x (M+2) extended array)
    %     Left  ghost (col 1):   h mirrored,  q negated
    %     Right ghost (col M+2): h mirrored,  q negated
    %----------------------------------------------------------------------
    U_ext = [ U(1,1),  U(1,:), U(1,M);
             -U(2,1),  U(2,:), -U(2,M) ];

    %----------------------------------------------------------------------
    % 3c. Numerical fluxes at all M+1 interfaces
    %     F_all(:,k) = flux between U_ext(:,k) and U_ext(:,k+1)
    %     k=1   -> left boundary interface  (x_{1/2})
    %     k=M+1 -> right boundary interface (x_{M+1/2})
    %----------------------------------------------------------------------
    F_all = zeros(2, M+1);

    for k = 1 : M+1
        UL = U_ext(:, k);
        UR = U_ext(:, k+1);

        switch scheme
            case 'godunov'
                F_all(:, k) = GodunovFlux(UL, UR, g);
            case 'laxwendroff'
                F_all(:, k) = LaxWendroffFlux(UL, UR, g);
            otherwise
                error('MainFVM: unknown scheme "%s". Use ''godunov'' or ''laxwendroff''.', scheme);
        end
    end

    %----------------------------------------------------------------------
    % 3d. Conservative explicit update
    %     U(:,i)^{n+1} = U(:,i)^n - (dt/dx) * ( F_{i+1/2} - F_{i-1/2} )
    %     F_{i-1/2} = F_all(:, i)     (left  face of cell i)
    %     F_{i+1/2} = F_all(:, i+1)   (right face of cell i)
    %----------------------------------------------------------------------
    U = U - (dt / dx) * (F_all(:, 2:M+1) - F_all(:, 1:M));

    %----------------------------------------------------------------------
    % 3e. Advance time
    %----------------------------------------------------------------------
    t = t + dt;

end

fprintf('Reached t = %.6f  (target T = %.6f)  |  %d steps\n', ...
        t, Data.T, iter);
fprintf('============================================================\n');

%==========================================================================
% 4. PACK OUTPUT
%==========================================================================
Solutions.x = x;
Solutions.h = U(1,:)';
Solutions.q = U(2,:)';
Solutions.u = U(2,:)' ./ max(U(1,:)', eps);
Solutions.t = t;

end