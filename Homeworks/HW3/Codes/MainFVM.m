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
x  = (x0 + dx/2 : dx : xL - dx/2)';   % cell centers, M×1

%==========================================================================
% 2. INITIAL CONDITIONS  (state matrix U is 2×M)
%==========================================================================
U = zeros(2, M);
U(1,:) = Data.h0(x)';   % row 1: h
U(2,:) = Data.q0(x)';   % row 2: q

%==========================================================================
% 3. TIME LOOP — explicit, adaptive dt
%==========================================================================
t      = 0;
iter   = 0;
max_iter = 1e6;   % safety cap — prevents silent infinite loops

while t < Data.T

    iter = iter + 1;
    if iter > max_iter
        warning('MainFVM: max iterations (%d) reached. Stopping early at t = %.6f.', max_iter, t);
        break
    end

    %----------------------------------------------------------------------
    % 3a. Compute adaptive dt from CFL condition
    %----------------------------------------------------------------------

    % FIX: positivity clamp — Lax-Wendroff reconstruction can produce h < 0
    % near shocks. Clamp here before computing wave speeds.
    U(1,:) = max(U(1,:), 0);

    h_cells = max(U(1,:), eps);              % guard h > 0 for divisions
    u_cells = U(2,:) ./ h_cells;            % velocity at all cells
    c_cells = sqrt(g .* h_cells);           % celerity at all cells

    c_max = max(abs(u_cells) + c_cells);    % maximum wave speed

    % FIX: guard against NaN/Inf in c_max (triggered by NaN in U)
    if ~isfinite(c_max) || c_max < eps
        warning('MainFVM: c_max = %g at t = %.6f. Stopping early.', c_max, t);
        break
    end

    dt = Data.CFL * dx / c_max;

    % Clip last step so we land exactly at T
    dt = min(dt, Data.T - t);

    %----------------------------------------------------------------------
    % 3b. Extend U with ghost cells for reflecting wall BC
    %     Ghost cell LEFT  (index 0): h_mirror, q_flipped
    %     Ghost cell RIGHT (index M+1): h_mirror, q_flipped
    %----------------------------------------------------------------------
    U_ext = zeros(2, M + 2);
    U_ext(:, 2:M+1) = U;                         % interior cells
    U_ext(1, 1)     =  U(1, 1);                  % left ghost: mirror h
    U_ext(2, 1)     = -U(2, 1);                  % left ghost: flip q
    U_ext(1, M+2)   =  U(1, M);                  % right ghost: mirror h
    U_ext(2, M+2)   = -U(2, M);                  % right ghost: flip q

    %----------------------------------------------------------------------
    % 3c. Compute numerical fluxes at all M+1 interfaces
    %     Interface i (in 1-based indexing) lies between cell i and i+1
    %     in the extended array (indices i and i+1 of U_ext, offset by 1)
    %
    %     F_all(:, i) = flux at interface x_{i-1/2}   (left face of cell i)
    %     We compute F at interfaces 0..M, stored as F_all(:,1..M+1)
    %----------------------------------------------------------------------
    F_all = zeros(2, M + 1);

    for i = 1 : M+1
        % U_ext column i   corresponds to left  state (cell i-1 extended)
        % U_ext column i+1 corresponds to right state (cell i   extended)
        UL = U_ext(:, i);
        UR = U_ext(:, i+1);

        if strcmp(scheme, 'godunov')
            F_all(:, i) = GodunovFlux(UL, UR, g);
        elseif strcmp(scheme, 'laxwendroff')
            F_all(:, i) = LaxWendroffFlux(UL, UR, g);
        else
            error('MainFVM: unknown scheme "%s". Use ''godunov'' or ''laxwendroff''.', scheme);
        end
    end

    %----------------------------------------------------------------------
    % 3d. Conservative FVM update
    %     U(:,i)^{n+1} = U(:,i)^n - (dt/dx) * ( F_{i+1/2} - F_{i-1/2} )
    %     F_{i-1/2} = F_all(:, i)
    %     F_{i+1/2} = F_all(:, i+1)
    %----------------------------------------------------------------------
    U = U - (dt / dx) * (F_all(:, 2:M+1) - F_all(:, 1:M));

    % FIX: clamp h >= 0 after update (2nd line of defense against NaN cascade)
    U(1,:) = max(U(1,:), 0);

    %----------------------------------------------------------------------
    % 3e. Advance time
    %----------------------------------------------------------------------
    t = t + dt;

end

fprintf('Reached t = %.6f  (target T = %.6f)  |  %d time steps\n', t, Data.T, iter);
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