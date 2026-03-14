% MainFVM.m
% Main FVM solver for the 1D nonlinear Shallow Water Equations.
%
% Solves:  U_t + F(U)_x = 0
%   U = [h; q],   F(U) = [q; q^2/h + 0.5*g*h^2]
%
% Boundary conditions: REFLECTING WALL enforced via ghost cells.
% Time integration:    explicit Forward Euler, adaptive CFL-based dt.
% Spatial scheme:      FVM on a uniform grid of M cells.
%
% LaxWendroffFlux requires 4 inputs (MUSCL-minmod): U_{i-1}, U_i, U_{i+1},
% U_{i+2}. MainFVM builds a 2-layer ghost-cell extension (M+4 columns) so
% that boundary interfaces can also access the second neighbor.
%
% INPUT:
%   Data    (struct)  test case parameters — see DataTest.m
%   scheme  (char)    'godunov' or 'laxwendroff'
%
% OUTPUT:
%   Solutions (struct):  .x  .h  .q  .u  .t

function Solutions = MainFVM(Data, scheme)

fprintf('============================================================\n');
fprintf('Test: %s  |  Scheme: %s  |  M = %d  |  CFL = %.2f\n', ...
        Data.name, scheme, Data.M, Data.CFL);

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
h_dry    = 1e-6;      % cells below this depth are treated as dry
t        = 0;
iter     = 0;
max_iter = 10000;     % ~700 steps expected for M=200,CFL=0.9,T=1

while t < Data.T * (1 - 1e-12)

    iter = iter + 1;
    if iter > max_iter
        warning('MainFVM: max iterations (%d) at t = %.6f. Stopping.', max_iter, t);
        break
    end

    %----------------------------------------------------------------------
    % 3a. Dry-cell treatment: zero BOTH h and q below h_dry
    %     (avoids u = q/h ~ q/0 -> Inf -> NaN in c_max)
    %----------------------------------------------------------------------
    dry       = U(1,:) < h_dry;
    U(1,dry)  = 0;
    U(2,dry)  = 0;

    %----------------------------------------------------------------------
    % 3b. Adaptive dt from CFL condition
    %----------------------------------------------------------------------
    h_safe = max(U(1,:), eps);
    c_max  = max(abs(U(2,:) ./ h_safe) + sqrt(g .* h_safe));

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
    % 3c. Build extended array with TWO ghost layers on each side (M+4)
    %
    %     Layout:  [ g2 | g1 | cell_1 ... cell_M | g1 | g2 ]
    %     Index:     1     2    3          M+2     M+3   M+4
    %
    %     Reflecting wall:  h_ghost = h_boundary,  q_ghost = -q_boundary
    %     Both ghost layers use the same boundary value (zero-gradient
    %     assumption for the second layer).
    %----------------------------------------------------------------------
    U_ext = zeros(2, M + 4);
    U_ext(:, 3:M+2) = U;                          % interior

    % Left ghosts (cols 2 and 1)
    U_ext(1,2) =  U(1,1);   U_ext(2,2) = -U(2,1);
    U_ext(1,1) =  U(1,1);   U_ext(2,1) = -U(2,1);

    % Right ghosts (cols M+3 and M+4)
    U_ext(1,M+3) =  U(1,M);  U_ext(2,M+3) = -U(2,M);
    U_ext(1,M+4) =  U(1,M);  U_ext(2,M+4) = -U(2,M);

    %----------------------------------------------------------------------
    % 3d. Numerical fluxes at all M+1 interfaces
    %
    %     Interface k (k = 1..M+1) lies between interior cells k-1 and k
    %     (with cells 0 and M+1 being the first ghost layers).
    %
    %     In the extended array (interior starts at col 3):
    %       U_LL = U_ext(:, k+1)   (i-1)
    %       U_L  = U_ext(:, k+2)   (i  )
    %       U_R  = U_ext(:, k+3)   (i+1)
    %       U_RR = U_ext(:, k+4)   (i+2)
    %----------------------------------------------------------------------
    F_all = zeros(2, M+1);

    for k = 1 : M+1
        ULL = U_ext(:, k);
        UL  = U_ext(:, k+1);
        UR  = U_ext(:, k+2);
        URR = U_ext(:, k+3);

        switch scheme
            case 'godunov'
                % Godunov only uses the two adjacent states
                F_all(:,k) = GodunovFlux(UL, UR, g);
            case 'laxwendroff'
                % MUSCL-minmod needs all four neighbours
                F_all(:,k) = LaxWendroffFlux(ULL, UL, UR, URR, g);
            otherwise
                error('MainFVM: unknown scheme "%s".', scheme);
        end
    end

    %----------------------------------------------------------------------
    % 3e. Conservative FVM update
    %----------------------------------------------------------------------
    U = U - (dt / dx) * (F_all(:, 2:M+1) - F_all(:, 1:M));

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
Solutions.x = x;
Solutions.h = U(1,:)';
Solutions.q = U(2,:)';
Solutions.u = U(2,:)' ./ max(U(1,:)', eps);
Solutions.t = t;

end