function [Errors, Solutions, Femregion, Data] = MainSWE(Data, nEl)
%% [Errors, Solutions, Femregion, Data] = MainSWE(Data, nEl)
%==========================================================================
%  Spectral Element solver for the 1-D Linearised Shallow Water Equations
%
%    eta_t + q_x        = 0        (mass conservation)
%    q_t   + g*H*eta_x  = - gamma*q   (momentum, optional friction)
%
%  on (0,L) x (0,T] with PERIODIC or WALL boundary conditions.
%
%  Spatial discretisation  : SEM (Legendre-Gauss-Lobatto nodes)
%  Time integration        : theta-method (theta=0.5 => Crank-Nicolson)
%
%  INPUT
%    Data  : struct built by DataTest.m
%    nEl   : number of spectral elements
%
%  OUTPUT
%    Errors    : struct  .L2   — L2 norm of eta error at t=T
%    Solutions : struct  .uh, .u_ex, .q, .x
%    Femregion : FE-region struct (from CreateFemregion)
%    Data      : (possibly enriched) Data struct
%==========================================================================

fprintf('============================================================\n');
fprintf(' MainSWE — test: %s   Ne = %d   p = %d   dt = %.2e\n', ...
        Data.name, nEl, Data.p, Data.dt);

% =========================================================================
% 1.  MESH AND FEM REGION
% =========================================================================
[Region]    = CreateMesh(Data, nEl);
[Femregion] = CreateFemregion(Data, Region);

ndof         = Femregion.ndof;
nln          = Femregion.nln;          % p+1 local nodes per element
ne           = Femregion.ne;
connectivity = Femregion.connectivity;
coord        = Femregion.coord;        % (ndof,1) physical coordinates

% =========================================================================
% 2.  ASSEMBLE GLOBAL MASS MATRIX  M  AND DERIVATIVE MATRIX  B
%
%  M_{ij} = int_0^L phi_i * phi_j dx       (diagonal, LGL mass)
%  B_{ij} = int_0^L phi_i * d(phi_j)/dx dx (transport operator)
%
%  At LGL nodes:  phi_i(xi_k) = delta_{ik}  =>
%    M_loc(i,j) = BJ * w_i * delta_{ij}
%    B_loc(i,j) = w_i * D_hat(i,j)           [BJ cancels exactly]
%  where D_hat is the (nln x nln) spectral derivative matrix on [-1,1].
% =========================================================================
fprintf('Assembling M and B ... \n');

[nodes_1D, w_1D] = xwlgl(nln);
[~, GradPhi]     = basis_and_der_at_lgl(nodes_1D, nln);
% GradPhi(k,1,j) = D_hat(k,j) = d phi_j / d xi  at  xi_k

% Extract the (nln x nln) reference-element derivative matrix
D_hat = squeeze(GradPhi(:, 1, :));   % (nln x nln)
W_diag = diag(w_1D);                 % diagonal weight matrix

% Local matrices (element-independent in the LGL case for B)
B_loc_ref = W_diag * D_hat;          % (nln x nln) local transport matrix

M = sparse(ndof, ndof);
B = sparse(ndof, ndof);

for ie = 1 : ne
    iglo       = connectivity(1:nln, ie);
    [BJ, ~]    = GetJacobian(coord(iglo, :), nodes_1D);
    % BJ = (x_b - x_a)/2  (Jacobian of the affine map)

    % Local mass: diagonal, scaled by Jacobian
    M_loc = BJ * W_diag;             % (nln x nln), diagonal

    % Local transport: BJ cancels between physical derivative and dx
    % B_loc(i,j) = int_{-1}^{1} phi_i(xi) * (1/BJ) * d phi_j/d xi * BJ d xi
    %            = int_{-1}^{1} phi_i(xi) * d phi_j/d xi  d xi  (BJ-independent)
    %            = w_i * D_hat(i,j)
    B_loc = B_loc_ref;               % (nln x nln)

    % Global assembly (standard FEM scatter-add)
    M(iglo, iglo) = M(iglo, iglo) + M_loc;
    B(iglo, iglo) = B(iglo, iglo) + B_loc;
end

% =========================================================================
% 3.  APPLY BOUNDARY CONDITIONS
% =========================================================================

if strcmp(Data.boundary, 'PP')
    % -------------------------------------------------------------------
    %  PERIODIC BC:  identify DOF 1 (x=0) with DOF ndof (x=L)
    %  "Fuse" last row/column into first row/column, then discard last.
    % -------------------------------------------------------------------

    % Fuse mass (diagonal => only (1,1) and (end,end) matter)
    M(1,1)     = M(1,1)   + M(end,end);
    M          = M(1:end-1, 1:end-1);

    % Fuse transport (full row and column)
    B(1, :)    = B(1, :)  + B(end, :);
    B(:, 1)    = B(:, 1)  + B(:, end);
    B          = B(1:end-1, 1:end-1);

    ndof_per   = ndof - 1;            % effective DOF count after periodicity
    x_per      = coord(1:end-1, 1);   % (ndof_per x 1) unique coordinates

    % -------------------------------------------------------------------
    % 4a.  BUILD GLOBAL BLOCK SYSTEM   M_glob * dU/dt + S_glob * U = 0
    %
    %   M_glob = blkdiag(M, M)
    %   S_glob = [ 0        ,   B    ]
    %            [ g*H * B  ,  gamma*M]     (gamma=0 without friction)
    % -------------------------------------------------------------------
    M_glob = blkdiag(M, M);

    Z = sparse(ndof_per, ndof_per);
    if Data.use_friction
        S_glob = [Z,                 B           ; ...
                  Data.g*Data.H*B,  Data.gamma*M];
    else
        S_glob = [Z,                 B           ; ...
                  Data.g*Data.H*B,  Z           ];
    end

    % -------------------------------------------------------------------
    % 5a.  INITIAL CONDITIONS
    % -------------------------------------------------------------------
    eta0_vec = Data.eta0(x_per);
    q0_vec   = Data.q0(x_per);
    U        = [eta0_vec; q0_vec];    % (2*ndof_per x 1)

    % -------------------------------------------------------------------
    % 6.  THETA-METHOD TIME LOOP
    %
    %  (M_glob + theta*dt*S_glob) U^{k+1} = (M_glob - (1-theta)*dt*S_glob) U^k
    % -------------------------------------------------------------------
    theta  = 0.5;                     % Crank-Nicolson
    dt     = Data.dt;
    Nsteps = round(Data.T / dt);

    LHS     = M_glob + theta       * dt * S_glob;
    RHS_mat = M_glob - (1 - theta) * dt * S_glob;

    % Pre-factorize LHS (constant throughout the simulation)
    [L_fac, U_fac, P_fac] = lu(LHS);

    fprintf('Time loop: %d steps, dt = %.2e ...\n', Nsteps, dt);
    prog = 0;
    fprintf('Progress: %3d%%\n', prog);

    for k = 1 : Nsteps
        rhs = RHS_mat * U;
        U   = U_fac \ (L_fac \ (P_fac * rhs));

        if mod(k, max(1, floor(Nsteps/10))) == 0
            prog = round(100 * k / Nsteps);
            fprintf('\b\b\b\b%3d%%', prog);
        end
    end
    fprintf('\n');

    % -------------------------------------------------------------------
    % 7.  POST-PROCESSING
    % -------------------------------------------------------------------
    eta    = U(1 : ndof_per);
    q      = U(ndof_per+1 : end);

    % Exact solution evaluated on the periodic DOFs at t = T
    eta_ex = Data.uex(x_per, Data.T);

    % L2 error via the periodic mass matrix (discrete inner product)
    diff    = eta - eta_ex;
    err_L2  = sqrt(diff' * M * diff);

    fprintf(' L2 error (eta):  %.6e\n', err_L2);
    fprintf('============================================================\n');

    Errors    = struct('L2', err_L2);
    Solutions = struct('uh',   eta,    ...
                       'u_ex', eta_ex, ...
                       'q',    q,      ...
                       'x',    x_per);

elseif strcmp(Data.boundary, 'WW')
    % -------------------------------------------------------------------
    %  WALL BC:  q(0,t) = q(L,t) = 0  (Dirichlet, no-flux)
    %            eta: free at walls    (natural / Neumann)
    %
    %  Implementation: eliminate the q DOFs at x=0 (DOF 1) and x=L
    %  (DOF ndof) from the q block of the system.
    %  Free indices in the global state vector [eta(1..ndof); q(1..ndof)]:
    %    free_eta : indices  1 .. ndof
    %    free_q   : indices  ndof+2 .. 2*ndof-1  (interior q DOFs)
    % -------------------------------------------------------------------
    x_all    = coord(:, 1);           % (ndof x 1)

    % Full block system (no BC applied yet)
    M_glob_f = blkdiag(M, M);
    Z_f      = sparse(ndof, ndof);
    S_glob_f = [Z_f,                B ; ...
                Data.g*Data.H*B,   Z_f];

    % Free-DOF mask: all eta, interior q (exclude q at x=0 and x=L)
    free_eta = (1 : ndof);
    free_q   = ndof + (2 : ndof-1);
    free     = [free_eta, free_q];    % total free DOFs

    % Reduced matrices
    M_r   = M_glob_f(free, free);
    S_r   = S_glob_f(free, free);

    % Initial condition (q=0 everywhere, including at walls)
    eta0_vec = Data.eta0(x_all);
    q0_all   = Data.q0(x_all);
    U_full   = [eta0_vec; q0_all];
    U        = U_full(free);

    % Theta-method matrices
    theta   = 0.5;
    dt      = Data.dt;
    Nsteps  = round(Data.T / dt);

    LHS_r   = M_r + theta       * dt * S_r;
    RHS_mat = M_r - (1 - theta) * dt * S_r;

    [L_fac, U_fac, P_fac] = lu(LHS_r);

    fprintf('Time loop (wall BC): %d steps, dt = %.2e ...\n', Nsteps, dt);
    prog = 0;
    fprintf('Progress: %3d%%\n', prog);

    for k = 1 : Nsteps
        rhs = RHS_mat * U;
        U   = U_fac \ (L_fac \ (P_fac * rhs));

        if mod(k, max(1, floor(Nsteps/10))) == 0
            prog = round(100 * k / Nsteps);
            fprintf('\b\b\b\b%3d%%', prog);
        end
    end
    fprintf('\n');

    % Reconstruct full solution (wall DOFs remain zero)
    U_full(free) = U;
    eta = U_full(1 : ndof);
    q   = U_full(ndof+1 : end);

    fprintf('============================================================\n');

    Errors    = struct('L2', NaN);    % no exact solution for wall BC
    Solutions = struct('uh',   eta,   ...
                       'u_ex', 0*eta, ...
                       'q',    q,     ...
                       'x',    x_all);

else
    error('MainSWE: boundary type ''%s'' not supported.', Data.boundary);
end

end % function MainSWE