% =========================================================================
% NUMERICAL MODELING AND SIMULATION FOR ACOUSTICS - HOMEWORK 2
% Problem 1: 1D Shallow Water Equations / Second-Order Wave Equation
% Method: Spectral Element Method (SEM) with Newmark Time Integration
% =========================================================================
clear; clc; close all;

% =========================================================================
%% 1. Problem Parameters
% =========================================================================
L = 1.0;            % Domain length 
H = 1.0;            % Mean depth 
g = 9.81;           % Gravity 
c = sqrt(g*H);      % Wave speed
T = 0.5;            % Final time     

% Discretization study parameters 
Ne_array = [10, 20, 40, 80]; % Number of elements
N = 4;                       % Polynomial deg
dt = 1e-5;                   % Time step

% Initialize arrays to store errors
L2_errors = zeros(length(Ne_array), 1);
% =========================================================================
%% 2. Convergence Study Loop
% =========================================================================
fprintf('Starting Convergence Study...\n');
fprintf('--------------------------------------------------\n');
fprintf('Ne \t N \t dt \t\t L2 Error\n');
fprintf('--------------------------------------------------\n');

for iter = 1:length(Ne_array)
    Ne = Ne_array(iter);
    
    % --- A. MESH GENERATION & PERIODIC MAPPING ---
    % Get Reference GLL points (xi) and weights (w), and derivative matrix (D_hat)
    [xi, w, D_hat] = gll_routine(N);
    
    % Create uniform elements
    x_edges = linspace(0, L, Ne+1);
    
    % Total degrees of freedom for periodic domain
    Ndof = Ne * N; 
    global_x = zeros(Ndof, 1);
    
    % Connectivity matrix
    % For periodic Bound.Conds, the last node of the last element maps back to node 1.
    LM = zeros(Ne, N+1);
    for e = 1:Ne
        for i = 1:(N+1)
            global_idx = (e-1)*N + i;
            if global_idx > Ndof
                LM(e, i) = 1; % Periodic wrap 
            else
                LM(e, i) = global_idx;
            end
        end
        % Map local GLL points to global physical coordinates
        Le = x_edges(e+1) - x_edges(e);
        x_local = x_edges(e) + (xi + 1) * Le / 2;
        global_x(LM(e, 1:N)) = x_local(1:N); % Exclude last node to avoid overwrite
    end
    % =========================================================================
    % --- B. MATRIX ASSEMBLY ---
    % =========================================================================
    % Initialize sparse global matrices
    M_global = sparse(Ndof, Ndof);
    K_global = sparse(Ndof, Ndof);
    
    for e = 1:Ne
        Le = x_edges(e+1) - x_edges(e);
        J = Le / 2; % Jacobian of the affine mapping
        
        % Local Mass Matrix : diagonal (GLL quadrature)
        M_loc = diag(w * J);
        
        % Local Stiffness Matrix
        % K_ij = gH * \int (phi_i' * phi_j') dx
        K_loc = g * H * (D_hat' * diag(w) * D_hat) / J;
        
        % Assemble into global matrices using Connectivity Matrix (LM)
        for i = 1:(N+1)
            I = LM(e, i);
            for j = 1:(N+1)
                J_idx = LM(e, j);
                M_global(I, J_idx) = M_global(I, J_idx) + M_loc(i, j);
                K_global(I, J_idx) = K_global(I, J_idx) + K_loc(i, j);
            end
        end
    end
    % =========================================================================
    % --- C. INITIAL CONDITIONS ---
    % =========================================================================
    % eta(x,0) = exp(-50*(x-0.5)^2) 
    eta0 = exp(-50 * (global_x - 0.5).^2);
    
    % q(x,0) = c*eta(x,0)  implies eta_t(x,0) = -q_x(x,0) = -c*eta_x(x,0)
    % Exact spatial derivative for the initial velocity
    eta_x_0 = -100 * (global_x - 0.5) .* exp(-50 * (global_x - 0.5).^2);
    eta_dot0 = -c * eta_x_0; 
    % =========================================================================
    % --- D. NEWMARK TIME INTEGRATION ---
    % =========================================================================
    % Average Acceleration Method 
    % beta = 0.25; gamma = 0.5;
    gamma = 0.505; 
    beta = 0.25 * (0.5 + gamma)^2;
    
    % Initial acceleration: M * \eta_{\ddot{0}} = -K * eta0
    eta_ddot0 = -M_global \ (K_global * eta0);
    
    % Pre-allocate state vectors
    eta = eta0;
    eta_dot = eta_dot0;
    eta_ddot = eta_ddot0;
    
    % Pre-factorize the iteration matrix for efficiency
    A = M_global + beta * dt^2 * K_global;
    
    Nt = ceil(T / dt);
    for n = 1:Nt
        % Predictor step
        eta_pred = eta + dt * eta_dot + (dt^2 / 2) * (1 - 2*beta) * eta_ddot;
        eta_dot_pred = eta_dot + dt * (1 - gamma) * eta_ddot;
        
        % Solve for new acceleration
        RHS = -K_global * eta_pred;
        eta_ddot_new = A \ RHS;
        
        % Corrector step
        eta = eta_pred + beta * dt^2 * eta_ddot_new;
        eta_dot = eta_dot_pred + gamma * dt * eta_ddot_new;
        eta_ddot = eta_ddot_new;
    end
    
    % --- E. ERROR COMPUTATION ---
    % Exact solution at T: \eta(x-cT, 0) accounting for periodicity 
    X_shifted = mod(global_x - c * T, L);
    eta_exact = exp(-50 * (X_shifted - 0.5).^2);
    
    % L2 Error norm using the global mass matrix 
    error_vec = eta_exact - eta;
    L2_errors(iter) = sqrt(error_vec' * M_global * error_vec);
    
    fprintf('%d \t %d \t %.1e \t %.4e\n', Ne, N, dt, L2_errors(iter));
end
% =========================================================================
%% 3. Plotting Final Results (Finest Mesh)
% =========================================================================
figure('Name', 'SEM Shallow Water Solver', 'Position', [100, 100, 800, 400]);

% Append the first node to the end for a seamless periodic plot
x_plot = [global_x; L];
eta_plot = [eta; eta(1)];
eta_exact_plot = [eta_exact; eta_exact(1)];

plot(x_plot, eta_exact_plot, 'k-', 'LineWidth', 2, 'DisplayName', 'Exact Solution');
hold on; grid on;
plot(x_plot, eta_plot, 'ro--', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', sprintf('SEM (Ne=%d, N=%d)', Ne, N));
xlabel('Domain x [m]', 'Interpreter', 'latex');
ylabel('Free Surface Elevation $\eta(x,T)$', 'Interpreter', 'latex');
title(sprintf('Wave Propagation at T = %.2f s', T), 'Interpreter', 'latex');
legend('Location', 'best', 'Interpreter', 'latex');
xlim([0 L]);

% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function [xi, w, D] = gll_routine(N)
    % Computes GLL nodes, weights, and derivative matrix
    % for polynomial degree N on the reference domain [-1, 1].
    
    % Nodes
    xi = zeros(N+1, 1);
    xi(1) = -1; xi(end) = 1;
    if N > 1
        % Roots of the derivative of the Legendre polynomial L_N
        % Calculated via eigenvalues of Jacobi matrix for L'_N
        diag_elements = zeros(N-1, 1);
        off_diag = sqrt((1:N-2).*((1:N-2)+2) ./ ( (2*(1:N-2)+1).*(2*(1:N-2)+3) ));
        Jmat = diag(off_diag, 1) + diag(off_diag, -1);
        xi(2:N) = sort(eig(Jmat));
    end
    
    % Weights
    L_N = zeros(N+1, 1);
    for i = 1:N+1
        [L_N(i), ~] = legendre_eval(N, xi(i));
    end
    w = 2 ./ (N * (N + 1) * L_N.^2);
    
    % Derivative Matrix
    D = zeros(N+1, N+1);
    for i = 1:N+1
        for j = 1:N+1
            if i ~= j
                [L_Ni, ~] = legendre_eval(N, xi(i));
                [L_Nj, ~] = legendre_eval(N, xi(j));
                D(i,j) = (L_Ni / L_Nj) / (xi(i) - xi(j));
            elseif i == j && i == 1
                D(i,j) = -N*(N+1)/4;
            elseif i == j && i == N+1
                D(i,j) = N*(N+1)/4;
            else
                D(i,j) = 0;
            end
        end
    end
end

function [L_val, L_prime_val] = legendre_eval(N, x)
    % Evaluates the Legendre polynomial of degree N and its derivative at x
    if N == 0
        L_val = 1; L_prime_val = 0; return;
    elseif N == 1
        L_val = x; L_prime_val = 1; return;
    end
    L_prev = 1; L_curr = x;
    dL_prev = 0; dL_curr = 1;
    for k = 1:N-1
        L_next = ((2*k + 1)*x*L_curr - k*L_prev) / (k + 1);
        dL_next = ((2*k + 1)*x*dL_curr + (2*k + 1)*L_curr - k*dL_prev) / (k + 1);
        L_prev = L_curr; L_curr = L_next;
        dL_prev = dL_curr; dL_curr = dL_next;
    end
    L_val = L_curr;
    L_prime_val = dL_curr;
end
% =========================================================================
%% 4. Convergence Analysis (Calculate Rates and Plot)
% =========================================================================
fprintf('\n--- Convergence Rates (r) ---\n');
r_values = zeros(length(Ne_array)-1, 1);
for i = 1:length(Ne_array)-1
    % Calculate r = ln(E_coarse / E_fine) / ln(Ne_fine / Ne_coarse)
    r_values(i) = log(L2_errors(i) / L2_errors(i+1)) / log(Ne_array(i+1) / Ne_array(i));
    fprintf('Between Ne=%d and Ne=%d: r = %.2f\n', Ne_array(i), Ne_array(i+1), r_values(i));
end

% Generate Log-Log Plot
figure('Name', 'L2 Error Convergence', 'Position', [200, 200, 600, 500]);
loglog(Ne_array, L2_errors, 'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'SEM $L^2$ Error');
hold on; grid on;

% Plot reference slope O(h^(N+1))
ref_slope = (Ne_array / Ne_array(1)).^(-(N+1)); 
ref_line = L2_errors(1) * ref_slope;
loglog(Ne_array, ref_line, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('$\\mathcal{O}(h^{%d})$ Reference', N+1));

% Formatting
xlabel('Number of Elements $N_e \propto h^{-1}$', 'Interpreter', 'latex');
ylabel('$L^2$ Error Norm $||\eta - \eta_h||_{L^2}$', 'Interpreter', 'latex');
title(sprintf('Spectral Element Convergence ($N=%d$)', N), 'Interpreter', 'latex');
legend('Location', 'southwest', 'Interpreter', 'latex');

% Improve axis ticks for log scale readability
set(gca, 'XTick', Ne_array);
xlim([min(Ne_array)*0.9, max(Ne_array)*1.1]);

%% 4. Point 5: Optional Studies
% Study A: Reflecting Walls ONLY (\gamma = 0)

% Study B: Friction ONLY with Periodic BCs (\gamma > 0)

fprintf('\nRunning Point 5: Isolated Studies on Boundaries and Friction...\n');

% Shared Parameters
gamma_fric = 2.0; % Positive friction coefficient for Study B
Ne = Ne_array(end); % Use finest mesh
beta = 0.25; gamma_nw = 0.5; % Newmark parameters
Nt = ceil(T / dt);

% =========================================================================
%% STUDY A: REFLECTING WALLS (NO FRICTION)
% =========================================================================
% 1. Non-Periodic Mesh
Ndof_wall = Ne * N + 1; 
global_x_wall = zeros(Ndof_wall, 1);
LM_wall = zeros(Ne, N+1);
x_edges = linspace(0, L, Ne+1);

for e = 1:Ne
    for i = 1:(N+1)
        LM_wall(e, i) = (e-1)*N + i; % Standard contiguous mapping
    end
    Le = x_edges(e+1) - x_edges(e);
    x_local = x_edges(e) + (xi + 1) * Le / 2;
    global_x_wall(LM_wall(e, 1:N+1)) = x_local;
end

% 2. Matrix Assembly
M_wall = sparse(Ndof_wall, Ndof_wall);
K_wall = sparse(Ndof_wall, Ndof_wall);

for e = 1:Ne
    Le = x_edges(e+1) - x_edges(e);
    J = Le / 2;
    M_loc = diag(w * J);
    K_loc = g * H * (D_hat' * diag(w) * D_hat) / J;
    
    for i = 1:(N+1)
        I = LM_wall(e, i);
        for j = 1:(N+1)
            J_idx = LM_wall(e, j);
            M_wall(I, J_idx) = M_wall(I, J_idx) + M_loc(i, j);
            K_wall(I, J_idx) = K_wall(I, J_idx) + K_loc(i, j);
        end
    end
end

% 3. Initial Conditions & Newmark Integration (Conservative)
eta0_w = exp(-50 * (global_x_wall - 0.5).^2);
eta_x_0_w = -100 * (global_x_wall - 0.5) .* exp(-50 * (global_x_wall - 0.5).^2);
eta_dot0_w = -c * eta_x_0_w; 

eta_ddot0_w = -M_wall \ (K_wall * eta0_w);
eta_w = eta0_w; eta_dot_w = eta_dot0_w; eta_ddot_w = eta_ddot0_w;

A_wall = M_wall + beta * dt^2 * K_wall; % No damping matrix

for n = 1:Nt
    eta_pred = eta_w + dt * eta_dot_w + (dt^2 / 2) * (1 - 2*beta) * eta_ddot_w;
    eta_dot_pred = eta_dot_w + dt * (1 - gamma_nw) * eta_ddot_w;
    
    RHS = -K_wall * eta_pred;
    eta_ddot_new = A_wall \ RHS;
    
    eta_w = eta_pred + beta * dt^2 * eta_ddot_new;
    eta_dot_w = eta_dot_pred + gamma_nw * dt * eta_ddot_new;
    eta_ddot_w = eta_ddot_new;
end

% =========================================================================
%% STUDY B: FRICTION WITH PERIODIC BOUNDARIES
% =========================================================================
% Note: We reuse the periodic matrices (M_global, K_global) from Point 4.
% If you run this standalone, ensure M_global and K_global are assembled as in Point 4.

% 1. Friction Matrix
C_per = gamma_fric * M_global;

% 2. Initial Conditions & Newmark Integration (Dissipative)
eta_ddot0_fric = -M_global \ (C_per * eta_dot0 + K_global * eta0);
eta_fric = eta0; eta_dot_fric = eta_dot0; eta_ddot_fric = eta_ddot0_fric;

A_fric = M_global + gamma_nw * dt * C_per + beta * dt^2 * K_global;

for n = 1:Nt
    eta_pred = eta_fric + dt * eta_dot_fric + (dt^2 / 2) * (1 - 2*beta) * eta_ddot_fric;
    eta_dot_pred = eta_dot_fric + dt * (1 - gamma_nw) * eta_ddot_fric;
    
    RHS = -K_global * eta_pred - C_per * eta_dot_pred;
    eta_ddot_new = A_fric \ RHS;
    
    eta_fric = eta_pred + beta * dt^2 * eta_ddot_new;
    eta_dot_fric = eta_dot_pred + gamma_nw * dt * eta_ddot_new;
    eta_ddot_fric = eta_ddot_new;
end

% =========================================================================
%% PLOTTING THE COMPARISON
% =========================================================================
figure('Name', 'Point 5 Isolated Studies', 'Position', [100, 100, 1000, 450]);

% Plot A: Reflecting Walls
subplot(1, 2, 1);
plot(global_x_wall, eta0_w, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Initial $\eta(x,0)$');
hold on; grid on;
plot(global_x_wall, eta_w, 'r-', 'LineWidth', 2, 'DisplayName', 'Reflecting Walls ($\gamma = 0$)');
yline(1, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Unitary Elevation ($\\eta = 1 $)'));
xlabel('Domain $x$ [m]', 'Interpreter', 'latex');
ylabel('Elevation $\eta(x,T)$', 'Interpreter', 'latex');
title('Isolated Effect of Reflecting Walls', 'Interpreter', 'latex');
legend('Location', 'best', 'Interpreter', 'latex');
xlim([0 L]);

% Plot B: Friction (Periodic)
subplot(1, 2, 2);
x_plot_per = [global_x; L]; % Append for periodic visual wrap
eta0_plot_per = [eta0; eta0(1)];
eta_fric_plot = [eta_fric; eta_fric(1)];
eta_exact_plot = [eta_exact; eta_exact(1)]; % From Point 4 exact solution

plot(x_plot_per, eta0_plot_per, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Initial $\eta(x,0)$');
hold on; grid on;
plot(x_plot_per, eta_exact_plot, 'k-', 'LineWidth', 1, 'DisplayName', 'Exact Undamped');
plot(x_plot_per, eta_fric_plot, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Friction ($\\gamma = %.1f$)', gamma_fric));
yline(1, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('Unitary Elevation ($\\eta = 1 $)'));
xlabel('Domain $x$ [m]', 'Interpreter', 'latex');
ylabel('Elevation $\eta(x,T)$', 'Interpreter', 'latex');
title('Isolated Effect of Friction (Periodic)', 'Interpreter', 'latex');
legend('Location', 'best', 'Interpreter', 'latex');
xlim([0 L]);
ylim([0, 1.8]);
