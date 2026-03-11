%% RunMainSWE.m
%==========================================================================
%  Top-level driver for Homework 2 — Linearised Shallow Water Equations
%  SEM spatial discretisation + Crank-Nicolson time integration
%
%  Sections:
%   POINT 4  — Snapshot at T=0.5: numerical vs exact solution
%   POINT 4  — h-convergence  (fix p=4, vary Ne = 10,20,40,80)
%   POINT 4  — p-convergence  (fix Ne=5, vary p = 1,2,3,4,5,6)
%   POINT 5a — Reflecting walls  (wall BC, energy split)
%   POINT 5b — Friction term  (compare periodic with/without gamma)
%==========================================================================

clear; close all; clc;

%% ── Paths ────────────────────────────────────────────────────────────────
addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing
addpath SemLib

%% ── Global style ─────────────────────────────────────────────────────────
fs  = 25;           % axis / tick / legend font size
fst = 50;           % title font size
lw  = 2;            % line width
ms  = 8;            % marker size
c1  = '#0072BD';    % blue
c2  = '#D95319';    % orange
c3  = '#77AC30';    % green
c4  = '#7E2F8E';    % purple

% =========================================================================
%  POINT 4 — Snapshot: numerical solution vs exact at t = T
% =========================================================================
fprintf('\n=== POINT 4: Snapshot at T = 0.5 ===\n');

Data_snap      = DataTest('HW2_P4');
Data_snap.p    = 4;
Data_snap.dt   = 1e-3;
Data_snap.calc_errors = false;

[~, Sol_snap, Fem_snap, ~] = MainSWE(Data_snap, 20);

figure('Name','P4.1_eta_snapshot','NumberTitle','off');
plot(Sol_snap.x, Sol_snap.u_ex, '-',  'Color', c1, 'LineWidth', lw); hold on;
plot(Sol_snap.x, Sol_snap.uh,   '--', 'Color', c2, 'LineWidth', lw, ...
     'MarkerSize', ms);
xlabel('x',           'FontSize', fs);
ylabel('\eta(x, T)',  'FontSize', fs);
title('SEM solution vs exact — T = 0.5', ...
      'FontSize', fst, 'FontWeight', 'bold');
legend('Exact \eta', 'Numerical \eta', ...
       'FontSize', fs, 'Location', 'best');
grid on; box on;
ax = gca; ax.FontSize = fs;

% ── Also plot the discharge q ─────────────────────────────────────────────
q_ex = Data_snap.c * Sol_snap.u_ex;     % exact q = c * eta (no friction)

figure('Name','P4.2_q_snapshot','NumberTitle','off');
plot(Sol_snap.x, q_ex,       '-',  'Color', c1, 'LineWidth', lw); hold on;
plot(Sol_snap.x, Sol_snap.q, '--', 'Color', c2, 'LineWidth', lw);
xlabel('x',           'FontSize', fs);
ylabel('q(x, T)',     'FontSize', fs);
title('Discharge q vs exact — T = 0.5', ...
      'FontSize', fst, 'FontWeight', 'bold');
legend('Exact q', 'Numerical q', 'FontSize', fs, 'Location', 'best');
grid on; box on;
ax = gca; ax.FontSize = fs;

% =========================================================================
%  POINT 4 — h-convergence  (fix p = 4, vary Ne)
%  dt is chosen small enough to isolate the spatial error.
%  For high-order SEM (p=4), use dt = h/c * 0.05 to keep temporal
%  error below spatial error on the coarser meshes.
% =========================================================================
fprintf('\n=== POINT 4: h-convergence (p=4, Ne=10,20,40) ===\n');

Data_conv   = DataTest('HW2_P4');
Data_conv.p = 4;
Data_conv.calc_errors = false;

% Ne=[10,20,40]: the Gaussian (effective width ~0.07) is well-resolved
% by Ne=10 (h=0.1, 10 elements cover the domain) so all three points
% are in the ASYMPTOTIC regime and rates should cluster near p+1=5.
%
% dt = alpha * h^{(p+1)/2} = alpha * h^2.5  balances O(dt^2) ~ O(h^5).
% alpha=0.05  =>  dt^2 / h^5 = alpha^2 = 0.0025  (temporal < 0.25% of spatial).
Ne_list  = [10, 20, 40];
nNe      = length(Ne_list);
alpha_dt = 0.05;

h_vec   = zeros(1, nNe);
eL2_vec = zeros(1, nNe);
dt_used = zeros(1, nNe);

for i = 1 : nNe
    Ne = Ne_list(i);
    h  = Data_conv.L / Ne;
    Data_conv.dt = alpha_dt * h^((Data_conv.p+1)/2);

    [Err, ~, Fem, ~] = MainSWE(Data_conv, Ne);

    h_vec(i)   = Fem.h;
    eL2_vec(i) = Err.L2;
    dt_used(i) = Data_conv.dt;

    fprintf('  Ne = %3d,  h = %.4f,  L2 err = %.4e,  dt = %.2e\n', ...
            Ne, Fem.h, Err.L2, Data_conv.dt);
end

% ── Convergence rates ──────────────────────────────────────────────────
rates = zeros(1, nNe-1);
for i = 1 : nNe-1
    rates(i) = log(eL2_vec(i)/eL2_vec(i+1)) / log(h_vec(i)/h_vec(i+1));
end
fprintf('\n  Convergence rates (L2, eta):\n');
for i = 1 : nNe-1
    fprintf('    Ne = %2d -> %2d :  rate = %.2f\n', ...
            Ne_list(i), Ne_list(i+1), rates(i));
end

% ── Reference slope O(h^{p+1}), anchored at COARSEST point (Ne=10) ────
p_ref     = Data_conv.p;
h_ref     = linspace(min(h_vec)*0.7, max(h_vec)*1.1, 200);
slope_ref = eL2_vec(1) * (h_ref / h_vec(1)).^(p_ref+1);

figure('Name','P4.3_h_convergence','NumberTitle','off');
loglog(h_vec, eL2_vec, '-o', 'Color', c2, 'LineWidth', lw, ...
       'MarkerSize', ms, 'MarkerFaceColor', c2); hold on;
loglog(h_ref, slope_ref, '--', 'Color', c1, 'LineWidth', lw-0.5);

% Rate labels above midpoints
for i = 1 : nNe-1
    xm = sqrt(h_vec(i) * h_vec(i+1));
    ym = sqrt(eL2_vec(i) * eL2_vec(i+1)) * 2.5;
    text(xm, ym, sprintf('rate = %.1f', rates(i)), ...
         'FontSize', fs-4, 'Color', c2, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center');
end

xlabel('h = L/N_e',           'FontSize', fs);
ylabel('L^2 error on \eta',   'FontSize', fs);
title(sprintf('h-convergence — SEM, p = %d', p_ref), ...
      'FontSize', fst, 'FontWeight', 'bold');
legend({'L^2 error', sprintf('O(h^{%d})', p_ref+1)}, ...
       'FontSize', fs, 'Location', 'northwest');
grid on; box on;
ax = gca; ax.FontSize = fs;

% =========================================================================
%  POINT 4 — p-convergence  (fix Ne = 5, vary p = 1 .. 6)
% =========================================================================
fprintf('\n=== POINT 4: p-convergence (Ne=5, p=1..6) ===\n');

Data_pconv    = DataTest('HW2_P4');
% dt fixed small: temporal error well below spatial for all p considered.
% At p=1, Ne=5 the spatial error can be O(h^2) ≈ 0.04, so dt=1e-3 is fine.
% At p=5, Ne=5 the spatial error drops to ~1e-8; dt=1e-4 gives dt^2=1e-8.
% Use dt=5e-5 as safe common choice.
Data_pconv.dt = 5e-5;
Data_pconv.calc_errors = false;

Ne_fixed = 5;
% CreateMesh supports only p = 1 .. 5
p_list   = [1, 2, 3, 4, 5];
nP       = length(p_list);

eL2_p  = zeros(1, nP);
DoF_p  = zeros(1, nP);           % total DOF (Ne * p)

for i = 1 : nP
    Data_pconv.p = p_list(i);
    [Err, ~, Fem, ~] = MainSWE(Data_pconv, Ne_fixed);
    eL2_p(i) = Err.L2;
    DoF_p(i) = Fem.ndof - 1;     % periodic DOF count
    fprintf('  p = %d,  DoF = %3d,  L2 err = %.4e\n', ...
            p_list(i), DoF_p(i), eL2_p(i));
end

figure('Name','P4.4_p_convergence','NumberTitle','off');
semilogy(p_list, eL2_p, '-s', 'Color', c3, 'LineWidth', lw, ...
         'MarkerSize', ms, 'MarkerFaceColor', c3);

% Annotate algebraic decay rate between consecutive p values
for i = 1 : nP-1
    xm = 0.5*(p_list(i) + p_list(i+1));
    ym = sqrt(eL2_p(i) * eL2_p(i+1));
    rate_p = log(eL2_p(i)/eL2_p(i+1));   % drop per unit p in log scale
    text(xm, ym * 0.4, sprintf('\\Delta = %.1f', rate_p), ...
         'FontSize', fs-5, 'Color', c3, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center');
end

xlabel('Polynomial degree p',    'FontSize', fs);
ylabel('L^2 error on \eta',      'FontSize', fs);
title(sprintf('p-convergence — SEM, N_e = %d', Ne_fixed), ...
      'FontSize', fst, 'FontWeight', 'bold');
legend({'L^2 error'}, 'FontSize', fs, 'Location', 'northeast');
grid on; box on;
ax = gca; ax.FontSize = fs;

% =========================================================================
%  POINT 5a — Reflecting walls  (wall BC)
% =========================================================================
fprintf('\n=== POINT 5a: Reflecting walls ===\n');

Data_wall   = DataTest('HW2_P5a');
Data_wall.p = 4;
Data_wall.dt = 1e-3;
Data_wall.calc_errors = false;

[~, Sol_wall, ~, ~] = MainSWE(Data_wall, 40);

figure('Name','P5a.1_eta_wall','NumberTitle','off');
plot(Sol_wall.x, Sol_wall.uh, '-', 'Color', c1, 'LineWidth', lw);
xlabel('x',           'FontSize', fs);
ylabel('\eta(x, T)',  'FontSize', fs);
title('Reflecting walls — \eta at T = 0.5', ...
      'FontSize', fst, 'FontWeight', 'bold');
grid on; box on;
ax = gca; ax.FontSize = fs;

figure('Name','P5a.2_q_wall','NumberTitle','off');
plot(Sol_wall.x, Sol_wall.q, '-', 'Color', c2, 'LineWidth', lw);
xlabel('x',      'FontSize', fs);
ylabel('q(x,T)', 'FontSize', fs);
title('Reflecting walls — q at T = 0.5', ...
      'FontSize', fst, 'FontWeight', 'bold');
grid on; box on;
ax = gca; ax.FontSize = fs;

% =========================================================================
%  POINT 5b — Friction term  (compare periodic run with / without gamma)
% =========================================================================
fprintf('\n=== POINT 5b: Friction (gamma=1, periodic BC) ===\n');

% Run WITHOUT friction (P4 for reference)
Data_nofrict   = DataTest('HW2_P4');
Data_nofrict.p = 4;
Data_nofrict.dt = 1e-3;
[~, Sol_nofrict, ~, ~] = MainSWE(Data_nofrict, 20);

% Run WITH friction (P5b)
Data_frict   = DataTest('HW2_P5b');
Data_frict.p = 4;
Data_frict.dt = 1e-3;
[~, Sol_frict, ~, ~] = MainSWE(Data_frict, 20);

figure('Name','P5b.1_eta_friction','NumberTitle','off');
plot(Sol_nofrict.x, Sol_nofrict.uh, '-',  'Color', c1, 'LineWidth', lw); hold on;
plot(Sol_frict.x,   Sol_frict.uh,   '--', 'Color', c4, 'LineWidth', lw);
plot(Sol_nofrict.x, Sol_nofrict.u_ex, ':', 'Color', c3, 'LineWidth', lw-0.5);
xlabel('x',           'FontSize', fs);
ylabel('\eta(x, T)',  'FontSize', fs);
title('Friction effect on \eta — T = 0.5', ...
      'FontSize', fst, 'FontWeight', 'bold');
legend({'No friction', sprintf('Friction \\gamma = %g', Data_frict.gamma), ...
        'Exact (no friction)'}, ...
       'FontSize', fs, 'Location', 'best');
grid on; box on;
ax = gca; ax.FontSize = fs;

figure('Name','P5b.2_q_friction','NumberTitle','off');
plot(Sol_nofrict.x, Sol_nofrict.q, '-',  'Color', c1, 'LineWidth', lw); hold on;
plot(Sol_frict.x,   Sol_frict.q,   '--', 'Color', c4, 'LineWidth', lw);
xlabel('x',      'FontSize', fs);
ylabel('q(x,T)', 'FontSize', fs);
title('Friction effect on q — T = 0.5', ...
      'FontSize', fst, 'FontWeight', 'bold');
legend({'No friction', sprintf('Friction \\gamma = %g', Data_frict.gamma)}, ...
       'FontSize', fs, 'Location', 'best');
grid on; box on;
ax = gca; ax.FontSize = fs;

% ── Energy decay plot for friction case ───────────────────────────────────
% Re-run P5b collecting energy at each step to show dissipation
fprintf('\n  Re-running P5b to track energy over time ...\n');

Data_E    = DataTest('HW2_P5b');
Data_E.p  = 4;
Data_E.dt = 1e-3;
Ne_E      = 20;

[Region_E]    = CreateMesh(Data_E, Ne_E);
[Fem_E]       = CreateFemregion(Data_E, Region_E);

ndof_E = Fem_E.ndof;
nln_E  = Fem_E.nln;

[nodes_E, w_E] = xwlgl(nln_E);
[~, GP_E]      = basis_and_der_at_lgl(nodes_E, nln_E);
D_E   = squeeze(GP_E(:,1,:));
W_E   = diag(w_E);
Bloc_E = W_E * D_E;

M_E = sparse(ndof_E, ndof_E);
B_E = sparse(ndof_E, ndof_E);
for ie = 1 : Fem_E.ne
    ig = Fem_E.connectivity(1:nln_E, ie);
    [BJ_E, ~] = GetJacobian(Fem_E.coord(ig,:), nodes_E);
    M_E(ig, ig) = M_E(ig, ig) + BJ_E * W_E;
    B_E(ig, ig) = B_E(ig, ig) + Bloc_E;
end

% Periodic fusion
M_E(1,1)  = M_E(1,1)  + M_E(end,end);
M_E       = M_E(1:end-1, 1:end-1);
B_E(1,:)  = B_E(1,:)  + B_E(end,:);
B_E(:,1)  = B_E(:,1)  + B_E(:,end);
B_E       = B_E(1:end-1, 1:end-1);

ndp = ndof_E - 1;
x_E = Fem_E.coord(1:end-1, 1);

Mg  = blkdiag(M_E, M_E);
Sg  = [sparse(ndp,ndp),        B_E          ; ...
       Data_E.g*Data_E.H*B_E,  Data_E.gamma*M_E];

eta0_E = Data_E.eta0(x_E);
q0_E   = Data_E.q0(x_E);
U_E    = [eta0_E; q0_E];

theta_E  = 0.5;
dt_E     = Data_E.dt;
Nst_E    = round(Data_E.T / dt_E);
LHS_E    = Mg + theta_E * dt_E * Sg;
RHSm_E   = Mg - (1 - theta_E) * dt_E * Sg;
[Lf_E, Uf_E, Pf_E] = lu(LHS_E);

% Reference energy without friction (P4)
Data_E0   = DataTest('HW2_P4');
Data_E0.p = 4;  Data_E0.dt = dt_E;
[Region_E0]  = CreateMesh(Data_E0, Ne_E);
[Fem_E0]     = CreateFemregion(Data_E0, Region_E0);
M_E0 = sparse(Fem_E0.ndof, Fem_E0.ndof);
B_E0 = sparse(Fem_E0.ndof, Fem_E0.ndof);
for ie = 1 : Fem_E0.ne
    ig = Fem_E0.connectivity(1:nln_E, ie);
    [BJ0,~] = GetJacobian(Fem_E0.coord(ig,:), nodes_E);
    M_E0(ig,ig) = M_E0(ig,ig) + BJ0 * W_E;
    B_E0(ig,ig) = B_E0(ig,ig) + Bloc_E;
end
M_E0(1,1)  = M_E0(1,1)  + M_E0(end,end);
M_E0       = M_E0(1:end-1, 1:end-1);
B_E0(1,:)  = B_E0(1,:)  + B_E0(end,:);
B_E0(:,1)  = B_E0(:,1)  + B_E0(:,end);
B_E0       = B_E0(1:end-1, 1:end-1);

Mg0  = blkdiag(M_E0, M_E0);
Sg0  = [sparse(ndp,ndp), B_E0; Data_E0.g*Data_E0.H*B_E0, sparse(ndp,ndp)];
U_E0 = [Data_E0.eta0(x_E); Data_E0.q0(x_E)];
LHS0 = Mg0 + theta_E*dt_E*Sg0;
RHSm0= Mg0 - (1-theta_E)*dt_E*Sg0;
[Lf0,Uf0,Pf0] = lu(LHS0);

% Energy:  E = 0.5 * ( g * eta^T M eta  +  (1/H) * q^T M q )
g_E = Data_E.g;  H_E = Data_E.H;
energy_fun = @(U_vec) 0.5 * (g_E * U_vec(1:ndp)'   * M_E  * U_vec(1:ndp) + ...
                              (1/H_E) * U_vec(ndp+1:end)' * M_E * U_vec(ndp+1:end));

nsave   = min(500, Nst_E);
save_every = max(1, floor(Nst_E / nsave));
t_save  = zeros(1, ceil(Nst_E/save_every)+1);
E_frict = zeros(1, ceil(Nst_E/save_every)+1);
E_nofri = zeros(1, ceil(Nst_E/save_every)+1);

cnt = 1;
t_save(cnt)  = 0;
E_frict(cnt) = energy_fun(U_E);
E_nofri(cnt) = energy_fun(U_E0);

for k = 1 : Nst_E
    U_E  = Uf_E  \ (Lf_E  \ (Pf_E  * (RHSm_E  * U_E)));
    U_E0 = Uf0   \ (Lf0   \ (Pf0   * (RHSm0   * U_E0)));
    if mod(k, save_every) == 0
        cnt = cnt + 1;
        t_save(cnt)  = k * dt_E;
        E_frict(cnt) = energy_fun(U_E);
        E_nofri(cnt) = energy_fun(U_E0);
    end
end
t_save  = t_save(1:cnt);
E_frict = E_frict(1:cnt);
E_nofri = E_nofri(1:cnt);

figure('Name','P5b.3_energy','NumberTitle','off');
plot(t_save, E_nofri / E_nofri(1), '-',  'Color', c1, 'LineWidth', lw); hold on;
plot(t_save, E_frict / E_frict(1), '--', 'Color', c4, 'LineWidth', lw);
xlabel('t',                       'FontSize', fs);
ylabel('E(t) / E(0)',             'FontSize', fs);
title('Energy conservation / dissipation', ...
      'FontSize', fst, 'FontWeight', 'bold');
legend({'No friction (CN, \theta=0.5)', ...
        sprintf('Friction \\gamma = %g', Data_E.gamma)}, ...
       'FontSize', fs, 'Location', 'best');
grid on; box on;
ax = gca; ax.FontSize = fs;

fprintf('\nAll plots generated.\n');