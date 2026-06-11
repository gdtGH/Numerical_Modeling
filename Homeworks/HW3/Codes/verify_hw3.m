% verify_hw3.m
% Verification driver for the HW3 finite-volume SWE solvers.
%
%   A. Exact Riemann star state (h*, u*) on the three initial data
%   B. Consistency:  GodunovFlux(U,U) == F(U)
%   C. Mass conservation of the 1st-order Godunov run (closed walls)
%   D. Godunov 1st-order: shock-front position / left plateau (dam break)
%   E. Pure Lax-Wendroff on the stationary-shock case: blow-up reporting
%
% Run from the Codes/ folder:  >> verify_hw3
%
% Tests A and B are pure unit tests on the flux (no time stepping).
% Tests C-E launch MainFVM and therefore take a little longer.

clear; clc;
addpath Solvers
g = 9.81;

fprintf('############################################################\n');
fprintf('# HW3 VERIFICATION  (exact Godunov + pure Lax-Wendroff)\n');
fprintf('############################################################\n');

%==========================================================================
% A. EXACT RIEMANN STAR STATE
%    Recompute (h*, u*) with the same wave functions used by GodunovFlux
%    and compare with the reference values from the assignment.
%==========================================================================
fprintf('\n=== A. Exact Riemann star state (u = 0 everywhere) ===\n');

cases = { 'dam break  (hL=1.0, hR=0.1)', 1.0, 0.1,  0.396,  2.321;
          'reverse    (hL=0.1, hR=1.0)', 0.1, 1.0,  0.396, -2.321;
          'stationary (hL=1.0, hR=0.5)', 1.0, 0.5,  0.727,  0.923 };

tolA = 1e-2;
allA = true;
for k = 1:size(cases,1)
    name = cases{k,1};
    hL = cases{k,2};  hR = cases{k,3};
    hstar_ref = cases{k,4};  ustar_ref = cases{k,5};

    [hstar, ustar] = star_state(hL, 0, hR, 0, g);

    okh = abs(hstar - hstar_ref) < tolA;
    oku = abs(ustar - ustar_ref) < tolA;
    allA = allA && okh && oku;

    fprintf(['  %-28s  h* = %7.4f (ref %6.3f) %s   ' ...
             'u* = %8.4f (ref %7.3f) %s\n'], ...
            name, hstar, hstar_ref, tick(okh), ustar, ustar_ref, tick(oku));
end
fprintf('  --> A %s\n', verdict(allA));

%==========================================================================
% B. CONSISTENCY  GodunovFlux(U,U) == F(U)
%==========================================================================
fprintf('\n=== B. Consistency  GodunovFlux(U,U) == F(U) ===\n');
h = 0.7;  q = 0.3;
U = [h; q];
F_exact = [q; q^2/h + 0.5*g*h^2];
F_num   = GodunovFlux(U, U, g);
errB    = norm(F_num - F_exact, inf);
okB     = errB < 1e-14;
fprintf('  U = [%.3f; %.3f]\n', h, q);
fprintf('  F(U)            = [%.10f; %.10f]\n', F_exact(1), F_exact(2));
fprintf('  GodunovFlux(U,U)= [%.10f; %.10f]\n', F_num(1),  F_num(2));
fprintf('  |error|_inf = %.3e\n', errB);
fprintf('  --> B %s\n', verdict(okB));

%==========================================================================
% C. MASS CONSERVATION (1st-order Godunov, reflecting walls = closed system)
%==========================================================================
fprintf('\n=== C. Mass conservation (Godunov 1st order) ===\n');
Tests = DataTest();
okC = true;
for k = 1:numel(Tests)
    Data = Tests{k};
    dx   = (Data.domain(2)-Data.domain(1)) / Data.M;
    xc   = (Data.domain(1)+dx/2 : dx : Data.domain(2)-dx/2)';
    m0   = sum(Data.h0(xc)) * dx;

    Sol  = MainFVM(Data, 'godunov');
    m1   = sum(Sol.h) * dx;

    err  = abs(m1 - m0);
    okk  = err < 1e-10;
    okC  = okC && okk;
    fprintf('  %-16s  m0 = %.12f  mT = %.12f  |dm| = %.3e  %s\n', ...
            Data.name, m0, m1, err, tick(okk));
end
fprintf('  --> C %s\n', verdict(okC));

%==========================================================================
% D. GODUNOV 1st ORDER — dam-break shock front and left plateau
%==========================================================================
fprintf('\n=== D. Godunov 1st order — dam-break diagnostics ===\n');
Data = Tests{1};                 % dam break, T = 1.0
Sol  = MainFVM(Data, 'godunov');
report_dambreak(Sol, 'T = 1.0 (front already reflected off the wall)');

Data012   = Data;  Data012.T = 0.12;     % early time, front still travelling
Sol012    = MainFVM(Data012, 'godunov');
report_dambreak(Sol012, 'T = 0.12 (free-running front, expect x ~ 0.87)');
% Shock speed estimate from the front position at T = 0.12
xf = front_position(Sol012);
fprintf('  estimated front position x_f(0.12) = %.4f  ->  sigma ~ %.3f m/s\n', ...
        xf, (xf - 0.5)/0.12);

%==========================================================================
% E. PURE LAX-WENDROFF — stationary-shock blow-up
%==========================================================================
fprintf('\n=== E. Pure Lax-Wendroff blow-up (stationary shock) ===\n');
SolLW = MainFVM(Tests{3}, 'laxwendroff');
if SolLW.diverged
    fprintf('  Blow-up DETECTED.\n');
    fprintf('  last stable time   t = %.6f\n', SolLW.t);
    fprintf('  blow-up step time  t = %.6f\n', SolLW.t_blowup);
else
    fprintf('  No blow-up flagged (ran to t = %.6f).\n', SolLW.t);
end
fprintf('  h range at last stable state: [%.4f, %.4f]\n', ...
        min(SolLW.h), max(SolLW.h));

fprintf('\n############################################################\n');
fprintf('# SUMMARY:  A %s | B %s | C %s\n', ...
        verdict(allA), verdict(okB), verdict(okC));
fprintf('############################################################\n');

%==========================================================================
% LOCAL HELPERS
%==========================================================================
function [hstar, ustar] = star_state(hL, uL, hR, uR, g)
    phi   = @(h) wf(h,hL,g) + wf(h,hR,g) + (uR - uL);
    h_hi  = max(hL, hR);
    while phi(h_hi) < 0, h_hi = 2*h_hi; end
    hstar = fzero(phi, [1e-6, h_hi]);
    ustar = 0.5*(uL+uR) + 0.5*(wf(hstar,hR,g) - wf(hstar,hL,g));
end

function f = wf(h, hk, g)
    if h <= hk
        f = 2*(sqrt(g*h) - sqrt(g*hk));
    else
        f = (h - hk)*sqrt(0.5*g*(1/h + 1/hk));
    end
end

function xf = front_position(Sol)
    % Shock front = location of the steepest depth jump.  The rarefaction
    % fan is spread over many cells (gentle slope) whereas the shock drops
    % over ~1 cell, so the global max of |dh| singles out the front.
    h  = Sol.h;  x = Sol.x;
    [~, k] = max(abs(diff(h)));
    xf = 0.5*(x(k) + x(k+1));           % interface between the two cells
end

function report_dambreak(Sol, tag)
    h = Sol.h;  x = Sol.x;
    % left plateau = mean depth over the left tenth of the domain
    plateau = mean(h(x < 0.1));
    xf      = front_position(Sol);
    fprintf('  [%s]\n', tag);
    fprintf('     reached t = %.4f | left plateau h ~ %.4f | front x_f ~ %.4f | h in [%.4f, %.4f]\n', ...
            Sol.t, plateau, xf, min(h), max(h));
end

function s = tick(ok),    if ok, s = '[OK]';   else, s = '[XX]'; end, end
function s = verdict(ok), if ok, s = 'PASS';   else, s = 'FAIL'; end, end
