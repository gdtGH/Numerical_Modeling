% RunMainFVM.m
% Top-level script for HW3 — Nonlinear SWE with FVM.
%
% Runs all three Riemann test cases (Dam Break, Reverse Dam Break,
% Stationary Shock) with both the Godunov and Lax-Wendroff schemes.
% Produces one figure per test case, each with two subplots:
%   Left panel:  water height h(x, T)
%   Right panel: discharge    q(x, T)

clear; close all; clc;

addpath Solvers

%==========================================================================
% GLOBAL STYLE — edit ONLY here, never inside solver files
%==========================================================================
fs  = 25;           % axis/tick/legend font size
fst = 30;           % title font size
lw  = 2;            % line width
c1  = '#0072BD';    % Godunov      (blue)
c2  = '#D95319';    % Lax-Wendroff (orange)

%==========================================================================
% LOAD TEST CASES
%==========================================================================
Tests = DataTest();

%==========================================================================
% TEST LABELS (for figure titles)
%==========================================================================
test_titles = {'Dam Break', 'Reverse Dam Break', 'Stationary Shock'};

%==========================================================================
% MAIN LOOP — one figure per test
%==========================================================================
for k = 1 : numel(Tests)

    Data = Tests{k};

    fprintf('\n========== Running test %d/%d: %s ==========\n', ...
            k, numel(Tests), Data.name);

    %----------------------------------------------------------------------
    % Run both schemes
    %----------------------------------------------------------------------
    Sol_G  = MainFVM(Data, 'godunov');
    Sol_LW = MainFVM(Data, 'laxwendroff');

    %----------------------------------------------------------------------
    % Figure: two subplots side by side
    %----------------------------------------------------------------------
    fig = figure('Units', 'normalized', 'Position', [0.05, 0.15, 0.88, 0.60]);

    % ---- Left panel: water height h(x, T) --------------------------------
    subplot(1, 2, 1);

    plot(Sol_G.x,  Sol_G.h,  '-',  'Color', c1, 'LineWidth', lw); hold on;
    plot(Sol_LW.x, Sol_LW.h, '--', 'Color', c2, 'LineWidth', lw);

    xlabel('x',          'FontSize', fs);
    ylabel('h(x, T)',    'FontSize', fs);
    title([test_titles{k}, ' — Water Height'], ...
          'FontSize', fst, 'FontWeight', 'bold');
    legend('Godunov (1st order)', 'Lax-Wendroff (2nd order)', ...
           'FontSize', fs, 'Location', 'best');
    grid on; box on;
    ax = gca; ax.FontSize = fs;

    % ---- Right panel: discharge q(x, T) ----------------------------------
    subplot(1, 2, 2);

    plot(Sol_G.x,  Sol_G.q,  '-',  'Color', c1, 'LineWidth', lw); hold on;
    plot(Sol_LW.x, Sol_LW.q, '--', 'Color', c2, 'LineWidth', lw);

    xlabel('x',          'FontSize', fs);
    ylabel('q(x, T)',    'FontSize', fs);
    title([test_titles{k}, ' — Discharge'], ...
          'FontSize', fst, 'FontWeight', 'bold');
    legend('Godunov (1st order)', 'Lax-Wendroff (2nd order)', ...
           'FontSize', fs, 'Location', 'best');
    grid on; box on;
    ax = gca; ax.FontSize = fs;

end
