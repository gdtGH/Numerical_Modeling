% RunMainFVM.m
% Top-level script for HW3 — Nonlinear SWE with FVM.
%
% Runs all three Riemann test cases (Dam Break, Reverse Dam Break,
% Stationary Shock) with both the Godunov and Lax-Wendroff schemes.
% Produces one figure per test case, each with two subplots:
%   Left panel:  water height h(x, T)
%   Right panel: discharge    q(x, T)
%
% Figures are saved as .png to out_dir following the naming convention:
%   X.Y_Description.png   (LaTeX report image convention)

clear; close all; clc;

addpath Solvers

%==========================================================================
% OUTPUT DIRECTORY
%==========================================================================
out_dir = 'Output';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%==========================================================================
% GLOBAL STYLE — edit ONLY here
%==========================================================================
fs  = 20;           % axis/tick/legend font size
fst = 28;           % sgtitle font size
fsb = 22;           % subplot title font size
lw  = 2;            % line width
c1  = '#0072BD';    % Godunov       (blue)
c2  = '#D95319';    % Lax-Wendroff  (orange)

%==========================================================================
% LOAD TEST CASES
%==========================================================================
Tests = DataTest();

%==========================================================================
% FIGURE NAMES AND TITLES
%==========================================================================
fig_names   = {'3.1_DamBreak', '3.2_RevDamBreak', '3.3_Shock'};
test_titles = {'Dam Break',    'Reverse Dam Break', 'Stationary Shock'};

%==========================================================================
% MAIN LOOP — one figure per test case
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
    % Figure
    %----------------------------------------------------------------------
    figure('Name', fig_names{k}, 'NumberTitle', 'off', ...
           'Units', 'normalized', 'Position', [0.05, 0.15, 0.88, 0.60]);

    % Shared super-title for the whole figure
    sgtitle(test_titles{k}, 'FontSize', fst, 'FontWeight', 'bold');

    % ---- Left panel: water height h(x, T) --------------------------------
    subplot(1, 2, 1);

    plot(Sol_G.x,  Sol_G.h,  '-',  'Color', c1, 'LineWidth', lw); hold on;
    plot(Sol_LW.x, Sol_LW.h, '--', 'Color', c2, 'LineWidth', lw);

    xlabel('x',         'FontSize', fs);
    ylabel('h(x, T)',   'FontSize', fs);
    title('Water height  h', 'FontSize', fsb, 'FontWeight', 'bold');
    legend('Godunov (1^{st} order)', 'Lax-Wendroff (2^{nd} order)', ...
           'FontSize', fs, 'Location', 'best');
    grid on; box on;
    ax = gca; ax.FontSize = fs;

    % ---- Right panel: discharge q(x, T) ----------------------------------
    subplot(1, 2, 2);

    plot(Sol_G.x,  Sol_G.q,  '-',  'Color', c1, 'LineWidth', lw); hold on;
    plot(Sol_LW.x, Sol_LW.q, '--', 'Color', c2, 'LineWidth', lw);

    xlabel('x',         'FontSize', fs);
    ylabel('q(x, T)',   'FontSize', fs);
    title('Discharge  q', 'FontSize', fsb, 'FontWeight', 'bold');
    legend('Godunov (1^{st} order)', 'Lax-Wendroff (2^{nd} order)', ...
           'FontSize', fs, 'Location', 'best');
    grid on; box on;
    ax = gca; ax.FontSize = fs;

    %----------------------------------------------------------------------
    % Save figure
    %----------------------------------------------------------------------
    exportgraphics(gcf, fullfile(out_dir, [get(gcf, 'Name') '.png']), ...
                   'Resolution', 150);
    fprintf('Saved: %s.png\n', fig_names{k});

end

fprintf('\nAll plots generated and saved to: %s\n', out_dir);