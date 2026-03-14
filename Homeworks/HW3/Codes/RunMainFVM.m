% RunMainFVM.m
% Top-level script for HW3 — Nonlinear SWE with FVM.
%
% PART 1 — Godunov vs Lax-Wendroff (standard flux)
%   All three Riemann tests: Dam Break, Reverse Dam Break, Stationary Shock.
%   One figure per test, two subplots: h(x,T) and q(x,T).
%   Figures: 3.1_DamBreak.png  3.2_RevDamBreak.png  3.3_Shock.png
%
% PART 2 — Standard Godunov vs Modified Godunov (logarithmic pressure flux)
%   Dam Break and Stationary Shock only.
%   One figure per test, two subplots: h(x,T) and q(x,T).
%   Figures: 3.4_Modified_DamBreak.png  3.5_Modified_Shock.png

clear; close all; clc;

%% ========================================================================
% Inizializing
% =========================================================================

addpath Solvers

% OUTPUT DIRECTORY
out_dir = 'Output';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% GLOBAL STYLE — edit ONLY here
fs  = 25;           % axis/tick/legend font size
fst = 50;           % sgtitle font size
fsb = 32;           % subplot title font size
lw  = 2;            % line width
c1  = '#0072BD';    % primary   — Godunov / Standard  (blue,   solid)
c2  = '#D95319';    % secondary — LW / Modified       (orange, dashed)

% LOAD ALL TEST CASES
Tests = DataTest();

%% =========================================================================
% Point 2 — Standard flux: Godunov vs Lax-Wendroff
% =========================================================================

fig_names_1   = {'3.1_DamBreak', '3.2_RevDamBreak', '3.3_Shock'};
test_titles_1 = {'Dam Break',    'Reverse Dam Break', 'Stationary Shock'};

fprintf('\n%s\n', repmat('=',1,60));
fprintf('  PART 1 — Godunov vs Lax-Wendroff (standard flux)\n');
fprintf('%s\n', repmat('=',1,60));

for k = 1 : numel(Tests)

    Data = Tests{k};

    fprintf('\n--- Test %d/%d: %s ---\n', k, numel(Tests), Data.name);

    Sol_G  = MainFVM(Data, 'godunov');
    Sol_LW = MainFVM(Data, 'laxwendroff');

    figure('Name', fig_names_1{k}, 'NumberTitle', 'off', ...
           'Units', 'normalized', 'Position', [0.05, 0.15, 0.88, 0.60]);

    sgtitle(test_titles_1{k}, 'FontSize', fst, 'FontWeight', 'bold');

    % Left panel — h
    subplot(1, 2, 1);
    plot(Sol_G.x,  Sol_G.h,  '-',  'Color', c1, 'LineWidth', lw); hold on;
    plot(Sol_LW.x, Sol_LW.h, '--', 'Color', c2, 'LineWidth', lw);
    xlabel('x',         'FontSize', fs);
    ylabel('h(x, T)',   'FontSize', fs);
    title('Water height  h', 'FontSize', fsb, 'FontWeight', 'bold');
    legend('Godunov (1^{st} order)', 'Lax-Wendroff (2^{nd} order)', ...
           'FontSize', fs, 'Location', 'best');
    grid on; box on; ax = gca; ax.FontSize = fs;

    % Right panel — q
    subplot(1, 2, 2);
    plot(Sol_G.x,  Sol_G.q,  '-',  'Color', c1, 'LineWidth', lw); hold on;
    plot(Sol_LW.x, Sol_LW.q, '--', 'Color', c2, 'LineWidth', lw);
    xlabel('x',         'FontSize', fs);
    ylabel('q(x, T)',   'FontSize', fs);
    title('Discharge  q', 'FontSize', fsb, 'FontWeight', 'bold');
    legend('Godunov (1^{st} order)', 'Lax-Wendroff (2^{nd} order)', ...
           'FontSize', fs, 'Location', 'best');
    grid on; box on; ax = gca; ax.FontSize = fs;

    exportgraphics(gcf, fullfile(out_dir, [get(gcf,'Name') '.png']), ...
                   'Resolution', 150);
    fprintf('Saved: %s.png\n', fig_names_1{k});

end

%% =========================================================================
% Point 3 — Modified pressure flux: Standard vs Modified Godunov
% =========================================================================

h_ref = 1;          % reference water height for logarithmic pressure

fig_names_2   = {'3.4_Modified_DamBreak', '3.5_Modified_Shock'};
test_titles_2 = {'Dam Break',             'Stationary Shock'};
test_idx      = [1, 3];   % indices in Tests cell array

fprintf('\n%s\n', repmat('=',1,60));
fprintf('  PART 2 — Standard vs Modified Godunov (h_ref = %g)\n', h_ref);
fprintf('%s\n', repmat('=',1,60));

for k = 1 : numel(test_idx)

    Data = Tests{test_idx(k)};

    fprintf('\n--- Test %d/%d: %s ---\n', k, numel(test_idx), Data.name);

    Sol_std = MainFVM(Data, 'godunov');
    Sol_mod = MainFVM(Data, 'modified_godunov', h_ref);

    figure('Name', fig_names_2{k}, 'NumberTitle', 'off', ...
           'Units', 'normalized', 'Position', [0.05, 0.15, 0.88, 0.60]);

    sgtitle([test_titles_2{k}, '  —  Standard vs Modified Pressure'], ...
            'FontSize', fst, 'FontWeight', 'bold');

    % Left panel — h
    subplot(1, 2, 1);
    plot(Sol_std.x, Sol_std.h, '-',  'Color', c1, 'LineWidth', lw); hold on;
    plot(Sol_mod.x, Sol_mod.h, '--', 'Color', c2, 'LineWidth', lw);
    xlabel('x',       'FontSize', fs);
    ylabel('h(x, T)', 'FontSize', fs);
    title('Water height  h', 'FontSize', fsb, 'FontWeight', 'bold');
    legend('Standard Godunov', ...
           sprintf('Modified Godunov  (h_{ref} = %g)', h_ref), ...
           'FontSize', fs, 'Location', 'best');
    grid on; box on; ax = gca; ax.FontSize = fs;

    % Right panel — q
    subplot(1, 2, 2);
    plot(Sol_std.x, Sol_std.q, '-',  'Color', c1, 'LineWidth', lw); hold on;
    plot(Sol_mod.x, Sol_mod.q, '--', 'Color', c2, 'LineWidth', lw);
    xlabel('x',       'FontSize', fs);
    ylabel('q(x, T)', 'FontSize', fs);
    title('Discharge  q', 'FontSize', fsb, 'FontWeight', 'bold');
    legend('Standard Godunov', ...
           sprintf('Modified Godunov  (h_{ref} = %g)', h_ref), ...
           'FontSize', fs, 'Location', 'best');
    grid on; box on; ax = gca; ax.FontSize = fs;

    exportgraphics(gcf, fullfile(out_dir, [get(gcf,'Name') '.png']), ...
                   'Resolution', 150);
    fprintf('Saved: %s.png\n', fig_names_2{k});

end

fprintf('\nAll plots generated and saved to: %s\n', out_dir);