% RunMainFVM.m — HW3, nonlinear SWE with FVM.
% PART 1: Godunov vs Lax-Wendroff (standard flux) — Dam Break, Reverse Dam Break, Stationary Shock.
% PART 2: Standard vs Modified Godunov (logarithmic pressure) — Dam Break, Stationary Shock.
% Each test produces two single-panel figures: <name>_h.png and <name>_q.png.

clear; close all; clc;

%% Setup
addpath Solvers

out_dir = 'Output';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

fs = 25;            % axis/legend font size
lw = 2;             % line width
c1 = '#0072BD';     % Godunov / Standard  (blue,   solid)
c2 = '#D95319';     % LW / Modified       (orange, dashed)

Tests = DataTest();

%% PART 1 — Godunov vs Lax-Wendroff
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

    lw_legend = 'Lax-Wendroff';
    if Sol_LW.diverged
        sub_txt = sprintf('First-order Godunov vs pure Lax-Wendroff (last stable t = %.3f s)', Sol_LW.t);
    else
        sub_txt = 'First-order Godunov vs pure Second-order Lax-Wendroff';
    end

    % h(x,T)
    figure('Name', [fig_names_1{k} '_h'], 'NumberTitle', 'off');
    set(gcf, 'Units', 'pixels', 'WindowState', 'maximized');
    plot(Sol_G.x,  Sol_G.h,  '-',  'Color', c1, 'LineWidth', lw); hold on;
    plot(Sol_LW.x, Sol_LW.h, '--', 'Color', c2, 'LineWidth', lw);
    xlabel('x',       'FontSize', fs);
    ylabel('h(x, T)', 'FontSize', fs);
    if exist('subtitle','file') || ~verLessThan('matlab','9.9')   % R2020b+
        title([test_titles_1{k} ' — Water height  h(x,T)'], 'FontSize', fs+3, 'FontWeight','bold');
        subtitle(sub_txt, 'FontSize', fs-3);
    else
        title({[test_titles_1{k} ' — Water height  h(x,T)'], ...
               ['\fontsize{' num2str(fs-3) '}' sub_txt]}, ...
              'FontSize', fs+3, 'FontWeight','bold');
    end
    legend('Godunov', lw_legend, 'FontSize', fs, 'Location', 'best');
    grid on; box on; ax = gca; ax.FontSize = fs;
    exportgraphics(gcf, fullfile(out_dir, [fig_names_1{k} '_h.png']), 'Resolution', 200);
    fprintf('Saved: %s_h.png\n', fig_names_1{k});

    % q(x,T)
    figure('Name', [fig_names_1{k} '_q'], 'NumberTitle', 'off');
    set(gcf, 'Units', 'pixels', 'WindowState', 'maximized');
    plot(Sol_G.x,  Sol_G.q,  '-',  'Color', c1, 'LineWidth', lw); hold on;
    plot(Sol_LW.x, Sol_LW.q, '--', 'Color', c2, 'LineWidth', lw);
    xlabel('x',       'FontSize', fs);
    ylabel('q(x, T)', 'FontSize', fs);
    if exist('subtitle','file') || ~verLessThan('matlab','9.9')   % R2020b+
        title([test_titles_1{k} ' — Discharge  q(x,T)'], 'FontSize', fs+3, 'FontWeight','bold');
        subtitle(sub_txt, 'FontSize', fs-3);
    else
        title({[test_titles_1{k} ' — Discharge  q(x,T)'], ...
               ['\fontsize{' num2str(fs-3) '}' sub_txt]}, ...
              'FontSize', fs+3, 'FontWeight','bold');
    end
    legend('Godunov', lw_legend, 'FontSize', fs, 'Location', 'best');
    grid on; box on; ax = gca; ax.FontSize = fs;
    exportgraphics(gcf, fullfile(out_dir, [fig_names_1{k} '_q.png']), 'Resolution', 200);
    fprintf('Saved: %s_q.png\n', fig_names_1{k});

end

%% PART 2 — Standard vs Modified Godunov
h_ref = 1;          % reference height for logarithmic pressure

fig_names_2   = {'3.4_Modified_DamBreak', '3.5_Modified_Shock'};
test_titles_2 = {'Dam Break',             'Stationary Shock'};
test_idx      = [1, 3];

fprintf('\n%s\n', repmat('=',1,60));
fprintf('  PART 2 — Standard vs Modified Godunov (h_ref = %g)\n', h_ref);
fprintf('%s\n', repmat('=',1,60));

for k = 1 : numel(test_idx)

    Data = Tests{test_idx(k)};
    fprintf('\n--- Test %d/%d: %s ---\n', k, numel(test_idx), Data.name);

    Sol_std = MainFVM(Data, 'godunov');
    Sol_mod = MainFVM(Data, 'modified_godunov', h_ref);

    sub_txt    = sprintf('Standard vs modified-pressure Godunov  (h_{ref} = %g)', h_ref);
    mod_legend = 'Modified';

    % h(x,T)
    figure('Name', [fig_names_2{k} '_h'], 'NumberTitle', 'off');
    set(gcf, 'Units', 'pixels', 'WindowState', 'maximized');
    plot(Sol_std.x, Sol_std.h, '-',  'Color', c1, 'LineWidth', lw); hold on;
    plot(Sol_mod.x, Sol_mod.h, '--', 'Color', c2, 'LineWidth', lw);
    xlabel('x',       'FontSize', fs);
    ylabel('h(x, T)', 'FontSize', fs);
    if exist('subtitle','file') || ~verLessThan('matlab','9.9')   % R2020b+
        title([test_titles_2{k} ' — Water height  h(x,T)'], 'FontSize', fs+3, 'FontWeight','bold');
        subtitle(sub_txt, 'FontSize', fs-3);
    else
        title({[test_titles_2{k} ' — Water height  h(x,T)'], ...
               ['\fontsize{' num2str(fs-3) '}' sub_txt]}, ...
              'FontSize', fs+3, 'FontWeight','bold');
    end
    legend('Standard', mod_legend, 'FontSize', fs, 'Location', 'best');
    grid on; box on; ax = gca; ax.FontSize = fs;
    exportgraphics(gcf, fullfile(out_dir, [fig_names_2{k} '_h.png']), 'Resolution', 200);
    fprintf('Saved: %s_h.png\n', fig_names_2{k});

    % q(x,T)
    figure('Name', [fig_names_2{k} '_q'], 'NumberTitle', 'off');
    set(gcf, 'Units', 'pixels', 'WindowState', 'maximized');
    plot(Sol_std.x, Sol_std.q, '-',  'Color', c1, 'LineWidth', lw); hold on;
    plot(Sol_mod.x, Sol_mod.q, '--', 'Color', c2, 'LineWidth', lw);
    xlabel('x',       'FontSize', fs);
    ylabel('q(x, T)', 'FontSize', fs);
    if exist('subtitle','file') || ~verLessThan('matlab','9.9')   % R2020b+
        title([test_titles_2{k} ' — Discharge  q(x,T)'], 'FontSize', fs+3, 'FontWeight','bold');
        subtitle(sub_txt, 'FontSize', fs-3);
    else
        title({[test_titles_2{k} ' — Discharge  q(x,T)'], ...
               ['\fontsize{' num2str(fs-3) '}' sub_txt]}, ...
              'FontSize', fs+3, 'FontWeight','bold');
    end
    legend('Standard', mod_legend, 'FontSize', fs, 'Location', 'best');
    grid on; box on; ax = gca; ax.FontSize = fs;
    exportgraphics(gcf, fullfile(out_dir, [fig_names_2{k} '_q.png']), 'Resolution', 200);
    fprintf('Saved: %s_q.png\n', fig_names_2{k});

end

fprintf('\nAll plots generated and saved to: %s\n', out_dir);
