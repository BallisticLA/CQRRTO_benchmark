%% plot_applications_results.m — Plot CQRRT_linop_composite_applications benchmark results
%  (Binary supports both sparse direct and composite FEM inputs; "applications" refers to
%   the three downstream tasks: generalized LS, generalized singular values, singular vectors.)
%
% Usage:
%   plot_applications_results(data_dir, results_csv, breakdown_csv)
%   plot_applications_results(data_dir, results_csv, breakdown_csv, plot_mode)
%   plot_applications_results(data_dir, results_csv, breakdown_csv, plot_mode, title_suffix)
%   plot_applications_results(data_dir, results_csv, breakdown_csv, plot_mode, title_suffix, main_tab, bd_tab)
%
% Arguments:
%   data_dir      — directory containing the CSV files
%   results_csv   — filename of *_gsvd_results.csv
%   breakdown_csv — filename of *_gsvd_breakdown.csv
%   plot_mode     — 'best_speed' (default), 'worst_ortho', or 'best_ortho'
%   title_suffix  — optional string appended to figure titles
%   main_tab      — optional uitab handle: plot main figure into this tab
%   bd_tab        — optional uitab handle: plot breakdown into this tab
%
% Subplot order: (1) Timing  (2) Orthogonality  (3) Error panel  (4) Memory
% Algorithm display order: CQRRT_linop, CQRRT_expl, CQRRT_linop_stb, CholQR, sCholQR3, sCholQR3_se
%   (CSV names: sCholQR3_basic -> displayed as sCholQR3; sCholQR3 -> displayed as sCholQR3_se)

function plot_applications_results(data_dir, results_csv, breakdown_csv, plot_mode, title_suffix, main_tab, bd_tab)

if nargin < 4, plot_mode = 'best_speed'; end
if nargin < 5, title_suffix = ''; end
if nargin < 6, main_tab = []; end
if nargin < 7, bd_tab = []; end

% =========================================================================
%  Wong colorblind-friendly palette
% =========================================================================
w_blue      = [  0 114 178] / 255;   % #0072B2
w_orange    = [230 159   0] / 255;   % #E69F00
w_skyblue   = [ 86 180 233] / 255;   % #56B4E9
w_green     = [  0 158 115] / 255;   % #009E73
w_vermilion = [213  94   0] / 255;   % #D55E00
w_purple    = [204 121 167] / 255;   % #CC79A7
w_gray      = [0.65 0.65 0.65];
w_ltgray    = [0.85 0.85 0.85];

% =========================================================================
%  Algorithm display order and name mapping
%  CSV name          -> display label (tex interpreter)
%  sCholQR3_basic    -> sCholQR3   (standard pseudocode)
%  sCholQR3          -> sCholQR3_se (storage-efficient)
% =========================================================================
alg_csv_order  = {'CQRRT_linop', 'CQRRT_expl', 'CQRRT_linop_stb', 'CQRRT_linop_stb_bqrrp', 'CholQR', 'sCholQR3_basic', 'sCholQR3'};
alg_disp_names = {'CQRRT\_linop', 'CQRRT\_expl', 'CQRRT\_linop\_stb', 'CQRRT\_linop\_stb\_bqrrp', 'CholQR', 'sCholQR3', 'sCholQR3\_se'};

% =========================================================================
%  Load results CSV
% =========================================================================
results_path   = fullfile(data_dir, results_csv);
breakdown_path = fullfile(data_dir, breakdown_csv);

n_comments = count_comment_lines(results_path);
opts = detectImportOptions(results_path, 'NumHeaderLines', n_comments);
T = readtable(results_path, opts);

algorithms = T.algorithm;
qr_time    = T.qr_time_us;
orth_error = T.orth_error;
total_a    = T.total_a_time_us;
total_b    = T.total_b_time_us;
total_c    = T.total_c_time_us;
peak_rss   = T.peak_rss_kb;
analytical = T.analytical_kb;
m_val      = T.m(1);
n_val      = T.n(1);


% =========================================================================
%  Sort algorithms by display order; build display labels
% =========================================================================
unique_algs_csv = unique(algorithms, 'stable');
unique_algs = {};
disp_labels = {};
for k = 1:numel(alg_csv_order)
    if any(strcmp(unique_algs_csv, alg_csv_order{k}))
        unique_algs{end+1} = alg_csv_order{k};
        disp_labels{end+1} = alg_disp_names{k};
    end
end
% Append any algorithm not in the ordered list (safety net)
for k = 1:numel(unique_algs_csv)
    if ~any(strcmp(unique_algs, unique_algs_csv{k}))
        unique_algs{end+1} = unique_algs_csv{k};
        disp_labels{end+1} = strrep(unique_algs_csv{k}, '_', '\_');
    end
end
n_algs = numel(unique_algs);

% Mode display string
switch plot_mode
    case 'best_speed',  mode_str = 'Best Speed';
    case 'worst_ortho', mode_str = 'Worst Ortho';
    otherwise,          mode_str = 'Best Ortho';
end

% =========================================================================
%  Select representative run per algorithm
% =========================================================================
sel_idx = zeros(n_algs, 1);
for a = 1:n_algs
    mask    = strcmp(algorithms, unique_algs{a});
    indices = find(mask);
    switch plot_mode
        case 'best_speed'
            [~, best] = min(total_a(mask));
            sel_idx(a) = indices(best);
        case 'worst_ortho'
            [~, worst] = max(orth_error(mask));
            sel_idx(a) = indices(worst);
        otherwise
            [~, best] = min(orth_error(mask));
            sel_idx(a) = indices(best);
    end
end

% =========================================================================
%  Main figure: timing | orthogonality | error panel | memory
% =========================================================================
if isempty(main_tab)
    figure('Position', [100, 100, 800, 900]);
    parent_main = gcf;
else
    parent_main = main_tab;
end
n_main_rows = 3;
tl_main = tiledlayout(parent_main, n_main_rows, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---- Subplot 1: Timing ----
nexttile(tl_main);
bar_data = zeros(n_algs, 4);
for a = 1:n_algs
    i = sel_idx(a);
    bar_data(a, :) = [qr_time(i), total_a(i), total_b(i), total_c(i)] / 1000;
end
b = bar(bar_data);
set(gca, 'XTick', 1:n_algs, 'XTickLabel', disp_labels, 'TickLabelInterpreter', 'tex', 'FontSize', 11);
ylabel('Time (ms)', 'FontSize', 12);
title(sprintf('Timing (%s)', mode_str), 'FontSize', 14, 'FontWeight', 'bold');
legend({'Q-less QR', 'App (a): Gen. LS', 'App (b): Gen. svals', 'App (c): Gen. svecs'}, ...
       'Location', 'northeastoutside', 'FontSize', 10);
grid on;
b(1).FaceColor = w_blue;
b(2).FaceColor = w_orange;
b(3).FaceColor = w_green;
b(4).FaceColor = w_purple;

% ---- Subplot 2: Orthogonality ----
ax_orth = nexttile(tl_main);
orth_vals = zeros(n_algs, 1);
for a = 1:n_algs
    orth_vals(a) = orth_error(sel_idx(a));
end
b_orth = bar(ax_orth, orth_vals);
b_orth.FaceColor = w_blue;
b_orth.FaceAlpha = 0.85;
y2_floor = 1e-16; if min(orth_vals) < 1e-16, y2_floor = min(orth_vals) / 10; end
set(ax_orth, 'YScale', 'log', 'YLim', [y2_floor, 1], 'YTick', [1e-15, 1e-10, 1e-5, 1], ...
    'XTick', 1:n_algs, 'XTickLabel', disp_labels, 'TickLabelInterpreter', 'tex', 'FontSize', 11);
ylabel(ax_orth, '$\|Q^TQ - I\|_F / \sqrt{n}$', 'Interpreter', 'latex', 'FontSize', 12);
title(ax_orth, sprintf('Orthogonality (%s)', mode_str), 'FontSize', 14, 'FontWeight', 'bold');
grid(ax_orth, 'on');

% ---- Subplot 3: Memory ----
nexttile(tl_main);
mem_data = zeros(n_algs, 2);
for a = 1:n_algs
    i = sel_idx(a);
    mem_data(a, :) = [peak_rss(i), analytical(i)] / 1024;
end
bm = bar(mem_data);
set(gca, 'XTick', 1:n_algs, 'XTickLabel', disp_labels, ...
    'TickLabelInterpreter', 'tex', 'FontSize', 11);
ylabel('Memory (MB)', 'FontSize', 12);
title('Peak Working Memory', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Peak RSS', 'Analytical'}, 'Location', 'northeastoutside', 'FontSize', 10);
grid on;
bm(1).FaceColor = w_green;
bm(2).FaceColor = w_orange;
ylim([0, max(mem_data(:)) * 1.15]);

% Supertitle
if isempty(title_suffix)
    main_title = sprintf('Applications Benchmark (%d x %d)', m_val, n_val);
else
    main_title = sprintf('Applications Benchmark (%d x %d) - %s', m_val, n_val, title_suffix);
end
if ~isempty(main_tab), main_tab.Title = main_title; end
sgtitle(tl_main, main_title, 'FontSize', 15, 'FontWeight', 'bold');

% =========================================================================
%  Breakdown figure
% =========================================================================
if isempty(breakdown_csv), return; end

n_comments_bd = count_comment_lines(breakdown_path);
bd_opts = detectImportOptions(breakdown_path, 'NumHeaderLines', n_comments_bd);
T_bd = readtable(breakdown_path, bd_opts);
bd_algs  = T_bd.algorithm;
% Name-based selection (2026-08-27): every t* column, in order. The old
% positional `T_bd{:, 5:end}` is the read pattern that silently grabbed the
% wrong column once before (07-31 awk incident class); a schema change ahead
% of the t-block would have shifted every slice here without an error.
t_names  = T_bd.Properties.VariableNames(startsWith(T_bd.Properties.VariableNames, 't'));
bd_times = T_bd{:, t_names};

% Operation-to-color mapping (same operation = same color across all subplots):
%
%   CQRRT variants (linop, expl, linop_stb):
%     gray       Alloc / Copy (M→dense)
%     blue       Sketch (SASO)
%     skyblue    QR on sketch + R_sk inversion
%     green      Fwd / Precond — apply R_sk^{-1} to M
%     orange     Gram — build preconditioned Gram (Adj+TRMM merged for linop; SYRK for expl)
%     purple     Cholesky of Gram
%     vermilion  Finalize — final TRMM R := R_chol * R_sk  (same op across all three)
%     ltgray     Other
%
%   sCholQR variants:
%     gray       Alloc
%     green      Fwd (linop forward passes)
%     orange     Adj (linop adjoint passes)
%     blue       Syrk / Gemm2 (rank-k Gram update, passes 2 & 3)
%     purple     Cholesky
%     vermilion  TRSM+Upd (sCholQR3_basic: trsm1 + upd; sCholQR3_se: updates)
%     skyblue    Q Mat (final Q materialization)
%     ltgray     Other

% Consistent color semantics across all CQRRT variants:
%   Alloc/Copy  gray       — buffer allocation or M copy
%   Sketch      blue       — SASO application
%   QR+inv      skyblue    — QR on sketch + R_sk inversion
%   Fwd/Precond green      — applying R_sk^{-1} to M (linop: forward passes; expl: TRSM in-place)
%   Gram        orange     — building preconditioned Gram matrix
%                            (linop: Adj passes + TRMM scaling merged; expl: single SYRK)
%   Chol        purple     — Cholesky of Gram
%   Finalize    vermilion  — final TRMM R := R_chol * R_sk (same op in all three)
%   Other        light gray — unaccounted overhead
bd_defs = struct( ...
    'CQRRT_linop', struct('total_col', 11, 'disp', 'CQRRT\_linop', ...
        'groups', {{
            {'Alloc',    [1],    w_gray}
            {'Sketch',   [2],    w_blue}
            {'QR+Tri.Inv',[3 4], w_skyblue}
            {'Fwd',      [5],    w_green}      % forward linop passes: apply R_sk^{-1} to M
            {'Gram',     [6 7],  w_orange}     % Adj passes + TRMM_gram = full preconditioned Gram
            {'Chol',     [8],    w_purple}
            {'Finalize', [9],    w_vermilion}  % TRMM: R := R_chol * R_sk (undo preconditioning)
            {'Other',     [10],   w_ltgray}
        }}), ...
    'CQRRT_expl', struct('total_col', 10, 'disp', 'CQRRT\_expl', ...
        'groups', {{
            {'Copy (M→dense)', [3],  w_gray}      % memcpy M_materialized (slot [2], repurposed trtri)
            {'Sketch',         [1],  w_blue}       % saso
            {'QR',             [2],  w_skyblue}    % QR on sketch
            {'Precond',        [4],  w_green}      % TRSM: M := M * R_sk^{-1} (apply R_sk^{-1} to M)
            {'Gram',           [5],  w_orange}     % SYRK: G = M^T M (preconditioned Gram)
            {'Chol',           [7],  w_purple}
            {'Finalize',       [8],  w_vermilion}  % TRMM: R := R_chol * R_sk (undo preconditioning)
            {'Other',           [9],  w_ltgray}
        }}), ...
    'CQRRT_linop_stb', struct('total_col', 11, 'disp', 'CQRRT\_linop\_stb', ...
        'groups', {{
            {'Alloc',    [1],    w_gray}
            {'Sketch',   [2],    w_blue}
            {'QR+Tri.Inv',[3 4], w_skyblue}
            {'Fwd',      [5],    w_green}      % forward linop passes: apply R_sk^{-1} to M
            {'Gram',     [6 7],  w_orange}     % Adj passes + GEMM_gram = full preconditioned Gram
            {'Chol',     [8],    w_purple}
            {'Finalize', [9],    w_vermilion}  % TRMM: R := R_chol * R_sk (undo preconditioning)
            {'Other',     [10],   w_ltgray}
        }}), ...
    'CQRRT_linop_stb_bqrrp', struct('total_col', 11, 'disp', 'CQRRT\_linop\_stb\_bqrrp', ...
        'groups', {{
            {'Alloc',    [1],    w_gray}
            {'Sketch',   [2],    w_blue}
            {'QR+BQRRP', [3 4],  w_skyblue}   % QR on sketch + BQRRP preconditioning
            {'Fwd',      [5],    w_green}
            {'Gram',     [6 7],  w_orange}
            {'Chol',     [8],    w_purple}
            {'Finalize', [9],    w_vermilion}
            {'Other',    [10],   w_ltgray}
        }}), ...
    'CholQR', struct('total_col', 6, 'disp', 'CholQR', ...
        'groups', {{
            {'Alloc', [1],   w_gray}
            {'Fwd',   [2],   w_green}
            {'Adj',   [3],   w_orange}
            {'Chol',  [4],   w_purple}
            {'Other',  [5],   w_ltgray}
        }}), ...
    'sCholQR3_basic', struct('total_col', 15, 'disp', 'sCholQR3', ...
        'groups', {{
            {'Alloc',    [1],       w_gray}
            {'Fwd',      [2 6],    w_green}      % fwd1, fwd_q
            {'Adj',      [3],      w_orange}     % adj1
            {'Syrk',     [7 10],   w_blue}       % syrk2, syrk3
            {'Chol',     [4 8 11], w_purple}     % chol1, chol2, chol3
            {'TRSM+Upd', [5 9 12], w_vermilion}  % trsm1 (Q from pass-1 chol), upd2, upd3
            {'Q Mat',    [13],     w_skyblue}
            {'Other',     [14],     w_ltgray}
        }}), ...
    'sCholQR3', struct('total_col', 18, 'disp', 'sCholQR3\_se', ...
        'groups', {{
            {'Alloc',  [1],          w_gray}
            {'Fwd',    [2 6 11],     w_green}
            {'Adj',    [3 7 12],     w_orange}
            {'Gemm',   [8 13],       w_blue}
            {'Chol',   [4 9 14],     w_purple}
            {'Update', [5 10 15],    w_vermilion}
            {'Q Mat',  [16],         w_skyblue}
            {'Other',   [17],         w_ltgray}
        }}) ...
);

% Build breakdown list in display order (only algorithms present in CSV)
bd_csv_order = {'CQRRT_linop', 'CQRRT_expl', 'CQRRT_linop_stb', 'CQRRT_linop_stb_bqrrp', 'CholQR', 'sCholQR3_basic', 'sCholQR3'};
bd_alg_list = {};
for k = 1:numel(bd_csv_order)
    alg = bd_csv_order{k};
    if isfield(bd_defs, alg) && any(strcmp(unique_algs, alg))
        bd_alg_list{end+1} = alg;
    end
end
n_bd_algs = numel(bd_alg_list);

% Compute unified y-axis limit
max_total_ms = 0;
for k = 1:n_bd_algs
    alg = bd_alg_list{k};
    def = bd_defs.(alg);
    a = find(strcmp(unique_algs, alg), 1);
    bi = find(strcmp(bd_algs, alg));
    sel_run  = T.run(sel_idx(a));
    bd_runs  = T_bd.run(bi);
    match    = find(bd_runs == sel_run, 1);
    if isempty(match), match = 1; end
    max_total_ms = max(max_total_ms, bd_times(bi(match), def.total_col) / 1000);
end
bd_ylim = [0, max_total_ms * 1.08];

% Split into row 1 (CQRRT algs) and row 2 (all others)
cqrrt_set = {'CQRRT_linop', 'CQRRT_expl', 'CQRRT_linop_stb', 'CQRRT_linop_stb_bqrrp'};
bd_row1 = {};
bd_row2 = {};
for k = 1:numel(bd_alg_list)
    if any(strcmp(cqrrt_set, bd_alg_list{k}))
        bd_row1{end+1} = bd_alg_list{k};
    else
        bd_row2{end+1} = bd_alg_list{k};
    end
end
n_row1 = numel(bd_row1);
n_row2 = numel(bd_row2);
n_cols = max(n_row1, n_row2);

if isempty(bd_tab)
    fig_bd    = figure('Position', [150, 150, 350 * n_cols, 800]);
    parent_bd = fig_bd;
else
    fig_bd    = [];
    parent_bd = bd_tab;
end

% Single grid so every tile is the same size.
% Shorter top row is centered by offsetting into blank leading tiles.
tl_bd       = tiledlayout(parent_bd, 2, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
row1_offset = floor((n_cols - n_row1) / 2);

% --- Row 1: CQRRT algorithms ---
for k = 1:n_row1
    alg  = bd_row1{k};
    def  = bd_defs.(alg);
    grps = def.groups;

    a  = find(strcmp(unique_algs, alg), 1);
    bi = find(strcmp(bd_algs, alg));
    sel_run = T.run(sel_idx(a));
    bd_runs = T_bd.run(bi);
    match = find(bd_runs == sel_run, 1);
    if isempty(match), match = 1; end
    row = bi(match);
    times_us = bd_times(row, :);

    n_grps     = numel(grps);
    grp_vals   = zeros(1, n_grps);
    grp_labels = cell(1, n_grps);
    grp_colors = cell(1, n_grps);
    for g = 1:n_grps
        grp_labels{g} = grps{g}{1};
        grp_colors{g} = grps{g}{3};
        grp_vals(g)   = sum(times_us(grps{g}{2})) / 1000;
    end

    nonzero   = grp_vals > 0;
    vals_nz   = grp_vals(nonzero);
    labels_nz = grp_labels(nonzero);
    colors_nz = grp_colors(nonzero);

    nexttile(tl_bd, row1_offset + k);
    bplot = bar(1, fliplr(vals_nz), 'stacked');
    for i = 1:numel(bplot)
        bplot(i).FaceColor = colors_nz{numel(bplot) + 1 - i};
        bplot(i).FaceAlpha = 0.9;
    end
    ylim(bd_ylim);
    if k == 1, ylabel('Time (ms)', 'FontSize', 10); end
    xlabel(def.disp, 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
    lgd = legend(flip(bplot), labels_nz);
    lgd.FontSize = 7;
    lgd.Location = 'northeastoutside';
    lgd.Interpreter = 'none';
    set(gca, 'XTick', []);
    grid on;
    set(gca, 'FontSize', 10);
end

% --- Row 2: non-CQRRT algorithms ---
for k = 1:n_row2
    alg  = bd_row2{k};
    def  = bd_defs.(alg);
    grps = def.groups;

    a  = find(strcmp(unique_algs, alg), 1);
    bi = find(strcmp(bd_algs, alg));
    sel_run = T.run(sel_idx(a));
    bd_runs = T_bd.run(bi);
    match = find(bd_runs == sel_run, 1);
    if isempty(match), match = 1; end
    row = bi(match);
    times_us = bd_times(row, :);

    n_grps     = numel(grps);
    grp_vals   = zeros(1, n_grps);
    grp_labels = cell(1, n_grps);
    grp_colors = cell(1, n_grps);
    for g = 1:n_grps
        grp_labels{g} = grps{g}{1};
        grp_colors{g} = grps{g}{3};
        grp_vals(g)   = sum(times_us(grps{g}{2})) / 1000;
    end

    nonzero   = grp_vals > 0;
    vals_nz   = grp_vals(nonzero);
    labels_nz = grp_labels(nonzero);
    colors_nz = grp_colors(nonzero);

    nexttile(tl_bd, n_cols + k);
    bplot = bar(1, fliplr(vals_nz), 'stacked');
    for i = 1:numel(bplot)
        bplot(i).FaceColor = colors_nz{numel(bplot) + 1 - i};
        bplot(i).FaceAlpha = 0.9;
    end
    ylim(bd_ylim);
    if k == 1, ylabel('Time (ms)', 'FontSize', 10); end
    xlabel(def.disp, 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
    lgd = legend(flip(bplot), labels_nz);
    lgd.FontSize = 7;
    lgd.Location = 'northeastoutside';
    lgd.Interpreter = 'none';
    set(gca, 'XTick', []);
    grid on;
    set(gca, 'FontSize', 10);
end

if isempty(title_suffix)
    bd_title = sprintf('QR Runtime Breakdown (%d x %d, %s)', m_val, n_val, mode_str);
else
    bd_title = sprintf('QR Runtime Breakdown (%d x %d, %s) - %s', m_val, n_val, mode_str, title_suffix);
end
if ~isempty(bd_tab)
    bd_tab.Title = bd_title;
end
sgtitle(tl_bd, bd_title, 'FontSize', 13, 'FontWeight', 'bold');

fprintf('Plots generated for %s\n', results_path);
fprintf('  Mode: %s, Algorithms: %s, Matrix: %d x %d\n', ...
        mode_str, strjoin(disp_labels, ', '), m_val, n_val);
end

