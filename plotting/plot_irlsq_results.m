%% plot_irlsq_results.m — Plot CQRRT_linop_irlsq benchmark results
%
% The sparse-input IR-LSQ benchmark runs every Q-less QR variant selected by
% method_mask through the IR-LSQ pipeline (Algorithm 1, Epperly–Meier–Nakatsukasa
% 2025) on a tall sparse matrix loaded from a .mtx file. Each (algorithm, run)
% does Q-less QR → restarted IR (x_0 = 0; up to ir_n_steps rounds, outer_tol and
% the LS-floor exit end it early) with inner CG. (The "2-step IR" wording that
% used to live here described the pre-2026-08-07 fixed-round scheme.)
%
% CSV schema (results; all reads are NAME-based, this comment is documentation):
%   algorithm, run, m, n, qr_status, qr_time_us, peak_rss_kb, analytical_kb,
%   orth_error, ir_total_us, ir_outer_iters, ir_inner_iters_total,
%   ls_residual_norm, ls_solution_error, [reg schema: kappa/mu/precision columns,]
%   ir_inner_capped, ir_inner_relres, ir_inner_best_relres, ir_inner_best_iter,
%   cond_precond, ir_setup_us, lsqr_iters, engine_status, stop_reason
% (the last three columns and a REAL ir_setup_us arrived 2026-08-27; a
%  *_rounds.csv sidecar carries per-round engine records.)
%
% ls_residual_norm uses the Higham normwise backward-error metric
%   ||A x - b|| / (||A||_2 * ||x|| + ||b||),
% drivable to machine epsilon for a backward-stable LS solver.
%
% ls_solution_error is -1 / NaN when no ground-truth x_true exists (FEM b = L^{-1} r).
%
% CSV schema (breakdown):
%   algorithm, run, phase, t0..t10
%   phase in {QR, IR}.  IR layout (6): outer_total, inner_cg_total, trsm, fwd, adj, other.
%   Since 2026-08-27 Blendenpik-family rows carry real IR entries too (refine rows:
%   engine split; published rows: the LSQR op split with a 0 inner_cg slot).
%
% Usage:
%   plot_irlsq_results(data_dir, results_csv, breakdown_csv)
%   plot_irlsq_results(data_dir, results_csv, breakdown_csv, title_suffix)
%   plot_irlsq_results(data_dir, results_csv, breakdown_csv, title_suffix, main_tab, bd_tab)
%   plot_irlsq_results(..., main_tab, bd_tab, timing_agg)
%
% timing_agg ('best' default | 'mean') controls how multi-run CSVs (num_runs > 1,
% 2026-08-05) collapse to one row per method; see aggregate_runs.m.

function plot_irlsq_results(data_dir, results_csv, breakdown_csv, title_suffix, main_tab, bd_tab, timing_agg)

if nargin < 4, title_suffix = ''; end
if nargin < 5, main_tab = []; end
if nargin < 6, bd_tab = []; end
if nargin < 7 || isempty(timing_agg), timing_agg = 'best'; end

% =========================================================================
%  Wong colorblind-friendly palette (matches plot_applications_results.m)
% =========================================================================
w_blue      = [  0 114 178] / 255;
w_orange    = [230 159   0] / 255;
w_skyblue   = [ 86 180 233] / 255;
w_green     = [  0 158 115] / 255;
w_vermilion = [213  94   0] / 255;
w_purple    = [204 121 167] / 255;
w_gray      = [0.65 0.65 0.65];
w_ltgray    = [0.85 0.85 0.85];

% =========================================================================
%  Algorithm display order.  GEQP3-stabilized variant (CQRRT_linop_stb) is
%  intentionally absent: the benchmark dropped that path.
% =========================================================================
% "Blendenpik_refine"/"Blendenpik_cold_refine" appear in [NEW 08-09]+ CSVs: Blendenpik's
% OWN preconditioner and answer handed to our refinement solver, so the published rows and
% the like-for-like comparison sit side by side.
% "Blendenpik_cold" appears in [NEW 08-05]+ CSVs (warm x_0 is Blendenpik-only and
% both variants run). "Blendenpik" keeps its bare label; the era note carries policy.
% "unpreconditioned" (2026-09-15, Oleg's Figure 6 question): the refinement
% engine on the raw operator, no factor (mask bit 128, run as the 0915_up1
% overlay era and merged per cell by merge_overlay_rows). Last in the order, as
% on the Toeplitz figure; it has no QR phase and no orthogonality entry.
alg_csv_order  = {'CQRRTO_linop','CQRRTO_linop_gemmL','CQRRT_linop','CQRRT_linop_gemmL', 'CholQR', 'CholQR2', 'sCholQR3_basic', 'sCholQR3', 'Blendenpik', 'Blendenpik_cold', 'Blendenpik_refine', 'Blendenpik_cold_refine', 'unpreconditioned'};
% 2026-09-04 (Max): the word "refinement" is banned from figure text. The two
% plotted Blendenpik rows are BOTH the refined (shared-engine) variants, so they
% take the bare labels; the single-shot rows are excluded from plots and keep
% distinguishing display names in case they are ever re-enabled.
alg_disp_names = {'CQRRTO','CQRRTO (GEMM left)','CQRRTO','CQRRTO (GEMM left)', 'CholQR', 'CholQR2', 'sCholQR3 (no blocking)', 'sCholQR3', 'Blendenpik (single-shot)', 'Blendenpik (single-shot, zero x_0)', 'Blendenpik (sketch-and-solve x_0)', 'Blendenpik', 'unpreconditioned'};

% =========================================================================
%  Load CSVs.  count_comment_lines is a helper in this directory; we rely on
%  it to skip the leading "#"-prefixed metadata block.
% =========================================================================
results_path   = fullfile(data_dir, results_csv);
breakdown_path = fullfile(data_dir, breakdown_csv);
if ~isfile(results_path)
    error('plot_irlsq_results: file not found: %s', results_path);
end
n_skip = count_comment_lines(results_path);
opts = detectImportOptions(results_path, 'NumHeaderLines', n_skip);
% With very few data rows (the mask-limited diag cells carry only 2),
% detectImportOptions can fail to promote the column line to variable names
% and falls back to Var1..VarN. Setting VariableNamesLine alone is NOT enough:
% setvartype validates against the STALE VariableNames property. Read the
% column line ourselves and assign the names directly (2026-07-31, hit on the
% warmir/coldbp 2x2 cells).
if ~ismember('algorithm', opts.VariableNames)
    fid = fopen(results_path, 'r');
    for s = 1:n_skip, fgetl(fid); end
    names = strtrim(strsplit(fgetl(fid), ','));
    fclose(fid);
    assert(numel(names) == numel(opts.VariableNames), ...
        'plot_irlsq_results: %d header names vs %d detected columns in %s', ...
        numel(names), numel(opts.VariableNames), results_path);
    opts.VariableNames = names;
    opts.DataLines     = [n_skip + 2, Inf];
end
opts = setvartype(opts, 'algorithm', 'char');
T = readtable(results_path, opts);

% 2026-09-14 (Max): the accuracy panel is the data error ||Ax-b||/||b||, read
% from the per-round sidecar (last cycle of each run), so the FEM and Toeplitz
% figures share one metric with its floor at the data noise level. Joined
% BEFORE aggregate_runs so the kept run carries its own value.
T = attach_data_error(T, strrep(results_path, '_results.csv', '_rounds.csv'));

% Multi-run CSVs (num_runs > 1, 2026-08-05): collapse to one row per algorithm
% up front, so every column copy below sees the aggregated table. agg_note
% ('best of 5 runs' / 'mean of 5 runs') lands in the figure title.
[T, agg_note] = aggregate_runs(T, timing_agg);

% 2026-08-24 (Oleg): paper-figure roster trim -- PARTIALLY OVERRIDDEN 2026-09-04
% (Max): the two cold Blendenpik rows are back in the plots. The reason given at
% the time was that on a NOISE-FREE benchmark the warm x_0 (sketch-and-solve on a
% consistent system) starts AT the u*kappa accuracy floor, so the warm refine row
% showed no iterative work at all and only the cold rows exercised the solve.
%
% 2026-09-14: THAT PREMISE NO LONGER HOLDS for FEM2. The plotted FEM2 era is now
% [0911 n11], which carries noise_level=1e-11 on the right-hand side, so the warm
% start enters at the noise level rather than the floor (measured x0_relres
% 1.33e-11, against 3.6e-15 in the noise-free eras). Both warm rows now do real
% work in all three cells: Blendenpik_refine runs 3 outer / 17 inner iterations,
% and the plain warm "Blendenpik" row runs 1 outer / 80-82 inner and exits on
% ne_floor -- it no longer exits at 0 LSQR iterations, which was the stated
% ground for excluding it.
%
% RULING 2026-09-15 (Oleg): the sketch-and-solve-started row is removed from the
% paper figures altogether; only the zero-start row remains, labelled
% 'Blendenpik'. The single-shot rows stay excluded.
% NOTE this comment applies to FEM2 only; the Toeplitz benchmark is still on the
% noise-free [0902 ns1] era, so plot_toeplitz_results.m keeps the old reasoning.
% sCholQR3_basic (no-blocking control) stays excluded regardless.
% 2026-09-04 (Max): CQRRT_linop_gemmL (RANDLAPACK_GRAM_LEFT=gemm arm) is plot-excluded:
% the paper's PCholQR model applies the inverse preconditioner by triangular
% solves, which is the library-default TRSM arm; the GEMM-left arm stays in the
% merged CSVs for re-enable. Display name is CQRRTO, the paper's algorithm name.
% 2026-09-15 (Oleg): the sketch-and-solve-started row (Blendenpik_refine) is OUT
% of the paper figures; the zero-start row is the published algorithm and takes
% the bare label 'Blendenpik'. The warm-start wall-time slice disappears with it
% (the legend entry is dropped automatically when no plotted row has one).
alg_plot_exclude = {'sCholQR3_basic', 'Blendenpik', 'Blendenpik_cold', 'Blendenpik_refine', 'CQRRT_linop_gemmL', 'CQRRTO_linop_gemmL'};  % single-shot rows + warm-start row + GEMM-left arm excluded
T(ismember(T.algorithm, alg_plot_exclude), :) = [];

algorithms       = T.algorithm;
qr_status        = T.qr_status;
qr_time          = T.qr_time_us;
peak_rss_kb      = T.peak_rss_kb;
analytical_kb    = T.analytical_kb;
ir_total_us      = T.ir_total_us;
ir_inner_total   = T.ir_inner_iters_total;
ls_residual_norm = T.ls_residual_norm;
ls_data_error    = T.ls_data_error;      % ||Ax-b||/||b||, attach_data_error (2026-09-14)
ls_solution_err  = T.ls_solution_error;
orth_error       = T.orth_error;
m_val            = T.m(1);
n_val            = T.n(1);
% (The adaptive-shift retries and canonical-GEQRF-rate PANELS were removed
%  2026-08-24 per Oleg; since 2026-08-31 chol_retries drives the asterisk on
%  the method labels instead of a panel of its own.)

% =========================================================================
%  Sort algorithms by display order
% =========================================================================
unique_algs_csv = unique(algorithms, 'stable');
unique_algs = {};
disp_labels = {};
for k = 1:numel(alg_csv_order)
    if any(strcmp(unique_algs_csv, alg_csv_order{k}))
        unique_algs{end+1} = alg_csv_order{k}; %#ok<AGROW>
        disp_labels{end+1} = alg_disp_names{k}; %#ok<AGROW>
    end
end
% Catch any extras (defensive)
for k = 1:numel(unique_algs_csv)
    if ~any(strcmp(unique_algs, unique_algs_csv{k}))
        unique_algs{end+1} = unique_algs_csv{k}; %#ok<AGROW>
        disp_labels{end+1} = strrep(unique_algs_csv{k}, '_', '\_'); %#ok<AGROW>
    end
end
n_algs = numel(unique_algs);

% =========================================================================
%  Per-algorithm row lookup + failure flags. Multi-run selection already
%  happened in aggregate_runs (one row per algorithm at this point); this
%  block is kept as a defensive no-op reduction in case a table ever reaches
%  here unaggregated.
% =========================================================================
sel_idx     = zeros(n_algs, 1);
sel_failed  = false(n_algs, 1);
for a = 1:n_algs
    mask     = strcmp(algorithms, unique_algs{a});
    indices  = find(mask);
    succ     = qr_status(mask) == 0;
    if any(succ)
        succ_idx = indices(succ);
        total = qr_time(succ_idx) + ir_total_us(succ_idx);
        [~, best] = min(total);
        sel_idx(a) = succ_idx(best);
    else
        sel_idx(a) = indices(1);
        sel_failed(a) = true;
    end
end

% =========================================================================
%  Mark adaptive-shift rescues on the method labels (2026-08-31, Max). A
%  nominally unshifted method whose Cholesky broke down and was rescued by the
%  adaptive shift measures the SHIFT-RESCUED variant of its algorithm (the
%  cholqr_primitive contract caveat), and on the boundary cells WHICH branch
%  runs can flip with the BLAS thread count (native_ill, 2026-08-31 report).
%  The label carries an asterisk and a greyed footnote explains it; pre-08-27
%  CSVs lack the column, hence the guard.
% =========================================================================
has_retries = ismember('chol_retries', T.Properties.VariableNames);
rescued = false(n_algs, 1);
if has_retries
    for a = 1:n_algs
        if ~sel_failed(a) && T.chol_retries(sel_idx(a)) > 0
            rescued(a) = true;
            disp_labels{a} = [disp_labels{a} '*'];
        end
    end
end

% Backward-error termination (2026-09-18): when the era ran with the engine's
% oracle on (header echoes be_tol > 0), each method is judged from its RECORDED
% final estimate (be_final <= be_tol), not from the exit code. A method that
% stopped short of the target is CENSORED: its bars are washed out in every
% panel via the tick label (a dagger), since the paper exports each tile on its
% own; a live row that never ran the oracle gets a double dagger. Eras without
% the knob keep the old rendering. See be_censoring.m.
BE = be_censoring(T, sel_idx, sel_failed, results_path);
for a = 1:n_algs
    if BE.censored(a),  disp_labels{a} = [disp_labels{a} '^{\dagger}']; end
    if BE.no_oracle(a), disp_labels{a} = [disp_labels{a} '^{\ddagger}']; end
end

% =========================================================================
%  Main figure: 1x4 layout (2026-09-07, per Oleg: one line, same footprint as
%  the synthetic-experiment panels; the peak-memory panel is dropped from the
%  paper figures, its numbers are stated in the text)
%    (1) Stacked timing QR+IR        (2) Inner CG iterations
%    (3) Q-factor orthogonality loss (4) Normwise backward error
%  Optional tiles, appended to the right when enabled: memory (SHOW_MEMORY),
%  LS forward error (SHOW_FORWARD). Per-panel paper export lives in run_all.m.
%  (QR-build canonical rate + adaptive-shift retries panels removed 2026-08-24
%   per Oleg's review.)
% =========================================================================
if isempty(main_tab)
    figure('Position', [100, 100, 1700, 430]);
    parent_main = gcf;
else
    parent_main = main_tab;
end
SHOW_MEMORY   = false;  % 2026-09-07 (Oleg): peak-memory panel out of the paper figures
SHOW_BACKWARD = true;   % accuracy panel on (4th tile); metric chosen below
SHOW_FORWARD  = false;  % 2026-09-04 (Max): forward-error panel dropped from the paper layout
% 2026-09-14 (Max): 'data' plots ||Ax-b||/||b|| (ls_data_error, from the rounds
% sidecar), the same metric as the Toeplitz data-error panel, so both benchmark
% figures read directly against the 1e-11 noise floor. 'backward' restores the
% Higham normwise backward error ||Ax-b||/(||A||||x||+||b||) used 2026-09-04 to
% 09-11. Under a noisy right-hand side that measure ranked a wrong large-norm
% solution BEST on native_ill (CholQR2, ||x|| in the denominator), which is why
% it was retired from the paper figures.
ACCURACY_METRIC = 'data';
n_tiles = 4 + SHOW_MEMORY + SHOW_FORWARD;
tl_main = tiledlayout(parent_main, 1, n_tiles, 'TileSpacing', 'compact', 'Padding', 'compact');

title_label = sprintf('Sparse IR-LSQ Benchmark — %d \\times %d', m_val, n_val);
if ~isempty(title_suffix)
    title_label = sprintf('%s — %s', title_label, title_suffix);
end
if ~isempty(agg_note)
    title_label = sprintf('%s — %s', title_label, agg_note);
end
title(tl_main, title_label, 'FontWeight', 'bold');

x_pos = 1:n_algs;

% ---- (1) Stacked timing: QR + warm start (x_0 build) + solve ----
% Since 2026-08-27: ir_setup_us is a REAL column (the Blendenpik x0 build; the
% dead "embedded warm start" annotation block that used to live here targeted a
% row the roster filter had removed and matched a header phrase no CSV carries).
% The solve splits into inner-CG work vs OUTER REFINEMENT when the breakdown CSV
% is available (Max's two-color request): the inner-CG color (orange) MATCHES
% the iterations panel, so a bar's count and its cost are the same color. The
% second segment is Algorithm 1 lines 5 and 7 (recompute the true residual,
% map it through R^-T A^T, apply the correction R^-1 dy), which costs about one
% inner iteration per ROUND, so the two segments scale with different counters.
nexttile(tl_main, 1);
qr_ms  = arrayfun(@(i) qr_time(i)/1000,     sel_idx);
ir_ms  = arrayfun(@(i) ir_total_us(i)/1000, sel_idx);
ws_ms  = zeros(n_algs, 1);
if ismember('ir_setup_us', T.Properties.VariableNames)
    ws_ms = arrayfun(@(i) max(T.ir_setup_us(i), 0)/1000, sel_idx);
end
% Inner-CG slice from the breakdown CSV (IR row, t1 = inner_cg_total).
inner_ms = ir_ms; ovh_ms = zeros(n_algs, 1);
if isfile(breakdown_path)
    n_skip_b0 = count_comment_lines(breakdown_path);
    Tb0 = readtable(breakdown_path, 'NumHeaderLines', n_skip_b0);
    for a = 1:n_algs
        if sel_failed(a), continue; end
        match = strcmp(Tb0.algorithm, unique_algs{a}) & strcmp(Tb0.phase, 'IR');
        run_i = T.run(sel_idx(a));
        if run_i >= 0, match = match & Tb0.run == run_i; end
        idx = find(match, 1);
        if ~isempty(idx) && Tb0.t1(idx) > 0
            inner_ms(a) = Tb0.t1(idx)/1000;
            ovh_ms(a)   = max(ir_ms(a) - inner_ms(a), 0);
        end
    end
end
qr_ms(sel_failed) = 0; inner_ms(sel_failed) = 0; ovh_ms(sel_failed) = 0; ws_ms(sel_failed) = 0;
b = bar(x_pos, [qr_ms, ws_ms, inner_ms, ovh_ms], 'stacked');
b(1).FaceColor = w_blue;      b(1).DisplayName = 'QR / sketch build';
b(2).FaceColor = w_vermilion; b(2).DisplayName = 'warm start (x_0 build)';
b(3).FaceColor = w_orange;    b(3).DisplayName = 'solve: inner CG';
b(4).FaceColor = w_gray;      b(4).DisplayName = 'solve: outer loop';
if ~any(ws_ms > 0),  delete(b(2)); end           % legend stays clean on old data
if ~any(ovh_ms > 0), b(3).DisplayName = 'solve'; delete(b(4)); end
ylabel('Time (ms)'); title('Wall-time per algorithm');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
legend('Location', 'northeast'); grid on; box on;   % 2026-09-04 (Max): legend top-right
tot_ms = qr_ms + ws_ms + inner_ms + ovh_ms;
ylim([0, 1.35*max(tot_ms)]);   % headroom so the top-right legend clears the tallest bar (2026-09-04)
for a = 1:n_algs
    if sel_failed(a)   % 2% of the axis, not y=1 ms, so the label is visible at any scale
        text(x_pos(a), 0.02*max(1, max(tot_ms)), 'FAIL', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'Color', w_vermilion);
    end
end

% ---- Accuracy: data error (default) or Higham normwise backward error ----
% 2026-09-04 (Max): paper layout = wall clock, iterations, storage, orthogonality
% loss, accuracy (5th). The forward-error panel is off (SHOW_FORWARD).
if SHOW_BACKWARD
nexttile(tl_main, 4);   % accuracy: 4th panel (2026-09-07 order)
switch ACCURACY_METRIC
    case 'data'
        resid = arrayfun(@(i) ls_data_error(i), sel_idx);
        acc_ylabel = '||Ax - b|| / ||b||';
        acc_title  = 'Data error';
    case 'backward'
        resid = arrayfun(@(i) ls_residual_norm(i), sel_idx);
        acc_ylabel = '||Ax - b|| / (||A||\cdot||x|| + ||b||)';
        acc_title  = 'Normwise backward error';
    otherwise
        error('plot_irlsq_results: ACCURACY_METRIC must be ''data'' or ''backward'', got ''%s''', ACCURACY_METRIC);
end
resid(sel_failed) = NaN;
resid(resid < 0) = NaN;
bar(x_pos, resid, 'FaceColor', w_skyblue); set(gca, 'YScale', 'log');
% Dynamic decade limits (2026-08-27, ported from the Toeplitz plotter's 08-24
% fix): the old hard floor at 1e-16 sat ABOVE the best methods' values on
% native_ill (5.3e-17 to 6.9e-17), and log bars draw from the axis bottom, so
% the three most accurate methods rendered as MISSING bars. One decade of
% margin below the finite minimum keeps every bar visible.
fin = resid(isfinite(resid) & resid > 0);
if ~isempty(fin)
    ylim([10^(floor(log10(min(fin))) - 1), 10^ceil(log10(max(fin)))]);
else
    ylim([1e-16, 1e0]);
end
ylabel(acc_ylabel); title(acc_title);
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
grid on; box on;
yl = ylim;
for a = 1:n_algs
    if sel_failed(a) || isnan(resid(a))
        text(x_pos(a), yl(2)*0.9, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    end
end

end  % SHOW_BACKWARD

% ---- (3) Memory: peak RSS vs analytical prediction ----
% analytical <= 0 = "no analytical model" (-1 sentinel since 2026-08-27; 0 in
% older CSVs): no bar, rather than a real-looking zero-MB bar.
if SHOW_MEMORY
nexttile(tl_main, 5);   % optional 5th tile (2026-09-07)
mem_peak = arrayfun(@(i) peak_rss_kb(i),   sel_idx) / 1024;     % MB
mem_pred = arrayfun(@(i) analytical_kb(i), sel_idx) / 1024;     % MB
mem_pred(mem_pred <= 0) = NaN;
mem_peak(sel_failed) = NaN;  mem_pred(sel_failed) = NaN;
b = bar(x_pos, [mem_peak, mem_pred], 'grouped');
b(1).FaceColor = w_purple;     b(1).DisplayName = 'Peak RSS';
b(2).FaceColor = w_ltgray;     b(2).DisplayName = 'Analytical';
ylabel('Memory (MB)'); title('Peak vs predicted working memory');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
legend('Location', 'northwest'); grid on; box on;
end  % SHOW_MEMORY

% ---- (4) Inner CG iterations (total across the outer IR rounds) ----
% Algorithmic signal: lower = R is a better preconditioner for A^T A. Outer
% rounds VARY per method (up to ir_n_steps, with outer_tol and the LS-floor
% exit ending runs early; the old "fixed at 2 by construction" note described
% the pre-2026-08-07 scheme); their per-round overhead is the gray segment in
% panel (1), so plotting inner totals here stays honest.
nexttile(tl_main, 2);   % iterations: 2nd panel (2026-09-04 order)
inner_iters = arrayfun(@(i) ir_inner_total(i), sel_idx);
inner_iters(sel_failed) = NaN;
% NO warm-start segment here (2026-08-28, Max). The sketch-and-solve x0 build
% costs TIME but ZERO iterations, so the "time-equivalent iterations" base
% segment that used to sit under each bar was a manufactured quantity in the
% wrong units. It belongs in the wall-time panel, where it is a real measured
% duration, and nowhere else. This panel shows exactly one thing: iterations
% actually performed.
bar(x_pos, inner_iters, 'FaceColor', w_orange);
ylabel('Inner CG iterations (total)');
% Backward-error termination (2026-09-18, see be_censoring.m and the tick-label
% markers above): censored methods (stopped short of the target) and rows that
% never ran the oracle are washed out here as well; the bar-top text stays the
% bare count so it fits at paper size, the dagger on the tick label carries the
% verdict into every exported tile.
if BE.active
    wash = (BE.censored | BE.no_oracle) & ~isnan(inner_iters(:));
    if any(wash)
        % NaN-mask the full vector so the overlay bars keep the base bars' width
        % (bar() sizes bars from the spacing of the x values it is given).
        iters_c = inner_iters; iters_c(~wash) = NaN;
        hold on;
        bar(x_pos, iters_c, 'FaceColor', 'w', 'FaceAlpha', 0.6, ...
            'EdgeColor', w_orange, 'LineStyle', '--', 'LineWidth', 1.2);
    end
    title(sprintf('Inner CG iterations to backward-error target (%.1e)', BE.be_tol));
else
    title('Inner CG iterations to convergence');
end
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
grid on; box on;
yl = ylim;
for a = 1:n_algs
    if sel_failed(a) || isnan(inner_iters(a))
        text(x_pos(a), yl(2)*0.9, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    else
        % Annotate the actual integer count on top of the bar.
        text(x_pos(a), inner_iters(a), sprintf('%d', inner_iters(a)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontWeight', 'bold');
    end
end

% ---- (4) Orthogonality loss in Q-factor: ||Q^T Q - I||_F / sqrt(n) ----
nexttile(tl_main, 3);   % orthogonality: 3rd panel (2026-09-07 order)
orth_vals = arrayfun(@(i) orth_error(i), sel_idx);
orth_vals(sel_failed) = NaN;
orth_vals(orth_vals < 0) = NaN;
% w_green, NOT vermilion: orange is reserved for the warm-start segments
% (2026-07-31, Max: one orange thing per figure).
bar(x_pos, orth_vals, 'FaceColor', w_green); set(gca, 'YScale', 'log');
ylim([1e-16, 1e0]);
ylabel('||Q^T Q - I||_F / \surd n'); title('Q-factor orthogonality loss');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
grid on; box on;
yl = ylim;
for a = 1:n_algs
    if sel_failed(a)   % FAIL is distinct from N/A: the build broke vs no R exists (as on the Toeplitz figure)
        text(x_pos(a), yl(2)*0.9, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    elseif strcmp(unique_algs{a}, 'unpreconditioned') || isnan(orth_vals(a))
        text(x_pos(a), yl(2)*0.9, 'N/A', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_gray);
    end
end

% ---- (6) LS solution error: ||x - x_true|| / ||x_true|| (forward error) ----
% Reinstated 2026-08-31 (Max). The IR stop test is on the BACKWARD error, so a
% weak preconditioner can exit early with a backward error that passes the test
% while its solution sits orders further from x_true (native_ill CholQR: 3
% inner iterations, worst forward error in the cell). Panels (2) and (5) alone
% cannot show that; this panel is the discriminating one. -1 sentinel = no
% ground truth for the row.
if SHOW_FORWARD
nexttile(tl_main, 5 + SHOW_MEMORY);   % optional trailing tile (2026-09-07)
fwd = arrayfun(@(i) ls_solution_err(i), sel_idx);
fwd(sel_failed) = NaN;
fwd(fwd < 0) = NaN;
bar(x_pos, fwd, 'FaceColor', w_blue); set(gca, 'YScale', 'log');
% Dynamic decade limits, same rationale as panel (2): log bars draw from the
% axis bottom, so a hard floor can hide the most accurate methods entirely.
fin = fwd(isfinite(fwd) & fwd > 0);
if ~isempty(fin)
    ylim([10^(floor(log10(min(fin))) - 1), 10^ceil(log10(max(fin)))]);
else
    ylim([1e-16, 1e0]);
end
ylabel('||x - x_{true}|| / ||x_{true}||'); title('LS solution error (forward)');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
grid on; box on;
yl = ylim;
for a = 1:n_algs
    if sel_failed(a)
        text(x_pos(a), yl(2)*0.9, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    elseif isnan(fwd(a))
        text(x_pos(a), yl(2)*0.9, 'N/A', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_gray);
    end
end
end  % SHOW_FORWARD

% Footnote for the asterisked labels REMOVED per Max (2026-09-04): the asterisk
% stays on the x labels; its meaning (adaptive-shift-rescued row) is documented
% in run_all.m's era note and the dev log rather than on the figure.

% =========================================================================
%  Breakdown figure (optional): IR-LSQ phase breakdown stacked bar
%  Layout (6 fields): outer_total, inner_cg_total, trsm, fwd, adj, other
%  Rendered only when the caller supplies bd_tab (2026-08-27): the breakdown
%  CSV itself is also consumed by panel (1)'s solve split above, and passing it
%  must not force this extra figure into a tab-export run.
% =========================================================================
if isfile(breakdown_path) && ~isempty(bd_tab)
    n_skip_b = count_comment_lines(breakdown_path);
    Tb = readtable(breakdown_path, 'NumHeaderLines', n_skip_b);

    if isempty(bd_tab)
        figure('Position', [200, 100, 900, 500]);
        parent_bd = gcf;
    else
        parent_bd = bd_tab;
    end
    tl_bd = tiledlayout(parent_bd, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl_bd, sprintf('IR-LSQ runtime breakdown — %s', title_label), 'FontWeight', 'bold');

    % Build a per-algorithm matrix [inner_cg | trsm | fwd | adj | other] in ms.
    % All 5 segments come from IterRefineLSQ's populate_times() and sum to
    % outer_total = t0. (x_0 = 0 now: there is no sketch-and-solve initial guess,
    % so ir_total_us matches outer_total up to loop overhead — no x0 segment.)
    ir_breakdown = zeros(n_algs, 5);
    for a = 1:n_algs
        if sel_failed(a), continue; end
        run_idx = T.run(sel_idx(a));
        match = strcmp(Tb.algorithm, unique_algs{a}) & strcmp(Tb.phase, 'IR');
        if run_idx >= 0
            match = match & Tb.run == run_idx;   % 'best': the selected run's rows
            idx = find(match, 1);
            if isempty(idx), continue; end
            ir_breakdown(a, 1) = Tb.t1(idx)  / 1000;  % inner_cg ex. fwd/adj/trsm
            ir_breakdown(a, 2) = Tb.t2(idx)  / 1000;  % trsm
            ir_breakdown(a, 3) = Tb.t3(idx)  / 1000;  % fwd
            ir_breakdown(a, 4) = Tb.t4(idx)  / 1000;  % adj
            ir_breakdown(a, 5) = Tb.t5(idx)  / 1000;  % other
        else
            % run == -1 sentinel from aggregate_runs 'mean': average the
            % breakdown across all runs so it matches the averaged totals.
            if ~any(match), continue; end
            ir_breakdown(a, 1) = mean(Tb.t1(match)) / 1000;
            ir_breakdown(a, 2) = mean(Tb.t2(match)) / 1000;
            ir_breakdown(a, 3) = mean(Tb.t3(match)) / 1000;
            ir_breakdown(a, 4) = mean(Tb.t4(match)) / 1000;
            ir_breakdown(a, 5) = mean(Tb.t5(match)) / 1000;
        end
    end

    nexttile(tl_bd);
    b = bar(x_pos, ir_breakdown, 'stacked');
    b(1).FaceColor = w_orange;    b(1).DisplayName = 'inner CG control (axpy/dot)';
    b(2).FaceColor = w_skyblue;   b(2).DisplayName = 'TRSM';
    b(3).FaceColor = w_green;     b(3).DisplayName = 'J fwd';
    b(4).FaceColor = w_vermilion; b(4).DisplayName = 'J^T adj';
    b(5).FaceColor = w_gray;      b(5).DisplayName = 'axpy/copy/nrm2';
    ylabel('Time (ms)'); title('IR-LSQ phase breakdown');
    xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
    legend('Location', 'northeastoutside'); grid on; box on;
    for a = 1:n_algs
        if sel_failed(a)
            text(x_pos(a), 1, 'FAIL', 'HorizontalAlignment', 'center', ...
                 'FontWeight', 'bold', 'Color', w_vermilion);
        end
    end
end

end
