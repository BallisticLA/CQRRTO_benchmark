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
alg_csv_order  = {'CQRRT_linop','CQRRT_linop_gemmL', 'CholQR', 'CholQR2', 'sCholQR3_basic', 'sCholQR3', 'Blendenpik', 'Blendenpik_cold', 'Blendenpik_refine', 'Blendenpik_cold_refine'};
alg_disp_names = {'CQRRT','CQRRT (GEMM left)', 'CholQR', 'CholQR2', 'sCholQR3 (no blocking)', 'sCholQR3', 'Blendenpik', 'Blendenpik (zero x_0)', 'Blendenpik + refinement', 'Blendenpik (zero x_0) + refinement'};

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

% Multi-run CSVs (num_runs > 1, 2026-08-05): collapse to one row per algorithm
% up front, so every column copy below sees the aggregated table. agg_note
% ('best of 5 runs' / 'mean of 5 runs') lands in the figure title.
[T, agg_note] = aggregate_runs(T, timing_agg);

% 2026-08-24 (Oleg): paper-figure roster trim. sCholQR3_basic (the no-blocking
% control) and every Blendenpik row except the stabilized "Blendenpik +
% refinement" are dropped from the PLOTS only; the CSVs keep all rows.
alg_plot_exclude = {'sCholQR3_basic', 'Blendenpik', 'Blendenpik_cold', 'Blendenpik_cold_refine'};
T(ismember(T.algorithm, alg_plot_exclude), :) = [];

algorithms       = T.algorithm;
qr_status        = T.qr_status;
qr_time          = T.qr_time_us;
peak_rss_kb      = T.peak_rss_kb;
analytical_kb    = T.analytical_kb;
ir_total_us      = T.ir_total_us;
ir_inner_total   = T.ir_inner_iters_total;
ls_residual_norm = T.ls_residual_norm;
orth_error       = T.orth_error;
m_val            = T.m(1);
n_val            = T.n(1);
% (The adaptive-shift retries and canonical-GEQRF-rate panels were removed
%  2026-08-24 per Oleg; chol_retries stays in the CSVs, just unplotted.)

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
%  Main figure: 2x3 layout (5 panels; bottom-right tile blank)
%    row 1: (1) Stacked timing QR+IR  (2) Normwise backward error
%           (3) Memory peak vs analytical
%    row 2: (4) Q-factor orthogonality loss
%           (5) Inner CG iterations (setup/work stacked when setup exists)
%  (QR-build canonical rate + adaptive-shift retries panels removed 2026-08-24
%   per Oleg's review.)
% =========================================================================
if isempty(main_tab)
    figure('Position', [100, 100, 1100, 900]);
    parent_main = gcf;
else
    parent_main = main_tab;
end
tl_main = tiledlayout(parent_main, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

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
% The solve splits into inner-CG work vs restart overhead when the breakdown CSV
% is available (Max's two-color request): the inner-CG color (orange) MATCHES
% the iterations panel, so a bar's count and its cost are the same color.
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
b(4).FaceColor = w_gray;      b(4).DisplayName = 'solve: restart overhead';
if ~any(ws_ms > 0),  delete(b(2)); end           % legend stays clean on old data
if ~any(ovh_ms > 0), b(3).DisplayName = 'solve'; delete(b(4)); end
ylabel('Time (ms)'); title('Wall-time per algorithm');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
legend('Location', 'northwest'); grid on; box on;
tot_ms = qr_ms + ws_ms + inner_ms + ovh_ms;
for a = 1:n_algs
    if sel_failed(a)   % 2% of the axis, not y=1 ms, so the label is visible at any scale
        text(x_pos(a), 0.02*max(1, max(tot_ms)), 'FAIL', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'Color', w_vermilion);
    end
end

% ---- (2) Residual: Higham normwise backward error ----
nexttile(tl_main, 2);
resid = arrayfun(@(i) ls_residual_norm(i), sel_idx);
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
ylabel('||Ax - b|| / (||A||\cdot||x|| + ||b||)'); title('Normwise backward error');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
grid on; box on;
yl = ylim;
for a = 1:n_algs
    if sel_failed(a) || isnan(resid(a))
        text(x_pos(a), yl(2)*0.9, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    end
end

% ---- (3) Memory: peak RSS vs analytical prediction ----
% analytical <= 0 = "no analytical model" (-1 sentinel since 2026-08-27; 0 in
% older CSVs): no bar, rather than a real-looking zero-MB bar.
nexttile(tl_main, 3);
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

% ---- (4) Inner CG iterations (total across the outer IR rounds) ----
% Algorithmic signal: lower = R is a better preconditioner for A^T A. Outer
% rounds VARY per method (up to ir_n_steps, with outer_tol and the LS-floor
% exit ending runs early; the old "fixed at 2 by construction" note described
% the pre-2026-08-07 scheme); their per-round overhead is the gray segment in
% panel (1), so plotting inner totals here stays honest.
nexttile(tl_main, 5);
inner_iters = arrayfun(@(i) ir_inner_total(i), sel_idx);
inner_iters(sel_failed) = NaN;
% Setup-vs-work split (2026-07-30, Max): when a method carries prep work
% (the sketch-and-solve warm start), show it as a base segment in
% TIME-EQUIVALENT iterations, so the stacked total stays proportional
% to the solve time in panel (1) and the two panels cannot disagree.
setup_iters = zeros(n_algs, 1);
if ismember('ir_setup_us', T.Properties.VariableNames)
    su  = max(T.ir_setup_us(sel_idx), 0);
    wus = ir_total_us(sel_idx);                % solve time (setup is its own slot now)
    setup_iters = su ./ wus .* inner_iters;    % warm start in time-equiv iterations
    setup_iters(~isfinite(setup_iters) | wus <= 0) = 0;   % degenerate denominator: suppress
end
hb = bar(x_pos, [setup_iters, inner_iters], 'stacked');
hb(1).FaceColor = w_vermilion; hb(1).DisplayName = 'warm start (time-equiv iters)';
hb(2).FaceColor = w_orange;    hb(2).DisplayName = 'inner CG';
if any(setup_iters > 0), legend('Location', 'northwest'); end
ylabel('Inner CG iterations (total)');
title('Inner CG iterations to convergence');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
grid on; box on;
yl = ylim;
for a = 1:n_algs
    if sel_failed(a) || isnan(inner_iters(a))
        text(x_pos(a), yl(2)*0.9, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    else
        % Annotate the actual integer count on top of the stacked bar.
        text(x_pos(a), setup_iters(a) + inner_iters(a), sprintf('%d', inner_iters(a)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontWeight', 'bold');
    end
end

% ---- (4) Orthogonality loss in Q-factor: ||Q^T Q - I||_F / sqrt(n) ----
nexttile(tl_main, 4);
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
    if sel_failed(a) || isnan(orth_vals(a))
        text(x_pos(a), yl(2)*0.9, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    end
end

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
