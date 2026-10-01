%% plot_toeplitz_results.m — plot the Toeplitz LS benchmark (Oleg's 2nd experiment).
%
% One row per method in the CSV (num_runs=1). 5-panel layout matching the App-1
% plotter: wall-time (build+solve), data error, memory (peak vs analytical),
% Q-orthogonality loss, solver iterations. (The recovery/forward-error bar, the
% canonical-GEQRF-rate panel and the Cholesky-retries panel were removed
% 2026-08-24 per Oleg's review.)
%
% CSV schema (see toeplitz_ls_benchmark.cc; all reads below are NAME-based, so
% this comment is documentation, not a parsing contract):
%   algorithm,run,m,n,qr_status,qr_time_us,solve_time_us,peak_rss_kb,analytical_kb,
%   orth_error,iterations,solver_flag,solver_relres,aug_relres,normal_relres,
%   data_relres,recovery_error,cond_estimate,chol_retries,
%   solve_fwd_us,solve_adj_us,solve_trsm_us,setup_us,solver,pcg_rounds,
%   lsqr_iters,stop_reason,t_inner_us,t_fwd_inner_us,t_adj_inner_us,
%   t_trsm_inner_us,t_overhead_us,total_row_us
% (the last block was added 2026-08-27: pcg_restarts renamed pcg_rounds, and the
%  inner-kernel vs restart-overhead solve split became recordable. Older CSVs
%  simply lack those columns and render with a single solve segment.)
%
% Usage:
%   plot_toeplitz_results(data_dir, results_csv, title_suffix, main_tab)
%   plot_toeplitz_results(data_dir, results_csv, title_suffix, main_tab, timing_agg)
%
% timing_agg ('best' default | 'mean') controls how multi-run CSVs (num_runs > 1,
% 2026-08-05) collapse to one row per method; see aggregate_runs.m.

function plot_toeplitz_results(data_dir, results_csv, title_suffix, main_tab, timing_agg)

if nargin < 3, title_suffix = ''; end
if nargin < 4, main_tab = []; end
if nargin < 5 || isempty(timing_agg), timing_agg = 'best'; end

% Wong colorblind-friendly palette (matches the App-1 plotter).
w_blue = [0 114 178]/255;  w_orange = [230 159 0]/255;  w_skyblue = [86 180 233]/255;
w_green = [0 158 115]/255; w_vermilion = [213 94 0]/255; w_purple = [204 121 167]/255;
w_gray = [0.65 0.65 0.65];

% "Blendenpik_refine"/"Blendenpik_cold_refine" appear in [NEW 08-09]+ CSVs: Blendenpik's
% OWN preconditioner and answer handed to our restarted-PCG refinement, so the published
% rows and the like-for-like comparison sit side by side (as published, Blendenpik has no
% refinement loop, which otherwise conflates preconditioner quality with solver structure).
% "Blendenpik_cold" appears in [NEW 08-05]+ CSVs (warm x_0 is Blendenpik-only and
% both variants run). "Blendenpik" keeps its bare label because its x_0 policy is
% era-dependent (old cold campaigns disabled its warm start); the era note carries it.
alg_csv_order  = {'CQRRTO_linop','CQRRTO_linop_gemmL','CQRRT_linop','CQRRT_linop_gemmL','CholQR','CholQR2','sCholQR3_basic','sCholQR3','Blendenpik','Blendenpik_cold','Blendenpik_refine','Blendenpik_cold_refine','unpreconditioned'};
% 2026-09-04 (Max): "refinement" banned from figure text; refined variants take
% the bare labels (single-shot rows are plot-excluded, names kept for re-enable).
alg_disp_names = {'CQRRTO','CQRRTO (GEMM left)','CQRRTO','CQRRTO (GEMM left)','CholQR','CholQR2','sCholQR3 (no blocking)','sCholQR3','Blendenpik (single-shot)','Blendenpik (single-shot, zero x_0)','Blendenpik (sketch-and-solve x_0)','Blendenpik','unprecond'};   % 'unprecond' (2026-09-27, Oleg): the CSV name stays 'unpreconditioned', only the label is short

results_path = fullfile(data_dir, results_csv);
if ~isfile(results_path), error('plot_toeplitz_results: file not found: %s', results_path); end
n_skip = count_comment_lines(results_path);
opts = detectImportOptions(results_path, 'NumHeaderLines', n_skip);
opts = setvartype(opts, 'algorithm', 'char');
T = readtable(results_path, opts);

% Multi-run CSVs (num_runs > 1, 2026-08-05): collapse to one row per method
% before anything indexes the table. agg_note lands in the figure title.
% Join the KW sidecar BEFORE aggregate_runs so the kept run carries its own
% value (same rule as attach_data_error in plot_irlsq_results.m). Joining after
% aggregation happens to work for timing_agg='best', which preserves the run
% index, but would silently empty the panel under 'mean', where run is set to -1.
T = attach_kw_backward_error(T, strrep(results_path, '_results.csv', '_backward_error.csv'));
[T, agg_note] = aggregate_runs(T, timing_agg);

% 2026-08-24 (Oleg): paper-figure roster trim -- PARTIALLY OVERRIDDEN 2026-09-04
% (Max): the two cold Blendenpik rows are back in the plots. On this noise-free
% benchmark the warm x_0 (sketch-and-solve on a consistent system) starts AT the
% u*kappa accuracy floor, so the warm refine row shows no iterative work at all;
% the cold rows are the ones that exercise the solve and are the honest
% comparator for published (zero-x_0) Blendenpik. sCholQR3_basic (no-blocking
% control) and the plain warm "Blendenpik" row (exits at 0 LSQR iterations via
% the rescaled stop test; ruling on that row still open) stay excluded.
% Flag the roster change to Oleg rather than letting the figure drift silently.
% 2026-09-04 (Max): CQRRT_linop_gemmL (RANDLAPACK_GRAM_LEFT=gemm arm) is plot-excluded:
% the paper's PCholQR model applies the inverse preconditioner by triangular
% solves, which is the library-default TRSM arm; the GEMM-left arm stays in the
% merged CSVs for re-enable. Display name is CQRRTO, the paper's algorithm name.
% 2026-09-15 (Oleg): the sketch-and-solve-started row (Blendenpik_refine) is OUT
% of the paper figures; the zero-start row is the published algorithm and takes
% the bare label 'Blendenpik'. The warm-start wall-time slice disappears with it.
alg_plot_exclude = {'sCholQR3_basic', 'Blendenpik', 'Blendenpik_cold', 'Blendenpik_refine', 'CQRRT_linop_gemmL', 'CQRRTO_linop_gemmL'};  % single-shot rows + warm-start row + GEMM-left arm excluded
T(ismember(T.algorithm, alg_plot_exclude), :) = [];

algs = T.algorithm;
m_val = T.m(1); n_val = T.n(1);

% Order rows by display order (skip methods not present).
unique_algs = {}; disp_labels = {}; sel = [];
for k = 1:numel(alg_csv_order)
    idx = find(strcmp(algs, alg_csv_order{k}), 1);
    if ~isempty(idx)
        unique_algs{end+1} = alg_csv_order{k}; %#ok<AGROW>
        disp_labels{end+1} = alg_disp_names{k}; %#ok<AGROW>
        sel(end+1) = idx; %#ok<AGROW>
    end
end
% Defensive extras (2026-08-27, mirroring plot_irlsq_results): any CSV method not
% in the display-order list is APPENDED under its raw name rather than silently
% dropped, so a renamed or new benchmark row cannot vanish from the figure.
extras = setdiff(unique(algs, 'stable'), alg_csv_order, 'stable');
for k = 1:numel(extras)
    idx = find(strcmp(algs, extras{k}), 1);
    unique_algs{end+1} = extras{k}; %#ok<AGROW>
    disp_labels{end+1} = strrep(extras{k}, '_', '\_'); %#ok<AGROW>
    sel(end+1) = idx; %#ok<AGROW>
    warning('plot_toeplitz_results: method "%s" not in display order; appended.', extras{k});
end
sel = sel(:);   % column, so metric vectors are Nx1 and [a,b] is Nx2 for grouped/stacked bars
na = numel(sel);
x = 1:na;
failed = (T.qr_status(sel) ~= 0);

% Mark adaptive-shift rescues on the method labels (2026-08-31, Max; mirrors
% plot_irlsq_results). chol_retries > 0 means the Cholesky broke down during
% the build and the adaptive shift rescued it, so the row measures the
% shift-rescued variant of its algorithm (cholqr_primitive contract caveat).
% On the prolate Toeplitz cells this is CholQR's normal state: a single
% unshifted pass at this conditioning cannot succeed. Pre-08-27 CSVs lack the
% column, hence the guard.
rescued = false(na, 1);
if ismember('chol_retries', T.Properties.VariableNames)
    for a = 1:na
        if ~failed(a) && T.chol_retries(sel(a)) > 0
            rescued(a) = true;   % no label marker since 2026-09-27: the paper's tables carry the retry counts
        end
    end
end

% Backward-error termination (2026-09-18): when the era ran with the engine's
% oracle on (header echoes be_tol > 0), each method is judged from its RECORDED
% final estimate (be_final <= be_tol), not from the exit code. A method that
% stopped short of the target is CENSORED: its bars are washed out in every
% panel via the tick label (a dagger), since the paper exports each tile on its
% own; a live row that never ran the oracle (LSQR rows) gets a double dagger.
% Eras without the knob keep the old rendering. See be_censoring.m.
BE = be_censoring(T, sel, failed, results_path);
% Dagger / double-dagger label suffixes REMOVED (2026-09-27, Max): the censored
% row is still washed out in the iteration panel below (BE.censored), and the
% paper's caption and Toeplitz table carry the verdict in words.


% --- Figure / tab ---
if isempty(main_tab), figure('Position',[100 100 2050 430]); parent = gcf; else parent = main_tab; end
% 1x5 layout (1x4 from 2026-09-07 per Oleg): wall-time, inner CG iterations,
% orthogonality loss, data error, KW backward error; one line, same footprint as
% the synthetic panels. The peak-memory panel was retired from the paper figures
% on 2026-09-07 (its numbers go in the text) and its disabled code was removed
% 2026-09-21; git history has it. Per-panel paper export lives in run_all.m,
% whose PAPER_SUFFIX list must match the tile order.
% 2026-09-21 (Max): 5th tile = sketched Karlson-Walden backward error, the
% quantity Epperly's stopping rule targets. Read from the *_backward_error.csv
% sidecar, which is written whether or not the oracle drove termination, so an
% era run with be_tol_mult=0 is comparable here to one run with it on.
SHOW_KW = true;
n_tiles = 4 + SHOW_KW;
tl = tiledlayout(parent, 1, n_tiles, 'TileSpacing','compact', 'Padding','compact');
ttl = sprintf('Toeplitz LS Benchmark — %d \\times %d', m_val, n_val);
if ~isempty(title_suffix), ttl = sprintf('%s — %s', ttl, title_suffix); end
if ~isempty(agg_note), ttl = sprintf('%s — %s', ttl, agg_note); end
title(tl, ttl, 'FontWeight','bold', 'FontSize', 13);

% (1) Wall-time: build + warm start (x0 build) + solve, stacked, ms. Since
% 2026-08-27 the solve splits into inner-CG work vs OUTER REFINEMENT (Algorithm 1
% lines 5 and 7: recompute the true residual b - Ax, map it to the NE space
% via R^-T A^T, and apply the correction R^-1 dy). That is not overhead in the
% sense of waste: it is the refinement step itself, and it costs about one
% inner iteration per round (same 1 fwd + 1 adj + 2 trsv), which is why the
% two segments scale with DIFFERENT counters (iterations vs rounds), Max's two-color request:
% the inner-CG color (orange) MATCHES the iterations panel, so the count a bar
% carries and the time it cost are the same color. Vermilion = warm-start x0
% build (real for the Blendenpik + refinement row since the 08-27 redesign).
% Old CSVs without t_inner_us render a single orange solve segment.
nexttile(tl, 1);
build_ms = T.qr_time_us(sel)/1000;  solve_ms = T.solve_time_us(sel)/1000;
ws_ms = zeros(na, 1);
if ismember('setup_us', T.Properties.VariableNames)
    ws_ms = max(T.setup_us(sel), 0)/1000;
end
have_split = ismember('t_inner_us', T.Properties.VariableNames);
if have_split
    inner_ms = T.t_inner_us(sel)/1000;
    ovh_ms   = T.t_overhead_us(sel)/1000;
    nosplit  = (T.t_inner_us(sel) < 0);            % -1 = no history (LSQR rows)
    inner_ms(nosplit) = solve_ms(nosplit);          % lumped into the inner color
    ovh_ms(nosplit)   = 0;
else
    inner_ms = solve_ms; ovh_ms = zeros(na, 1);
end
build_ms(failed) = 0; inner_ms(failed) = 0; ovh_ms(failed) = 0; ws_ms(failed) = 0;
% Plotted in SECONDS since 2026-09-27 (same reason as plot_irlsq_results: the
% exponent label of a millisecond axis does not survive the 110-pt paper panel).
% The vectors keep their _ms names; every use below is in seconds.
build_ms = build_ms/1e3; ws_ms = ws_ms/1e3; inner_ms = inner_ms/1e3; ovh_ms = ovh_ms/1e3;
b = bar(x, [build_ms, ws_ms, inner_ms, ovh_ms], 'stacked');
b(1).FaceColor = w_blue;      b(1).DisplayName = 'build';        % short legend entries for the four-across paper row (2026-09-27)
b(2).FaceColor = w_vermilion; b(2).DisplayName = 'warm start (x_0 build)';
b(3).FaceColor = w_orange;    b(3).DisplayName = 'inner CG';
b(4).FaceColor = w_gray;      b(4).DisplayName = 'outer loop';
if ~any(ws_ms > 0),  delete(b(2)); end
if ~any(ovh_ms > 0), b(3).DisplayName = 'solve'; delete(b(4)); end
ylabel('Time (s)'); title('Wall-time per algorithm');   % linear scale, starts at 0
tot_ms = build_ms + ws_ms + inner_ms + ovh_ms;
ylim([0 max(1, 1.15*max(tot_ms))]);
xticks(x); xticklabels(disp_labels); xtickangle(35); legend('Location','northeast'); grid on; box on;
for a = 1:na   % failed rows: explicit marker instead of a zero-height bar (2026-08-27)
    if failed(a)
        text(x(a), 0.02*max(1, max(tot_ms)), 'FAIL', 'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', 'FontWeight','bold', 'Color', w_vermilion);
    end
end

% (2) Accuracy: data rel error (||Tx-b||/||b||) only. The recovery/forward-error
% bar (||x-xtrue||/||xtrue||, green) was removed 2026-08-24 per Oleg's review;
% recovery_error stays in the CSVs, just unplotted.
nexttile(tl, 4);   % data error: 4th panel (2026-09-07 order)
dre = T.data_relres(sel);
dre(failed) = NaN; dre(dre<0) = NaN;
bar(x, dre, 'FaceColor', w_skyblue); set(gca,'YScale','log');
% Explicit y-limits: the values cluster at the data noise floor (~1e-11) with
% unpreconditioned ~1 decade above; MATLAB's auto limits put the axis bottom
% ABOVE the floor bars, rendering them invisible (log bars draw from the axis
% bottom). One decade of margin below the min keeps every bar visible.
fin = dre(isfinite(dre) & dre > 0);
if ~isempty(fin)
    ylim([10^(floor(log10(min(fin))) - 1), 10^ceil(log10(max(fin)))]);
end
ylabel('||Tx-b||/||b||'); title('Data error');
xticks(x); xticklabels(disp_labels); xtickangle(35); grid on; box on;

% (4b) Sketched Karlson-Walden backward error -- the criterion Epperly's rule
% stops on. Deliberately NOT the data error: a run can meet this target while its
% residual is still orders above the noise floor.
if SHOW_KW
nexttile(tl, 5);
kw = T.be_kw_measured(sel);
kw(failed) = NaN; kw(kw < 0) = NaN;
bar(x, kw, 'FaceColor', w_purple); set(gca,'YScale','log');
fin = kw(isfinite(kw) & kw > 0);
if ~isempty(fin)
    ylim([10^(floor(log10(min(fin))) - 1), 10^ceil(log10(max(fin)))]);
else
    ylim([1e-18 1e0]);
end
if isfield(BE, 'active') && BE.active && BE.be_tol > 0
    hold on; yline(BE.be_tol, '--', sprintf('target %.1e', BE.be_tol), ...
                   'Color', w_vermilion, 'LineWidth', 1.1, ...
                   'LabelHorizontalAlignment','left'); hold off;
end
ylabel('sketched KW backward error / ||A||_F');
title('Backward error (Karlson-Walden)');
xticks(x); xticklabels(disp_labels); xtickangle(35); grid on; box on;
yl = ylim;
for a = 1:numel(x)
    if isnan(kw(a))
        text(x(a), yl(2)*0.9, 'n/a', 'HorizontalAlignment','center', ...
             'FontWeight','bold', 'Color', w_vermilion);
    end
end
end  % SHOW_KW

% (3) Memory: peak RSS vs analytical (MB). analytical <= 0 means "no analytical
% model for this row" (-1 sentinel since 2026-08-27; 0 in older CSVs): render as
% no bar plus an N/A tag rather than a real zero-MB bar.
% (4) Solver work count: inner CG iterations of the shared PCG-NE engine (the
% recorded quantity for every plotted row since the 08-27 refine redesign; the
% pre-08-27 "LSQR" titles misnamed it). Orange matches the wall-time panel's
% inner-CG segment.
%
% NO warm-start segment here (2026-08-28, Max). The x0 build costs TIME but ZERO
% iterations, so the "time-equivalent iterations" base segment that used to sit
% under each bar was a manufactured quantity in the wrong units. It belongs in
% the wall-time panel, where it is a real measured duration, and nowhere else.
% This panel now shows exactly one thing: iterations actually performed.
nexttile(tl, 2);   % iterations: 2nd panel (2026-09-04 order)
iters = T.iterations(sel); iters(failed) = NaN;
bar(x, iters, 'FaceColor', w_orange);
% Backward-error termination (2026-09-18, see be_censoring.m and the tick-label
% markers above): censored methods (stopped short of the target) and rows that
% never ran the oracle are washed out here as well; the bar-top text stays the
% bare count so it fits at paper size, the dagger on the tick label carries the
% verdict into every exported tile.
if BE.active
    wash = BE.censored | BE.no_oracle;
    if any(wash)
        % NaN-mask the full vector so the overlay bars keep the base bars' width
        % (bar() sizes bars from the spacing of the x values it is given).
        iters_c = iters; iters_c(~wash) = NaN;
        hold on;
        bar(x, iters_c, 'FaceColor', 'w', 'FaceAlpha', 0.6, ...
            'EdgeColor', w_orange, 'LineStyle', '--', 'LineWidth', 1.2);
    end
    if strcmp(BE.mode, 'floor')
        title(sprintf('Inner CG iterations to the data-error floor (%.1e)', BE.floor));
    else
        title(sprintf('Inner CG iterations to backward-error target (%.1e)', BE.be_tol));
    end
else
    title('Inner CG iterations (PCG-NE) to convergence');
end
ylabel('Inner CG iterations');   % same wording as the FEM panel; capital I (Max, 2026-09-27)
xticks(x); xticklabels(disp_labels); xtickangle(35); grid on; box on;
for a = 1:na
    if failed(a)
        text(x(a), 0, 'FAIL', 'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
             'FontWeight','bold', 'Color', w_vermilion);
        continue;
    end
    text(x(a), iters(a), sprintf('%d', iters(a)), 'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
end

% (5) Q-orthogonality loss ||R^-T A'A R^-1 - I||_F/sqrt(n). -1 => no preconditioner.
% w_green so it cannot be confused with the reserved segment colors (vermilion =
% warm start, orange = inner CG; comment corrected 2026-08-27, it used to name
% the wrong color). Dynamic top limit so a kappa(R)^2-regression value above 1
% pins visibly instead of clipping silently (2026-08-27); floor stays at 1e-16.
nexttile(tl, 3);   % orthogonality: 3rd panel (2026-09-07 order)
orth = T.orth_error(sel); orth(failed) = NaN; orth(orth<0) = NaN;
ymax = 1e1;
finite_orth = orth(isfinite(orth));
if ~isempty(finite_orth), ymax = max(1e1, 10^ceil(log10(max(finite_orth)))); end
bar(x, orth, 'FaceColor', w_green); set(gca,'YScale','log'); ylim([1e-16 ymax]);
ylabel('||R^{-T}A^TA R^{-1}-I||_F/\surd n'); title('Preconditioner orthogonality loss');
xticks(x); xticklabels(disp_labels); xtickangle(35); grid on; box on;
yl = ylim;
for a = 1:na
    if failed(a)   % FAIL is distinct from N/A: the build broke vs no R exists (2026-08-27)
        text(x(a), yl(2)*0.5, 'FAIL', 'HorizontalAlignment','center','FontWeight','bold','Color',w_vermilion);
    elseif strcmp(unique_algs{a}, 'unpreconditioned') || isnan(orth(a))
        text(x(a), yl(2)*0.5, 'N/A', 'HorizontalAlignment','center','FontWeight','bold','Color',w_gray);
    end
end

% Footnote for the asterisked labels REMOVED per Max (2026-09-04): the asterisk
% stays on the x labels; its meaning (adaptive-shift-rescued row) is documented
% in run_all.m's era note and the dev log rather than on the figure.

end
