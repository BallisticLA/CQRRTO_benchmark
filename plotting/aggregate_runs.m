%% aggregate_runs.m — collapse multi-run benchmark CSVs to one row per algorithm.
%
% As of 2026-08-05 the Toeplitz and IR-LSQ benchmarks record num_runs repetitions
% per method (a 'run' column), because a single run at 4-7 solver iterations is
% at the wall-clock noise floor. This helper reduces the table back to the
% one-row-per-algorithm shape every plotter expects.
%
%   [Tagg, agg_note] = aggregate_runs(T)            % mode = 'best' (default)
%   [Tagg, agg_note] = aggregate_runs(T, 'mean')
%
% Modes:
%   'best'  Keep, per algorithm, the WHOLE row of the successful run with the
%           smallest total wall time (qr_time_us + solve/ir time). Whole-row
%           keeps timings, iteration counts, and accuracy mutually consistent
%           (they all come from the same run). Conventional "true cost" timing.
%   'mean'  Per algorithm, average every numeric column over the successful
%           runs. Count-like columns are rounded so sprintf('%d', ...) labels
%           in the plotters stay valid. The 'run' column is set to -1 as a
%           sentinel meaning "aggregate of several runs" (plotters that pair
%           rows with a breakdown CSV must special-case it).
%
% Failed runs (qr_status ~= 0) are excluded from aggregation; if every run of
% an algorithm failed, its first row is kept as-is so failures stay visible.
%
% Tables without a 'run' column (all campaigns before 2026-08-05) pass through
% untouched, as do single-run tables; agg_note is then ''. Otherwise agg_note
% is e.g. 'best of 5 runs' for figure titles, so provenance lands in the PDF.

function [Tagg, agg_note] = aggregate_runs(T, mode)

if nargin < 2 || isempty(mode), mode = 'best'; end
if ~any(strcmp(mode, {'best', 'mean'}))
    error('aggregate_runs: mode must be ''best'' or ''mean'', got ''%s''', mode);
end

agg_note = '';
if ~ismember('run', T.Properties.VariableNames), Tagg = T; return; end

algs = unique(T.algorithm, 'stable');
max_runs = 0;

% The total-time column pair differs between the two schemas.
if ismember('ir_total_us', T.Properties.VariableNames)
    solve_col = 'ir_total_us';        % IR-LSQ (irlsq_reg) results
    setup_col = 'ir_setup_us';        % warm-start x0 build (real since 2026-08-27)
elseif ismember('solve_time_us', T.Properties.VariableNames)
    solve_col = 'solve_time_us';      % Toeplitz LS results
    setup_col = 'setup_us';
else
    error('aggregate_runs: no ir_total_us or solve_time_us column; unknown schema');
end

% Columns whose mean must stay integer-valued for the plotters' %d labels.
count_cols = {'iterations', 'ir_inner_iters_total', 'ir_outer_iters', ...
              'chol_retries', 'solver_flag', 'ir_inner_capped', ...
              'ir_inner_best_iter', 'qr_status'};

keep = zeros(numel(algs), 1);
mean_rows = cell(numel(algs), 1);
for a = 1:numel(algs)
    rows = find(strcmp(T.algorithm, algs{a}));
    max_runs = max(max_runs, numel(rows));
    succ = rows(T.qr_status(rows) == 0);
    if isempty(succ)
        keep(a) = rows(1);            % all failed: keep one row, visibly failed
        continue;
    end
    % Backward-error termination (2026-09-18): the min-total-time pick below is
    % ~95% QR-build jitter, so if the runs of one method disagree on HOW they
    % ended (e.g. 'be' in two runs, 'floor' in one, for a method sitting near the
    % target) the displayed verdict would be a coin flip. Restrict the candidates
    % to the modal stop_reason first (preferring 'be' on a tie) and say so.
    if ismember('stop_reason', T.Properties.VariableNames) && numel(succ) > 1
        rs = string(T.stop_reason(succ)); rs = rs(:);
        [u, ~, j] = unique(rs);
        if numel(u) > 1
            cnt = accumarray(j, 1);
            modal = u(cnt == max(cnt));
            if any(modal == "be"), modal = "be"; else, modal = modal(1); end
            warning('aggregate_runs: %s runs disagree on stop_reason (%s); keeping the %s runs', ...
                    algs{a}, strjoin(cellstr(rs'), '/'), modal);
            if strcmp(mode, 'mean')
                error('aggregate_runs: mode ''mean'' cannot average runs with different stop_reason (%s)', algs{a});
            end
            succ = succ(rs == modal);
        end
        if strcmp(mode, 'mean') && any(rs == "be")
            error('aggregate_runs: mode ''mean'' is not defined for an oracle-terminated era (stop_reason ''be''); use ''best''');
        end
    end
    % Best-total includes the warm-start setup slice (2026-08-27; before, a warm
    % row's x0 build did not count against its run selection).
    total = T.qr_time_us(succ) + T.(solve_col)(succ);
    if ismember(setup_col, T.Properties.VariableNames)
        total = total + max(T.(setup_col)(succ), 0);
    end
    [~, ibest] = min(total);
    keep(a) = succ(ibest);
    if strcmp(mode, 'mean')
        row = T(succ(ibest), :);      % template keeps non-numeric columns
        for v = T.Properties.VariableNames
            col = v{1};
            if ~isnumeric(T.(col)), continue; end
            vals = T.(col)(succ);
            % Sentinel masking (2026-08-27): -1 means "not measured" in several
            % numeric columns; averaging it as data biases or corrupts the mean
            % (mean(-1, 1e-14) ~ -0.5, which the plotters then NaN away).
            vals(vals == -1) = NaN;
            mval = mean(vals, 'omitnan');
            if isnan(mval), mval = -1; end   % all-sentinel column stays sentinel
            if any(strcmp(col, count_cols)), mval = round(mval); end
            row.(col) = mval;
        end
        row.run = -1;                 % sentinel: not a single run
        mean_rows{a} = row;
    end
end

Tagg = T(keep, :);
if strcmp(mode, 'mean')
    for a = 1:numel(algs)
        if ~isempty(mean_rows{a}), Tagg(a, :) = mean_rows{a}; end
    end
end

if max_runs > 1
    agg_note = sprintf('%s of %d runs', mode, max_runs);
end

end
