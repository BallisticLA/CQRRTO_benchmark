%% merge_overlay_rows.m — append rows from a second era's cell to a base cell's CSV set.
%
% Some rows are cheaper to run as their own era than to fold into a queued one.
% The first case (2026-09-15) is the unpreconditioned reference row (mask bit
% 128) that Oleg asked to see on the FEM figures: it needed a new binary, and
% rebuilding the campaign install would have replaced the binary under the
% still-queued 0914_be1 cells, so it ran as era 0915_up1 on its own. This helper
% merges such rows into the base cell deterministically at plot time, the way
% merge_gram_arms merges the two Gram arms:
%
%   [out_dir, csv_name] = merge_overlay_rows(base_dir, base_csv, overlay_cell_dir, rows)
%
% reads base_dir/base_csv, takes from the NEWEST non-empty results CSV in
% overlay_cell_dir only the DATA rows whose algorithm is in `rows`, and writes
% base_dir/merged_overlay/<base_csv> = base file verbatim + one comment line
% naming the overlay file and its build + the kept rows. The same is done for
% the same-stem _breakdown, _rounds and _backward_error sidecars when both sides
% have them. Returns the merged directory and the (unchanged) CSV name, so the
% caller swaps them in for the plotter; the merge is re-done on every call.
%
% Provenance: the overlay's RANDLAPACK_GIT_COMMIT is written as
% "# OVERLAY_COMMIT=<sha> rows=<...> from=<path>", deliberately NOT as a second
% RANDLAPACK_GIT_COMMIT line: assert_era_commit reads the first stamp per file
% and the merged file is still the base era's file. The overlay era gets its own
% assert_era_commit call in run_all, and both builds go into the figure note.
% Column lines must match exactly, else this errors rather than misalign.

function [out_dir, csv_name] = merge_overlay_rows(base_dir, base_csv, overlay_cell_dir, rows)

csv_name = base_csv;
out_dir  = fullfile(base_dir, 'merged_overlay');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

ov_results = newest_csv(overlay_cell_dir, '*_results.csv');
if isempty(ov_results)
    error('merge_overlay_rows: no non-empty *_results.csv under %s', overlay_cell_dir);
end
ov_stem = regexprep(ov_results, '_results\.csv$', '');
ov_commit = header_commit(fullfile(overlay_cell_dir, ov_results));

kinds = {'results', 'breakdown', 'rounds', 'backward_error'};
for k = 1:numel(kinds)
    base_path = fullfile(base_dir, regexprep(base_csv, '_results\.csv$', ['_' kinds{k} '.csv']));
    ov_path   = fullfile(overlay_cell_dir, [ov_stem '_' kinds{k} '.csv']);
    if ~isfile(base_path)
        continue;   % the base cell has no such sidecar; nothing to merge into
    end
    if ~isfile(ov_path)
        % Keep the plotter's strrep-derived sidecar lookups working: copy the base
        % sidecar through unchanged, and say so.
        copyfile(base_path, fullfile(out_dir, regexprep(base_csv, '_results\.csv$', ['_' kinds{k} '.csv'])));
        fprintf('  merge_overlay_rows: overlay has no %s sidecar; base copied unchanged\n', kinds{k});
        continue;
    end
    merge_pair(base_path, ov_path, out_dir, rows, ov_commit);
end
fprintf('  merge_overlay_rows: rows {%s} from %s (build %s) -> %s\n', ...
    strjoin(rows, ','), ov_results, ov_commit(1:min(7, numel(ov_commit))), out_dir);
end

function name = newest_csv(d, pattern)
    f = dir(fullfile(d, pattern));
    f = f([f.bytes] > 0);
    name = '';
    if isempty(f), return; end
    [~, ord] = sort({f.name});     % timestamped names -> newest last
    name = f(ord(end)).name;
end

function c = header_commit(p)
    c = 'unstamped';
    lines = readlines(p);
    for i = 1:numel(lines)
        L = char(lines(i));
        if isempty(L) || L(1) ~= '#', break; end
        tok = regexp(L, 'RANDLAPACK_GIT_COMMIT=([0-9a-f]{7,40})', 'tokens', 'once');
        if ~isempty(tok), c = tok{1}; return; end
    end
end

function merge_pair(base_path, ov_path, out_dir, rows, ov_commit)
    [~, stem, ext] = fileparts(base_path);
    b = readlines(base_path); b = b(strlength(b) > 0);
    o = readlines(ov_path);   o = o(strlength(o) > 0);
    % Column line = first non-comment line on each side; they must agree.
    bc = find(~startsWith(b, '#'), 1); oc = find(~startsWith(o, '#'), 1);
    if isempty(bc) || isempty(oc) || ~strcmp(b(bc), o(oc))
        error('merge_overlay_rows: column mismatch between %s and %s', base_path, ov_path);
    end
    keep = strings(0, 1);
    for i = oc+1:numel(o)
        alg = extractBefore(o(i), ',');
        if any(strcmp(alg, rows)), keep(end+1, 1) = o(i); end %#ok<AGROW>
    end
    % Rows already present in the base under the same names are dropped in favour
    % of the overlay, so a rerun of the overlay era never doubles them.
    base_rows = b(bc+1:end);
    base_alg  = extractBefore(base_rows, ',');
    base_rows = base_rows(~ismember(base_alg, rows));
    fid = fopen(fullfile(out_dir, [stem ext]), 'w');
    fprintf(fid, '%s\n', b(1:bc-1));
    fprintf(fid, '# OVERLAY_COMMIT=%s rows=%s from=%s (%d rows appended by merge_overlay_rows)\n', ...
        ov_commit, strjoin(rows, ','), ov_path, numel(keep));
    fprintf(fid, '%s\n', b(bc));
    if ~isempty(base_rows), fprintf(fid, '%s\n', base_rows); end
    if ~isempty(keep),      fprintf(fid, '%s\n', keep); end
    fclose(fid);
end
