%% merge_gram_arms.m — build a merged two-arm CSV set for one campaign cell.
%
% Since the 0812_d2 era, campaigns run each cell TWICE: the trsm arm (default
% Gram left factor) and the gemm arm (RANDLAPACK_GRAM_LEFT=gemm), in
% <cell_dir>/trsm_left and <cell_dir>/gemm_left. The figures show both CQRRT
% variants side by side, which used to rely on a HAND-merged `*_both` tree
% (gemm-arm CQRRT_linop rows renamed CQRRT_linop_gemmL) — an unaudited manual
% step the 2026-08-27 audit flagged. This helper performs that merge
% deterministically at plot time:
%
%   merged_csv = merge_gram_arms(cell_dir, pattern)
%
% reads the NEWEST non-empty CSV matching `pattern` in each arm, renames the
% gemm arm's exact-token CQRRT_linop rows to CQRRT_linop_gemmL, appends them to
% the trsm arm's rows, and writes <cell_dir>/merged/<trsm csv name>. Returns ''
% when either arm is missing (caller falls back to a flat single-arm cell).
% Sidecar files (same stem, _breakdown/_rounds) are merged the same way when
% present in both arms. The merge is re-done on every call, so a re-pulled arm
% can never leave a stale merged file behind.

function merged_csv = merge_gram_arms(cell_dir, pattern)

merged_csv = '';
arm_t = fullfile(cell_dir, 'trsm_left');
arm_g = fullfile(cell_dir, 'gemm_left');
if ~exist(arm_t, 'dir') || ~exist(arm_g, 'dir'), return; end

ft = newest_csv(arm_t, pattern);
fg = newest_csv(arm_g, pattern);
if isempty(ft) || isempty(fg), return; end

out_dir = fullfile(cell_dir, 'merged');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

merged_csv = merge_pair(fullfile(arm_t, ft), fullfile(arm_g, fg), out_dir);

% Sidecars: same-stem *_breakdown.csv and *_rounds.csv, when both arms have them.
for suffix = {'breakdown', 'rounds'}
    pt = regexprep(pattern, 'results', suffix{1});
    st = newest_csv(arm_t, pt); sg = newest_csv(arm_g, pt);
    if ~isempty(st) && ~isempty(sg)
        merge_pair(fullfile(arm_t, st), fullfile(arm_g, sg), out_dir);
    end
end

end

function name = newest_csv(d, pattern)
    f = dir(fullfile(d, pattern));
    f = f([f.bytes] > 0);
    name = '';
    if isempty(f), return; end
    [~, ord] = sort({f.name});     % timestamped names -> newest last
    name = f(ord(end)).name;
end

function out_name = merge_pair(trsm_path, gemm_path, out_dir)
    [~, stem, ext] = fileparts(trsm_path);
    out_name = [stem, ext];
    t_lines = readlines(trsm_path);
    g_lines = readlines(gemm_path);
    % Drop empty lines (readlines yields one for the trailing newline): a blank
    % line INSIDE the file would end count_comment_lines' header scan early.
    t_lines = t_lines(strlength(t_lines) > 0);
    g_lines = g_lines(strlength(g_lines) > 0);
    % gemm arm: keep only DATA rows whose first field is exactly CQRRT_linop,
    % renamed. Exact-token match so CQRRT_linop_gemmL can never be re-renamed.
    keep = strings(0, 1);
    for i = 1:numel(g_lines)
        L = g_lines(i);
        if startsWith(L, 'CQRRT_linop,')
            keep(end+1, 1) = "CQRRT_linop_gemmL," + extractAfter(L, 'CQRRT_linop,'); %#ok<AGROW>
        end
    end
    fid = fopen(fullfile(out_dir, out_name), 'w');
    fprintf(fid, '%s\n', t_lines{:});
    if ~isempty(keep), fprintf(fid, '%s\n', keep{:}); end
    fclose(fid);
end
