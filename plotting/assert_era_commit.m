function commit = assert_era_commit(camp_dir, era)
%ASSERT_ERA_COMMIT  Refuse to plot an era whose CSVs straddle more than one build.
%
%   commit = assert_era_commit(camp_dir, era) scans every *_results.csv under
%   camp_dir (all subdirectories, merged/ included), extracts the
%   RANDLAPACK_GIT_COMMIT each writer stamped into its comment header, and
%   errors unless every stamped file carries the SAME commit. Returns that
%   commit (full SHA as char; 'unstamped' for a uniformly pre-stamp era).
%
%   Why this exists (2026-08-29): per-file provenance makes a mixed-build era
%   DETECTABLE, not harmless. An era whose cells ran against two builds is no
%   longer internally comparable, yet it looks complete -- right cell count,
%   right names, plausible numbers, every row honestly stamped. Nothing stops
%   it being plotted as one thing unless something compares the stamp ACROSS
%   files and refuses. This is that check. Precedent: the 2026-07-31
%   FunNystromPP figures rendered clean from a superseded campaign and were
%   only caught on 08-03 -- "renders clean" is not "renders the right data".
%
%   Stamp policy: writers stamp RANDLAPACK_GIT_COMMIT since 2026-08-27.
%     * all files stamped, one commit  -> pass, return it
%     * stamped files disagree         -> ERROR (mixed-build era)
%     * stamped and unstamped mixed    -> ERROR (era straddles the stamping
%                                         boundary; needs a human look)
%     * all files unstamped            -> WARN and return 'unstamped' (a
%                                         retired pre-08-27 era re-enabled for
%                                         comparison; single-build cannot be
%                                         verified, only tolerated)

f = dir(fullfile(camp_dir, '**', '*_results.csv'));
f = f([f.bytes] > 0);
if isempty(f)
    error('assert_era_commit:noData', ...
        '%s: no non-empty *_results.csv under %s', era, camp_dir);
end

commits = strings(numel(f), 1);
for k = 1:numel(f)
    p = fullfile(f(k).folder, f(k).name);
    fid = fopen(p, 'r');
    if fid < 0
        error('assert_era_commit:unreadable', '%s: cannot open %s', era, p);
    end
    cleaner = onCleanup(@() fclose(fid));
    while true
        ln = fgetl(fid);
        if ~ischar(ln) || isempty(ln) || ln(1) ~= '#'
            break;   % header block ended without a stamp
        end
        tok = regexp(ln, 'RANDLAPACK_GIT_COMMIT=([0-9a-f]{7,40})', 'tokens', 'once');
        if ~isempty(tok)
            commits(k) = tok{1};
            break;
        end
    end
    clear cleaner;
end

stamped = commits ~= "";
uniq = unique(commits(stamped));

if ~any(stamped)
    warning('assert_era_commit:unstamped', ...
        ['%s: none of the %d CSVs carry a RANDLAPACK_GIT_COMMIT stamp ' ...
         '(pre-2026-08-27 era). Single-build consistency CANNOT be verified.'], ...
        era, numel(f));
    commit = 'unstamped';
    return;
end

if any(~stamped)
    bad = f(~stamped);
    list = join(string(fullfile({bad.folder}, {bad.name})), newline);
    error('assert_era_commit:partialStamp', ...
        ['%s: %d of %d CSVs have no commit stamp while the rest do -- the era ' ...
         'straddles the stamping boundary and may mix builds. Unstamped files:\n%s'], ...
        era, sum(~stamped), numel(f), list);
end

if numel(uniq) ~= 1
    lines = strings(numel(f), 1);
    for k = 1:numel(f)
        lines(k) = sprintf('  %s -> %s', fullfile(f(k).folder, f(k).name), commits(k));
    end
    error('assert_era_commit:mixedBuilds', ...
        ['%s: CSVs straddle %d different builds -- the era is not internally ' ...
         'comparable and will not be plotted as one thing. Re-run the odd cells ' ...
         'or split the eras. File -> commit:\n%s'], ...
        era, numel(uniq), join(lines, newline));
end

commit = char(uniq);
fprintf('  [%s] build provenance OK: %d CSVs, single commit %s\n', ...
    era, numel(f), commit(1:min(7, numel(commit))));
end
