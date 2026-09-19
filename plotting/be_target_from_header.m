function tol = be_target_from_header(results_path)
%BE_TARGET_FROM_HEADER Backward-error termination target echoed in a results CSV.
%   tol = be_target_from_header(results_path) scans the leading '#' provenance
%   lines of a benchmark results CSV for the engine knob echo 'be_tol=<value>'
%   (written by both drivers next to the other solver knobs since 2026-09-18)
%   and returns the value. Returns -1 when the field is absent (eras before the
%   knob existed) or negative (the oracle was off), so callers can test
%   `be_target_from_header(p) > 0` to know whether a run's iteration count means
%   "iterations to the backward-error target" or "iterations to the LS-floor exit".
    tol = -1;
    fid = fopen(results_path, 'r');
    if fid == -1, return; end
    closer = onCleanup(@() fclose(fid));
    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        if isempty(line) || line(1) ~= '#', break; end
        % 'be_tol=' does not match inside 'be_tol_mult=' (different next character).
        tok = regexp(line, 'be_tol=([-+0-9.eE]+)', 'tokens', 'once');
        if ~isempty(tok)
            tol = str2double(tok{1});
            return;
        end
    end
end
