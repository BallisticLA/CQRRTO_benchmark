function C = be_censoring(T, rows, failed, results_path)
%BE_CENSORING Per-method backward-error verdicts for an era run with the oracle on.
%   C = be_censoring(T, rows, failed, results_path) reads the era's target from the
%   CSV header (be_target_from_header) and, when the oracle was on, classifies the
%   selected rows from the RECORDED final estimate rather than from the exit code:
%     C.active     scalar: the header echoes be_tol > 0
%     C.be_tol     the target (relative to ||A||_F), -1 when off
%     C.ran        the row ran the engine with the oracle (t_be_us >= 0)
%     C.converged  ran and be_final <= be_tol and the exit was not a breakdown
%     C.censored   ran and not converged: it stopped without reaching the target
%     C.no_oracle  live row that never ran the oracle (LSQR rows, rows merged in
%                  from an era without the knob); a count, not a verdict
%     C.be_final   the recorded final estimate (NaN where absent)
%   All logical fields are columns aligned with `rows`. `failed` marks rows whose
%   build failed (qr_status ~= 0); those are neither censored nor no_oracle.
%
%   The engine's own verdict (stop_reason 'be') is cross-checked against the
%   recorded estimate and a disagreement is an error, so a figure can never show
%   a 'be' bar whose estimate sits above the target or vice versa. A 'tol' exit
%   (LS tolerance met first) is judged by be_final alone.
    n = numel(rows);
    C.be_tol    = be_target_from_header(results_path);
    C.active    = C.be_tol > 0;
    C.ran       = false(n, 1);
    C.converged = false(n, 1);
    C.censored  = false(n, 1);
    C.no_oracle = false(n, 1);
    C.be_final  = nan(n, 1);
    if ~C.active, return; end

    need = {'t_be_us', 'be_final', 'stop_reason'};
    missing = need(~ismember(need, T.Properties.VariableNames));
    if ~isempty(missing)
        error('be_censoring: %s echoes be_tol=%g but lacks column(s) %s; the CSV comes from an older binary or the era is mislabelled', ...
              results_path, C.be_tol, strjoin(missing, ', '));
    end
    live        = ~failed(:);
    C.ran       = (T.t_be_us(rows) >= 0); C.ran = C.ran(:) & live;
    C.be_final  = T.be_final(rows);       C.be_final = C.be_final(:);
    reasons     = string(T.stop_reason(rows)); reasons = reasons(:);
    % A negative be_final is the "no measurement" sentinel (older binaries write -1
    % for a cold row that exited before any round); it must never read as "below
    % the target". Newer binaries write +inf for that case, which fails the test
    % on its own. Only a measured value can convict or acquit.
    measured    = C.ran & (C.be_final >= 0);
    C.converged = measured & (C.be_final <= C.be_tol) & reasons ~= "breakdown";
    C.censored  = measured & ~C.converged;
    C.no_oracle = live & ~measured;

    % Consistency: for measured engine exits other than 'tol', 'be' must mean
    % converged and anything else must mean not converged. A 'breakdown' row whose
    % returned iterate met the target is deliberately counted as censored (the
    % factor was unusable), so it is excluded here.
    judged = measured & reasons ~= "tol" & reasons ~= "breakdown";
    bad = judged & ((reasons == "be") ~= C.converged);
    if any(bad)
        names = cellstr(string(T.algorithm(rows(bad))));
        error('be_censoring: stop_reason and be_final disagree for %s in %s', strjoin(names, ', '), results_path);
    end
end
