%% attach_kw_backward_error.m — add be_kw_measured to a results table.
%
% The *_backward_error.csv sidecar records the sketched Karlson-Walden backward
% error (Epperly-Meier-Nakatsukasa 2024, eq. 4.2, relative to ||A||_F) for every
% (algorithm, run). It is a MEASUREMENT, written whenever the reference is built,
% independently of whether the backward-error oracle was driving termination.
% That is what makes it the one accuracy metric comparable across an era that ran
% with be_tol_mult > 0 and one that ran with the oracle off: the results column
% be_final is -1 whenever the oracle was off, but this sidecar is populated in
% both.
%
%   T = attach_kw_backward_error(T, be_path)
%
% Adds T.be_kw_measured. Rows with no sidecar record get NaN, never a silent
% zero, so a missing era renders an EMPTY panel rather than a wrong one.
%
% 2026-09-21 (Max): added so the stopping-rule comparison can show the quantity
% Epperly's rule actually targets, alongside the data error it does not control.

function T = attach_kw_backward_error(T, be_path)

n = height(T);
T.be_kw_measured = nan(n, 1);
if ~isfile(be_path)
    warning('attach_kw_backward_error: no sidecar %s; KW panel will be empty', be_path);
    return;
end

n_skip = count_comment_lines(be_path);
Tb = readtable(be_path, 'NumHeaderLines', n_skip, 'TextType', 'char');
if ~ismember('be_kw_theta', Tb.Properties.VariableNames)
    warning('attach_kw_backward_error: %s has no be_kw_theta column; KW panel will be empty', be_path);
    return;
end

alg_T  = cellstr(T.algorithm);
alg_Tb = cellstr(Tb.algorithm);
join_on_run = ismember('run', T.Properties.VariableNames) && ...
              ismember('run', Tb.Properties.VariableNames);

for i = 1:n
    sel = strcmp(alg_Tb, alg_T{i});
    if join_on_run, sel = sel & (Tb.run == T.run(i)); end
    if any(sel)
        v = Tb.be_kw_theta(sel);
        v = v(isfinite(v) & v >= 0);      % -1 marks "not measured"
        if ~isempty(v), T.be_kw_measured(i) = v(1); end
    end
end

end
