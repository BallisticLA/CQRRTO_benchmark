%% attach_data_error.m — add ls_data_error (||Ax-b||/||b||) to an IR-LSQ results table.
%
% The results CSV carries the Higham normwise backward error (ls_residual_norm)
% but not the plain relative residual. The per-round sidecar
% (*_irlsq_reg_rounds.csv, present since 2026-08-27) records ls_relres after
% every refinement cycle, so the LAST round of each (algorithm, run) is that
% run's final data error. This helper joins it onto T BEFORE aggregate_runs,
% so whichever run the aggregation keeps carries its own value.
%
%   T = attach_data_error(T, rounds_path)
%
% Rows with no rounds record are the single-shot LSQR rows (Blendenpik and
% Blendenpik_cold, engine_status == -1); for those LSQR's own final relative
% residual, ir_inner_relres, is the same quantity and is used instead. Anything
% else without a record gets NaN, never a silent zero. A missing sidecar file
% yields all-NaN plus a warning, so an old era renders an EMPTY accuracy panel
% rather than a wrong one.
%
% Why not derive it from ls_residual_norm: that needs ||A|| and ||x||, which the
% CSV does not carry.
%
% 2026-09-14 (Max): introduced when the FEM accuracy panel switched from the
% normwise backward error to the data error, so the FEM and Toeplitz figures
% share one metric that reads directly against the 1e-11 noise floor.

function T = attach_data_error(T, rounds_path)

n = height(T);
T.ls_data_error = nan(n, 1);
if ~isfile(rounds_path)
    warning('attach_data_error: no rounds sidecar %s; data-error panel will be empty', rounds_path);
    return;
end

n_skip = count_comment_lines(rounds_path);
Tr = readtable(rounds_path, 'NumHeaderLines', n_skip, 'TextType', 'char');
alg_T  = cellstr(T.algorithm);
alg_Tr = cellstr(Tr.algorithm);
join_on_run = ismember('run', T.Properties.VariableNames) && ...
              ismember('run', Tr.Properties.VariableNames);
lsqr_fallback = ismember('engine_status', T.Properties.VariableNames) && ...
                ismember('ir_inner_relres', T.Properties.VariableNames);

for i = 1:n
    sel = strcmp(alg_Tr, alg_T{i});
    if join_on_run, sel = sel & (Tr.run == T.run(i)); end
    if any(sel)
        rounds = Tr.round(sel);
        relres = Tr.ls_relres(sel);
        [~, j] = max(rounds);           % final refinement cycle of this run
        T.ls_data_error(i) = relres(j);
    elseif lsqr_fallback && T.engine_status(i) == -1
        T.ls_data_error(i) = T.ir_inner_relres(i);   % single-shot LSQR row
    elseif ismember('engine_status', T.Properties.VariableNames) && T.engine_status(i) == 5 ...
            && ismember('x0_relres', T.Properties.VariableNames) && T.x0_relres(i) >= 0
        % Zero-round oracle exit (2026-09-18): the warm start already met the
        % backward-error target, no round ran, and the returned x IS x0.
        T.ls_data_error(i) = T.x0_relres(i);
    end
end

end
