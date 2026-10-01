function [fl, ok, col] = data_error_floor_from_header(results_path, T, rows)
% Attainable data error of a noisy consistent LS problem: the noise component orthogonal to
% range(A), E||e_perp|| = noise * sqrt(1 - n/m). noise comes from the CSV '#' header
% (noise_level= for FEM, relnoise= for Toeplitz); m, n from the table. ok = false when any
% ingredient is missing, so callers fall back to no marking. col names the table's data-error
% column: ls_data_error (FEM, joined by attach_data_error) or data_relres (Toeplitz results CSV).
fl = -1; ok = false; col = '';
if ismember('ls_data_error', T.Properties.VariableNames), col = 'ls_data_error';
elseif ismember('data_relres', T.Properties.VariableNames), col = 'data_relres';
else, return; end
fid = fopen(results_path, 'r'); if fid < 0, return; end
noise = -1;
while true
    line = fgetl(fid);
    if ~ischar(line) || isempty(line) || line(1) ~= '#', break; end
    tok = regexp(line, '(?:noise_level|relnoise)=([-+0-9.eE]+)', 'tokens', 'once');
    if ~isempty(tok), noise = str2double(tok{1}); break; end
end
fclose(fid);
if ~(noise > 0), return; end
if ~all(ismember({'m','n'}, T.Properties.VariableNames)), return; end
m = T.m(rows(1)); n = T.n(rows(1));
if ~(m > n && n > 0), return; end
fl = noise * sqrt(1 - n / m); ok = true;
end
