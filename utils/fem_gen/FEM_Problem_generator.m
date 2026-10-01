function FEM_Problem_generator(out_dir, basename)
%FEM_Problem_generator  Generate 2D FEM (K, M, P) test matrices and export
% as Matrix Market files for the IR-LSQ benchmark.
%
% Per Oleg's correction (Slack, 2026-05-25), the saved Problem struct contains:
%   Problem.K = data.K          % stiffness  (pure, NOT the regularized data.Kmetric)
%   Problem.M = data.M          % mass
%   Problem.V = data.P          % sparse prolongation (trial basis)
%   Problem.L = chol(data.M)    % factor of MASS (with fill-reducing permutation)
%
% This produces the operator family
%   J = L^{-1} * K * P
% with L*L' = data.M (so B^{-1} = L^{-1}, A = K, U_h = V).  Because L is the
% factor of M (not K), the composite does NOT algebraically collapse:
%   J' * J = P' * K' * L^{-T} * L^{-1} * K * P = P' * K * M^{-1} * K * P
% which is the standard Petrov-Galerkin coarse-grid quadratic form.
%
% Outputs:
%   <out_dir>/<basename>.mat
%   <out_dir>/<basename>_<m>x<n>_K.mtx     (stiffness,    symmetric)
%   <out_dir>/<basename>_<m>x<n>_M.mtx     (mass,         symmetric)
%   <out_dir>/<basename>_<m>x<n>_V.mtx     (prolongation, general)
%
% Usage:
%   FEM_Problem_generator;                                       % defaults
%   FEM_Problem_generator(out_dir, basename);                    % custom paths

if nargin < 1
    % 2026-09-21: was an absolute machine-specific path to the retired
    % input_matrices/ tree. Resolve relative to this file instead: utils/fem_gen
    % -> utils -> repo root -> matrices/.
    out_dir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'matrices');
end
if nargin < 2
    basename = 'FEM_Problem_2_corrected';
end

% --- Mesh / regularization knobs (unchanged from prior generator) ---------
opts = struct('HmaxFine',   0.01, ...
              'HmaxCoarse', 0.03, ...
              'AlphaReg',   1e-14, ...   % not used here (we no longer use Kmetric)
              'Plot',       false);

fprintf('Building 2D FEM benchmark via build_gsvd_benchmark_2d...\n');
data = build_gsvd_benchmark_2d(opts);

% --- Per Oleg's correction --------------------------------------------------
Problem.K = data.K;          % pure stiffness (NOT data.Kmetric)
Problem.M = data.M;          % mass
Problem.V = data.P;          % prolongation (trial basis)

% --- Cholesky of MASS (not stiffness), with fill-reducing permutation ------
fprintf('Computing chol(data.M) with fill-reducing permutation...\n');
[R, flag, perm] = chol(data.M, 'vector');
if flag ~= 0
    error('Cholesky factorization of data.M failed (flag=%d). Check that M is SPD.', flag);
end
Pc = zeros(1, length(perm));
Pc(perm) = 1:length(perm);
L = R(:, Pc)';
rel_err = norm(L*L' - data.M, 'fro') / norm(data.M, 'fro');
fprintf('  chol density: nnz(R)/n = %.3f\n', nnz(R)/size(R,1));
fprintf('  ||L*L'' - data.M||_F / ||data.M||_F = %.3e\n', rel_err);
Problem.L       = L;
Problem.options = opts;

[m, n] = size(Problem.V);

% --- Save .mat ---------------------------------------------------------------
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
mat_path = fullfile(out_dir, [basename, '.mat']);
fprintf('\nSaving %s\n', mat_path);
save(mat_path, 'Problem', '-v7.3');

% --- Export Matrix Market files --------------------------------------------
suffix = sprintf('%dx%d', m, n);
K_mtx  = fullfile(out_dir, sprintf('%s_%s_K.mtx', basename, suffix));
M_mtx  = fullfile(out_dir, sprintf('%s_%s_M.mtx', basename, suffix));
V_mtx  = fullfile(out_dir, sprintf('%s_%s_V.mtx', basename, suffix));

fprintf('Writing %s\n', K_mtx); write_mtx_sym(K_mtx, Problem.K);
fprintf('Writing %s\n', M_mtx); write_mtx_sym(M_mtx, Problem.M);
fprintf('Writing %s\n', V_mtx); write_mtx_gen(V_mtx, Problem.V);
fprintf('\nDone.\n');

end


% ===========================================================================
% Local helpers (Matrix Market writers)
% ===========================================================================

function write_mtx_sym(fname, A)
% Write a sparse symmetric matrix in Matrix Market format (lower-triangle only).
    [rows, cols, vals] = find(tril(A));
    fid = fopen(fname, 'w');
    if fid < 0, error('Cannot open %s for writing', fname); end
    fprintf(fid, '%%%%MatrixMarket matrix coordinate real symmetric\n%%\n');
    fprintf(fid, '%d %d %d\n', size(A,1), size(A,2), numel(vals));
    fprintf(fid, '%d %d %.15E\n', [rows(:)'; cols(:)'; vals(:)']);
    fclose(fid);
end


function write_mtx_gen(fname, A)
% Write a sparse general (non-symmetric) matrix in Matrix Market format.
    [rows, cols, vals] = find(A);
    fid = fopen(fname, 'w');
    if fid < 0, error('Cannot open %s for writing', fname); end
    fprintf(fid, '%%%%MatrixMarket matrix coordinate real general\n%%\n');
    fprintf(fid, '%d %d %d\n', size(A,1), size(A,2), numel(vals));
    fprintf(fid, '%d %d %.15E\n', [rows(:)'; cols(:)'; vals(:)']);
    fclose(fid);
end
