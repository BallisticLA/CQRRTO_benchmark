%% resize_fem_triple.m — Generate a resized (K, M, V) FEM matrix triple
%
% Companion to resize_fem_pair.m (which only resizes K + V).  Use this for
% the post-Oleg-correction FEM2 setup where the benchmark operator is
% J = L^{-1} * K * V with L = chol(M), so we need to resize M alongside K
% and V consistently.
%
% Kron-Cholesky identity:  chol(kron(I_k, M)) = kron(I_k, chol(M)).
% Hence for the 'expand' mode (K_big = kron(I_k, K), M_big = kron(I_k, M),
% V_big = kron(I_k, V)), we get  J_big = kron(I_k, J)  — block-diagonal,
% preserves kappa(J_big) = kappa(J) exactly.
% For 'expand_rows' (V_big = repmat(V, k, 1)), J_big = repmat(J, k, 1) and
% the singular values scale by sqrt(k) with ratio unchanged.
%
% Modes that scale both m and n (m/n ratio preserved):
%
%   resize_fem_triple(K, M, V, output_dir, 'shrink', p)
%       Extract the p×p principal submatrices of K and M; first p rows of V.
%       Drops V columns that become entirely zero after the row cut.
%       Output: <base>_<m_out>x<n_out>_{K,M,V}.mtx
%
%   resize_fem_triple(K, M, V, output_dir, 'expand', k)
%       K_big = kron(I_k, K),  M_big = kron(I_k, M),  V_big = kron(I_k, V).
%       kappa(L^{-1} K V) exactly preserved.
%
%   resize_fem_triple(K, M, V, output_dir, 'target', m_target)
%       Expand+trim to exactly m_target rows; n scales with k.
%       kappa exactly preserved if m_target is a multiple of m.
%
% Modes that scale only m (n held fixed, kappa exactly preserved):
%
%   resize_fem_triple(K, M, V, output_dir, 'expand_rows', k)
%       K_big = kron(I_k, K),  M_big = kron(I_k, M),  V_big = repmat(V, k, 1).
%       Singular values of L_big^{-1} K_big V_big scale by sqrt(k); ratio unchanged.
%
%   resize_fem_triple(K, M, V, output_dir, 'target_rows', m_target)
%       Expand+trim to exactly m_target rows; n always held fixed.
%       kappa exactly preserved if m_target is a multiple of m.
%
%   <base> is derived from K_file by stripping the trailing _K.mtx suffix and
%   any trailing _<digits>x<digits>.  All output filenames encode the actual
%   output dimensions as <m>x<n>.
%
% Example:
%   dir = fullfile(repo_root, 'matrices');   % (was output_matrices/, retired 2026-07-27)
%   resize_fem_triple( ...
%       fullfile(dir, 'FEM_Problem_2_corrected_75824x8304_K.mtx'), ...
%       fullfile(dir, 'FEM_Problem_2_corrected_75824x8304_M.mtx'), ...
%       fullfile(dir, 'FEM_Problem_2_corrected_75824x8304_V.mtx'), ...
%       dir, 'expand', 6);    % -> ...454944x49824_{K,M,V}.mtx

function resize_fem_triple(K_file, M_file, V_file, output_dir, mode, param)

if ~ismember(mode, {'shrink', 'expand', 'expand_rows', 'target', 'target_rows'})
    error('resize_fem_triple: unknown mode ''%s''', mode);
end
if param ~= floor(param) || param <= 0
    error('resize_fem_triple: param must be a positive integer');
end

% --- read inputs ---
fprintf('Reading %s ... ', K_file);
K = read_mtx_sym(K_file);
fprintf('(%d x %d, nnz=%d)\n', size(K,1), size(K,2), nnz(K));

fprintf('Reading %s ... ', M_file);
M = read_mtx_sym(M_file);
fprintf('(%d x %d, nnz=%d)\n', size(M,1), size(M,2), nnz(M));

fprintf('Reading %s ... ', V_file);
V = read_mtx_gen(V_file);
fprintf('(%d x %d, nnz=%d)\n', size(V,1), size(V,2), nnz(V));

[m, n] = size(V);
assert(size(K,1) == m, 'K row count (%d) does not match V row count (%d)', size(K,1), m);
assert(size(M,1) == m, 'M row count (%d) does not match V row count (%d)', size(M,1), m);

% --- apply operation ---
switch mode
    case 'shrink'
        p = param;
        if p >= m
            error('shrink size p=%d must be less than m=%d', p, m);
        end
        fprintf('\nShrinking to %d x %d ...\n', p, p);
        K_new = K(1:p, 1:p);
        M_new = M(1:p, 1:p);
        V_new = V(1:p, :);
        keep  = any(V_new, 1);
        V_new = V_new(:, keep);
        n_new = full(sum(keep));
        fprintf('  K: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, m, nnz(K), p, p, nnz(K_new));
        fprintf('  M: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, m, nnz(M), p, p, nnz(M_new));
        fprintf('  V: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d] (%d zero cols dropped)\n', ...
            m, n, nnz(V), p, n_new, nnz(V_new), n - n_new);
        m_out = p; n_out = n_new;

    case 'expand'
        k = param;
        fprintf('\nExpanding by factor %d ...\n', k);
        Ik = speye(k);
        K_new = kron(Ik, K);
        M_new = kron(Ik, M);
        V_new = kron(Ik, V);
        fprintf('  K: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, m, nnz(K), k*m, k*m, nnz(K_new));
        fprintf('  M: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, m, nnz(M), k*m, k*m, nnz(M_new));
        fprintf('  V: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, n, nnz(V), k*m, k*n, nnz(V_new));
        fprintf('  Note: kappa(L^{-1} K V) is exactly preserved.\n');
        m_out = k*m; n_out = k*n;

    case 'expand_rows'
        k = param;
        fprintf('\nExpanding rows by factor %d (n held fixed) ...\n', k);
        Ik = speye(k);
        K_new = kron(Ik, K);
        M_new = kron(Ik, M);
        V_new = repmat(V, k, 1);
        fprintf('  K: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, m, nnz(K), k*m, k*m, nnz(K_new));
        fprintf('  M: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, m, nnz(M), k*m, k*m, nnz(M_new));
        fprintf('  V: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, n, nnz(V), k*m, n, nnz(V_new));
        fprintf('  Note: kappa(L^{-1} K V) is exactly preserved (singular values scale by sqrt(k)).\n');
        m_out = k*m; n_out = n;

    case 'target'
        m_target = param;
        if m_target == m
            error('target size equals current m=%d; nothing to do', m);
        elseif m_target < m
            % pure shrink
            fprintf('\nShrinking to target %d x %d ...\n', m_target, m_target);
            K_new = K(1:m_target, 1:m_target);
            M_new = M(1:m_target, 1:m_target);
            V_new = V(1:m_target, :);
            keep  = any(V_new, 1);
            V_new = V_new(:, keep);
            n_new = full(sum(keep));
        else
            % expand then trim
            k = ceil(m_target / m);
            fprintf('\nTargeting m=%d: expanding by k=%d (to %d), then trimming ...\n', ...
                m_target, k, k*m);
            Ik = speye(k);
            K_tmp = kron(Ik, K);
            M_tmp = kron(Ik, M);
            V_tmp = kron(Ik, V);
            if k * m == m_target
                K_new = K_tmp;  M_new = M_tmp;  V_new = V_tmp;
                n_new = k * n;
                fprintf('  Exact multiple — no trim needed. kappa(L^{-1} K V) exactly preserved.\n');
            else
                K_new = K_tmp(1:m_target, 1:m_target);
                M_new = M_tmp(1:m_target, 1:m_target);
                V_new = V_tmp(1:m_target, :);
                keep  = any(V_new, 1);
                V_new = V_new(:, keep);
                n_new = full(sum(keep));
                fprintf('  Trimmed %d rows (%.1f%% of last copy). kappa may differ slightly.\n', ...
                    k*m - m_target, 100*(k*m - m_target)/m);
            end
        end
        fprintf('  K: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, m, nnz(K), m_target, m_target, nnz(K_new));
        fprintf('  M: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, m, nnz(M), m_target, m_target, nnz(M_new));
        fprintf('  V: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, n, nnz(V), m_target, n_new, nnz(V_new));
        m_out = m_target; n_out = size(V_new, 2);

    case 'target_rows'
        m_target = param;
        if m_target == m
            error('target size equals current m=%d; nothing to do', m);
        elseif m_target < m
            fprintf('\nShrinking to target %d x %d (n held fixed) ...\n', m_target, m_target);
            K_new = K(1:m_target, 1:m_target);
            M_new = M(1:m_target, 1:m_target);
            V_new = V(1:m_target, :);
            keep  = any(V_new, 1);
            V_new = V_new(:, keep);
            n_new = full(sum(keep));
        else
            k = ceil(m_target / m);
            fprintf('\nTargeting m=%d (n held fixed): expanding by k=%d (to %d), then trimming ...\n', ...
                m_target, k, k*m);
            Ik = speye(k);
            K_tmp = kron(Ik, K);
            M_tmp = kron(Ik, M);
            V_tmp = repmat(V, k, 1);
            if k * m == m_target
                K_new = K_tmp;  M_new = M_tmp;  V_new = V_tmp;
                fprintf('  Exact multiple — no trim needed. kappa(L^{-1} K V) exactly preserved.\n');
            else
                K_new = K_tmp(1:m_target, 1:m_target);
                M_new = M_tmp(1:m_target, 1:m_target);
                V_new = V_tmp(1:m_target, :);
                keep  = any(V_new, 1);
                V_new = V_new(:, keep);
                fprintf('  Trimmed %d rows (%.1f%% of last copy). kappa may differ slightly.\n', ...
                    k*m - m_target, 100*(k*m - m_target)/m);
            end
            n_new = size(V_new, 2);
        end
        fprintf('  K: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, m, nnz(K), m_target, m_target, nnz(K_new));
        fprintf('  M: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, m, nnz(M), m_target, m_target, nnz(M_new));
        fprintf('  V: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, n, nnz(V), m_target, n_new, nnz(V_new));
        m_out = m_target; n_out = size(V_new, 2);
end

% --- derive output filenames ---
% Extract a short problem name: strip _K.mtx and any trailing _<m>x<n>.
% e.g. "FEM_Problem_2_corrected_75824x8304_K" -> "FEM_Problem_2_corrected"
[~, K_name, ~] = fileparts(K_file);
base = K_name;
if length(base) > 2 && strcmp(base(end-1:end), '_K')
    base = base(1:end-2);
end
base = regexprep(base, '_\d+x\d+$', '');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
K_out = fullfile(output_dir, sprintf('%s_%dx%d_K.mtx', base, m_out, n_out));
M_out = fullfile(output_dir, sprintf('%s_%dx%d_M.mtx', base, m_out, n_out));
V_out = fullfile(output_dir, sprintf('%s_%dx%d_V.mtx', base, m_out, n_out));

% --- write outputs ---
fprintf('\nWriting %s ...\n', K_out); write_mtx_sym(K_out, K_new);
fprintf('Writing %s ...\n', M_out); write_mtx_sym(M_out, M_new);
fprintf('Writing %s ...\n', V_out); write_mtx_gen(V_out, V_new);
fprintf('Done.\n');

end


% =========================================================================
% Local helpers (Matrix Market readers/writers; duplicated from resize_fem_pair)
% =========================================================================

function A = read_mtx_sym(fname)
% Read a symmetric Matrix Market coordinate file.
% Reconstructs the full symmetric matrix from the stored lower-triangle entries.
    [rows, cols, vals, m, n] = read_mtx_coords(fname);
    A = sparse(rows, cols, vals, m, n);
    mask = rows ~= cols;
    A = A + sparse(cols(mask), rows(mask), vals(mask), m, n);
end

function A = read_mtx_gen(fname)
    [rows, cols, vals, m, n] = read_mtx_coords(fname);
    A = sparse(rows, cols, vals, m, n);
end

function [rows, cols, vals, m, n] = read_mtx_coords(fname)
    fid = fopen(fname, 'r');
    if fid < 0, error('Cannot open: %s', fname); end
    line = fgetl(fid);
    while ~isempty(line) && line(1) == '%'
        line = fgetl(fid);
    end
    dims       = sscanf(line, '%d %d %d');
    m          = dims(1);
    n          = dims(2);
    nnz_count  = dims(3);
    data       = fscanf(fid, '%d %d %lf', [3, nnz_count]);
    fclose(fid);
    rows = data(1,:)';
    cols = data(2,:)';
    vals = data(3,:)';
end

function write_mtx_sym(fname, A)
    [rows, cols, vals] = find(tril(A));
    fid = fopen(fname, 'w');
    fprintf(fid, '%%%%MatrixMarket matrix coordinate real symmetric\n%%\n');
    fprintf(fid, '%d %d %d\n', size(A,1), size(A,2), numel(vals));
    fprintf(fid, '%d %d %.15E\n', [rows(:)'; cols(:)'; vals(:)']);
    fclose(fid);
end

function write_mtx_gen(fname, A)
    [rows, cols, vals] = find(A);
    fid = fopen(fname, 'w');
    fprintf(fid, '%%%%MatrixMarket matrix coordinate real general\n%%\n');
    fprintf(fid, '%d %d %d\n', size(A,1), size(A,2), numel(vals));
    fprintf(fid, '%d %d %.15E\n', [rows(:)'; cols(:)'; vals(:)']);
    fclose(fid);
end
