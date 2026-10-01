function gen_fem2_kappa_variants(out_dir, targets, cal_points, do_verify)
%GEN_FEM2_KAPPA_VARIANTS  Generate full-rank FEM2 IR-LSQ matrices at target kappas.
%
% Produces (K, M, V) Matrix Market triples for the operator J = L^{-1} K V
% (L = chol(M)) at a set of TARGET condition numbers, using the smooth
% coefficient cCoeffContrast (Coeff='contrast' in build_gsvd_benchmark_2d).
% Unlike Oleg's cCoeffIllConditioned (near-disconnecting ligament, kappa~2e13
% with a spurious near-null mode), this gives a CLEAN, FULL-RANK operator whose
% kappa is dialed to the target.
%
% CALIBRATION is a power-law fit, NOT slope-1. Measured 2026-06-18, kappa(J)
% scales as ~ rho^0.81 (the smooth field's *effective* contrast is sub-linear in
% rho), so each target's rho comes from a log-log fit over known (rho,kappa)
% points -- default the two we measured. Pass more points to refine.
%
% VERIFICATION is OFF by default and deliberately so: the dense values-only SVD
% of J (75824x8304, ~5 GB) is the slow step (minutes-to-hours under memory
% pressure), and it is REDUNDANT -- the ISAAC irlsq_reg benchmark reports the
% achieved condition number as `kappa_measured` (= kappa(R)) in its CSV (read the
% `dd` cells). So we generate fast (3 builds + writes, no SVD) and read the true
% kappa back from the run. Pass do_verify=true to confirm locally (adds one SVD
% per target; uses a fill-reduced Cholesky to keep it as fast as possible).
%
% Inputs:
%   out_dir    : where the .mtx land (default matrices/contrast/)
%   targets    : target condition numbers (default [1e7 1e9 1e11])
%   cal_points : Nx2 [rho kappa] measurements that fix the power law.
%                Default [1e2 1.711e4; 5.8445e4 3.006e6] (measured 2026-06-18).
%   do_verify  : logical, run the dense SVD check per target (default false).
%
% Usage:
%   gen_fem2_kappa_variants;                                   % fast, no SVD
%   gen_fem2_kappa_variants([], [1e7 1e9 1e11], [], true);     % + local verify
%   gen_fem2_kappa_variants([], 1e11);                         % single target

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(out_dir)
    % 2026-09-21: matrices/ consolidation (was input_matrices next to this script)
    out_dir = fullfile(fileparts(here), 'matrices', 'contrast');
end
if nargin < 2 || isempty(targets),    targets    = [1e7 1e9 1e11]; end
if nargin < 3 || isempty(cal_points), cal_points = [1e2 1.711e4; 5.8445e4 3.006e6]; end
if nargin < 4 || isempty(do_verify),  do_verify  = false; end

addpath(fullfile(here, 'fem_gen'));    % build_gsvd_benchmark_2d
addpath(fullfile(here, 'utils'));      % resize_fem_triple
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% --- Power-law fit: log10(kappa) = a*log10(rho) + b -------------------------
pc    = polyfit(log10(cal_points(:,1)), log10(cal_points(:,2)), 1);
a = pc(1);  b = pc(2);
rho_for = @(k) 10.^((log10(k) - b) / a);
fprintf('=== Calibration over %d point(s): log10(kappa) = %.4f*log10(rho) + %.4f ===\n', ...
        size(cal_points,1), a, b);
fprintf('    points:'); fprintf(' (rho=%.3g,kappa=%.3g)', cal_points.'); fprintf('\n');
fprintf('    verify=%d  (kappa is also reported by the ISAAC benchmark as kappa_measured)\n\n', do_verify);

for i = 1:numel(targets)
    kap_t = targets(i);
    rho_t = rho_for(kap_t);
    label = sprintf('1e%d', round(log10(kap_t)));
    base  = sprintf('FEM_Problem_2_kappa%s', label);
    extrap = rho_t < min(cal_points(:,1)) || rho_t > max(cal_points(:,1));
    fprintf('=== Target kappa = %.0e  ->  rho = %.4e%s   (basename %s) ===\n', ...
            kap_t, rho_t, tern(extrap,'  [EXTRAPOLATED]',''), base);

    data = build_fem2_contrast(rho_t);
    K = data.K;  M = data.M;  V = data.P;
    [m, n] = size(V);

    if do_verify
        [kap_ach, smin, smax] = verify_kappa(K, M, V);
        fprintf('  achieved kappa = %.3e (target %.0e), sigma_max=%.3e sigma_min=%.3e\n', ...
                kap_ach, kap_t, smax, smin);
    end

    suffix = sprintf('%dx%d', m, n);
    Kf = fullfile(out_dir, sprintf('%s_%s_K.mtx', base, suffix));
    Mf = fullfile(out_dir, sprintf('%s_%s_M.mtx', base, suffix));
    Vf = fullfile(out_dir, sprintf('%s_%s_V.mtx', base, suffix));
    fprintf('  writing small triple (%s) ...\n', suffix);
    write_mtx_sym(Kf, K);  write_mtx_sym(Mf, M);  write_mtx_gen(Vf, V);

    resize_fem_triple(Kf, Mf, Vf, out_dir, 'expand', 2);   % medium 151648x16608
    resize_fem_triple(Kf, Mf, Vf, out_dir, 'expand', 4);   % large  303296x33216
    fprintf('\n');
end

fprintf('Done: %d kappa-variant(s) x 3 sizes written to %s\n', numel(targets), out_dir);
fprintf('Next: rsync to ISAAC $DATA; submit. Read achieved kappa from the dd-cell CSV\n');
fprintf('  column kappa_measured. If a target is off, add its (rho,kappa) to cal_points\n');
fprintf('  and re-run just that target.\n');
end


% ===========================================================================
% Local helpers
% ===========================================================================

function data = build_fem2_contrast(rho)
% build_gsvd_benchmark_2d with the smooth-contrast coefficient at the production
% mesh resolution (matches FEM_Problem_generator.m struct-call form).
    opts = struct('HmaxFine',   0.01, ...
                  'HmaxCoarse', 0.03, ...
                  'AlphaReg',   1e-14, ...   % unused (we don't form Kmetric)
                  'Plot',       false, ...
                  'Coeff',      'contrast', ...
                  'Rho',        rho);
    data = build_gsvd_benchmark_2d(opts);
end

function [kap, smin, smax] = verify_kappa(K, M, V)
% Dense values-only SVD of J = L^{-1} K V. Uses a fill-reduced Cholesky of M and
% a row permutation (singular values are permutation-invariant) to cut the
% Cholesky/solve fill, and forces a dense RHS so the triangular solve is BLAS-3
% rather than a near-full sparse-fill solve. This is the ONLY accurate sigma_min
% route here (the Gram/eig route squares kappa and dies above ~1e8).
    p  = symamd(M);
    Lp = chol(M(p,p), 'lower');             % fill-reduced factor
    Jp = Lp \ full(K(p,:) * V);             % dense -> BLAS-3 solve; same sing. values
    s  = svd(Jp);                           % one output => values only (gesdd 'N')
    smax = s(1);  smin = s(end);  kap = smax / smin;
end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end

function write_mtx_sym(fname, A)
% Sparse symmetric -> Matrix Market (lower-triangle only). Matches
% FEM_Problem_generator.m / resize_fem_triple.m exactly.
    [rows, cols, vals] = find(tril(A));
    fid = fopen(fname, 'w');
    if fid < 0, error('Cannot open %s for writing', fname); end
    fprintf(fid, '%%%%MatrixMarket matrix coordinate real symmetric\n%%\n');
    fprintf(fid, '%d %d %d\n', size(A,1), size(A,2), numel(vals));
    fprintf(fid, '%d %d %.15E\n', [rows(:)'; cols(:)'; vals(:)']);
    fclose(fid);
end

function write_mtx_gen(fname, A)
% Sparse general -> Matrix Market.
    [rows, cols, vals] = find(A);
    fid = fopen(fname, 'w');
    if fid < 0, error('Cannot open %s for writing', fname); end
    fprintf(fid, '%%%%MatrixMarket matrix coordinate real general\n%%\n');
    fprintf(fid, '%d %d %d\n', size(A,1), size(A,2), numel(vals));
    fprintf(fid, '%d %d %.15E\n', [rows(:)'; cols(:)'; vals(:)']);
    fclose(fid);
end
