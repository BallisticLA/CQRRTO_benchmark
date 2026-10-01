function gen_fem2_ill(out_dir, do_verify)
%GEN_FEM2_ILL  Generate Oleg's ORIGINAL ill-conditioned FEM2 operator (Coeff='ill').
%
% Writes (K, M, V) Matrix Market triples for J = L^{-1} K V (L = chol(M)) using
% build_gsvd_benchmark_2d's 'ill' coefficient -- Oleg's original
% cCoeffIllConditioned: a near-disconnecting ligament on the complex plate
% geometry, documented there as kappa ~ 2e13 by design (the GSVD stress test).
% This is the closest reproducible stand-in for "his unmodified original"; the
% literal FEM_Problem_2.mat he sent is not on disk and was only kappa 1.4e3.
%
% Produces EXACTLY ONE operator -- whatever the generator natively outputs at its
% production mesh -- and nothing else:
%   FEM_Problem_2_ill_75824x8304_{K,M,V}.mtx
%
% Deliberately NO size ladder. gen_fem2_kappa_variants.m follows its small triple
% with resize_fem_triple 'expand' 2 and 4, but those larger sizes are OUR Kron
% enlargements, not anything Oleg's generator produces. The point of this dataset
% is to run his unmodified original, so expanding it would defeat the purpose.
% (If a bigger version is ever wanted, call resize_fem_triple on the output
% explicitly and label it as ours, not as his.)
%
% REQUIRES the MATLAB PDE Toolbox (createpde / generateMesh / decsg / ...).
% Installed 2026-07-27 via:
%   mpm install --release=R2026a --destination=/usr/local/MATLAB/R2026a \
%               --products Partial_Differential_Equation_Toolbox
%
% ---------------------------------------------------------------------------
% EXPECTATION, measured 2026-07-27 on coarse proxy meshes BEFORE running this:
%
%   Hmax(F/C)     m x n      kappa      kappa^colnorm   colspread   sigma_min
%   0.10/0.30    805x107    2.42e9         9.62e4        5.72e8      6.8e-2
%   0.06/0.18   2156x284    9.61e9         1.19e6        3.11e9      4.8e-2
%   0.04/0.12   4811x578    1.10e13        2.03e3        8.74e12     9.2e-5
%
% The nominal kappa is genuinely huge, but it is almost entirely COLUMN SCALING
% (colspread ~ kappa), so kappa^colnorm -- what CholeskyQR actually feels, since
% it factors diagonal column scaling out -- stays at 1e3..1e6 and is NOT monotone
% in mesh refinement. Expect this dataset to be EASY for CholQR despite its
% headline condition number. It is worth running as Oleg's original reference,
% but the genuinely hard test is the butterfly family at kappa^colnorm = 1e10
% (gen_fem2_hard_variants), where colspread ~ 1.
%
% Also expect sigma_min to fall with refinement (9.2e-5 at the finest proxy):
% the ligament produces a near-null mode, so the production-mesh operator may be
% effectively RANK-DEFICIENT. That is the known reason we moved off this
% generator in the first place; breakdowns in the benchmark are informative here,
% not a bug -- read qr_status and the new ir_inner_capped column.
% ---------------------------------------------------------------------------
%
% Inputs:
%   out_dir   : where the .mtx land (default matrices/ill/)
%   do_verify : logical, dense values-only SVD check of the small triple
%               (default false). At 75824x8304 that is a ~5 GB dense solve --
%               do NOT enable it on a memory-constrained machine.
%
% Usage:
%   gen_fem2_ill;                      % generate, no local SVD
%   gen_fem2_ill([], true);            % + verify kappa locally (heavy)

here = fileparts(mfilename('fullpath'));
repo_root = fileparts(here);   % 2026-09-21: matrices/ consolidation (was input_matrices_ill)
if nargin < 1 || isempty(out_dir),   out_dir   = fullfile(repo_root, 'matrices', 'ill'); end
if nargin < 2 || isempty(do_verify), do_verify = false; end
addpath(fullfile(here, 'fem_gen'));    % build_gsvd_benchmark_2d
addpath(fullfile(here, 'utils'));      % (not needed for the single triple; kept so a
                                       %  caller can reach resize_fem_triple by hand)
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

if exist('createpde', 'file') ~= 2
    error('gen_fem2_ill:noPDEToolbox', ...
          ['This generator needs the PDE Toolbox (createpde not found). Install with:\n' ...
           '  mpm install --release=R2026a --destination=/usr/local/MATLAB/R2026a ' ...
           '--products Partial_Differential_Equation_Toolbox']);
end

base = 'FEM_Problem_2_ill';
fprintf('=== Oleg original ill-conditioned FEM2 (Coeff=''ill'') -> %s ===\n', out_dir);

t0 = tic;
opts = struct('HmaxFine',   0.01, ...     % production mesh, same as the contrast family
              'HmaxCoarse', 0.03, ...
              'AlphaReg',   1e-14, ...    % unused (Kmetric not formed)
              'Plot',       false, ...
              'Coeff',      'ill');
data = build_gsvd_benchmark_2d(opts);
fprintf('  assembled in %.1f min\n', toc(t0)/60);

K = data.K;  M = data.M;  V = data.P;
[m, n] = size(V);
fprintf('  small triple: m=%d n=%d nnz(V)=%d\n', m, n, nnz(V));

if do_verify
    p  = symamd(M);
    Lp = chol(M(p,p), 'lower');
    Jp = Lp \ full(K(p,:) * V);
    s  = svd(Jp);
    cn = vecnorm(Jp);
    sn = svd(Jp ./ cn);
    fprintf('  kappa=%.3e  kappa^colnorm=%.3e  colspread=%.3e  sigma_min=%.3e\n', ...
            s(1)/s(end), sn(1)/sn(end), max(cn)/min(cn), s(end));
end

suffix = sprintf('%dx%d', m, n);
Kf = fullfile(out_dir, sprintf('%s_%s_K.mtx', base, suffix));
Mf = fullfile(out_dir, sprintf('%s_%s_M.mtx', base, suffix));
Vf = fullfile(out_dir, sprintf('%s_%s_V.mtx', base, suffix));
fprintf('  writing triple (%s) ...\n', suffix);
write_mtx_sym(Kf, K);  write_mtx_sym(Mf, M);  write_mtx_gen(Vf, V);

fprintf('\nDone: 1 operator (3 files) written to %s\n', out_dir);
fprintf('Next: copy to ISAAC $DATA, then run with kappa_target=1.\n');
fprintf('Read the achieved kappa back from the dd-cell CSV column kappa_measured.\n');
end


% ===========================================================================
% Writers -- byte-identical to gen_fem2_kappa_variants.m / resize_fem_triple.m
% ===========================================================================

function write_mtx_sym(fname, A)
% Sparse symmetric -> Matrix Market (lower triangle only).
    A = tril(sparse(A));
    [i, j, v] = find(A);
    fid = fopen(fname, 'w');
    fprintf(fid, '%%%%MatrixMarket matrix coordinate real symmetric\n');
    fprintf(fid, '%d %d %d\n', size(A,1), size(A,2), numel(v));
    fprintf(fid, '%d %d %.17g\n', [i, j, v].');
    fclose(fid);
end

function write_mtx_gen(fname, A)
% Sparse general -> Matrix Market.
    A = sparse(A);
    [i, j, v] = find(A);
    fid = fopen(fname, 'w');
    fprintf(fid, '%%%%MatrixMarket matrix coordinate real general\n');
    fprintf(fid, '%d %d %d\n', size(A,1), size(A,2), numel(v));
    fprintf(fid, '%d %d %.17g\n', [i, j, v].');
    fclose(fid);
end
