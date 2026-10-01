function gen_fem2_hard_variants(out_dir, targets, nxf, nyf, r, seed, nstages)
%GEN_FEM2_HARD_VARIANTS  Hard, full-rank IR-LSQ operators at target kappa^colnorm.
%
% Fully local (NO PDE Toolbox). Builds J = L^{-1} K V from a UNIFORM-coefficient P1
% FEM (genuine, well-conditioned K, M) whose coarse basis carries the difficulty:
%
%     V_hard = (V .* D) * W,   D = diag(g.^linspace(0,1,n)),   W = sparse butterfly.
%
% D injects the target conditioning ALGEBRAICALLY (a geometric column scaling -- the
% same axis the benchmark's own geometric_colscale uses), and the butterfly rotation
% W spreads it across all columns so it SURVIVES column equilibration:
% kappa^colnorm(J) ~ g. This construction is MESH-INDEPENDENT (g sets kappa^colnorm
% directly, no h-dependent calibration), so g = target and a kappa^colnorm measured at
% a small proxy mesh transfers to any size -- unlike a coefficient-contrast knob,
% whose kappa drifts ~2 orders with mesh refinement.
%
% Why not a dense random rotation: a full orthogonal W hardens equally but DENSIFIES V
% (5 GB at original size). A 2-3 stage butterfly (random-paired Givens layers) hardens
% just as well while keeping V*W sparse (fill ~8x, n-independent) -- ~55 MB .mtx at the
% original 75824x8304, writable and shippable. See utils/apply_butterfly.m.
%
% Sizes: assemble the base at (nxf,nyf,r), then Kron-expand x{1,2,4} (block-diagonal,
% preserves kappa AND kappa^colnorm exactly; columns don't mix across blocks).
%
% IMPORTANT: run the benchmark with kappa_target=1 (the default) on these files. The
% conditioning is already baked into V; a kappa_target>1 would re-apply geometric_colscale
% on top and you'd benchmark the wrong operator (the easy axis CholeskyQR ignores).
%
% Defaults produce the near-original campaign sizes {75466x8256, 150932x16512,
% 301864x33024} -- the structured mesh cannot hit the exact FEM_Problem_2 dims
% (75824x8304) because gcd forces unreachable prime coarse-DOF counts; these are the
% closest same-aspect (m/n ~ 9.14) sizes.
%
% Usage:
%   gen_fem2_hard_variants;                                  % near-original defaults
%   gen_fem2_hard_variants(out_dir, [1e7 1e9 1e11], 48,24,3) % small, fast (for tests)

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'utils'));
% 2026-09-21: the three sibling input_matrices*/ trees were consolidated into
% matrices/ on 2026-07-27; this default still pointed at the retired name, so a
% documented regeneration (gen_fem2_hard_variants([], target)) wrote to a stale
% path under utils/. repo_root = parent of utils/.
repo_root = fileparts(here);
if nargin < 1 || isempty(out_dir), out_dir = fullfile(repo_root, 'matrices', 'hard-kcolnorm'); end
if nargin < 2 || isempty(targets), targets = [1e7 1e9 1e11]; end   % target kappa^colnorm
if nargin < 3 || isempty(nxf), nxf = 390; end   % near-original: coarse n=8256, fine m=75466
if nargin < 4 || isempty(nyf), nyf = 195; end
if nargin < 5 || isempty(r),   r   = 3;   end
if nargin < 6 || isempty(seed), seed = 1;  end
if nargin < 7 || isempty(nstages), nstages = 3; end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

svd_ok = size_ok_for_svd(nxf, nyf);   % self-verify kappa^colnorm only when affordable

% --- assemble the well-conditioned base ONCE (uniform coefficient rho=1) ---
fprintf('assembling base FEM (uniform coeff) at %dx%d r%d ...\n', nxf, nyf, r);
[K, M, V] = gen_fem2_small_local([], 1, nxf, nyf, r, 'smooth', false);
n = size(V, 2);
fprintf('base: m=%d n=%d nnz(V)=%d\n\n', size(K,1), n, nnz(V));

for kt = targets
    g  = kt;                                   % kappa^colnorm ~ g (mesh-independent)
    D  = g .^ linspace(0, 1, n);               % geometric column scaling, span [1, g]
    Vh = apply_butterfly(V .* D, nstages, seed);
    thr = 1e-14 * max(abs(nonzeros(Vh)));      % drop denormal fill from the rotation
    Vh(abs(Vh) < thr) = 0;

    base = sprintf('FEM2_hard_kcolnorm1e%d', round(log10(kt)));
    if svd_ok
        [kap, kapc, smin] = compute_kappa_variants(K, M, Vh);
        fprintf('target=%.0e -> kappa^colnorm=%.3e (kappa=%.3e sigma_min=%.2e, full-rank) nnz/col=%.1f\n', ...
                kt, kapc, kap, smin, nnz(Vh)/n);
    else
        fprintf('target=%.0e -> generated (n=%d too large for local SVD; kappa^colnorm~g by the mesh-independent\n', kt, n);
        fprintf('             construction, verify at a small proxy mesh) nnz(V)=%d nnz/col=%.1f\n', nnz(Vh), nnz(Vh)/n);
    end

    for k = [1 2 4]
        Ik = speye(k);  Kk = kron(Ik, K);  Mk = kron(Ik, M);  Vk = kron(Ik, Vh);
        suf = sprintf('%dx%d', size(Vk, 1), size(Vk, 2));
        write_mtx_sym(fullfile(out_dir, sprintf('%s_%s_K.mtx', base, suf)), Kk);
        write_mtx_sym(fullfile(out_dir, sprintf('%s_%s_M.mtx', base, suf)), Mk);
        write_mtx_gen(fullfile(out_dir, sprintf('%s_%s_V.mtx', base, suf)), Vk);
    end
    fprintf('  wrote %s_{base,x2,x4}_{K,M,V}.mtx\n', base);
end
fprintf('\nDone: %d kappa^colnorm-variant(s) x 3 sizes -> %s\n', numel(targets), out_dir);
fprintf('REMINDER: run the benchmark with kappa_target=1 on these files.\n');
end


function ok = size_ok_for_svd(nxf, nyf)
    % dense SVD of J (m x n) is affordable up to a few thousand coarse DOFs
    ok = (nxf * nyf <= 96*48*4);
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
