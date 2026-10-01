function [K, M, V] = gen_fem2_small_local(out_dir, rho, nxf, nyf, r, coeff, rotate)
%GEN_FEM2_SMALL_LOCAL  Toolbox-free small FEM2-style (K, M, V) generator.
%
% Pure base-MATLAB P1 FEM on a structured triangulation of [0,4]x[0,2] -- needs
% NO PDE Toolbox, so it runs in WSL MATLAB. Builds the same operator family as
% FEM_Problem_2:  J = L^{-1} K V,  L = chol(M),  V = coarse->fine prolongation,
% with the smooth log-uniform conductivity c = rho^s (s in [-1/2,1/2]) so kappa(J)
% is dialed by the contrast rho -- exactly like cCoeffContrast in
% build_gsvd_benchmark_2d. Sized for LOCAL runs (a few thousand fine DOFs) so the
% whole generate -> C++ benchmark -> plot loop stays on this machine.
%
% It is NOT Oleg's slit/hole geometry -- it's a clean structured mesh, intended
% as a fast, full-rank, kappa-controlled testbed (e.g. to verify the sCholQR3
% shift fix) without burning ISAAC time.
%
% Writes K/M/V .mtx (same format as FEM_Problem_generator) and prints kappa(J),
% which it can afford to measure with a dense SVD at this size.
%
% Usage:
%   gen_fem2_small_local;                              % defaults below
%   gen_fem2_small_local(out_dir, rho)                 % tune contrast (kappa knob)
%   gen_fem2_small_local(out_dir, rho, nxf, nyf, r)    % full control

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'utils'));                 % compute_kappa_variants
% 2026-09-21: matrices/ consolidation (was input_matrices_local under utils/).
if nargin < 1, out_dir = fullfile(fileparts(here), 'matrices', 'local'); end
if nargin < 2 || isempty(rho), rho = 1e5; end
if nargin < 3 || isempty(nxf), nxf = 96; end     % fine cells in x
if nargin < 4 || isempty(nyf), nyf = 48; end     % fine cells in y  ([0,4]x[0,2] -> ~square cells)
if nargin < 5 || isempty(r),   r   = 3;  end      % coarsening factor (V: coarse->fine)
if nargin < 6 || isempty(coeff), coeff = 'smooth'; end
if nargin < 7 || isempty(rotate), rotate = false; end   % false | true | seed(scalar)
assert(mod(nxf, r) == 0 && mod(nyf, r) == 0, 'r must divide both nxf and nyf');

% Coefficient c(x,y):
%   'smooth'   = log-uniform gradient -> pure column scaling (easy: kappa^colnorm ~ tens)
%   'regional' = sharp piecewise blocks (candidate for a HARD kappa^colnorm, i.e.
%                ill-conditioning that survives column normalization)
if strcmpi(coeff, 'regional')
    coeffFn = @(cx, cy) cCoeffRegional(cx, cy, rho);
else
    coeffFn = @(cx, cy) cCoeffContrast(cx, cy, rho);
end

Lx = 4; Ly = 2;
[K, M, freeF] = assemble_fem(nxf, nyf, Lx, Ly, coeffFn);
[Vfull, freeC] = prolongation(nxf, nyf, r);
V = Vfull(freeF, freeC);
V = V(:, any(V, 1));                              % drop any all-zero coarse columns (defensive)

% Rotate the coarse basis V -> V*W (W random orthogonal, reproducibly seeded).
% Leaves kappa(J)=kappa(J*W) unchanged (orthogonal), but smears the column-norm
% spread across all columns so kappa^colnorm(J*W) ~ kappa(J) -- i.e. it converts the
% cheap (column-scaling) ill-conditioning of the coefficient field into GENUINE,
% equilibration-surviving hardness. Verified: coefficient contrast alone cannot do
% this (the coarse projection averages it out); the rotation is what makes the test
% hard. V*W densifies, so keep the base small and Kron-expand for larger sizes.
if ~isequal(rotate, false)
    seed = 0;  if isnumeric(rotate), seed = rotate; end
    rs = RandStream('threefry', 'Seed', seed);
    [W, ~] = qr(randn(rs, size(V, 2)));
    V = full(V) * W;
end

m = size(K, 1);  n = size(V, 2);
fprintf('fine free DOFs m=%d, coarse free DOFs n=%d, rho=%.3e, coeff=%s, rotate=%d\n', ...
        m, n, rho, coeff, ~isequal(rotate, false));

% Report BOTH kappa and the column-equilibrated kappa^colnorm (dense SVD, cheap here).
[kap, kapc, smin, smax, spread] = compute_kappa_variants(K, M, V);
fprintf('kappa(J)=%.3e  kappa^colnorm=%.3e  colspread=%.3e  smin=%.3e smax=%.3e\n', ...
        kap, kapc, spread, smin, smax);

% Write .mtx only when called as a command (no output captured); return-only otherwise.
if nargout == 0
    base   = sprintf('FEM2_local_%s_rho%g', coeff, rho);
    suffix = sprintf('%dx%d', m, n);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    Kf = fullfile(out_dir, sprintf('%s_%s_K.mtx', base, suffix));
    Mf = fullfile(out_dir, sprintf('%s_%s_M.mtx', base, suffix));
    Vf = fullfile(out_dir, sprintf('%s_%s_V.mtx', base, suffix));
    write_mtx_sym(Kf, K);  write_mtx_sym(Mf, M);  write_mtx_gen(Vf, V);
    fprintf('wrote %s_%s_{K,M,V}.mtx to %s\n', base, suffix, out_dir);
end
end


% ===========================================================================
function [K, M, freeF] = assemble_fem(nx, ny, Lx, Ly, coeffFn)
% P1 FEM stiffness + consistent mass on a structured triangulation, homogeneous
% Dirichlet (interior DOFs kept). coeffFn(cx,cy) gives conductivity at element centroids.
nnodes = (nx+1) * (ny+1);
xs = linspace(0, Lx, nx+1);  ys = linspace(0, Ly, ny+1);
nid = @(i, j) j*(nx+1) + i + 1;            % i=0..nx, j=0..ny -> 1-based id

% node coordinates
X = zeros(nnodes, 1);  Y = zeros(nnodes, 1);
for j = 0:ny
    for i = 0:nx
        k = nid(i, j);  X(k) = xs(i+1);  Y(k) = ys(j+1);
    end
end

% triangles: 2 per structured cell
ntri = 2 * nx * ny;
Tri = zeros(ntri, 3);  e = 0;
for j = 0:ny-1
    for i = 0:nx-1
        n00 = nid(i, j);  n10 = nid(i+1, j);  n01 = nid(i, j+1);  n11 = nid(i+1, j+1);
        e = e+1;  Tri(e, :) = [n00 n10 n01];
        e = e+1;  Tri(e, :) = [n10 n11 n01];
    end
end

% assemble
II = zeros(9*ntri, 1);  JJ = II;  SK = II;  SM = II;  idx = 0;
Mloc = [2 1 1; 1 2 1; 1 1 2];
for e = 1:ntri
    v  = Tri(e, :);
    x1 = X(v(1)); y1 = Y(v(1));  x2 = X(v(2)); y2 = Y(v(2));  x3 = X(v(3)); y3 = Y(v(3));
    A  = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
    b  = [y2-y3; y3-y1; y1-y2];                 % d phi / dx numerators
    c  = [x3-x2; x1-x3; x2-x1];                 % d phi / dy numerators
    ce = coeffFn(mean([x1 x2 x3]), mean([y1 y2 y3]));
    Ke = (ce / (4*A)) * (b*b' + c*c');
    Me = (A / 12) * Mloc;
    for a = 1:3
        for bb = 1:3
            idx = idx+1;  II(idx) = v(a);  JJ(idx) = v(bb);
            SK(idx) = Ke(a, bb);  SM(idx) = Me(a, bb);
        end
    end
end
Kf = sparse(II, JJ, SK, nnodes, nnodes);
Mf = sparse(II, JJ, SM, nnodes, nnodes);

% homogeneous Dirichlet -> keep interior nodes
free = true(nnodes, 1);
for j = 0:ny
    for i = 0:nx
        if i==0 || i==nx || j==0 || j==ny, free(nid(i, j)) = false; end
    end
end
freeF = find(free);
K = Kf(freeF, freeF);  M = Mf(freeF, freeF);
K = (K + K') / 2;  M = (M + M') / 2;            % clean tiny asymmetry from roundoff
end


% ===========================================================================
function [Vfull, freeC] = prolongation(nxf, nyf, r)
% Bilinear coarse->fine interpolation on nested structured grids.
nxc = nxf / r;  nyc = nyf / r;
nfine   = (nxf+1) * (nyf+1);
ncoarse = (nxc+1) * (nyc+1);
fid = @(i, j) j*(nxf+1) + i + 1;
cid = @(i, j) j*(nxc+1) + i + 1;

II = zeros(4*nfine, 1);  JJ = II;  S = II;  idx = 0;
for jf = 0:nyf
    for iff = 0:nxf
        ic = floor(iff / r);  jc = floor(jf / r);
        if ic == nxc, ic = nxc - 1; end
        if jc == nyc, jc = nyc - 1; end
        s = (iff - ic*r) / r;  t = (jf - jc*r) / r;     % local coords in [0,1]
        fn = fid(iff, jf);
        cc = [cid(ic, jc),     (1-s)*(1-t);
              cid(ic+1, jc),   s*(1-t);
              cid(ic, jc+1),   (1-s)*t;
              cid(ic+1, jc+1), s*t];
        for q = 1:4
            idx = idx+1;  II(idx) = fn;  JJ(idx) = cc(q, 1);  S(idx) = cc(q, 2);
        end
    end
end
Vfull = sparse(II, JJ, S, nfine, ncoarse);

free = true(ncoarse, 1);
for j = 0:nyc
    for i = 0:nxc
        if i==0 || i==nxc || j==0 || j==nyc, free(cid(i, j)) = false; end
    end
end
freeC = find(free);
end


% ===========================================================================
function c = cCoeffContrast(x, y, rho)
% Smooth log-uniform conductivity, contrast c_max/c_min = rho, no near-disconnect.
s = 0.5 * cos(pi*x/4) .* cos(pi*y/2);
c = rho .^ s;
end

function c = cCoeffRegional(x, y, rho)
% Sharp piecewise-constant conductivity, contrast c_max/c_min = rho
% (c in [rho^-1/2, rho^1/2]). Isolated SOFT islands in a stiff sea + a few stiff
% channels -> many sharp material interfaces and weakly-coupled island modes.
% NO near-disconnecting ligament, so the operator stays full rank. This is the
% Oleg-style regional contrast (his cCoeffIllConditioned minus the x~=1 ligament),
% tested as a candidate for a hard, equilibration-surviving kappa^colnorm.
% Domain [0,4] x [0,2]. Vectorized (accepts scalar or array x,y).
hi = sqrt(rho);  lo = 1/sqrt(rho);
c = hi * ones(size(x));                              % stiff sea
soft = (x>0.40 & x<0.95 & y>0.45 & y<1.05) | ...     % isolated soft islands
       (x>1.60 & x<2.15 & y>0.90 & y<1.55) | ...
       (x>2.85 & x<3.40 & y>0.25 & y<0.85) | ...
       (x>3.10 & x<3.60 & y>1.20 & y<1.75) | ...
       (x>0.90 & x<1.35 & y>1.35 & y<1.85);
c(soft) = lo;
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
