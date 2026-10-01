% Resolve paths relative to this script's location (works on any machine)
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);                                  % repo root
addpath(fullfile(script_dir, 'utils'));
addpath(fullfile(script_dir, 'plotting'));                            % plot drivers

tab_groups = {};   % {tabgroup handle, export prefix}

%% ========================================================================
%  CAMPAIGN ERAS (rebased 2026-08-07, Max)
%
%  The era ladder rolled forward: every pre-08-06 campaign (the 07-31
%  warm/cold pairs, the 08-05 lsqr rows, the 3kappa/1e6/hard_kcolnorm sets,
%  irlsq_diag, the June isaac-* dirs) is RETIRED and its local data DELETED
%  (2026-08-07 cleanup; the raw output survives on ISAAC under
%  ~/slurm/CQRRT/benchmark-out/). The 08-06 campaign is now [OLD 08-06];
%  [NEW 08-07] is the rerun on the unified restarted_pcg_ne engine.
%
%  Every era is marked in THREE places:
%    1. the figure Name,
%    2. the tab title (prefixed [OLD ...] / [NEW ...]),
%    3. the exported PDF filename prefix.
%  2026-07-30 (Max): era tags REMOVED from the in-plot super-titles -- the tab
%  strip and PDF filename carry the era; the panel title stays clean.
%
%  WHAT CHANGED in the [NEW 08-10] era (RandLAPACK commits 4c32aae -> 60f366a
%  -> d148e73 -> b5a2886 -> 6d026d5 -> a73d03d; read before comparing to OLD):
%    * ONE solver engine: both benchmarks now run restarted_pcg_ne.
%      IterRefineLSQ is a thin adapter over it (bitwise-identical solutions
%      pinned by test); inner_restarts is GONE -- every round is a restart
%      from the returned iterate against the TRUE residual.
%    * Per-round restart pacing: restart_drop default 1e-2 -> 1e-4 (Oleg's
%      pacing), plus an inner_abs_tol guard so rounds stop once below the
%      absolute target. FEM2 CLI slot 13 is ir_round_drop now (>= 1 rejected).
%    * Round caps raised to 20 in BOTH benchmarks (FEM2 ir_n_steps, Toeplitz
%      pcg_max_restarts). The OLD 3-4 round caps truncated the
%      unpreconditioned baseline ~5 orders above tol -- its OLD-era accuracy
%      bars are cap artifacts, not solver quality.
%    * Outer stagnation exit: 2 rounds without true-residual improvement
%      stop the method, so noise-floor methods no longer grind the cap.
%    * BLAS-2 THREAD GUARD (d148e73), the reason 08-07 was rerun as 08-08:
%      MKL threads the preconditioner's triangular solves badly at 64 threads.
%      Calibrated on the benchmark node (Gold 6430): dtrsv is fastest at 8-16
%      threads and degrades past that, so the solves are now capped (8 for
%      n <= 4000, 16 above). The preconditioned methods apply the factor twice
%      per inner iteration and the unpreconditioned baseline not at all, so
%      the overhead taxed exactly the methods that converge fastest. Measured
%      0807 -> 0808 on solve wall-clock: preconditioned methods 1.3-2.0x
%      faster, unpreconditioned ~1.0x (unchanged, as expected).
%      => 08-07 SOLVE TIMINGS ARE RETIRED as thread-contaminated; its accuracy
%      and iteration counts were never affected.
%    * FFT THREAD CAP (b5a2886) -- THE ACTUAL FIX for the wall-clock inversion,
%      and the reason 08-08 was rerun as 08-09. MKL's threaded FFT
%      INTERMITTENTLY STALLS at 64 threads: single DftiComputeForward calls
%      were caught taking 16-32 ms instead of ~0.1 ms, ~99% of the stall inside
%      the FFT call (fill/multiply/backward stay at 80-100 us). A method making
%      276 applies averages those away; one making 11 cannot, so the artifact
%      scaled INVERSELY with preconditioner quality and inverted the ranking.
%      Thread sweep, 3 repeats: at 8 threads CQRRT solves in 19 ms against
%      unpreconditioned's 200 ms (10.6x, matching 5-vs-267 iterations); at 64
%      both read ~344 ms. The transform is now capped (default 16). 08-08 SOLVE
%      TIMINGS ARE THEREFORE ALSO RETIRED; 08-09 is the first trustworthy
%      wall-clock. Accuracy and iteration counts were never affected in any era.
%    * FOUR BLENDENPIK ROWS (6d026d5). As published Blendenpik is sketch + QR +
%      LSQR with no refinement, so comparing it against Q-less methods that all
%      run through refinement conflated preconditioner quality with solver
%      structure. The published warm/cold rows are unchanged; two rows now hand
%      Blendenpik's OWN R and answer to our restarted-PCG engine. Measured:
%      cold Blendenpik's recovery error 3.008e+03 becomes 4.107e-04 with
%      refinement at IDENTICAL preconditioner conditioning, so its failure is
%      solver structure, not a bad preconditioner.
%      (a73d03d, why 08-09 was rerun as 08-10: 6d026d5 added those rows only to
%      the sparse-mode selector, but run_irlsq_reg has its OWN selector, so the
%      08-09 FEM2 cells produced 7 rows instead of 9 despite method_mask=127.
%      The 08-09 Toeplitz data was complete and correct; both families were
%      rerun together so one era means one commit.)
%    * Wall-clock is NOT comparable across eras (round policy changed).
% =========================================================================

% Plotted eras (2026-09-14): FEM2 = [0911 n11], Toeplitz = [0902 ns1]. ONE set of
% runs per benchmark. Both ran RandLAPACK dabde67 (92061e2 plus the
% RANDLAPACK_CHOL_SYMMETRIZE knob), EXCLUSIVE 64 threads, every cell on srm1526
% (the srm1527 twin; same Gold 6430 SKU -- cross-era wall-clock against the
% srm1527 eras carries the twin-node caveat from the 08-31 provenance note).
% The two benchmarks are on DIFFERENT era tags because only the FEM2 family was
% rerun with noise: the Toeplitz cells were never affected by the CLI slot bug
% (they do not take a noise argument on that path) and are unchanged since 09-02.
% Older eras are RETIRED from the plots; their data is retained under results/
% and on ISAAC, so re-adding a row below re-renders them.
%
%   [0911 n11] THE NOISY-RHS FEM2 ERA -- the plotted FEM2 set. Identical to
%              [0902 ns1] in every respect (same commit dabde67, same
%              RANDLAPACK_CHOL_SYMMETRIZE=0 arm, adaptive rescue ACTIVE, d=2n,
%              3 runs, EXCLUSIVE 64 threads, srm1526) EXCEPT that the right-hand
%              side now carries noise: b = A*x_true + noise with
%              noise_level = 1e-11, where every prior FEM2 era ran noise-free.
%              WHY THIS ERA EXISTS: with a noise-free b the least-squares problem
%              is consistent, so Blendenpik's sketch-and-solve initial guess is
%              already accurate to u*kappa and its warm row did no measurable
%              work -- it entered at the floor and stopped. That made the warm
%              Blendenpik row uninformative and did not match the published
%              method (Oleg, 09-07). With noise at 1e-11 the initial guess enters
%              at the noise level instead: measured x0_relres is 1.33e-11 in all
%              three cells, against 3.6e-15 in the noise-free eras, and the warm
%              refined row now runs 3 outer / 17 inner iterations.
%              PROVENANCE WARNING: era [0907 n11] was the FIRST attempt at this
%              and is NOT a noisy era -- gen_irlsq_reg_jobs.sh wrote the noise
%              level into the omega slot (argv[14]) and left the noise slot
%              (argv[13]) at 0, so every 0907 cell ran noise-free. It is retained
%              as a Platinum 8462Y+ hardware control, NOT as a noise era. Any
%              FEM2 era added below must be checked with
%              irlsq_reg/check_irlsq_reg_args.py before it is plotted.
%              Era-wide acceptance: qr_status=0 on all 72 rows (8 methods x 3
%              runs x 3 cells), zero failed factorizations, CQRRT 0 retries in
%              every cell, CholQR rescued (retries=1) in every cell.
%
%   [0902 ns1] THE NO-SYMMETRIZATION ARM (pre-B5 arithmetic on the audited,
%              fully instrumented binary). RANDLAPACK_CHOL_SYMMETRIZE=0: the
%              Gram is factorized with its upper triangle as computed, no
%              (G+G')/2 -- the pre-audit behavior, isolated after the 09-02
%              finding that the symmetrization's ULP-level re-rounding of two
%              knife-edge pivots (native_ill cols 8066/8303) was the sole cause
%              of the 0827->0829 CholQR flip (graded-pivot mechanism; see
%              dev-logs/2026-09-02-b-qless-qr-cholqr-graded-pivot-reinvestigation.md).
%              RANDLAPACK_CHOL_MAX_RETRIES unset: the adaptive shift rescue is
%              ACTIVE (pre-audit policy), so rescued rows carry chol_retries>0
%              and the chol_shift_* columns instead of FAIL stubs. Era-wide:
%              zero failed factorizations; CQRRT 0 retries in all 42 rows;
%              CholQR retries=1 in all 42 (unshifted potrf fails everywhere
%              with symmetrization off, native_ill included); FEM2 CQRRT rows
%              bit-identical to [0827 a1]. NOTE: the paper PRESCRIBES the
%              symmetrization -- this is the diagnostics/robustness arm, not a
%              candidate figure set for the paper without saying so in prose.
%
%   [0831 a1]  THE FIXED-ALGORITHM ERA. Every published Q-less row runs with
%              RANDLAPACK_CHOL_MAX_RETRIES=0, so a Cholesky breakdown reports as
%              a FAIL row (qr_status != 0) instead of silently switching the row
%              to the adaptive-shift-rescued variant of its method. That switch
%              is why earlier eras could not be read at face value: on the
%              FEM2 native_ill cell the POTRF verdict is decided by BLAS
%              summation order, so 64 vs 32 threads changed WHICH ALGORITHM a row
%              measured while the label stayed the same (see
%              randnla/reports/2026-08-31-iterations-vs-walltime-canonical.md).
%              Acceptance-checked by plotting/verify_0831_fixed_algorithm.py,
%              which pre-registers the per-row prediction from [0829 a1]
%              (chol_retries > 0 there => FAIL row here): 30/30 held.
%              CQRRT is unchanged from [0829 a1] (identical iteration counts,
%              QR times within run-to-run noise), which is the control: the knob
%              only ever touches the breakdown path.
%   [0829 a1]  NOT PLOTTED (row commented out below). The last era in which a
%              rescued row could masquerade as its unshifted method; kept for the
%              cross-era comparison the verifier relies on.
%   [0829 b1]  NOT PLOTTED. Same commit and node as 0829 a1 but NON-EXCLUSIVE on
%              32 CPUs: accuracy-grade only, WALL-CLOCK NOT CITABLE. Keep it as
%              the replication record -- it is the era that exposed the
%              thread-count-dependent branch in the first place.
%
% WHAT CHANGED vs [0827 a1] (read before comparing):
%   * cond_precond is REAL for the first time. The old code computed it
%     unconditionally despite the CLI flag; that bug is fixed and the SLURM
%     scripts now pass compute_cond=1 (they passed 0 before, which would now
%     honestly report -1).
%   * The warm Blendenpik row's LSQR stop test was rescaled to a TRUE relative
%     residual (audit item A7), so its iteration count and stop_reason move.
%     That is a correction, not a regression: the old row was not
%     tolerance-comparable to the other methods.
%   * cholqr_primitive now symmetrizes the Gram, (G+G')/2, before potrf, so
%     ULP-level drift against 0827_a1 is expected everywhere.
%   * The paper shift 11*n*eps*trace(G) is now the library DEFAULT; the
%     RANDLAPACK_SCHOLQR3_SHIFT=theory export the scripts still carry is a no-op.
%   * sCholQR3 breakdown CSVs are complete (t0..t17 + t_total); they were
%     truncated mid-iteration-3 before.
% Two-arm cells: each cell dir holds trsm_left/ + gemm_left/; the merge into one
% table (gemm-arm CQRRT_linop rows renamed CQRRT_linop_gemmL) happens HERE at
% plot time via plotting/merge_gram_arms.m (2026-08-27; the old hand-built
% *_both trees were an unaudited manual step). Flat single-arm cells still work.
% Retired era tags kept (unused) so re-enabling a comparison is one line in
% TOEP_CAMPAIGNS / FEM2_CAMPAIGNS rather than an archaeology exercise.
ERA_D2 = '[0812 d2]';  %#ok<NASGU>
NOTE_D2 = 'd=2n; theory shift; CQRRT: TRSM vs GEMM left factor; 1 run, shared node';  %#ok<NASGU>
ERA_A1  = '[0827 a1]';  %#ok<NASGU>
NOTE_A1 = 'audit-remediation rerun; d=2n; theory shift; 3 runs, exclusive node';  %#ok<NASGU>

% The two plotted eras. The b1 note leads with the timing caveat on purpose: it
% lands in the tab title and inside the exported PDF, so a wall-time panel can
% never be read as citable once the figure is separated from this file.
ERA_A2  = '[0829 a1]';  %#ok<NASGU>
NOTE_A2 = 'post-audit + paper-consistency rerun; d=2n; 3 runs; EXCLUSIVE 64 threads; timings citable';  %#ok<NASGU>
ERA_B1  = '[0829 b1]';  %#ok<NASGU>
NOTE_B1 = 'ACCURACY-GRADE (timings NOT citable: shared node, 32 threads); post-audit + paper-consistency rerun; d=2n; 3 runs';  %#ok<NASGU>

% [0831 a1] NOT PLOTTED (rows commented out below): fixed-algorithm rerun with
% RANDLAPACK_CHOL_MAX_RETRIES=0 on every published Q-less row, so a Cholesky
% breakdown reports as a FAIL row instead of silently measuring the
% shift-rescued variant (the 2026-08-31 native_ill thread-count finding).
% Same protocol otherwise: exclusive 64 threads, d=2n, 3 runs. Ran 08-31/09-01
% (6 cells srm1527; toep_large pulled separately); retired from the plots
% 2026-09-04 in favor of [0902 ns1].
ERA_A3  = '[0831 a1]';  %#ok<NASGU>
NOTE_A3 = 'fixed-algorithm rerun (chol_max_retries=0: breakdown = FAIL row); d=2n; 3 runs; EXCLUSIVE 64 threads; timings citable';  %#ok<NASGU>

% [0902 ns1] (PLOTTED): the no-symmetrization arm -- see the era block at the
% top of this section. The note leads with the arm identity on purpose: it lands
% in the tab title and the figure Name. 2026-09-04 (Max): the note is NO LONGER
% written into the in-plot super-title (paper figures carry only the problem
% label); the exported PDF filename prefix still carries the era tag.
ERA_NS1  = '[0902 ns1]';
NOTE_NS1 = 'NO-SYMMETRIZATION ARM (pre-B5: Gram as computed; adaptive shift rescue ACTIVE, rescued rows carry chol_retries>0); d=2n; 3 runs; EXCLUSIVE 64 threads (srm1526); timings citable within-era';

% [0911 n11] (PLOTTED, FEM2 ONLY): the noisy-RHS FEM2 era -- see the era block at
% the top of this section. Same arm and protocol as [0902 ns1]; the ONLY change
% is noise_level=1e-11 on the right-hand side, which is what makes the warm
% Blendenpik row measure anything. Ran 09-12 on srm1526 via long-bigmem, 3 whole
% cells (NOT the per-method split the ai-tenn attempt used).
% [0907 n11] is deliberately absent: despite the tag it ran noise-free (CLI slot
% bug, see the era block), and is kept only as a hardware control.
ERA_N11  = '[0911 n11]';
NOTE_N11 = 'NOISY RHS b = A*x_true + noise, noise_level=1e-11 (all earlier FEM2 eras were noise-free); no-symmetrization arm, adaptive shift rescue ACTIVE; d=2n; 3 runs; EXCLUSIVE 64 threads (srm1526); timings citable within-era';

% [0914 g1] (PLOTTED, Toeplitz): the [0902 ns1] configuration rerun from commit
% 3363980, which (a) passes the inner absolute guard to the Toeplitz driver, so the
% stagnation-probe cycles no longer run full inner solves and the per-method inner
% counts are the productive ones (CholQR 71 -> 13, Blendenpik 78 -> 25 locally),
% and (b) writes the sketched Karlson-Walden backward-error sidecar. Ran 09-16 on
% srm1526 via campus-bigmem, 64 CPUs non-exclusive, both Gram arms
% (trsm_left/gemm_left merged at plot time). [0902 ns1] stays plotted for the
% side-by-side; drop it once the paper's Section 5.2.2 numbers have moved.
ERA_G1  = '[0914 g1]';
% RETIRED 2026-09-21: the 0917 probe eras [0917 f16] (floor lowered to 1e-16)
% and [0917 f0] (floor off) were generated but NEVER RAN -- their cells were
% cancelled before execution and no results tree exists for either. Both
% questions have since been answered by eras that did run: [0921 f16] is the
% floor-1e-16 configuration with the outer tolerance matched across benchmarks,
% and [0919 kw1] runs with the floor off. Their declarations, notes and campaign
% rows are removed rather than left to look like plottable history.
% [0919 kw1] (BACKWARD-ERROR TERMINATION, both benchmarks): the 0914_g1 / 0914_be1
% configuration with the engine's outer success test replaced by Epperly's
% step-two criterion: a run ends once the sketched Karlson-Walden backward error
% of the iterate is <= sqrt(n)*u*||A||_F (be_tol_mult=1, checked once per round,
% oracle time excluded from every solve time), and the inner absolute floor is
% OFF (inner_abs_tol=0 / ir_inner_tol=0). Rows whose stop_reason is not 'be' did
% not reach the target and are drawn censored by the plotters. The FEM cells run
% the unpreconditioned row themselves (mask 247) instead of overlaying 0915_up1,
% whose CSVs predate the new be_kw / t_be_us columns. Decided with Oleg
% 2026-09-18 (report randnla/reports/2026-09-18-emn24-solver-stopping-rules.md).
ERA_KW1 = '[0919 kw1 Epperly-BE]';
NOTE_KW1_T = 'BACKWARD-ERROR TERMINATION (sketched Karlson-Walden <= sqrt(n) u ||A||_F, be_tol_mult=1; inner floor OFF; 0914_g1 otherwise); no-symmetrization arm, adaptive shift rescue ACTIVE; d=2n; 3 runs; 64 threads (campus-bigmem, non-exclusive); rows not marked be are censored; timings citable within-era';
NOTE_KW1_F = 'BACKWARD-ERROR TERMINATION (sketched Karlson-Walden <= sqrt(n) u ||A||_F, be_tol_mult=1; inner floor OFF; 0914_be1 otherwise); NOISY RHS noise_level=1e-11; no-symmetrization arm, adaptive shift rescue ACTIVE; d=2n; 3 runs; 64 threads (bigmem SPR, non-exclusive); rows not marked be are censored; timings citable within-era';
% [0921 f16] (OLD RULE, MATCHED FLOOR, both benchmarks; 2026-09-21): the direct
% companion to [0919 kw1]. Same binary (commit fed2f9a, the kw1 install root), same
% cells, same node; the ONLY differences are the three stopping knobs. The
% backward-error oracle is OFF (be_tol_mult=0, which is the pre-BE code path
% bit-for-bit), the inner absolute floor is 1e-16 (effectively zero) instead of
% eps^0.85, and the outer tolerance is 1e-12 in BOTH benchmarks instead of
% Toeplitz 1e-12 / FEM 10*eps. Purpose (Max, 2026-09-20): [0919 kw1] collapses the
% three preconditioned Q-less methods to 1-2 iterations on FEM and was judged to
% change the story too aggressively; this era shows what the old rule says when
% both benchmarks are configured identically. NOTE the outer tolerance is INERT in
% both benchmarks -- it is tested against the true LS relative residual, which
% floors at ~1e-10 against the 1e-11 noise, so no row has ever met it; the floor
% and the stagnation exit are what actually govern the counts here.
ERA_F16B = '[0921 f16 old-rule]';
NOTE_F16B_T = 'OLD STOPPING RULE, MATCHED (inner absolute floor 1e-16, outer tol 1e-12, be_tol_mult=0 = oracle off); companion to [0919 kw1], same binary/cells/node; no-symmetrization arm, adaptive shift rescue ACTIVE; d=2n; 3 runs; 64 threads (campus-bigmem SPR, non-exclusive); timings citable within-era';
NOTE_F16B_F = 'OLD STOPPING RULE, MATCHED (inner absolute floor 1e-16, outer tol 1e-12, be_tol_mult=0 = oracle off); companion to [0919 kw1], same binary/cells/node; NOISY RHS noise_level=1e-11; unpreconditioned row runs IN-era (mask 247); no-symmetrization arm, adaptive shift rescue ACTIVE; d=2n; 3 runs; 64 threads (long-bigmem SPR, non-exclusive); timings citable within-era';
ERA_P5 = '[1002 p5 in-loop-BE]';
NOTE_P5_T = 'kw1 rule, sketched KW target also polled every 5 inner iterations';
NOTE_P5_F = 'f16 rule, no mu augmentation';

NOTE_G1 = 'INNER GUARD FIX (inner_abs_tol passed to the Toeplitz driver; probe cycles no longer inflate inner counts); no-symmetrization arm, adaptive shift rescue ACTIVE; d=2n; 3 runs; 64 threads (campus-bigmem srm1526, non-exclusive); Karlson-Walden BE sidecar; timings citable within-era';

% ONE ERA PER BENCHMARK (2026-08-31, Max). The b1 rows are commented out, not
% deleted: b1 is the accuracy-grade replication and its data is retained under
% results/, so re-enabling the comparison is uncommenting one line.
% 2026-09-21 (Max): STOPPING-RULE COMPARISON. Only the two eras under comparison
% are active, so this run opens two Toeplitz windows and two FEM windows rather
% than eight. The other rows are commented out, not deleted -- uncomment to
% restore the previous side-by-side.
TOEP_CAMPAIGNS = {
    'toeplitz_ls_0919_kw1_pcg_ne', ERA_KW1,  NOTE_KW1_T,  {'small','fixedm','middle','large'}, false
    'toeplitz_ls_0921_f16_pcg_ne', ERA_F16B, NOTE_F16B_T, {'small','fixedm','middle','large'}, false
    'toeplitz_ls_1002_p5_pcg_ne',  ERA_P5,   NOTE_P5_T,   {'small','fixedm','middle','large'}, false
%   'toeplitz_ls_0914_g1_pcg_ne',  ERA_G1,  NOTE_G1,  {'small','fixedm','middle','large'},  false
%   'toeplitz_ls_0902_ns1_pcg_ne', ERA_NS1, NOTE_NS1, {'small','fixedm','middle','large'},  false
%   'toeplitz_ls_0831_a1_pcg_ne',  ERA_A3, NOTE_A3, {'small','fixedm','middle','large'},  false
%   'toeplitz_ls_0829_a1_pcg_ne',  ERA_A2, NOTE_A2, {'small','fixedm','middle','large'},  false
%   'toeplitz_ls_0829_b1_pcg_ne',  ERA_B1, NOTE_B1, {'small','fixedm','middle','large'},  false
};

% [0914 be1] (QUEUED 2026-09-14 on long-bigmem, FEM2): the 0911_n11 configuration
% rerun from commit 3363980 so every cell also writes the sketched Karlson-Walden
% backward-error sidecar. Becomes the plotted FEM2 era once pull_0914_be1.sh has
% run (the loop below skips an era whose results dir is absent).
ERA_BE1  = '[0914 be1]';
NOTE_BE1 = 'NOISY RHS b = A*x_true + noise, noise_level=1e-11; no-symmetrization arm, adaptive shift rescue ACTIVE; d=2n; 3 runs; 64 threads (long-bigmem SPR, non-exclusive); Karlson-Walden BE sidecar; timings citable within-era';
% OVERLAY ERA [0915 up1] (2026-09-15, Oleg's Figure 6 question): the
% unpreconditioned reference row (mask bit 128 = the refinement engine on the raw
% operator, no factor), run as its OWN era so the queued 0914_be1 cells and their
% pinned install were not touched. merge_overlay_rows appends its rows to each
% base cell at plot time; the figure note names both builds. Same partition,
% threads and knobs as 0914_be1; its wall-time bar is comparable only to the
% extent the two eras landed on the same node SKU (check the PROVENANCE lines).
OVERLAY_UP1 = 'irlsq_reg_0915_up1';

% FEM2 IR-LSQ campaigns: {data subdir, era tag, note, combos to look for, overlay}
%   overlay = '' or the results subdir of a second era whose rows are appended
%   per cell at plot time (see OVERLAY_UP1 above and plotting/merge_overlay_rows.m).
% 2026-09-21 (Max): STOPPING-RULE COMPARISON -- see the note on TOEP_CAMPAIGNS.
% Both active eras ran mask 247, so the unpreconditioned row is in-era and neither
% takes the 0915_up1 overlay.
FEM2_CAMPAIGNS = {
    'irlsq_reg_0919_kw1',      ERA_KW1,  NOTE_KW1_F,  {'dd'}, ''
    'irlsq_reg_0921_f16',      ERA_F16B, NOTE_F16B_F, {'dd'}, ''
    'irlsq_reg_1002_p5',       ERA_P5,   NOTE_P5_F,   {'dd'}, ''
%   'irlsq_reg_0914_be1',      ERA_BE1, NOTE_BE1, {'dd'}, OVERLAY_UP1
%   'irlsq_reg_0911_n11',      ERA_N11, NOTE_N11, {'dd'}, OVERLAY_UP1
%   'irlsq_reg_0902_ns1',      ERA_NS1, NOTE_NS1, {'dd'}, ''
%   'irlsq_reg_0831_a1',       ERA_A3, NOTE_A3, {'dd'}, ''
%   'irlsq_reg_0829_a1',       ERA_A2, NOTE_A2, {'dd'}, ''
%   'irlsq_reg_0829_b1',       ERA_B1, NOTE_B1, {'dd'}, ''
};

% 2026-08-05 (Max): the benchmarks now record num_runs repetitions per method
% (single-run timings at 4-7 solver iterations are at the noise floor).
% TIMING_AGG picks how multi-run CSVs collapse: 'best' (min-total run, default)
% or 'mean'. Single-run CSVs (every campaign before 08-05) are unaffected, and
% the aggregation is stamped into each figure title by aggregate_runs.
TIMING_AGG = 'best';

%% ========================================================================
%  Toeplitz least-squares benchmark (Oleg's 2nd experiment)
%
%  Prolate-Toeplitz regularized LS  min ||T x - b||^2 + lambda ||x||^2, augmented
%  A = [T; sqrt(lambda) I], solved by restarted PCG-NE with Q-less QR right
%  preconditioners (published Blendenpik rows solve by LSQR; unplotted).
%  5-panel plot_toeplitz_results layout (trimmed 2026-08-24 per Oleg: recovery/
%  forward-error bar, canonical-rate and shift-retries panels removed): wall-time
%  (build + x0 + inner CG + overhead), data error, peak-vs-analytical memory,
%  preconditioner orthogonality loss, inner CG iterations. Plotted methods: CQRRT (both
%  Gram arms) / CholQR / CholQR2 / sCholQR3 / Blendenpik + refinement /
%  unpreconditioned; sCholQR3_basic and the other Blendenpik rows stay in the
%  CSVs unplotted.
% =========================================================================
for cc = 1:size(TOEP_CAMPAIGNS, 1)
    [sub, era, note, order, want_sweep] = deal(TOEP_CAMPAIGNS{cc, :});
    toep_dir = fullfile(script_dir, 'results', sub);
    if ~exist(toep_dir, 'dir')
        fprintf('(skipped Toeplitz %s %s -- no data dir %s)\n', era, sub, toep_dir);
        continue;
    end
    % Refuse a mixed-build era before anything renders (2026-08-29): per-row
    % provenance only makes mixing DETECTABLE; this makes it FATAL. The build
    % goes into the note so every figure names the binary it came from.
    commit = assert_era_commit(toep_dir, era);
    note = sprintf('%s; build %s', note, commit(1:min(7, numel(commit))));
    if isempty(order)
        % Discover size subfolders, sorted by name (old dirs are zero-padded
        % n%05d, so name order == size order).
        d = dir(toep_dir);
        size_dirs = sort({d([d.isdir] & ~ismember({d.name}, {'.','..'})).name});
    else
        % Explicit order: the new tags (small/middle/large) do NOT sort into size
        % order alphabetically, and a mis-ordered x-axis reads as a real trend.
        size_dirs = order(cellfun(@(s) exist(fullfile(toep_dir, s), 'dir') == 7, order));
    end
    if isempty(size_dirs)
        fprintf('(skipped Toeplitz %s -- no size subfolders)\n', era); continue;
    end

    fig = figure('Name', sprintf('%s Toeplitz LS -- Q-less QR right preconditioners  (%s)', era, note), ...
                 'Position', [80 + 25*cc, 80, 1550, 900]);   % 1550: 5 tiles since 2026-09-21 (was 1250 for 4)
    tg = uitabgroup(fig);
    if want_sweep
        % GATED (2026-08-27 audit): plot_toeplitz_sweep predates the 08-24 roster
        % trim (no gemmL/refine rows, stale display names, off-scheme colors) and
        % would produce a wrong paper figure. Update it before re-enabling.
        error(['run_all: want_sweep is set but plot_toeplitz_sweep is stale ' ...
               '(pre-08-24 roster); update it before re-enabling the sweep tab.']);
    end
    any_size = false;
    for s = 1:numel(size_dirs)
        size_dir = fullfile(toep_dir, size_dirs{s});
        % Two-arm cells merge at plot time (merge_gram_arms); flat cells read as before.
        merged = merge_gram_arms(size_dir, '*_toeplitz_ls_results.csv');
        if ~isempty(merged)
            size_dir = fullfile(size_dir, 'merged');
            res_csv  = merged;
        else
            f = dir(fullfile(size_dir, '*_toeplitz_ls_results.csv'));
            f = f([f.bytes] > 0);                % empty CSV = job died before writing
            if isempty(f)
                fprintf('(skipped Toeplitz %s %s -- no non-empty CSV)\n', era, size_dirs{s});
                continue;
            end
            [~, ord] = sort({f.name});           % timestamped names -> newest last
            res_csv = f(ord(end)).name;
        end
        % Label the tab with the FULL m x n, not the folder name. The old folders are
        % named after n only (n%05d), so a tab reading "n16000" sits above a plot
        % titled "32000 x 16000" -- an easy misread of the tab as the matrix shape
        % (hit 2026-07-27). Read the shape out of the CSV so the two always agree.
        shape = toeplitz_shape_label(fullfile(size_dir, res_csv), size_dirs{s});
        mt = uitab(tg, 'Title', sprintf('%s %s', era, shape));
        plot_toeplitz_results(size_dir, res_csv, ...
            sprintf('prolate Toeplitz, PCG-NE: %s', shape), mt, TIMING_AGG);   % era note NOT in the super-title (2026-09-04, Max)
        any_size = true;
    end
    if any_size
        tab_groups(end+1, :) = {tg, sprintf('toeplitz_%s', era_slug(era))}; %#ok<SAGROW>
    else
        close(fig);
    end
end

%% ========================================================================
%  FEM_Problem_2 App-1 IR-LSQ (regularized) + Blendenpik
%
%  min ||A x - b||^2 solved by IterRefineLSQ (restarted rounds, up to ir_n_steps,
%  with outer_tol and LS-floor early exits; inner CG on the preconditioned
%  normal equations) with Q-less QR right preconditioners, plus
%  the Blendenpik competitor (SASO sketch + Householder QR + LSQR). One 5-panel
%  results tab per cell (2026-08-24 trim, per Oleg: canonical-rate and
%  shift-retries panels removed; sCholQR3_basic + non-refined Blendenpik rows
%  unplotted). Runtime-breakdown tabs intentionally suppressed (empty
%  breakdown CSV name).
%
%  Cells are DISCOVERED from the directory rather than built from a
%  sizes x kappas x combos product: the 1e10 campaign includes `native_ill_dd`,
%  which has no kappa label and no size ladder (it is Oleg's unmodified generator
%  at exactly one size), so the old triple loop silently skipped it.
% =========================================================================
for cc = 1:size(FEM2_CAMPAIGNS, 1)
    [sub, era, note, combos, overlay_sub] = deal(FEM2_CAMPAIGNS{cc, :});
    camp_dir = fullfile(script_dir, 'results', sub);
    if ~exist(camp_dir, 'dir')
        fprintf('(skipped FEM2 %s %s -- no data dir)\n', era, sub); continue;
    end
    % Same mixed-build refusal as the Toeplitz loop (2026-08-29).
    commit = assert_era_commit(camp_dir, era);
    note = sprintf('%s; build %s', note, commit(1:min(7, numel(commit))));
    % Overlay era (2026-09-15): checked for a single build on its own, named in
    % the note, merged per cell below. Absent = plot without the overlay rows.
    overlay_dir = '';
    if ~isempty(overlay_sub)
        overlay_dir = fullfile(script_dir, 'results', overlay_sub);
        if exist(overlay_dir, 'dir')
            ov_commit = assert_era_commit(overlay_dir, sprintf('%s overlay %s', era, overlay_sub));
            note = sprintf('%s; unpreconditioned row from %s build %s', note, overlay_sub, ov_commit(1:min(7, numel(ov_commit))));
        else
            fprintf('(FEM2 %s: overlay %s not pulled yet -- plotting WITHOUT the unpreconditioned row)\n', era, overlay_sub);
            overlay_dir = '';
        end
    end
    d = dir(camp_dir);
    cells = sort({d([d.isdir] & ~ismember({d.name}, {'.','..'})).name});
    for c = 1:numel(combos)
        combo = combos{c};
        % Cells for this precision combo, i.e. dirs ending _<combo>.
        mine = cells(endsWith(cells, ['_' combo]));
        if isempty(mine), continue; end
        fig = figure('Name', sprintf('%s FEM2 App 1 -- IR-LSQ (reg) -- %s  (%s)', era, combo, note), ...
                     'Position', [60 + 25*c + 40*cc, 60, 1450, 900]);   % 1450: 5 tiles since 2026-09-21 (was 1150 for 4)
        tg = uitabgroup(fig);
        any_cell = false;
        for k = 1:numel(mine)
            cell_dir = fullfile(camp_dir, mine{k});
            % Two-arm cells merge at plot time (merge_gram_arms); flat cells read as before.
            merged = merge_gram_arms(cell_dir, '*_irlsq_reg_results.csv');
            if ~isempty(merged)
                cell_dir = fullfile(cell_dir, 'merged');
                res_csv  = merged;
            else
                f = dir(fullfile(cell_dir, '*_irlsq_reg_results.csv'));
                f = f([f.bytes] > 0);
                if isempty(f), fprintf('(missing cell: %s/%s)\n', sub, mine{k}); continue; end
                [~, ord] = sort({f.name});
                % NEWEST, consistently. This used to take f(1) (the OLDEST) for FEM2
                % while the Toeplitz path took the newest -- so a re-run silently did
                % not show up on one of the two plots.
                res_csv = f(ord(end)).name;
            end
            if ~isempty(overlay_dir)
                ov_cell = fullfile(overlay_dir, mine{k});
                if exist(ov_cell, 'dir')
                    [cell_dir, res_csv] = merge_overlay_rows(cell_dir, res_csv, ov_cell, {'unpreconditioned'});
                else
                    fprintf('(FEM2 %s %s: no overlay cell in %s -- unpreconditioned row missing here)\n', era, mine{k}, overlay_sub);
                end
            end
            label = strrep(mine{k}, '_', '\_');
            mt = uitab(tg, 'Title', sprintf('%s %s', era, mine{k}));
            % Pass the breakdown CSV so the wall-time panel can split the solve
            % into inner CG vs outer refinement (2026-08-27); bd_tab stays []
            % which SUPPRESSES the separate breakdown figure.
            bd_csv = strrep(res_csv, '_results.csv', '_breakdown.csv');
            plot_irlsq_results(cell_dir, res_csv, bd_csv, ...
                sprintf('%s   [%s]', label, combo), mt, [], TIMING_AGG);   % era note NOT in the super-title (2026-09-04, Max)
            any_cell = true;
        end
        if any_cell
            tab_groups(end+1, :) = {tg, sprintf('fem2_irlsq_%s_%s', era_slug(era), combo)}; %#ok<SAGROW>
        else
            close(fig);
        end
    end
end

% (The [DIAG 07-29] inner-CG diagnostic section was removed 2026-08-05 (Max):
%  the stagnation-detection fix it motivated has landed and its question is
%  answered in the 07-29/07-30 session logs. Its CSVs were deleted locally in
%  the 2026-08-07 cleanup; ISAAC benchmark-out still has them if ever needed.)

%% ========================================================================
%  Export -- one vector PDF per tab. The era tag is already in the tab title, so
%  it survives into the slug; the prefix carries it too, so OLD and NEW PDFs can
%  never collide or be confused once separated from the figure window.
% =========================================================================
export_dir = fullfile(script_dir, 'figures');
if ~exist(export_dir, 'dir'), mkdir(export_dir); end
% Wipe ALL prior PDFs so figures/ holds only the current run's tabs.
old = dir(fullfile(export_dir, '*.pdf'));
for i = 1:numel(old), delete(fullfile(export_dir, old(i).name)); end
expfig = @(h, name) exportgraphics(h, fullfile(export_dir, [name, '.pdf']), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');

for g = 1:size(tab_groups, 1)
    tg     = tab_groups{g, 1};
    prefix = tab_groups{g, 2};
    for k = 1:numel(tg.Children)
        tab  = tg.Children(k);
        % Strip the leading [ERA ...] tag before slugging: the prefix already
        % carries the era, and keeping both gave names like
        % toeplitz_OLD_old_07_11_07_15_sweep.pdf.
        title_bare = regexprep(tab.Title, '^\s*\[[^\]]*\]\s*', '');
        slug = lower(regexprep(title_bare, '[^a-zA-Z0-9]+', '_'));
        slug = regexprep(slug, '^_+|_+$', '');
        fname = [prefix, '_', slug];
        expfig(tab, fname);
        fprintf('  Exported: %s.pdf\n', fname);
    end
end
fprintf('All figures exported to %s\n', export_dir);

%% ========================================================================
%  Paper panels (2026-09-07, per Oleg): every tile of every tab is ALSO exported
%  on its own at the footprint of the synthetic-experiment panels (150 x 116 pt,
%  placed at 0.24\textwidth in the paper), so the FEM2 / Toeplitz figures sit on
%  ONE line at the same size as Figure 3. Suffix order = tile order in the
%  plotters (time, iters, orth, err, bwd); the list must be extended whenever a
%  tile is added, or the new panel exports under a generic p<N> name. It assumes
%  the optional SHOW_MEMORY / SHOW_FORWARD tiles are off, which they are.
%  Panel titles are stripped (the caption
%  names the panels); the whole-tab PDFs above remain the provenance record
%  (super-title with era, cell, aggregation). Font size and tick angle are the
%  knobs to turn if the six x labels collide at this size.
% =========================================================================
PAPER_W_PT = 110; PAPER_H_PT = 170; PAPER_FONT_PT = 10; PAPER_XTICK_ANGLE = 90;   % 2026-09-27: the body figure is one row of four panels (Oleg), each printed at 0.235\textwidth = 87 pt in the SIAM class, so 10 pt text in a 110-pt panel prints at about 8 pt (Max: larger figure text); vertical tick labels keep the six method names apart. Was 150 x 210 pt, 7 pt, 60 degrees (2026-09-07).
PAPER_SUFFIX = {'time', 'iters', 'orth', 'err', 'bwd'};   % 'bwd' = KW backward error (2026-09-21)
paper_dir = fullfile(export_dir, 'paper');
if ~exist(paper_dir, 'dir'), mkdir(paper_dir); end
old = dir(fullfile(paper_dir, '*.pdf'));
for i = 1:numel(old), delete(fullfile(paper_dir, old(i).name)); end
for g = 1:size(tab_groups, 1)
    tg     = tab_groups{g, 1};
    prefix = tab_groups{g, 2};
    for k = 1:numel(tg.Children)
        tab = tg.Children(k);
        title_bare = regexprep(tab.Title, '^\s*\[[^\]]*\]\s*', '');
        slug = lower(regexprep(title_bare, '[^a-zA-Z0-9]+', '_'));
        slug = regexprep(slug, '^_+|_+$', '');
        tl = findobj(tab, 'Type', 'tiledlayout');
        if isempty(tl), continue; end
        axs = findobj(tl(1), 'Type', 'axes', '-depth', 1);
        % findobj lists children newest-first; order the tiles left to right.
        [~, ord] = sort(arrayfun(@(a) a.Layout.Tile, axs));
        axs = axs(ord);
        for p = 1:numel(axs)
            if p <= numel(PAPER_SUFFIX), sfx = PAPER_SUFFIX{p}; else, sfx = sprintf('p%d', p); end
            export_paper_panel(axs(p), fullfile(paper_dir, sprintf('%s_%s_%s.pdf', prefix, slug, sfx)), ...
                PAPER_W_PT, PAPER_H_PT, PAPER_FONT_PT, PAPER_XTICK_ANGLE);
        end
        fprintf('  Paper panels: %s_%s_{%s}.pdf\n', prefix, slug, strjoin(PAPER_SUFFIX, ','));
    end
end
fprintf('Paper panels exported to %s\n', paper_dir);

% ---- helpers -------------------------------------------------------------
function export_paper_panel(ax, out_pdf, w_pt, h_pt, font_pt, xtick_angle)
% Copy ONE tile (plus its legend, if any) into a standalone w x h point figure
% and export it as a vector PDF. The tab figure itself is left untouched.
    fig = figure('Visible', 'off', 'Color', 'white', 'Units', 'points', ...
                 'Position', [50 50 w_pt h_pt], ...
                 'PaperUnits', 'points', 'PaperSize', [w_pt h_pt], ...
                 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 w_pt h_pt]);
    cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
    % No legend in the paper panels (2026-09-27, Max): at 10 pt a three-entry
    % legend is as wide as the 70-pt plot box and covers the y axis. The
    % figure caption names the stacked segments bottom to top instead. The
    % interactive tab figures keep their legends.
    h = copyobj(ax, fig);
    ax2 = h(1);
    % Fixed plot box (2026-09-27, Max): every panel gets the same axes rectangle
    % inside the same w x h page, so the exported PDFs share one size and the
    % panels in the paper share one height and one baseline. OuterPosition
    % [0 0 1 1] let the plot box shrink with the length of the tick labels.
    % Margins: left = y tick labels + y label, bottom = vertical x tick labels.
    lm = 36; bm = 60; rm = 4; tm = 6;
    set(ax2, 'Units', 'points', 'Position', [lm, bm, w_pt - lm - rm, h_pt - bm - tm]);
    ax2.Title.String = '';                  % the figure caption names the panel
    % Bar-top iteration counts REMOVED from the paper panels (2026-09-27, Max):
    % the paper's tables carry the counts. FAIL and N/A tags are kept.
    tx = findobj(ax2, 'Type', 'text');
    for t = reshape(tx, 1, [])
        str = t.String; if iscell(str), str = strjoin(str, ' '); end
        if ~isempty(regexp(strtrim(char(str)), '^\d+$', 'once')), delete(t); end
    end
    ax2.FontSize = font_pt;                 % ticks + labels
    set([ax2.XLabel, ax2.YLabel], 'FontSize', font_pt - 1);   % one point below the ticks so the y label clears the plot-box height (Max, 2026-09-27)
    ax2.XAxis.TickLabelRotation = xtick_angle;
    % The wall-time panel had 1.6x headroom for its legend; without the legend
    % it needs none, so the linear axes are re-tightened to the bars. The
    % plotters set ylim with legend headroom (1.35x / 1.15x); undo it here.
    if strcmp(ax2.YScale, 'linear') && isempty(findobj(ax2, 'Type', 'text'))
        % Wall-time panel: fit the axis to the tallest stacked bar plus 8%.
        % YLimMode 'auto' ended the axis at a round tick below the bar top.
        bars = findobj(ax2, 'Type', 'bar');
        if ~isempty(bars)
            tops = arrayfun(@(b) max(b.YEndPoints(:)), bars);
            ax2.YLim = [0, 1.08*max(tops)];
        end
    elseif strcmp(ax2.YScale, 'linear') && ~isempty(findobj(ax2, 'Type', 'text'))
        % Bar-top labels (iteration counts) need headroom too, or the tallest
        % label sits on the frame at paper size (2026-09-18).
        yl = ax2.YLim; ax2.YLim = [yl(1), yl(1) + 1.15*(yl(2) - yl(1))];
    end
    set(findobj(ax2, 'Type', 'text'), 'FontSize', font_pt);   % bar-top counts, FAIL tags
    % print, not exportgraphics (2026-09-27): exportgraphics crops each panel to
    % its own content, so panels with longer labels came out at different sizes
    % and LaTeX scaled them differently. print honours PaperPosition: one
    % w x h page per panel, uncropped.
    print(fig, out_pdf, '-dpdf', '-painters');
end

function s = era_slug(era)
% '[OLD 07-11/07-15]' -> 'OLD'.  Keeps PDF prefixes short but unambiguous.
% warm/cold eras keep their qualifier so the two figure sets export to
% DISTINCT filenames instead of overwriting each other (2026-07-31).
    % Slug the WHOLE tag (2026-08-27 fix): the old leading-uppercase match
    % returned empty on digit-led tags like '[0812 d2]', so those exports
    % silently degraded to the generic 'ERA' prefix, and a future digit-led
    % era would have COLLIDED with them.
    t = regexprep(era, '[\[\]]', '');
    s = regexprep(strtrim(t), '[^A-Za-z0-9]+', '_');
    s = regexprep(s, '^_+|_+$', '');
    if isempty(s), s = 'ERA'; end
end
