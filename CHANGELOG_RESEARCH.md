# Changelog Research

## Scope

- Repo: `Dicklesworthstone/frankensim`
- Requested task: update the changelog using `changelog-md-workmanship`.
- Scope window: project inception on 2026-07-05 through
  `main@291c7dbfdc6fd366f6ee55a6dbc39f137e05afd6` on 2026-07-10.
- Public remote state when researched: local `HEAD` resolved to
  `291c7dbfdc6fd366f6ee55a6dbc39f137e05afd6`; `origin/main` resolved to
  `291c7dbfdc6fd366f6ee55a6dbc39f137e05afd6` immediately before the
  changelog commit.
- Working-tree policy: the changelog covers committed history only. Local
  uncommitted edits and scratch files were not used as changelog evidence,
  including the live Beads export, README, manifests and locks, crate source,
  contracts, tests, unsafe-capsule registry, and untracked probe/helper files.

## Evidence Sources

- Git history:
  - `git rev-list --count 291c7db` -> 999 commits.
  - `git rev-list --count 43d52f2..291c7db` -> 666 commits since the prior
    changelog baseline.
  - `git rev-list --count fb08842..291c7db` -> 377 commits since the previous
    changelog pass endpoint.
  - `git rev-list --count 319cb64..291c7db` -> 366 commits since the previous
    changelog endpoint.
  - `git rev-list --count d5873bf..291c7db` -> 318 commits since the 2026-07-09
    flagship expansion checkpoint.
  - `git rev-list --count e08e302..291c7db` -> 303 commits since the previous
    changelog endpoint.
  - `git rev-list --count 6725739..291c7db` -> 261 commits since that public
    mainline snapshot.
  - `git rev-list --count e993e76..291c7db` -> 142 commits in this update window.
  - `git rev-list --count fb34ef6..291c7db` -> 7 commits in the final follow-up.
  - `git log --reverse --no-merges --pretty=format:'%h %ad %s' --date=short`.
  - `git log --all --no-merges --pretty=format:'%H %h %ad %s' --date=short`.
  - `git diff --shortstat 43d52f2..291c7db` -> 774 files changed,
    142,929 insertions, 2,639 deletions.
  - `git diff --shortstat fb08842..291c7db` -> 475 files changed,
    63,442 insertions, 3,068 deletions.
  - `git diff --shortstat 319cb64..291c7db` -> 450 files changed,
    58,132 insertions, 3,092 deletions.
  - `git diff --shortstat d5873bf..291c7db` -> 388 files changed,
    43,989 insertions, 3,064 deletions.
  - `git diff --shortstat e08e302..291c7db` -> 384 files changed,
    42,270 insertions, 3,024 deletions.
  - `git diff --shortstat 6725739..291c7db` -> 347 files changed,
    36,242 insertions, 2,972 deletions.
  - `git diff --shortstat e993e76..291c7db` -> 250 files changed,
    24,976 insertions, 2,356 deletions.
  - `git diff --shortstat fb34ef6..291c7db` -> 19 files changed,
    1,043 insertions, 104 deletions.
- Version metadata:
  - `git for-each-ref refs/tags ...` -> no tags.
  - `gh release list --limit 100` -> no GitHub Releases.
- Tracker:
  - `git show HEAD:.beads/issues.jsonl | jq ...` for closed workstreams and
    close reasons at committed `HEAD`.
- Project docs:
  - `README.md`.
  - `AGENTS.md`.
  - crate-level `CONTRACT.md` files referenced through README/workstream scope.
  - Workspace count checks on 2026-07-10:
    - `git ls-tree -d --name-only HEAD:crates | rg '^fs-' | wc -l`
      -> 126 tracked `fs-*` crate directories.
    - `git ls-tree -r --name-only HEAD | rg '^crates/fs-[^/]+/CONTRACT\.md$' | wc -l`
      -> 126 tracked contracts.
    - `git ls-tree -r --name-only HEAD | rg '^crates/fs-[^/]+/tests/.*\.rs$' | wc -l`
      -> 251 tracked crate-level test files.

## Version Spine

| Node | Kind | Date | Notes |
|------|------|------|-------|
| `8e4c0a5` | Inception commit | 2026-07-05 | Initial FrankenSim plan. |
| `43d52f2` | Prior changelog baseline | 2026-07-07 | Latest state covered by the first changelog reconstruction. |
| `941a67e` | Earlier public checkpoint | 2026-07-08 | Value-of-information query planning; 621 commits. |
| `fb08842` | Prior changelog checkpoint | 2026-07-08 | Proof-robust, schedule, and flutter e2e campaigns; 622 commits. |
| `319cb64` | Prior public mainline snapshot | 2026-07-08 | Previous changelog endpoint; 633 commits. |
| `d5873bf` | Prior public mainline snapshot | 2026-07-09 | Ornithoid flagship contract, CutFEM-octree topopt lane, sparse-GP/MLMC/MOO decision loops, browser campaign tiers, proof-hardening fixes; 681 commits. |
| `9dc7417` | Public mainline snapshot | 2026-07-09 | Exact e2e certificate payloads, topopt marquee evidence hardening, and in-progress flagship replay-suite scaffold; 684 commits. |
| `1fe4ef5` | Public mainline snapshot | 2026-07-09 | Flagship replay-suite lint cleanup plus Beads tracker state for the fs-mesh v2 follow-up; 686 commits. |
| `20e9825` | Public mainline snapshot | 2026-07-09 | Frozen-golden flagship e2e suite; 687 commits. |
| `438128d` | Public mainline snapshot | 2026-07-09 | Mesh hull-facet encroachment protection and faster topology-optimization marquee witnesses; 688 commits. |
| `70e1d24` | Public mainline snapshot | 2026-07-09 | Boundary-splitting semantics and topopt move-application cleanup; 692 commits. |
| `d93ca59` | Public mainline snapshot | 2026-07-09 | Gated SME2 exploratory capsule and unsafe registration; 695 commits. |
| `e08e302` | Public mainline snapshot | 2026-07-09 | SME2 exploratory battery; 696 commits. |
| `6725739` | Public mainline snapshot | 2026-07-09 | Live browser flagship pipelines, mesh v2/v3 closure, sparse roofline/NUMA lanes, rand/FFT perf work, fail-closed IO/risk/probe hardening; 738 commits. |
| `e993e76` | Public mainline snapshot | 2026-07-10 | Deterministic `.powi` policy, declared RNG run identity, caller-owned cancellation gates, solver snapshot envelopes, GEMM perf evidence lanes, legal PLY face-first imports, standard empirical CVaR, and adjoint fail-closed regressions; 857 commits. |
| `f93c9de` | Public mainline snapshot | 2026-07-10 | Evidence-package schema and BLAKE3 roots, failure compounding, golden/claim/closure policy gates, machine-adaptive execution, fail-closed numerical fixes, and the first production GEMM autotune/roofline integration; 988 commits. |
| `fb34ef6` | Public mainline snapshot | 2026-07-10 | Trusted historical roofline baselines, an explicit rejection of global x86 FMA contraction, and a FrankenNumpy re-pin past the case-colliding corpus; 992 commits. |
| `2568c82` | Public mainline snapshot | 2026-07-10 | NEON interleaved-complex transpose, register-accumulator FEEC contractions, and a standalone zero-dependency clean-machine bootstrap; 998 commits. |
| `291c7db` | Public mainline snapshot | 2026-07-10 | Const-bound FEEC accumulation loops preserve the frozen golden while exposing a friendlier x86 SLP shape; 999 commits. |

No tags or GitHub Releases existed when researched.

## Coverage Ledger

| Chunk | Range | Status | Major themes |
|-------|-------|--------|--------------|
| 01 | 2026-07-05 | distilled | Plan, README, Beads, license, agent workflow. |
| 02 | 2026-07-06 foundations | distilled | Workspace scaffold, constellation, DSR, substrate, execution, bedrock numerics, evidence, ledger, geometry core. |
| 03 | 2026-07-06 representations | distilled | SDF, mesh, F-rep, meshing, transforms, Chebyshev, planning, optimization/image scaffolds. |
| 04 | 2026-07-07 core expansion | distilled | Constraint, GA, scenario/regime/material, query/time, FEEC, TileLang, topology, NURBS, operator DSL. |
| 05 | 2026-07-07 addendum flywheel | distilled | Three-color schema, falsifiers, recompute, physics VCS, semantic diff, speculation verifier/proposer/economics, governance, Phase 0 spine, assume-guarantee contracts. |
| 06 | latest solver slice | distilled | `fs-solver` mixed-precision Krylov refinement. |
| 07 | 43d52f2..15eb757 | distilled | NURBS/SDF conversion, CutFEM, matrix-free adjoint scaffold, physics-VCS bisect, `fs-solid` elasticity. |
| 08 | 15eb757..3f4e343 | distilled | Solid stability/structural elements, FEEC cohomology, ledger-DAG transposition, sheaf merge, topopt/topols, UQ, BO, anytime planner, whole-loop flywheel gate. |
| 09 | 3f4e343..6df2c03 | distilled | FMM/BEM, LBM, IGA, surrogate ladders, truss, neural reps, rendering, domain decomposition, Navier-Stokes, marquee runner. |
| 10 | 6df2c03..7049ca3 | distilled | LBM extensions, seismic frame, mesh coloring/recovery, contact, topology persistence, conformal hardening, lattice, e-racing, time slabs, Cheb variants, Payne-Hanek. |
| 11 | 7049ca3..fb08842 | distilled | Vortex-thruster QD campaign, DRO oracle, value-of-information queries, and three certified e2e capstones. |
| 12 | fb08842..319cb64 | distilled | Neural-shape and grammar campaigns, SensorForge, vessel flagship, metamaterial/truss/AnytimeBO/FlowCert e2e crates, inverse-trig AD, `fs-ad` bridge/Revolve/IFT integrations, vertex-patch Schwarz p-MG smoothing. |
| 13 | 319cb64..e08e302 | distilled | Self-knowledge and flywheel gates, Stokes block preconditioners, NSGA-III/MOEA/D/steering, sparse GP, adaptive MLMC/chance constraints, constrained polish engines, proposal-only generation, browser campaign tiers, differentiable rendering, ornithoid and CutFEM-octree flagships, exact e2e certificate payloads, topopt evidence hardening, frozen-golden flagship replay suite, mesh refinement protection, SME2 exploratory capsule, tracker state, proof-hygiene fixes. |
| 14 | e08e302..6725739 | distilled | Live `fs-wasm` flagship browser pipelines, `fs-mesh` v2/v3 closure, hex-dominant meshing scaffold, compact sparse roofline/NUMA work, deterministic rand/FFT perf lanes, PLY integer validation, non-finite risk/probe hardening, and packaging panic-surface cleanup. |
| 15 | 6725739..e993e76 | distilled | Deterministic integer-power policy, run-identity RNG replay, race/session cancellation ownership, versioned solver snapshots, f32/mixed/transposed/strided GEMM evidence lanes, standard empirical CVaR boundary weighting, legal PLY face-first import, and adjoint certificate fail-closed regressions. |
| 16 | e993e76..81c11df | distilled | Golden-coupling and claim-state policy, certificate/package hardening, canonical replay identity, deterministic math repairs, topology fixes, and failure-compounding families. |
| 17 | 81c11df..cdb62ee | distilled | Lock-driven bootstrap, canonical replay adoption, package schema v3, machine-adaptive execution, honest FFT/batched-GEMM measurements, and fail-closed roofline admission. |
| 18 | cdb62ee..fb34ef6 | distilled | Roofline receipt and historical-baseline admission, typed GEMM autotuning, transactional parallel GEMM, package schema v4/BLAKE3 ownership, closure-evidence enforcement, FFT/FEEC measured outcomes, production registry wiring, and the constellation case-collision re-pin. |
| 19 | fb34ef6..291c7db | distilled | Measured NEON complex transpose, FEEC register accumulators and const-bound loops with unchanged goldens, standalone pre-Cargo constellation bootstrap, and tracker handoffs preserving unresolved performance and admission boundaries. |

## Representative Commit Clusters

- Foundation:
  - `8e4c0a5` initial plan.
  - `a7e4d54` Rust workspace scaffold.
  - `4b31ce8` DSR-first policy.
- Substrate/runtime:
  - `47c7719` substrate probes.
  - `741979d` SIMD tiers.
  - `59b85cb` arenas and pools.
  - `39bd2f8` execution context and tile pools.
- Bedrock:
  - `089cc72` Philox streams.
  - `7456302` FFT core.
  - `44a883b` deterministic CSR core.
  - `e7dc872` interval arithmetic.
  - `a433365` GEMM/factorization core.
- Morph:
  - `397e325` region/chart crate.
  - `ebf97fb` SDF crate.
  - `124ed81` mesh crate.
  - `a295b6b` F-rep DAGs.
  - `4a17114` shape parameterizations.
  - `3948a82` topology certificates.
  - `f8a3cf7` NURBS algebra and trims.
- Flux:
  - `a968527` FEEC exterior calculus.
  - `89c1f82` operator DSL compiler.
  - `7aec6a5` solver battery.
  - `43d52f2` mixed-precision Krylov refinement.
- Geometry/solids/topology:
  - `279d611` certified NURBS-to-SDF conversion.
  - `c9e1a4b` SDF-to-NURBS refit.
  - `f781085` CutFEM on SDFs.
  - `15eb757` solid elasticity core.
  - `e396160` solid stability and Koiter scaffolding.
  - `bcc71d4` structural elements.
  - `4b4c24e` level-set topology optimization.
- Differentiation/planning/flywheel:
  - `41cadcc` ledger-DAG transposition.
  - `3fab970` gradient certificates.
  - `9ee7227` fidelity-ladder planner.
  - `1ade44a` anytime/refusal semantics.
  - `9ab427e` whole-loop flywheel harness.
  - `3f9a714` Phase 2 leverage gate.
  - `941a67e` value-of-information query ranking.
- Physics/numerics:
  - `7d95552` FEEC cohomology.
  - `57d775e` BDDC domain decomposition.
  - `8f031ad` pressure-robust Navier-Stokes.
  - `8c5b0e5` FMM and BEM.
  - `ded5b78` LBM extensions.
  - `7049ca3` Payne-Hanek trig reduction.
- Optimization/e2e:
  - `a76cbea` value-of-information crate.
  - `d947001` Wasserstein DRO oracle.
  - `c0f1a5c` SOS proof-carrying optimization.
  - `3c56418` quality-diversity archives.
  - `6df2c03` marquee study runner.
  - `8d27622` seismic frame flagship.
  - `dc1bf7f` certified vortex-thruster QD campaign.
  - `fb08842` three certified e2e campaigns.
  - `5cbdd90` neural-shape and grammar e2e campaigns.
  - `94404c4` SensorForge OED campaign.
  - `b95e00f` laminar-pour vessel flagship.
  - `8028bbc` metamaterial stiffness-density frontier.
  - `4c573d6` truss critical-load-path campaign.
  - `76e9c89` AnytimeBO.
  - `3eb480c` FlowCert CFD credibility map.
- Latest numerics/solver:
  - `922c835` deterministic inverse trig and AD `Real` operations.
  - `7575cdd` `fs-ad` bridge, Revolve, spill, and IFT integrations.
  - `319cb64` vertex-patch Schwarz p-MG smoothing.
- 2026-07-09 mainline expansion:
  - `7624964` self-knowledge e2e battery.
  - `3dcbfaa` Stokes block preconditioners and PMINRES.
  - `8c20bab` NSGA-III reference directions.
  - `61e9348` edge-aware differentiable rendering.
  - `fe8a2c2` `fs-plan` budget allocator.
  - `6d074b8` seismic UQ, anytime stopping, and adaptive MLMC.
  - `e2b3dff` sparse Gaussian-process lane.
  - `73c1c6b` proposal-only `fs-gen`.
  - `51adc98` steering and chance constraints.
  - `2f29772` constrained polish engines and Problem-IR runner.
  - `2199248` certified `fs-wasm` campaign tier.
  - `ffcb4a9` `fs-ornith` flagship crate.
  - `8c06c93` CutFEM-octree topology marquee lane.
  - `d5873bf` `fs-ornith` contract and flagship bead close.
  - `fbba704` flutter and neuroshape exact certificate payloads.
  - `1092e94` CutFEM-octree marquee evidence hardening.
  - `9dc7417` `fs-flagship-e2e` staged replay-suite scaffold.
  - `1337058` `fs-flagship-e2e` replay-suite lint posture.
  - `1fe4ef5` Beads tracker state for `frankensim-uee3`.
  - `20e9825` frozen-golden flagship e2e suite.
  - `438128d` hull-facet encroachment protection and topopt marquee witness tuning.
  - `70e1d24` boundary splitting and topopt move-application cleanup.
  - `a8e98b2` hull encroachment regression bound.
  - `d93ca59` gated SME2 exploratory capsule.
  - `e08e302` SME2 exploratory battery.
- Post-SME2 follow-up:
  - `e8c496e` CutFEM-octree topology marquee zero-remeshing contract.
  - `29fb1ea` SME2 streaming-mode GEMM prototype contract.
  - `f759916` deterministic ziggurat normal fast path.
  - `cb6dbde` Philox statistical stream battery.
  - `8254164` bitwise bulk Philox fills.
  - `648a891` c2r inverse FFT and 2D/3D pencil transforms.
  - `9a6d303` conforming interior facet recovery.
  - `6510b02` non-convex recovery facet triangulation.
  - `fa3d93d` 10-million-point `fs-mesh` perf lane closeout.
  - `151a643` hex-dominant meshing module.
  - `0cd9de2` compact/sharded SpMV and sparse roofline lane.
  - `dccc078` NUMA first-touch and real shared STREAM slices.
  - `b4f849f` invalid topology-optimization density-shape rejection.
  - `220d28c` fractional PLY face-list rejection.
  - `3842369` non-convex recovery and boundary-layer gates.
  - `4b77149` measured `fs-mesh` v3 boundary-layer decision.
  - `ff9ef32` CVaR/objective non-finite rejection.
  - `38f871f` non-finite probe budget rejection.
  - `fdffa73` provenance in evidence package roots.
  - `91478e2` package string-write panic cleanup.
  - `462f0f6` panic-arm removal from `fs-io` conformance checks.
  - `6725739` live browser flagship pipelines in `fs-wasm`.
- 2026-07-10 hardening and perf evidence:
  - `8c6db27` deterministic integer powers and policy checks.
  - `921c486` RNG streams bound to declared run identity.
  - `d31573b` race panic-total and cancellation wiring requirements.
  - `98ea5db` caller-owned race cancellation gates.
  - `361bb36` session pressure pauses bound to owned gates.
  - `b2cb2c2` versioned self-authenticating solver snapshot envelope.
  - `af0339e` packed f32 and mixed-precision GEMM paths.
  - `dbbffa8` transposed/strided op-form GEMM.
  - `5b8aeb7` GEMM and batched perf evidence lanes.
  - `7f6420f` standard empirical CVaR boundary weighting.
  - `4fbdefc` legal face-before-vertex PLY element order.
  - `e993e76` adjoint certificate fail-closed regressions.
- 2026-07-10 evidence-bound performance and package integrity:
  - `2ecb0cd` bind citable roofline runs to exact, re-derivable timed receipts.
  - `dc859e2` introduce bounded GEMM plans and the session autotune scaffold.
  - `666967d` make tuner evidence unit-safe, canonical, bounded, and explicitly
    non-statistical.
  - `f91b775` add private-staging transactional cancellation to parallel GEMM.
  - `3096ee5` bind tune rows and decisions to the complete execution identity.
  - `da8f438` wire the session-backed GEMM tuner into the shipped roofline
    registry, retaining its pre-admission persistence no-claim.
  - `f93c9de` derive target coverage from the production registry while keeping
    target attainment a separate evidence question.
  - `2d6d3f9` add governed fingerprint-specific historical axis baselines so
    sustained contention cannot self-authorize citable roofline evidence.
  - `005f0e2` document why globally enabling x86 FMA would break G5 bit
    determinism and must be replaced by narrow per-kernel capsules.
  - `fb34ef6` re-pin FrankenNumpy after upstream renamed case-colliding corpus
    files, restoring representable clean checkouts on default macOS.
  - `9711705` extract the dependency-free safe-Rust BLAKE3 owner.
  - `10d1e2c` migrate package/checker roots to domain-separated 32-byte BLAKE3.
  - `504d2a8` require closed bug beads to cite regression evidence or an
    explicit disposition.
  - `d4b9c04` retain the correct-but-slower six-step FFT behind its frontier
    gate instead of promoting it.
  - `acedb1b` remove the FEEC x86 libm-FMA call path while leaving the packed
    vectorization target open.
- 2026-07-10 post-checkpoint kernel and bootstrap follow-up:
  - `f7e9787` add the NEON interleaved-complex transpose capsule, scalar twin,
    and measured six-step relative lane while keeping the feature default-off.
  - `ea03f01` record the partial 27d3 result and release the open workstream.
  - `2168810` replace FEEC destination-memory accumulation with local register
    arrays while preserving the frozen golden and leaving x86 proof pending.
  - `a3edf0e` add the isolated zero-dependency clean-machine constellation
    bootstrap and its hermetic replay harness.
  - `2568c82` close 1t8i with explicit evidence and a residual literal
    blank-machine public-remote no-claim.
  - `291c7db` replace FEEC iterator/zip accumulator loops with explicit
    const-bound indexing while preserving the golden and pending ISA proofs.
- 2026-07-10 golden governance and counterexample compounding:
  - `10bdac0` add the golden-coupling registry and justified-bump protocol.
  - `a2a2b10` record four-quadrant deterministic-power evidence.
  - `2ca177e` couple README hashes, crate names, and sentinels to live source.
  - `3712602` refuse overflowing dimensional exponent arithmetic.
  - `d0c1e09` make certified evidence an opaque validated trust boundary.
  - `049538f` separate SOS algebra verification from rigorous value bounds.
  - `914ff55` introduce typed, length-prefixed canonical replay identities.
  - `f882987` land strict schema-v2 packages and solver-free checking.
  - `b426193` repair phantom topology tunnels from inconsistent connectivity.
  - `81c11df` close the failure-compounding workflow with four-quadrant
    manifest evidence.
- 2026-07-10 reproducible bootstrap, canonical replay, and honest performance:
  - `6fb8513` add lock-driven constellation bootstrap and repair a layer
    inversion exposed by policy checks.
  - `f96156d` derive feature-gated DSR lanes, visible skips, standalone WASM
    checks, and repository inventories.
  - `3a263c2` remove deterministic-trig drift and migrate capability identity.
  - `b24b0dc` finish deterministic BEM routing and restore full flagship metric
    streams.
  - `e00858c` add package composition receipts, falsifiers, and anchors in
    schema v3.
  - `cd185a6` align in-memory, serialized, and standalone package verification
    while bounding hostile inputs.
  - `fb3bfc8` add measured L3 topology, the affinity capsule, and CCD/first-touch
    A/B evidence.
  - `ceccb00` add advisory TilePool CCD pinning with bit-equivalence gates.
  - `ad80d5c` make autotuner replay and persistence refuse foreign or ambiguous
    records.
  - `535aa83` improve FFT raw throughput while recording that the 40% roofline
    target remains unmet.
  - `6d55c1b` add deterministic batched f32 and mixed GEMM without a performance
    claim.
  - `cdb62ee` refuse vacuous roofline passes when contention crushes both axes
    and the measured kernel together.
- Addendum:
  - `e43e3b1` three-color schema.
  - `39fd1a5` falsifier pairing.
  - `ea102b5` recompute store.
  - `772d975` physics VCS.
  - `2f2fe56` certified-speculation verifier.
  - `92636d0` speculation economics and ledger v3 telemetry.
  - `544eaee` Phase 0 spine gate.
  - `b33496e` assume-guarantee contracts.

## Open Questions

- The repo has no formal versioning scheme yet. Future changelog updates should
  split entries by tag/release once the first release is created.
- Links to Beads currently target the committed `.beads/issues.jsonl` file as a
  durable tracker source; a future issue viewer could provide more precise
  per-record URLs.
- `frankensim-yqug` was closed at `f93c9de` and reopened by `fb34ef6` with four
  active blockers: citable receipts must bind target/decision identity,
  production scheduling must traverse `TilePool`/`Cx`, durable keys must bind
  exact build/effective tier, and tune rows must not persist before roofline
  admission. Changelog prose must preserve that transition.
- `frankensim-7uq9` and `frankensim-t7x3` remain in progress even though the
  shared BLAKE3 owner and schema-v4 package migration landed materially; future
  updates should reconcile tracker state before calling those workstreams
  closed.
- FEEC's narrow baseline-x86 libm `fma` hole is closed, but the measured packed
  vectorization target remains red and `frankensim-cwjn` stays open. Do not
  infer target attainment from capsule registration.
- `frankensim-dfh3` closed at `2d6d3f9` with governed fingerprint-specific
  historical baselines. First-run measurements remain candidate evidence and
  cannot authorize themselves.
- `frankensim-1t8i` and `frankensim-7n2n` are closed at `2568c82`: the
  standalone bootstrap can run before Cargo resolves sibling paths, and the
  current FrankenNumpy pin checks out cleanly on both filesystem classes.
  Follow-up review still found documentation drift, no post-checkout
  cleanliness check for a newly cloned sibling, and silent fallback for a bare
  `--root` or `--from`; do not broaden the closure beyond the recorded replay
  and real-lock evidence.
- The committed NEON transpose facade uses unchecked `2 * n1 * n1` size
  arithmetic before entering its unsafe pointer loop. The measured performance
  result is real, but a safe-facade overflow regression remains required.
- The closure-evidence lint validates the vocabulary of close reasons; it does
  not resolve cited artifacts and prove that those external artifacts exist or
  pass.
- Uncommitted local work should only be added after it is committed and proven.
