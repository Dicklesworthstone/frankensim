# Changelog Research

## Scope

- Repo: `Dicklesworthstone/frankensim`
- Requested task: update the changelog using `changelog-md-workmanship`.
- Scope window: project inception on 2026-07-05 through
  `main@319cb64f052d76e15882ee53ace41092881c7fa8` on 2026-07-08.
- Public remote state when researched: `HEAD` and `origin/main` both resolved to
  `319cb64f052d76e15882ee53ace41092881c7fa8`.
- Working-tree policy: ignore uncommitted in-progress work outside
  `CHANGELOG.md` and `CHANGELOG_RESEARCH.md` (including current tracker,
  manifest, README, `fs-wasm`, e2e, recovery, and cache changes). This
  changelog covers committed history only.

## Evidence Sources

- Git history:
  - `git rev-list --count HEAD` -> 633 commits.
  - `git rev-list --count origin/main` -> 633 commits.
  - `git rev-list --count 43d52f2..HEAD` -> 300 commits since the prior
    changelog baseline.
  - `git rev-list --count fb08842..HEAD` -> 11 commits since the previous
    changelog pass endpoint.
  - `git log --reverse --no-merges --pretty=format:'%h %ad %s' --date=short`.
  - `git log --all --no-merges --pretty=format:'%H %h %ad %s' --date=short`.
  - `git diff --stat --compact-summary 43d52f2..HEAD` -> 492 files changed,
    85,660 insertions, 410 deletions.
  - `git diff --stat --compact-summary fb08842..HEAD` -> 53 files changed,
    5,432 insertions, 98 deletions.
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

## Version Spine

| Node | Kind | Date | Notes |
|------|------|------|-------|
| `8e4c0a5` | Inception commit | 2026-07-05 | Initial FrankenSim plan. |
| `43d52f2` | Prior changelog baseline | 2026-07-07 | Latest state covered by the first changelog reconstruction. |
| `941a67e` | Earlier public checkpoint | 2026-07-08 | Value-of-information query planning; 621 commits. |
| `fb08842` | Prior changelog checkpoint | 2026-07-08 | Proof-robust, schedule, and flutter e2e campaigns; 622 commits. |
| `319cb64` | Public mainline snapshot | 2026-07-08 | Latest `origin/main` researched; 633 commits. |

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
- Uncommitted local work outside the changelog files should only be added after
  it is committed and proven.
