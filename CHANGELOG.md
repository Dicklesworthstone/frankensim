# Changelog

This is a synthesized, agent-facing changelog for FrankenSim.

Scope window: project inception on 2026-07-05 through
[`main@438128d`](https://github.com/Dicklesworthstone/frankensim/commit/438128d988082f2183f406d75883b766fd6b7324)
on 2026-07-09.

This document was rebuilt from git history, tag/release metadata, the checked-in
Beads tracker, and the current README/contract surface. It is organized by
landed capabilities rather than raw diff order.

## Version Timeline

There are no git tags and no GitHub Releases as of
[`main@438128d`](https://github.com/Dicklesworthstone/frankensim/commit/438128d988082f2183f406d75883b766fd6b7324).

| Version | Kind | Date | Summary |
|---------|------|------|---------|
| [`main@438128d`](https://github.com/Dicklesworthstone/frankensim/commit/438128d988082f2183f406d75883b766fd6b7324) | Public mainline snapshot | 2026-07-09 | 688 commits, adding hull-facet encroachment protection to mesh refinement and tightening the topology-optimization marquee evidence fixtures. |
| [`main@1fe4ef5`](https://github.com/Dicklesworthstone/frankensim/commit/1fe4ef5c21737a92aa58dd02f901052de5b0a88b) | Public mainline snapshot | 2026-07-09 | 686 commits, adding the flagship replay-suite scaffold, topopt marquee evidence hardening, exact certificate payload refinements, and tracker state for the fs-mesh v2 follow-up. |
| [`main@9dc7417`](https://github.com/Dicklesworthstone/frankensim/commit/9dc7417090e62e7ab2756719034e2ab0d3b32875) | Mainline checkpoint | 2026-07-09 | 684 commits, registering the staged flagship replay suite and its contract. |
| [`main@d5873bf`](https://github.com/Dicklesworthstone/frankensim/commit/d5873bfd82a3c2dbe359c11aef4947a5def8cdba) | Public mainline snapshot | 2026-07-09 | 681 commits, adding the ornithoid flagship contract, CutFEM-octree topology lane, sparse-GP and MLMC decision loops, browser campaign tiers, and proof-hardening fixes. |
| [`main@319cb64`](https://github.com/Dicklesworthstone/frankensim/commit/319cb64f052d76e15882ee53ace41092881c7fa8) | Public mainline snapshot | 2026-07-08 | 633 commits, extending the campaign suite, inverse-trig/AD surface, and p-MG smoother evidence. |
| [`main@fb08842`](https://github.com/Dicklesworthstone/frankensim/commit/fb088428ae810ee3ddd893712588c1e64f6cc0c4) | Prior changelog checkpoint | 2026-07-08 | 622 commits, including the proof-robust, schedule, and flutter end-to-end campaign capstones. |
| [`origin/main@941a67e`](https://github.com/Dicklesworthstone/frankensim/commit/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f) | Earlier public checkpoint | 2026-07-08 | 621 commits, covering the large 2026-07-08 implementation wave through value-of-information query planning. |
| [`main@43d52f2`](https://github.com/Dicklesworthstone/frankensim/commit/43d52f2ed960ede2067eec353359063c86d9dbb3) | Prior changelog baseline | 2026-07-07 | Working Rust workspace with 333 commits, no formal release tag, and the Phase 0 spine closed. |
| [`8e4c0a5`](https://github.com/Dicklesworthstone/frankensim/commit/8e4c0a5c4f18aa7d0bd0add47e407b835c7a3b86) | Inception commit | 2026-07-05 | Initial FrankenSim plan. |

## 1. Project Foundation And Policy Spine

FrankenSim started as a plan-first repository and quickly became a real Rust
workspace with project-specific agent policy, Beads issue tracking, licensing,
README/illustration assets, DSR-first validation policy, and workspace policy
checks.

### Delivered capability

- Initial architecture plan and README explaining the deterministic,
  evidence-oriented simulation substrate.
- Beads tracker and `bv` workflow documentation for project workstreams.
- MIT license with the Anthropic/OpenAI rider.
- Rust workspace scaffold with flat `fs-*` crates.
- `xtask` policy gates for layer direction, dependency policy, contracts,
  unsafe capsules, and constellation locks.
- DSR documented as the preferred validation lane over GitHub Actions.

### Closed workstreams

- [`frankensim-epic-foundations-huq.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - workspace scaffold, layer enforcement, lint/format rails.
- [`frankensim-epic-foundations-huq.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - pinned Franken constellation integration.
- [`frankensim-epic-foundations-huq.4`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - CI/Gauntlet gate policy.

### Representative commits

- [`8e4c0a5`](https://github.com/Dicklesworthstone/frankensim/commit/8e4c0a5c4f18aa7d0bd0add47e407b835c7a3b86) - initial FrankenSim plan.
- [`48e45c7`](https://github.com/Dicklesworthstone/frankensim/commit/48e45c7e2965db98c29e85b7bf06cf2f021eb22e) - initialize Beads project tracker.
- [`7554599`](https://github.com/Dicklesworthstone/frankensim/commit/755459986a056368685250c932407da019af7a8b) - add the MIT license with AI rider.
- [`a7e4d54`](https://github.com/Dicklesworthstone/frankensim/commit/a7e4d54d904d131850ef238f4b8e32b366c60736) - add the initial Rust workspace scaffold.
- [`4b31ce8`](https://github.com/Dicklesworthstone/frankensim/commit/4b31ce842848f40c9548a21a6513bbc050f38aa6) - document DSR as the primary CI runner.

## 2. Substrate, Execution, Determinism, And Performance

The first implementation wave established the lower layers that later crates
depend on: hardware capability discovery, explicit execution contexts,
deterministic reductions, aligned allocation, tile identity, SIMD capsules,
autotuning, and roofline measurement.

### Delivered capability

- `fs-substrate` capability probes, machine fingerprints, tiled fields, Morton
  identities, bandwidth probes, and CCD-aware topology helpers.
- `fs-simd` scalar and capsule-backed SIMD facades.
- `fs-alloc` scoped arenas, aligned allocation, hugepage policy, arena pools,
  and site tracking.
- `fs-exec` execution contexts, cancellation gates, tile pools, deterministic
  reductions, speculative races, resumable solver state, and autotuner hooks.
- `fs-roofline` machine-axis probing, kernel registry, ledgered measurements,
  and staleness checks.

### Closed workstreams

- [`frankensim-epic-substrate-wf9.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - substrate capability probes and machine fingerprints.
- [`frankensim-epic-substrate-wf9.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - SIMD facades and capsules.
- [`frankensim-epic-substrate-wf9.4`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - scoped arenas and allocation pools.
- [`frankensim-epic-substrate-wf9.7`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - two-lane executor and `Cx` contract.
- [`frankensim-epic-perf-fz2.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - roofline harness.

### Representative commits

- [`47c7719`](https://github.com/Dicklesworthstone/frankensim/commit/47c771994b3294d1dda1949640e4c46fc920203a) - implement `fs-substrate` capability probes.
- [`741979d`](https://github.com/Dicklesworthstone/frankensim/commit/741979d7cc8f11b9ca1e0f556d28fc7095476ca0) - implement `fs-simd` scalar and capsule tiers.
- [`59b85cb`](https://github.com/Dicklesworthstone/frankensim/commit/59b85cb47f69f2fffbf2258ff119288de05c58a8) - implement `fs-alloc` scope arenas and pools.
- [`39bd2f8`](https://github.com/Dicklesworthstone/frankensim/commit/39bd2f8e686cb97623d322fd0bd8ab96e2e8db80) - add `fs-exec` tile context and pool modules.
- [`05e3c3d`](https://github.com/Dicklesworthstone/frankensim/commit/05e3c3d6554ca094f7a8d526834f0142107ff4fb) - add deterministic reduction support.

## 3. Bedrock Numerics And Certified Arithmetic

The BEDROCK layer filled in deterministic math, random streams, dense and sparse
linear algebra, FFTs, interval/affine/Taylor models, Chebyshev tools,
automatic differentiation, e-process inference, and mixed-precision ladders.

### Delivered capability

- `fs-math` deterministic elementary functions, strict complex arithmetic, EFT,
  and double-double helpers.
- `fs-rand` logical Philox streams, deterministic distributions, QMC/Sobol
  components, and replayable sampling.
- `fs-la` GEMM/factorization kernels, eigensolver scaffolding, mixed-precision
  refinement ladders, randomized NLA, and batched small-dense kernels.
- `fs-sparse` deterministic COO/CSR/BSR/SELL assembly and sparse operations.
- `fs-fft` FFT/DCT transform kernels.
- `fs-ivl` intervals, affine arithmetic, exact predicates, Taylor models, and
  Newton/Krawczyk certification helpers.
- `fs-cheb` Chebyshev function objects and Orr-Sommerfeld spectral probes.
- `fs-ad` forward duals and adjoint infrastructure.
- `fs-eproc` anytime-valid e-process machinery.

### Closed workstreams

- [`frankensim-epic-bedrock-6ys.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - logical Philox streams.
- [`frankensim-epic-bedrock-6ys.3`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - GEMM and dense kernels.
- [`frankensim-epic-bedrock-6ys.9`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - sparse formats and deterministic assembly.
- [`frankensim-epic-bedrock-6ys.12`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - interval and affine arithmetic.
- [`frankensim-epic-bedrock-6ys.13`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - Taylor models and certified roots.
- [`frankensim-epic-ascent-7tv.6`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - e-process inference core.

### Representative commits

- [`089cc72`](https://github.com/Dicklesworthstone/frankensim/commit/089cc72cd6d00db69e272fcfae13af42b1629f19) - implement logical Philox streams.
- [`7456302`](https://github.com/Dicklesworthstone/frankensim/commit/7456302d9ae1b6aaa96be2f7f04a8aeb8a12209f) - add FFT correctness-first transform core.
- [`44a883b`](https://github.com/Dicklesworthstone/frankensim/commit/44a883b17167dc0f848b13710208f9a606c35d67) - add deterministic CSR core.
- [`e7dc872`](https://github.com/Dicklesworthstone/frankensim/commit/e7dc872f211fa66bc3d41969ebe21c4c886dc97b) - add interval arithmetic core.
- [`a433365`](https://github.com/Dicklesworthstone/frankensim/commit/a433365ab46ce18726c20cb80ff25f527507e432) - add GEMM and factorization core.
- [`a8f88f8`](https://github.com/Dicklesworthstone/frankensim/commit/a8f88f82c3aaebc0864926e9b6125d7d27f1c14c) - add evidence and model-form certification core.

## 4. Geometry, Representation, Meshing, And Shape Control

The MORPH layer became a working representation stack rather than a plan:
regions/charts, SDFs, meshes, F-reps, NURBS, voxel fields, exact topology,
representation routing, geometric constraints, IO quarantine, and
parameterization levers all landed as crate-level surfaces with contracts.

### Delivered capability

- `fs-geom` chart/region contracts, conversion records, representation router,
  semantic diff, sheaf watertightness certificates, and sheaf repair.
- `fs-rep-sdf`, `fs-rep-mesh`, `fs-rep-frep`, `fs-rep-nurbs`, and
  `fs-rep-voxel` representation crates.
- `fs-mesh` Delaunay/refinement/remeshing scaffolding.
- `fs-xform` FFD, RBF, level-set, density, and composition parameterizations
  with Jacobian actions and foldover checks.
- `fs-topo` validity and topology certificates.
- `fs-query` geometry interrogation layer and `fs-geocon` geometric constraint
  primitives.
- `fs-io` import/export with quarantine receipts.

### Closed workstreams

- [`frankensim-epic-morph-wqd.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - region/chart core.
- [`frankensim-epic-morph-wqd.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - SDF representation crate.
- [`frankensim-epic-morph-wqd.4`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - mesh representation crate.
- [`frankensim-epic-morph-wqd.5`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - NURBS spline algebra and trims.
- [`frankensim-epic-morph-wqd.19`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - design parameterizations.
- [`frankensim-epic-morph-wqd.23`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - topology certificates.
- [`frankensim-epic-morph-wqd.25`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - import/export quarantine.

### Representative commits

- [`397e325`](https://github.com/Dicklesworthstone/frankensim/commit/397e3255e28599bbd1ebc99dee9fd12923b26837) - add `fs-geom` region and chart crate.
- [`ebf97fb`](https://github.com/Dicklesworthstone/frankensim/commit/ebf97fb63dc4b3c921b6b08e2d7db9f8ad0c1834) - add SDF representation crate.
- [`124ed81`](https://github.com/Dicklesworthstone/frankensim/commit/124ed81ceb7a8b89283d504177b9c5739e997208) - add mesh representation crate.
- [`a295b6b`](https://github.com/Dicklesworthstone/frankensim/commit/a295b6b3f6a380e79dd4cf4ec4f9fb6e4a47b4df) - add F-rep CSG DAGs and evaluators.
- [`4a17114`](https://github.com/Dicklesworthstone/frankensim/commit/4a17114c15002faf5d400706ef8971187da26614) - add design parameterizations.
- [`3948a82`](https://github.com/Dicklesworthstone/frankensim/commit/3948a82fb0793064c3428052335c59a46d70f8fe) - add validity and topology certificates.
- [`f8a3cf7`](https://github.com/Dicklesworthstone/frankensim/commit/f8a3cf7d60c56afdc121f53f85a9cbbb9c10d9eb) - implement NURBS algebra, trims, and honest Boolean policy.
- [`3432e6e`](https://github.com/Dicklesworthstone/frankensim/commit/3432e6e6d5f6a847cab52d00eaa47afff8eba2b0) - implement quarantined import/export.

## 5. Physics, Operators, Solvers, And Time Integration

The FLUX layer gained typed scenario/regime/material foundations, exterior
calculus infrastructure, operator generation, time integration, solver
state/reporting, p-multigrid, and mixed-precision Krylov refinement.

### Delivered capability

- `fs-scenario` typed boundary/load-case algebra.
- `fs-regime` Buckingham-Pi scaling, nondimensionalization, and admission
  predicates.
- `fs-material` constitutive kernels with consistent tangents.
- `fs-feec` exterior calculus complexes and high-order tensor/simplex families.
- `fs-opdsl` typed operator IR that generates primal/JVP/VJP/adjoint/DWR hooks.
- `fs-time` structure-preserving and resumable integrator batteries.
- `fs-solver` Krylov solvers, p-multigrid preconditioning, structured stall
  diagnosis, and mixed-precision refinement.
- `fs-iface` coupling-graph static checker.

### Closed workstreams

- [`frankensim-epic-flux-tfz.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - scenario algebra.
- [`frankensim-epic-flux-tfz.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - material kernels.
- [`frankensim-epic-flux-tfz.3`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - regime/admission gates.
- [`frankensim-epic-flux-tfz.4`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - operator DSL/compiler.
- [`frankensim-epic-flux-tfz.5`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - FEEC complexes.
- [`frankensim-epic-flux-tfz.6`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - high-order FEEC.
- [`frankensim-epic-flux-tfz.12`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - time integrators.

### Representative commits

- [`88c1b46`](https://github.com/Dicklesworthstone/frankensim/commit/88c1b46c85f0895f25cc84f8f85a60e1e28c95b5) - implement regime/admission gates.
- [`052d38b`](https://github.com/Dicklesworthstone/frankensim/commit/052d38bc5524e47f8ec931fcb74f9b852bd937e8) - implement material kernels.
- [`a968527`](https://github.com/Dicklesworthstone/frankensim/commit/a96852737b948c61233e6f8c6541233d3c9d10e2) - complete FEEC exterior calculus contract.
- [`89c1f82`](https://github.com/Dicklesworthstone/frankensim/commit/89c1f827e72a6e22a2a19b3d9a9b3737c5bf224c) - add typed operator IR compiler.
- [`de11ad1`](https://github.com/Dicklesworthstone/frankensim/commit/de11ad1caf20be601d165916685903ba45d0dfe8) - lock time-integrator evidence.
- [`7aec6a5`](https://github.com/Dicklesworthstone/frankensim/commit/7aec6a5b1d603d8fe96040dabdb752fb0baac976) - add Krylov and p-multigrid battery.
- [`43d52f2`](https://github.com/Dicklesworthstone/frankensim/commit/43d52f2ed960ede2067eec353359063c86d9dbb3) - add mixed-precision Krylov refinement.

## 6. ASCENT Optimization And Constraint Machinery

The optimization layer moved from problem sketches to typed optimization data,
constraints, geometric constraints, derivative-free engines, Goodhart checks,
and problem serialization.

### Delivered capability

- `fs-opt` typed objective/constraint graphs, manifold variables, PDE and
  stochastic nodes, differentiability routing, canonical serialization, and
  the Goodhart guard.
- `fs-constraint` typed hard/soft/chance/robust/cert/fabrication/code
  constraints, interval certificates, unsat cores, and ranked repairs.
- `fs-dfo` CMA-ES/IGO core with deterministic restart behavior and DFO
  no-claim boundaries.
- `fs-geocon` geometric constraint primitives for min-thickness, draft,
  symmetry, envelopes, and volume/mass.

### Closed workstreams

- [`frankensim-epic-ascent-7tv.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - optimization problem IR.
- [`frankensim-epic-ascent-7tv.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - constraint calculus.
- [`frankensim-epic-ascent-7tv.4`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - CMA-ES/DFO engines.
- [`frankensim-epic-epistype-qmao.5`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - Goodhart guard.
- [`frankensim-epic-morph-wqd.21`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - geometric constraint primitives.

### Representative commits

- [`5b86049`](https://github.com/Dicklesworthstone/frankensim/commit/5b860492632414045326fa7a260701eae515cc16) - scaffold `fs-opt`.
- [`64dbee1`](https://github.com/Dicklesworthstone/frankensim/commit/64dbee1cf45d28669e36e9eb658e08dab5827f78) - add the typed constraint calculus.
- [`f536e98`](https://github.com/Dicklesworthstone/frankensim/commit/f536e98f8254a443988830b2365e6a29d865dd3a) - add the Goodhart guard.
- [`b5f1210`](https://github.com/Dicklesworthstone/frankensim/commit/b5f12106da8c28f0a97e12ba2df968790537c8ac) - harden DFO/CMA behavior.
- [`19a4193`](https://github.com/Dicklesworthstone/frankensim/commit/19a4193b5356f8c910f3b3166ae52cdc8a068584) - add geometric constraint primitives.

## 7. HELM, Ledger, IR, Sessions, And Evidence

The orchestration layer gained evidence-carrying values, ledger schema,
time-travel queries, FrankenScript IR, admission checks, sessions/governor
support, planning models, image artifacts, and a narrow PV vertical skeleton.

### Delivered capability

- `fs-evidence` `Evidence<T>`/`Certified<T>` wrappers, model cards,
  discrepancy/bracketing records, falsifier pairing, and evidence gates.
- `fs-ledger` content-addressed artifacts, operations, event streams, tune rows,
  lineage, time travel, physics-VCS primitives, tombstone memory, and telemetry.
- `fs-ir` typed FrankenScript AST, JSON/s-expression syntax, admission checks,
  declarative query IR, and structured errors.
- `fs-plan` cost/error ledgers and planning models.
- `fs-session` capability tokens, governor checks, idempotency, and dry-run
  estimates.
- `fs-vskeleton` a photovoltaic vertical skeleton tying geometry, PDE,
  objective, adjoint, optimization, and ledger concepts together.
- `fs-img` deterministic PNG/EXR-style artifact plumbing and bias-labeled
  denoising.

### Closed workstreams

- [`frankensim-epic-helm-gp3.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - ledger core.
- [`frankensim-epic-helm-gp3.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - evidence wrappers.
- [`frankensim-epic-helm-gp3.3`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - ledger time travel.
- [`frankensim-epic-helm-gp3.4`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - IR core.
- [`frankensim-epic-helm-gp3.5`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - IR admission.
- [`frankensim-epic-helm-gp3.7`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - session/governor.
- [`frankensim-epic-foundations-huq.8`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - PV vertical skeleton milestone.

### Representative commits

- [`1e20bda`](https://github.com/Dicklesworthstone/frankensim/commit/1e20bdac6b18431dae1507c874721709cefc1fee) - add Design Ledger v0.
- [`a8f88f8`](https://github.com/Dicklesworthstone/frankensim/commit/a8f88f82c3aaebc0864926e9b6125d7d27f1c14c) - add evidence and certification core.
- [`33c1ed4`](https://github.com/Dicklesworthstone/frankensim/commit/33c1ed4c6ba49af09740e4b3c1cc9ec57449a30d) - complete ledger time travel.
- [`21369ba`](https://github.com/Dicklesworthstone/frankensim/commit/21369ba96892ba136bf8630a221005107f698be9) - add FrankenScript typed AST.
- [`df6e8dc`](https://github.com/Dicklesworthstone/frankensim/commit/df6e8dc93478dcf151de77b773ef1b3fd3429763) - implement IR static admission.
- [`c8e5ca9`](https://github.com/Dicklesworthstone/frankensim/commit/c8e5ca9523f7586d219799cbefbaf2d4e1b9dbfe) - implement session governor.
- [`8f71102`](https://github.com/Dicklesworthstone/frankensim/commit/8f71102ca9f82c68d69848dc5d9d493e164c435d) - complete `fs-img`.

## 8. Addendum Flywheel, Certified Speculation, And Governance

The addendum wave added a machine-readable governance layer and a speculative
execution/evidence flywheel: three-color epistemic schema, falsifier gates,
incremental recompute, physics VCS, semantic diff, proposer zoo, speculation
economics, risk registers, contract registries, and an executable Phase 0 spine
gate.

### Delivered capability

- Three-color ledger schema for verified/validated/estimated evidence.
- Falsifier pairing, consequence-by-doubt budgets, and no-falsifier-no-ship
  checks.
- Slack-bearing recompute store and tolerance-aware invalidation.
- Physics-VCS commit/checkout primitives and semantic diff.
- Certified speculation verifier, proposer zoo, and accept/reject economics.
- `fs-spececo`, `fs-verify`, `fs-probe`, `fs-ladder`, `fs-recompute`,
  `fs-govern`, and `fs-contract` addendum crates.
- Phase 0 spine gate: laundering refusal, falsifier pairing, tombstone protocol,
  interface-type checks, Goodhart guard, and budget-pie reporting.

### Closed workstreams

- [`frankensim-epic-epistype-qmao.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - three-color schema.
- [`frankensim-epic-epistype-qmao.4`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - falsifier pairing.
- [`frankensim-epic-flywheel-lmp4.1`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - certified speculation verifier.
- [`frankensim-epic-flywheel-lmp4.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - proposer zoo.
- [`frankensim-epic-flywheel-lmp4.3`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - accept/reject economics.
- [`frankensim-epic-flywheel-lmp4.6`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - recompute store.
- [`frankensim-epic-flywheel-lmp4.10`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - semantic diff.
- [`frankensim-epic-addendum-xpck.2`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - Phase 0 spine complete.
- [`frankensim-epic-flywheel-lmp4.14`](https://github.com/Dicklesworthstone/frankensim/blob/43d52f2ed960ede2067eec353359063c86d9dbb3/.beads/issues.jsonl) - assume-guarantee contracts.

### Representative commits

- [`e43e3b1`](https://github.com/Dicklesworthstone/frankensim/commit/e43e3b14c8db88bbc5263a98cfc9cfde19fcc2e3) - add the three-color epistemic schema and ledger gate.
- [`39fd1a5`](https://github.com/Dicklesworthstone/frankensim/commit/39fd1a5144ae211073fbf8100774f1b4366304b3) - add falsifier pairing.
- [`ea102b5`](https://github.com/Dicklesworthstone/frankensim/commit/ea102b596b725fe77808ee6a6c16473a8d2e5f11) - add the slack-bearing recompute store.
- [`772d975`](https://github.com/Dicklesworthstone/frankensim/commit/772d97508f1e19bf5d551c1147bd83b65f4e470f) - add physics-VCS base verbs.
- [`2f2fe56`](https://github.com/Dicklesworthstone/frankensim/commit/2f2fe56339aa811649d1ebd0702d0f1cf43e09b3) - add certified-speculation verifier v0.
- [`6c6f633`](https://github.com/Dicklesworthstone/frankensim/commit/6c6f6338e04fd412ba29c630f4d4bfd86a56d1d3) - add proposer zoo.
- [`92636d0`](https://github.com/Dicklesworthstone/frankensim/commit/92636d0bb5071a503ed872bb8337cd90b2676350) - land integrated speculation economics and ledger v3 telemetry.
- [`544eaee`](https://github.com/Dicklesworthstone/frankensim/commit/544eaee3c0fd27ab672c9f6cf36fbcbd7a99204a) - close the Phase 0 spine milestone.
- [`b33496e`](https://github.com/Dicklesworthstone/frankensim/commit/b33496e291b928f5113bd65c651b2d5b3b97b044) - add assume-guarantee component contracts.

## 9. Certified Geometry Conversion, CutFEM, Solids, And Topology Optimization

The post-baseline wave turned geometry and structural mechanics into a more
direct design loop: exact and certified representation conversion, CutFEM on
SDFs, elasticity and contact on cut geometry, structural stability, and both
density and level-set topology optimization landed as working crate surfaces.

### Delivered capability

- `fs-rep-nurbs` gained a certified NURBS-to-SDF converter and an SDF-to-NURBS
  refit path, giving the workspace a round-trip bridge between spline CAD,
  fields, F-reps, and sheaf watertightness checks.
- `fs-cutfem` added certified cut classification, depth-controlled cut
  quadrature, ghost penalties, Nitsche embedded boundary handling, quadtree
  adaptivity, and zero-meshing FEM on SDFs.
- `fs-solid` landed linear elasticity, finite-strain hyperelasticity through
  material cards, B-bar near-incompressibility mitigation, CutFEM frontends,
  stability continuation, buckling pencils, Koiter scaffolding, Cosserat rods,
  fiber sections, force-based beams, and SDF-native contact.
- `fs-topopt` added density topology optimization with robust three-field
  filtering, stress aggregation, eigenfrequency objectives, and golden
  batteries.
- `fs-topols` added level-set topology optimization over the CutFEM stack with
  WENO advection, fast-marching redistancing, velocity extension, topological
  derivatives, and deterministic snapshots.

### Closed workstreams

- [`frankensim-epic-morph-wqd.11`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - certified NURBS-to-SDF conversion.
- [`frankensim-epic-morph-wqd.12`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - SDF-to-NURBS refit.
- [`frankensim-epic-flux-tfz.13`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - `fs-solid` elasticity core.
- [`frankensim-epic-flux-tfz.14`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - advanced structural element kernels.
- [`frankensim-epic-flux-tfz.15`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - stability, continuation, and Koiter path.
- [`frankensim-epic-flux-tfz.16`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - SDF-native IPC contact.

### Representative commits

- [`279d611`](https://github.com/Dicklesworthstone/frankensim/commit/279d611817f854930a8dc20cd7605702a0f14c2f) - add the certified NURBS-to-SDF converter.
- [`c9e1a4b`](https://github.com/Dicklesworthstone/frankensim/commit/c9e1a4b984773dfed3fcd953fd251394567867ec) - close the certified SDF-to-NURBS round trip.
- [`f781085`](https://github.com/Dicklesworthstone/frankensim/commit/f781085a538063c1d6711777030efa3ca5817d19) - ship CutFEM on SDFs.
- [`15eb757`](https://github.com/Dicklesworthstone/frankensim/commit/15eb75700d69ca73d6fec1cc72b18f3f65666a01) - add the `fs-solid` elasticity core.
- [`e396160`](https://github.com/Dicklesworthstone/frankensim/commit/e3961608169ab16ea19beb69b368857ba5636221) - add stability, continuation, and Koiter scaffolding.
- [`bcc71d4`](https://github.com/Dicklesworthstone/frankensim/commit/bcc71d484b001d7188b920308971b7800c4cf543) - add structural rods, fiber sections, and force-based beams.
- [`4b4c24e`](https://github.com/Dicklesworthstone/frankensim/commit/4b4c24ef48d93be96cb95c4b9b362bb414d6c8d7) - add the `fs-topols` level-set topology optimization crate.
- [`5015e57`](https://github.com/Dicklesworthstone/frankensim/commit/5015e5795a805e37fda6f58b3f46e736c0428557) - add stress aggregate objective plumbing.
- [`0c66ce7`](https://github.com/Dicklesworthstone/frankensim/commit/0c66ce7741e3e48f9289f4a8e975d2ad12e74084) - add eigenfrequency objective support.
- [`7616efe`](https://github.com/Dicklesworthstone/frankensim/commit/7616efedd5e76fb5394a284dd7dab4bc867462b1) - close the SDF-native contact lane.

## 10. Differentiation, Adaptivity, Planning, And Self-Knowledge

The second addendum wave made derivatives, goal-oriented accuracy, and economic
planning first-class rather than incidental: gradients now carry evidence,
planners refuse when budgets are inadequate, and value-of-information work
turns ignorance into a priced queue.

### Delivered capability

- `fs-adjoint` gained a ledger-DAG VJP registry, fail-loud missing-gradient
  behavior, differentiability mitigations, interval residual gradient
  certificates, explanation objects, DWR accept gates, and exact Hessian-vector
  products for density misfit.
- `fs-dwr` added goal-oriented error estimation, anisotropic metric synthesis,
  marking, tile-level adaptivity scaffolding, and contract-level accuracy
  boundaries.
- `fs-ir` added a greedy fidelity-ladder planner, anytime query semantics,
  colored intervals, priced tightening hints, and explicit refusal reports.
- `fs-plan` added value-of-information query planning that ranks surrogate rung
  climbs and physical validation probes by decision impact per dollar.
- `fs-flywheel-e2e` added whole-loop, phase-1, and phase-2 gates that measure
  the compound effect of speculation, recompute, sheaf merges, tombstones,
  adjoint gradients, and planner ladders.

### Closed workstreams

- [`frankensim-epic-coupling-bk0o.1`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - ledger-DAG transposition.
- [`frankensim-epic-coupling-bk0o.3`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - gradient certificates and merge gates.
- [`frankensim-epic-flywheel-lmp4.16`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - fidelity-ladder planner.
- [`frankensim-epic-flywheel-lmp4.17`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - anytime and refusal semantics.
- [`frankensim-epic-flywheel-lmp4.18`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - whole-loop flywheel proof harness.
- [`frankensim-epic-addendum-xpck.4`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - Phase 2 leverage gate.

### Representative commits

- [`41cadcc`](https://github.com/Dicklesworthstone/frankensim/commit/41cadccb11349a1ab5b48374b5e5f875906cdfd7) - add ledger-DAG transposition.
- [`3fab970`](https://github.com/Dicklesworthstone/frankensim/commit/3fab970ad002fb16223020cd8dfde362c117e004) - add gradient certificates.
- [`e907bc6`](https://github.com/Dicklesworthstone/frankensim/commit/e907bc6796b5e8b840b52bc618eae1657a969bf2) - add the DWR goal-oriented accept test.
- [`9ee7227`](https://github.com/Dicklesworthstone/frankensim/commit/9ee7227858c7f09c86bd694b600a1bb263a013b1) - add the greedy fidelity-ladder planner.
- [`1ade44a`](https://github.com/Dicklesworthstone/frankensim/commit/1ade44a99cf9a67a6be8f6fa10b914209f9622c6) - add anytime and refusal semantics.
- [`9ab427e`](https://github.com/Dicklesworthstone/frankensim/commit/9ab427e3ca95b3e0e0556547c84234ce4ae37321) - land the whole-loop flywheel harness.
- [`3f9a714`](https://github.com/Dicklesworthstone/frankensim/commit/3f9a714868c2de4a8ab62359f903a7990c272b4c) - add the Phase 2 leverage gate.
- [`941a67e`](https://github.com/Dicklesworthstone/frankensim/commit/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f) - add value-of-information queries to `fs-plan`.

## 11. Extended Physics, Numerics, And Solver Kernels

The FLUX and BEDROCK layers broadened from core abstractions into usable
simulation kernels: harmonic FEEC analysis, domain decomposition, pressure-
robust flow, FMM/BEM, LBM, vortex particles, time slabs, and large-domain
deterministic trigonometry.

### Delivered capability

- `fs-feec` added cohomology, harmonic cochains, Hodge splitting, circulation
  extraction, and high-order vector-family spaces.
- `fs-dd` added BDDC domain decomposition with sheaf-framed edge coarse spaces
  and measured conditioning tables.
- `fs-flux` added a BDM1-P0 pressure-robust Navier-Stokes scaffold with
  boundary moments, conformance batteries, Picard stepping, and discrete
  adjoint hooks.
- `fs-fmm` and `fs-bem` added black-box Chebyshev FMM, Laplace panel methods,
  Hess-Smith airfoil screening, wake roll-up, and transpose/FMM validation.
- `fs-lbm` added D2Q9 core dynamics, rheology, thermal coupling, free-surface
  VOF, contact-line bracketing, and refinement.
- `fs-vpm`, `fs-iga`, and `fs-couple` added vortex-particle, isogeometric, and
  passive port-Hamiltonian coupling surfaces.
- `fs-time` gained time slabs as ledger cells, while `fs-math` extended
  deterministic trig with Payne-Hanek reduction across all finite `f64`.

### Closed workstreams

- [`frankensim-epic-flux-tfz.7`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - FEEC cohomology and harmonic decomposition.
- [`frankensim-epic-flux-tfz.11`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - domain decomposition.
- [`frankensim-epic-flux-tfz.17`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - pressure-robust Navier-Stokes.
- [`frankensim-epic-flux-tfz.19`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - LBM extensions.
- [`frankensim-r6r5`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - Payne-Hanek trig reduction.

### Representative commits

- [`7d95552`](https://github.com/Dicklesworthstone/frankensim/commit/7d95552cf6bbb1055296991fb43f119b27916635) - add FEEC cohomology.
- [`57d775e`](https://github.com/Dicklesworthstone/frankensim/commit/57d775e47f3d38f94b7c1d017c3c74888e551ec0) - add BDDC domain decomposition.
- [`8f031ad`](https://github.com/Dicklesworthstone/frankensim/commit/8f031ad48a56a3ac5e91023cd99782138918bcd4) - complete pressure-robust Navier-Stokes.
- [`8c5b0e5`](https://github.com/Dicklesworthstone/frankensim/commit/8c5b0e555e9d3d0f1a7840e941fed43890519008) - add `fs-fmm` and `fs-bem`.
- [`eeb919a`](https://github.com/Dicklesworthstone/frankensim/commit/eeb919a98ddb7dc8acb8aff87fb6fa85138eeb3b) - implement the FMM transpose panel operator.
- [`714e9e6`](https://github.com/Dicklesworthstone/frankensim/commit/714e9e6ef39d590ddfd4fa42254fec89a136cc4b) - add the D2Q9 LBM core.
- [`ded5b78`](https://github.com/Dicklesworthstone/frankensim/commit/ded5b78909a40b81b493289bea2cd6742aa99eb1) - complete LBM rheology, thermal, free-surface, and refinement extensions.
- [`43ea2dd`](https://github.com/Dicklesworthstone/frankensim/commit/43ea2dd57b29c33e221250cec35b5d8db814a6d5) - add the vortex-particle method core.
- [`8a90167`](https://github.com/Dicklesworthstone/frankensim/commit/8a9016744727cdeb0b9788df6999e86e9bdc38ab) - add time slabs as cells.
- [`7049ca3`](https://github.com/Dicklesworthstone/frankensim/commit/7049ca38995f62664d0715580ae7b8fe987f22ae) - ship Payne-Hanek trig reduction for all finite `f64`.

## 12. Optimization, UQ, Robust Design, And Archives

The ASCENT layer became much broader than the initial DFO core. It now includes
gradient optimization, Bayesian optimization, uncertainty quantification,
surrogates, proof-carrying optimization, value-of-information, fabrication
constraints, quality-diversity archives, robust design, and multiple explicit
design-space representations.

### Delivered capability

- `fs-ascent` added L-BFGS, trust-region Newton, Riemannian routing, augmented
  Lagrangians, Wolfe search, stopping policy, and gradient Pareto tracing.
- `fs-bo` added Gaussian-process Bayesian optimization, TuRBO trust regions,
  and multi-fidelity BO.
- `fs-dfo` expanded with log-domain Sinkhorn optimal transport, deterministic
  multi-objective utilities, and a finite-support Wasserstein DRO oracle.
- `fs-uq`, `fs-surrogate`, `fs-voi`, `fs-race`, and `fs-sos` added UQ,
  certify-or-escalate ROMs, active validation, anytime-valid e-racing, and SOS
  certificates.
- `fs-archive` added MAP-Elites/CVT quality-diversity archives, while
  `fs-lattice`, `fs-truss`, and `fs-fab` added homogenized lattice design,
  ground-structure sizing, and manufacturing/code-compliance constraints.
- `fs-shapeprog`, `fs-rep-neural`, `fs-toleralloc`, `fs-assimilate`, and
  `fs-robust` added program synthesis, Lipschitz-certified neural implicits,
  tolerance allocation, as-built/data assimilation, and objective epistemics.

### Closed workstreams

- [`frankensim-epic-ascent-7tv.18`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - value-of-information and active validation.
- [`frankensim-epic-ascent-7tv.20`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - fabrication and code-compliance constraints.
- [`frankensim-epic-ascent-7tv.7`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - e-racing with anytime-valid cancellation.
- [`frankensim-epic-ascent-7tv.10`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - proof-carrying SOS optimization.
- [`frankensim-epic-ascent-7tv.14`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - lattice homogenization and graded optimization.
- [`frankensim-vcia`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - Sinkhorn OT and DRO evidence lane.

### Representative commits

- [`f1d2028`](https://github.com/Dicklesworthstone/frankensim/commit/f1d202841d6cb4c29fb60a1fdbfd8624df1c8c69) - add the gradient optimizer stack.
- [`b93bd21`](https://github.com/Dicklesworthstone/frankensim/commit/b93bd21afce726eb3fddbd3b8c20db3d428a796d) - add gradient Pareto tracing.
- [`a76cbea`](https://github.com/Dicklesworthstone/frankensim/commit/a76cbea571cb0cb9b324a4734a200e40a63a8ca2) - add the value-of-information crate.
- [`dcfc2ac`](https://github.com/Dicklesworthstone/frankensim/commit/dcfc2ac3503004c5cf227e05282d47fdbb659ac2) - add TuRBO trust-region BO.
- [`2d3a3c4`](https://github.com/Dicklesworthstone/frankensim/commit/2d3a3c41de507156d9a22eac5385f649dc217f2c) - add multi-fidelity BO.
- [`961e2d6`](https://github.com/Dicklesworthstone/frankensim/commit/961e2d60cad948a0b48ac52256dbc89ee154bf8c) - add log-domain Sinkhorn optimal transport.
- [`382e616`](https://github.com/Dicklesworthstone/frankensim/commit/382e6169cfcfb951006ba53dfe9cd6b175454ee1) - add deterministic multi-objective utilities.
- [`d947001`](https://github.com/Dicklesworthstone/frankensim/commit/d9470015316f0e587cc797ff53b4a0dbf359d7fc) - add the discrete Wasserstein DRO oracle.
- [`c0f1a5c`](https://github.com/Dicklesworthstone/frankensim/commit/c0f1a5c2ef63582e173efc87be255e3af1242e6c) - add SOS proof-carrying optimization.
- [`3c56418`](https://github.com/Dicklesworthstone/frankensim/commit/3c56418c549420e1f62886af9cd1c31a78c38a38) - add MAP-Elites/CVT quality-diversity archives.
- [`4c3ec84`](https://github.com/Dicklesworthstone/frankensim/commit/4c3ec84cb5cdf5dcbf53e301accb3fc9e908be70) - ship lattice homogenization and graded optimization.
- [`f9fad53`](https://github.com/Dicklesworthstone/frankensim/commit/f9fad5309f23a18de2d9356cecb01d21938e8aa7) - add certified ground-structure truss sizing.

## 13. Rendering, Reporting, Browser Surface, And End-To-End Campaigns

The project now has several capstone paths where many crates compose into a
visible workflow: scientific rendering, browser demos, lab notebooks, flagship
optimization studies, and cross-layer end-to-end certification harnesses.

### Delivered capability

- `fs-render` added unbiased spectral path tracing, certified chart tracing,
  mixed-scene backends, and Woodcock volume rendering over live LBM fields.
- `fs-img`, `fs-viz`, and `fs-report` hardened PNG validation, added analytic
  scientific visualization primitives, and generated automatic lab notebooks
  with semantic diffs.
- `fs-wasm` added browser-oriented numerical kernel surfaces and tiered demo
  modules.
- `fs-marquee` added the raw-SDF, zero-meshing marquee study runner.
- `fs-frame` added the seismic-minimal frame flagship with truss layout,
  fiber-hinge time history, anytime-valid fragility, and CVaR sizing.
- `fs-thrust-e2e` added a certified quality-diversity vortex-thruster campaign
  composing VPM, evidence, surrogates, archives, and reports.
- [`fb08842`](https://github.com/Dicklesworthstone/frankensim/commit/fb088428ae810ee3ddd893712588c1e64f6cc0c4)
  added three additional certified end-to-end campaigns:
  proof-robust optimization, campaign scheduling, and flutter certification.

### Closed workstreams

- [`frankensim-epic-lumen-qfx.1`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - spectral path tracing.
- [`frankensim-epic-lumen-qfx.2`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - chart render backends.
- [`frankensim-epic-lumen-qfx.3`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - volume rendering.
- [`frankensim-epic-flagships-mye.1`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - marquee raw-SDF optimization study.
- [`frankensim-epic-flagships-mye.3`](https://github.com/Dicklesworthstone/frankensim/blob/941a67e3cfd3a2bbc1fc3209bf0a044bc073188f/.beads/issues.jsonl) - seismic-minimal frame flagship.

### Representative commits

- [`b089421`](https://github.com/Dicklesworthstone/frankensim/commit/b089421c5494e231d3623973251e91a30169efb0) - add the spectral path-tracing core.
- [`f8a927e`](https://github.com/Dicklesworthstone/frankensim/commit/f8a927e601372affc84484e041f8f99d513fe88a) - add certified chart render backends.
- [`0ccee78`](https://github.com/Dicklesworthstone/frankensim/commit/0ccee78a87ec4e81ca38554142e7506aea356bf5) - ship Woodcock volume rendering over live LBM fields.
- [`8aaecbc`](https://github.com/Dicklesworthstone/frankensim/commit/8aaecbcd1881349484bdabd3f71a5f38d250ee28) - add analytically verified scientific visualization primitives.
- [`6e5ae66`](https://github.com/Dicklesworthstone/frankensim/commit/6e5ae66563ec1267711c591ef23bbf47a4568596) - add automatic lab notebooks and semantic diffs.
- [`68768ba`](https://github.com/Dicklesworthstone/frankensim/commit/68768ba496da9ab6698041429f3b9c23a5d291f6) - add the standalone browser numerical-kernel surface.
- [`6df2c03`](https://github.com/Dicklesworthstone/frankensim/commit/6df2c03eb3a2ea5c8c731faad7cabda19475aabb) - ship the raw-SDF marquee study runner.
- [`8d27622`](https://github.com/Dicklesworthstone/frankensim/commit/8d27622657d0f4ace6bdc40991f295e9b958e34c) - ship the seismic-minimal frame flagship.
- [`dc1bf7f`](https://github.com/Dicklesworthstone/frankensim/commit/dc1bf7faf7ed6feafcda0f7c3c3eb78fe6d1d0e0) - add the certified quality-diversity vortex-thruster campaign.
- [`fb08842`](https://github.com/Dicklesworthstone/frankensim/commit/fb088428ae810ee3ddd893712588c1e64f6cc0c4) - adds proof-robust, schedule, and flutter certified end-to-end campaigns.

## 14. Latest Campaigns, AD Bridges, And p-Multigrid Smoothing

The latest 2026-07-08 wave adds more proof-carrying end-to-end examples and
tightens two core numerical surfaces: inverse-trig AD support and p-independent
multigrid smoothing.

### Delivered capability

- `fs-neuroshape-e2e` and `fs-grammar-e2e` add certified neural-shape topology
  and geometric-grammar campaigns, composing neural implicit charts,
  visualization/evidence, shape programs, archives, and fabrication checks.
- `fs-oed-e2e` adds SensorForge, a value-of-information sensor-placement
  campaign that stops when the posterior design choice is robust.
- `fs-vessel` adds the laminar-pour vessel flagship, combining Chebyshev
  stability objectives, free-surface LBM, CVaR robustification, e-racing, and
  spectral volume rendering.
- `fs-metamat-e2e`, `fs-truss-e2e`, `fs-adaptbo-e2e`, and `fs-flowcert-e2e`
  add certified metamaterial frontiers, optimal truss load paths, anytime BO,
  and CFD credibility maps.
- `fs-math` adds deterministic `asin`/`acos`, and `fs-ad` extends the `Real`
  surface with inverse-trig operations and chain-rule tests.
- `fs-ad` also gains the FrankenTorch tape bridge, binomial Revolve schedule,
  snapshot-store seam, and matrix-free IFT tangent route.
- `fs-solver` replaces p-MG Jacobi scaling with PU-symmetrized
  vertex-centered additive Schwarz smoothing and an exact r=1 Pavarino coarse
  term.

### Closed workstreams

- [`frankensim-epic-flagships-mye.4`](https://github.com/Dicklesworthstone/frankensim/blob/319cb64f052d76e15882ee53ace41092881c7fa8/.beads/issues.jsonl) - laminar-pour vessel flagship.
- [`frankensim-t88x`](https://github.com/Dicklesworthstone/frankensim/blob/319cb64f052d76e15882ee53ace41092881c7fa8/.beads/issues.jsonl) - inverse-trig `Real` operations for AD consumers.
- [`frankensim-o3ui`](https://github.com/Dicklesworthstone/frankensim/blob/319cb64f052d76e15882ee53ace41092881c7fa8/.beads/issues.jsonl) - `fs-ad` bridge, Revolve, spill, and IFT integrations.
- [`frankensim-x08j`](https://github.com/Dicklesworthstone/frankensim/blob/319cb64f052d76e15882ee53ace41092881c7fa8/.beads/issues.jsonl) - p-independent p-MG smoothing.

### Representative commits

- [`5cbdd90`](https://github.com/Dicklesworthstone/frankensim/commit/5cbdd903e8caf389e98eafe5ebc2436ac7a0093c) - add neural-shape topology and geometric grammar e2e campaigns.
- [`94404c4`](https://github.com/Dicklesworthstone/frankensim/commit/94404c4d8ef5f4fb7a33d9e46fddb32f51cb0aaa) - add SensorForge OED e2e.
- [`b95e00f`](https://github.com/Dicklesworthstone/frankensim/commit/b95e00ff7ee8140d0c1427824a23953af1353963) - ship the laminar-pour vessel flagship.
- [`8028bbc`](https://github.com/Dicklesworthstone/frankensim/commit/8028bbc29de12aff9fcb46849e688d0cd28bbab5) - add the metamaterial stiffness-density frontier campaign.
- [`4c573d6`](https://github.com/Dicklesworthstone/frankensim/commit/4c573d648bafd6d779def16e4ca7303fdb99f870) - add the optimal truss critical-load-path campaign.
- [`76e9c89`](https://github.com/Dicklesworthstone/frankensim/commit/76e9c89f62e91612b8886115c06ba8260c76cc94) - add AnytimeBO.
- [`922c835`](https://github.com/Dicklesworthstone/frankensim/commit/922c83597d2a79c1a862d9b8e37f5240b8ec9c8e) - add deterministic inverse trig and AD `Real` operations.
- [`3eb480c`](https://github.com/Dicklesworthstone/frankensim/commit/3eb480c8f4c1faa0d2b92a3b31e0d9db46256f0e) - add the FlowCert CFD credibility map.
- [`7575cdd`](https://github.com/Dicklesworthstone/frankensim/commit/7575cdd6bb0728794ef0fce4eae3512f9667b93d) - ship `fs-ad` adjoint integrations.
- [`319cb64`](https://github.com/Dicklesworthstone/frankensim/commit/319cb64f052d76e15882ee53ace41092881c7fa8) - add vertex-patch Schwarz smoothing to p-MG.

## 15. 2026-07-09 Mainline Expansion: Flagships, Decision Loops, And Proof Hygiene

The window from
[`319cb64`](https://github.com/Dicklesworthstone/frankensim/commit/319cb64f052d76e15882ee53ace41092881c7fa8)
through
[`438128d`](https://github.com/Dicklesworthstone/frankensim/commit/438128d988082f2183f406d75883b766fd6b7324)
is another implementation-heavy mainline slice, not a formal release. It adds
55 commits across 108 files, with 15,482 insertions and 323 deletions. The work
broadens flagship coverage, strengthens decision-making loops for optimization
and uncertainty, and fixes campaign checks that had been too correlated or too
underspecified to support the claims around them.

### Delivered capability

- `fs-ornith` lands as a smoke-tier ornithoid aircraft flagship with parameter
  decoding, BEM/VPM screening, e-raced candidate elimination, LBM refinement,
  SOS-style trim certificates, conformal lift/drag bands, Pareto-atlas rows,
  replay evidence, and a checked-in contract.
- `fs-topopt` gains a feature-gated `cutfem-marquee` lane for topology
  optimization over density-as-SDF fields, DWR-guided octree refinement,
  zero-remesh iteration logs, and topology/thickness witnesses.
- `fs-solver` adds Stokes-class saddle-point machinery: Schur approximations,
  block preconditioners, and PMINRES coverage over the FEEC mixed-space surface.
- `fs-dfo` expands from many-objective hypervolume accounting into NSGA-III,
  MOEA/D, contribution archiving, and world-fork steering. `fs-uq` adds chance
  constraints, seismic uncertainty helpers, anytime stopping, and adaptive MLMC
  admission checks.
- `fs-bo` adds inducing-point sparse Gaussian-process support and then hardens
  the sparse/exact comparison so approximate covariance loss is visible instead
  of smoothed over.
- `fs-ascent` adds interior-point and SQP constrained polish engines, a
  FrankenScipy oracle lane, Rosenbrock basin checks, and a Problem-IR study
  runner.
- `fs-render` adds edge-aware differentiable rendering and an inverse-rendering
  objective surface, while `fs-wasm` exposes deeper upper-stack kernels and a
  certified campaign browser tier.
- `fs-gen` lands as a proposal-only generation crate with compiler-enforced
  epistemic boundaries, making generated candidates explicit about what they
  are not allowed to certify.
- `fs-flagship-e2e` lands as an in-progress L6 replay-suite scaffold for the
  flagship family: smoke/mid/full stage wiring, metric-only content hashes,
  cross-flagship LBM and e-race audits, structured failure drills, forensic
  JSON rows, and a contract that keeps unfrozen smoke goldens explicit.
- `fs-mesh` hardens full-Ruppert refinement by checking candidate Steiner
  points against hull-facet diametral balls before insertion, while `fs-topopt`
  tightens the marquee tests around faster deterministic topology, benchmark,
  and thickness witnesses.
- End-to-end proof hygiene improves across `fs-flutter-e2e`,
  `fs-neuroshape-e2e`, `fs-flowcert-e2e`, `fs-truss-e2e`, and
  `fs-metamat-e2e`, including independent flutter stability checks, closed
  neural-shape boundary frames, steady-state credibility maps, corrected
  duality-gap enclosure direction, and PSD certificate symmetrization.
- HELM and flywheel surfaces gain a budget allocator, moonshot allocation
  tests, self-knowledge e2e battery, epistemic-engine acceptance gate, and
  Phase 3 horizon gate.

### Closed and active workstreams

- [`frankensim-epic-flagships-mye.2`](https://github.com/Dicklesworthstone/frankensim/blob/d5873bfd82a3c2dbe359c11aef4947a5def8cdba/.beads/issues.jsonl) - ornithoid multi-inlet aircraft flagship.
- [`frankensim-b7d0`](https://github.com/Dicklesworthstone/frankensim/blob/d5873bfd82a3c2dbe359c11aef4947a5def8cdba/.beads/issues.jsonl) - active CutFEM-octree topology-optimization marquee lane.
- [`frankensim-avuw`](https://github.com/Dicklesworthstone/frankensim/blob/d5873bfd82a3c2dbe359c11aef4947a5def8cdba/.beads/issues.jsonl) - Stokes-class saddle-point block preconditioners.
- [`frankensim-vcia`](https://github.com/Dicklesworthstone/frankensim/blob/d5873bfd82a3c2dbe359c11aef4947a5def8cdba/.beads/issues.jsonl) - many-objective, steering, and Wasserstein/DRO decision surfaces.
- [`frankensim-ijil`](https://github.com/Dicklesworthstone/frankensim/blob/d5873bfd82a3c2dbe359c11aef4947a5def8cdba/.beads/issues.jsonl) - constrained polish engines, FrankenScipy oracle, and Problem-IR runner.
- [`frankensim-epic-helm-gp3.9`](https://github.com/Dicklesworthstone/frankensim/blob/d5873bfd82a3c2dbe359c11aef4947a5def8cdba/.beads/issues.jsonl) - `fs-plan` allocation and replanning machinery.
- [`frankensim-epic-lumen-qfx.5`](https://github.com/Dicklesworthstone/frankensim/blob/d5873bfd82a3c2dbe359c11aef4947a5def8cdba/.beads/issues.jsonl) - differentiable rendering objective lane.
- [`frankensim-epic-flagships-mye.5`](https://github.com/Dicklesworthstone/frankensim/blob/9dc7417090e62e7ab2756719034e2ab0d3b32875/.beads/issues.jsonl) - active flagship e2e suite scaffold.
- [`frankensim-uee3`](https://github.com/Dicklesworthstone/frankensim/blob/1fe4ef5c21737a92aa58dd02f901052de5b0a88b/.beads/issues.jsonl) - active fs-mesh v2 follow-up.

### Representative commits

- [`7624964`](https://github.com/Dicklesworthstone/frankensim/commit/7624964781185c5d82879ba0e63aba9da304021d) - add the self-knowledge e2e battery.
- [`3dcbfaa`](https://github.com/Dicklesworthstone/frankensim/commit/3dcbfaaf3d3e90496719bb6a879500ca825e61b7) - add Stokes saddle-point block preconditioners and PMINRES.
- [`8c20bab`](https://github.com/Dicklesworthstone/frankensim/commit/8c20bab6da5ef694d4243ded896ac540d87502c8) - add the NSGA-III reference-direction lane.
- [`61e9348`](https://github.com/Dicklesworthstone/frankensim/commit/61e93489985b6cdbcea768535e6a12d12885f471) - add edge-aware differentiable rendering.
- [`fe8a2c2`](https://github.com/Dicklesworthstone/frankensim/commit/fe8a2c2061377e47293a3c638eb152337a97adcd) - add the `fs-plan` budget allocator.
- [`6d074b8`](https://github.com/Dicklesworthstone/frankensim/commit/6d074b8c92c396ac320815074640c63410f6df95) - add seismic uncertainty, anytime stopping, and adaptive MLMC.
- [`e2b3dff`](https://github.com/Dicklesworthstone/frankensim/commit/e2b3dffb1b4b58b7378e337e3134f4e4d95fb7ac) - add the sparse Gaussian-process lane in `fs-bo`.
- [`73c1c6b`](https://github.com/Dicklesworthstone/frankensim/commit/73c1c6b46f0f36673cc2bbda70164c0b11eee23c) - add `fs-gen` with proposal-only epistemics.
- [`51adc98`](https://github.com/Dicklesworthstone/frankensim/commit/51adc987414491b75ad4543fae7925e83804471a) - add steering and chance-constraint integration.
- [`2f29772`](https://github.com/Dicklesworthstone/frankensim/commit/2f29772b8536171929f2b4bae4ece8aa849caa3c) - add interior-point, SQP, oracle, and IR runner machinery.
- [`2199248`](https://github.com/Dicklesworthstone/frankensim/commit/219924842e0e0ebc71ec975b7d17bddf193db403) - expose the certified campaign tier in `fs-wasm`.
- [`ffcb4a9`](https://github.com/Dicklesworthstone/frankensim/commit/ffcb4a91e047b5f660d28083194da6daeb47b190) - add the smoke-tier `fs-ornith` crate.
- [`8c06c93`](https://github.com/Dicklesworthstone/frankensim/commit/8c06c93e7404b7bbd528eebc14d4af43914a8906) - add the feature-gated CutFEM-octree topology lane.
- [`d5873bf`](https://github.com/Dicklesworthstone/frankensim/commit/d5873bfd82a3c2dbe359c11aef4947a5def8cdba) - add the `fs-ornith` contract and close the ornithoid flagship bead.
- [`fbba704`](https://github.com/Dicklesworthstone/frankensim/commit/fbba7049d3f92714185f26ed1473503b78903281) - refine flutter and neuroshape certificate payloads.
- [`1092e94`](https://github.com/Dicklesworthstone/frankensim/commit/1092e946be32b211bf20f4748be2f44683feaeaf) - harden CutFEM-octree marquee boundary evidence.
- [`9dc7417`](https://github.com/Dicklesworthstone/frankensim/commit/9dc7417090e62e7ab2756719034e2ab0d3b32875) - add the staged flagship e2e replay-suite scaffold.
- [`1337058`](https://github.com/Dicklesworthstone/frankensim/commit/1337058ac6bc4b8f1dd9b7cdeb85db2fd0713bfd) - tighten the flagship replay-suite battery lint posture.
- [`1fe4ef5`](https://github.com/Dicklesworthstone/frankensim/commit/1fe4ef5c21737a92aa58dd02f901052de5b0a88b) - mark the fs-mesh v2 follow-up in progress.
- [`438128d`](https://github.com/Dicklesworthstone/frankensim/commit/438128d988082f2183f406d75883b766fd6b7324) - protect hull refinement and tighten marquee witnesses.

## Current Non-Release Status

- No package or crate release has been tagged yet.
- The canonical project state is `main`, not a versioned artifact; the latest
  implementation and tracker snapshot covered here is
  [`438128d`](https://github.com/Dicklesworthstone/frankensim/commit/438128d988082f2183f406d75883b766fd6b7324).
- The repository is actively changing; use crate `CONTRACT.md` files and Beads
  close reasons for detailed no-claim boundaries.
- The README describes the implemented workspace; the long-form plan remains
  `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md`.
- Uncommitted work in the local working tree is not covered by this changelog.

## Notes For Agents

- Start with the version timeline for chronology.
- Use thematic sections for architectural orientation.
- Use commit links for implementation evidence and Beads links for intent and
  acceptance criteria.
- Do not cite GitHub Releases for this repo until a real release exists.
- Treat DSR as the project validation source of truth unless `AGENTS.md`
  changes.
