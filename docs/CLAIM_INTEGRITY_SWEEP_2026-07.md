# Claim-Integrity Audit Sweep — 2026-07

Bead `frankensim-extreal-program-f85xj.2.2`. One full recorded pass over the
claim-emitting surfaces using the decision rules, gap shapes, and audit method
in [`docs/CLAIM_INTEGRITY.md`](CLAIM_INTEGRITY.md).

This was an **inventory**, not a repair: nothing was fixed in the sweep itself
(not even the doc-only findings) so the gate's open-P0 count stayed an accurate
measure of exposure rather than a number someone quietly worked down while
filing it.

> **Burn-down status (2026-07-22).** The repair wave that followed closed **35**
> of the 39 findings outright and left 2 partial and 2 deferred. The live
> inventory is now 47 defects, **40 closed / 7 open**, every open one owned, and
> `scripts/ci/claim_integrity_inventory.sh` exits 0. The per-surface verdict
> tables below describe the state **as swept**; consult the beads for current
> state. Still open: `.24` (fs-solver residual provenance — needs a cross-crate
> slice through `fs-bem`), `.31` (the vessel half of the e-race audit — needs a
> public race wrapper in `fs-vessel`), `.34` and `.39` (both blocked on files
> carrying another agent's in-flight work), plus the pre-existing `0547u`,
> `hukmw` and `2gs5h`.

## Result

| Measure | Value |
| --- | --- |
| New defects filed | 39 (`f85xj.2.5` … `f85xj.2.43`) |
| `severity:default-path` (P0) | 29 |
| `severity:gated` (P1) | 1 |
| `severity:doc-only` (P2) | 9 |
| Retro-tagged known instances (bead `.2.1`) | 8 (3 open, 5 closed) |
| Total live inventory after the sweep | 47 defects — 42 open, 5 closed |
| Open gating (P0) defects | 32, of which **29 are unowned** |

The 29 unowned P0s are the burn-down's starting value and are why
`scripts/ci/claim_integrity_inventory.sh` currently exits 1. That red is
correct: it is the honest statement that 29 verified over-claiming defects
exist and nobody is assigned to any of them. Triage-assignment, not
re-definition, is what turns it green.

## Method and its limits

Six independent auditors worked disjoint surface groups. Each was given the
definition verbatim and instructed to (1) write the strongest permitted
inference from the type and CONTRACT **before** reading the implementation, so
the implementation could not anchor the reading, (2) construct adversarial
instances rather than trust the tests, and (3) record CLEAN and UNAUDITED
verdicts as distinct outcomes.

Every finding carries a reachable repro and the honest claim the surface should
make instead; findings without both were not filed. Confidence is recorded per
bead: 24 HIGH, 13 MEDIUM, 2 LOW. The two LOW findings
(`f85xj.2.27`, `f85xj.2.40`) are filed with their reachability explicitly
unproven and may legitimately close as caller-error or not-a-defect.

**Recall check (honest answer: the sweep did not re-find the known-answer set).**
The two seeded known instances — `frankensim-flutter-boundary-witness-hukmw`
and `frankensim-0547u` — were both already remediated in the working tree or at
HEAD by the time the auditors read the code, so they were *verified as fixed*
rather than rediscovered. One auditor independently flagged the truncation
hazard in `lyapunov_boundary`/`eigen_boundary` from the field types alone
(hukmw's exact mechanism) before consulting the diff, which is partial credit;
neither auditor can claim the same for `0547u`, having approached that surface
through the already-typed `ComponentCountEvidence` API. The sweep's recall
against a *live* instance of these two shapes is therefore **unmeasured**, not
demonstrated. Both in-flight fixes were reviewed and judged sound.

## Surface verdicts

Bead IDs below are `f85xj.2.<n>`, shortened to `.<n>`.

### fs-evidence — the colour spine

| Surface | Verdict | Beads |
| --- | --- | --- |
| `color::color_of` / `compose` dispersion accounting | DEFECT | `.5` |
| `color::verified_from` / `validate_color_payload` (degenerate infinite interval) | DEFECT | `.6` |
| `DiscrepancyModel::{query,evidence_at}` / `DiscrepancyBand::max_rel` | DEFECT | `.7` |
| `StatisticalCertificate::rel_width` (EValue arm) | DEFECT | `.8` |
| `Evidence::assess` / `escalation_advice` | DEFECT | `.9` |
| `ModelBracket::evidence` (`in_domain`) | DEFECT | `.10` |
| `IntervalOp::Hull` rustdoc | DEFECT (doc) | `.11` |
| `Evidence::certified` / `CertifyError`; `Certified<T>` opacity | CLEAN | — |
| `color::compose` rank law, disjoint regimes, outward rounding | CLEAN | — |
| `admitted.rs` receipt gating, deny-all default | CLEAN (inherits `.6`) | — |
| `cards.rs`, `balance.rs` ontology/intervals/defect composition, `action.rs` cost states, `falsify.rs`, `vv/model.rs` schemas | CLEAN | — |
| `vv/model.rs` remaining ~4k lines of schemas; `vv/codec.rs` (~60 decoders); `identity.rs` (~5k lines, ~40 helpers); `balance.rs`/`action.rs` canonical encode/decode | **UNAUDITED** | — |

Confirmed still holding: the `wa8i` NaN-propagation fix, the `f5pr`
disjoint-regime demotion, and the `6pf9` declared/admitted split. Findings `.5`
and `.6` are siblings of those classes at sites the original fixes did not reach.

### fs-package / fs-checker — admission and release gates

| Surface | Verdict | Beads |
| --- | --- | --- |
| `coverage::rank_presence` (VerifiedColor / EstimatedColor rows) | DEFECT | `.12` |
| Anchor-subject binding described in both CONTRACTs | DEFECT (doc) | `.13` |
| Degenerate infinite enclosure admitted at `verify`/`verify_with` | DEFECT | `.6` (shared) |
| `Claim` sealed constructors, origin verifiers, `deny_all()` defaults | CLEAN | — |
| `try_merkle_root`, header/leaf/node domains, signature exclusion | CLEAN | — |
| `verify_with` capability transcript; composition-receipt re-derivation | CLEAN | — |
| Waiver MAC/expiry/taint cone; `ColorBreakdown`; `MagnitudeBudget` | CLEAN | — |
| `VerificationReceipt`, `AuthenticatedSignature`, `VerifiedPackage` binding | CLEAN | — |
| `fs_checker::check*` staging, release preflight, decision hash, `render_pie` | CLEAN | — |
| `semantic::verify_portable_semantics`, `exact-interval@1`, `bounded-linf-residual@1` | CLEAN | — |
| `to_json` claim emitters; `parse_falsifiers`; `receipt_catalog` byte primitives | **UNAUDITED** | — |

The self-consistency-versus-authenticity distinction is enforced where it
matters: `receipt_hash` documents itself as "unkeyed integrity, not independent
authenticity", and no callback runs before an expected-root/structural pass.

### fs-geom / fs-topo — geometry and topology claims

| Surface | Verdict | Beads |
| --- | --- | --- |
| `fs_topo::self_intersection_certificate` (vertex-sharing pairs skipped) | DEFECT | `.14` |
| Non-finite coordinates pass both certificates | DEFECT | `.15` |
| `manifold_certificate` `MisorientedEdge` sentinel | DEFECT | `.16` |
| `ManifoldReport::oriented` without a probe | DEFECT | `.17` |
| `fs_geom::fixtures::*::topology_hint` (`BettiBounds::exact`) | DEFECT | `.18` |
| `convert::SampledSdf::bound` / `nominal_field_bound` | DEFECT | `.19` |
| **Sheaf seam result** (`AdmittedSheafComplex::watertightness`, raw-vs-admitted split) | **CLEAN both directions** | — |
| **Outside-to-outside ray routine** (`validate_outside_ray_samples`) | **CLEAN both directions** | — |
| `Convert` promotion to `Certified`; `ConvertDiag` refusals; `check_agreement` | CLEAN | — |
| `fs_topo::cubical::{voxelize,betti,persistence0,verify_topology}`; `penalty` | CLEAN | — |
| `fs_geom::router` planner half (~3k lines); `sheaf_repair` numerics (~6k gated lines); `sheaf_merge`; `derived`/`derived_morphism` (~28k gated lines); `exit_path`; `diff` | **UNAUDITED** | — |

The two surfaces the README calls out by name — the seam result and the ray
routine — were audited specifically for both over- and under-enforcement and
are clean: `AdmittedSheafComplex` has a single private construction site, raw
complexes are pinned to `Unknown`/`NoClaim`, a non-empty-interface guard blocks
vacuous PASS, and the ray routine's toggles are typed as replay telemetry with
no promotion authority. No `0547u`-shaped defect exists inside this scope.

### fs-adjoint / fs-solver — certificates and solver reports

| Surface | Verdict | Beads |
| --- | --- | --- |
| `verify::verify_gradient` (zero-vs-zero probe) | DEFECT | `.20` |
| `transpose::fd_falsifier` absolute acceptance band | DEFECT | `.21` |
| `transpose::check_transpose` (zero probes / empty dims) | DEFECT (gated) | `.22` |
| `transpose::Vjp` cotangent lengths; `Tape::transpose` | DEFECT | `.23` |
| `krylov::SolveReport::{rel_residual,converged}`; PMINRES semidefinite guard | DEFECT | `.24` |
| `pmg::PMultigrid` coarse `PcgReport` discarded | DEFECT | `.25` |
| `ift::ift_gradient_matfree` zip truncation | DEFECT | `.26` |
| `hessian::density_misfit_hvp` vacuous tolerance | DEFECT (LOW) | `.27` |
| `LinearOp::apply_transpose` CONTRACT sentence | DEFECT (doc) | `.28` |
| **`transpose::VjpRegistry` missing/declared-seam refusal** | **CLEAN — claim verified** | — |
| **`dwr_accept` (all branches Estimated, sealed binding, non-promoting bracket)** | **CLEAN** | — |
| `ift::AdjointReport`, `sobolev`, `hadamard`, `timedep`, `certs::{adjoint_residual_bound,certify}`, `mitigate`, `explain::finalize` | CLEAN | — |
| `explain` attribution engines and fingerprint machinery; `dwr_accept` work/cancellation accounting; `pmg` smoother numerics; `block` operator bodies | **UNAUDITED** | — |

The specific README claim that ledger transposition "refuses missing or
declared non-differentiable seams loudly instead of returning a silent zero"
**holds as stated**. The silent-zero risk lives one level down, in the *shapes*
of registered VJP cotangents (`.23`) and in the falsifier meant to catch a zero
that slips through (`.21`).

### Campaign crates

| Surface | Verdict | Beads |
| --- | --- | --- |
| `fs_thrust_e2e` surrogate validity domain (descriptor hull misses `d`) | DEFECT | `.29` |
| `fs_thrust_e2e::full_color` Verified from impulse residual | DEFECT | `.30` |
| `fs-flagship-e2e` fe2e-006 "cross-flagship" e-race audit | DEFECT | `.31` |
| `fs-flagship-e2e` fe2e-007 budget-exhaustion + escalation drills | DEFECT | `.32` |
| `fs-flagship-e2e` fe2e-007 ledger crash-recovery drill | DEFECT | `.33` |
| `fs_flutter_e2e::FlutterReport::witness_color` | DEFECT | `.34` |
| `fs_flagship_e2e::lbm_core_roll_hash` coverage prose | DEFECT (doc) | `.36` |
| `fs-thrust-e2e` CONTRACT step-savings / equal-quality claims | DEFECT (doc) | `.37` |
| `FlutterReport::boundaries_agree` | ALREADY FILED — fix on disk reviewed sound | `hukmw` |
| `NeuroShapeReport` component evidence; WASM component fields | ALREADY FILED — fix reviewed sound | `0547u` |
| `campaign_rank` weakest-elite rule; `full_color` malformed-input handling; conformal `alpha` admissibility | CLEAN | — |
| `fs_ornith::certify`, `LdSurrogate`, `screen_generation`, `refine`; fe2e-001/002/003 goldens; fe2e-007(a) cancellation storm; fe2e-008 forensics | CLEAN | — |
| `fs-flagship-e2e/tests/production_scale.rs` | **UNAUDITED** | — |

Four sub-threshold observations were examined and deliberately **not** filed
(recorded so they are not mistaken for unaudited): `fs_ornith::Atlas::rows`
band/selection-bias framing; the fs-ornith CONTRACT "sign + order-band"
sentence; `NeuroShapeReport::nearest_surface_radius`'s stale field doc; and
`CampaignReport::best_drift` being uncoloured.

### fs-wasm and repository documentation

| Surface | Verdict | Beads |
| --- | --- | --- |
| `campaigns::flowcert` (fixed budget → `accurate` / Verified rank) | DEFECT | `.38` |
| `campaigns::neuroshape` `safe_radius` | DEFECT | `.39` |
| `dynamics::ga_motor_orbit` Err-arm fallback | DEFECT (LOW) | `.40` |
| `campaigns::fluttercert` rustdoc independence claim | DEFECT (doc) | `.35` |
| `flagships::run_vessel` slot [10] "certified min-max" | DEFECT (doc) | `.41` |
| README "Evidence, Packages, and Standards" Phase-0B sentence | DEFECT (doc) | `.42` |
| README `fs-mesh` "measured 10-million-point perf lane" | DEFECT (doc) | `.43` |
| `certified.rs::mandelbrot_certified`; `deep.rs::cutfem_quadtree`; `campaigns::{proofrobust,metamatcert,schedule_campaign,trusspath,sensorforge,grammarforge}`; `flagships::{run_ornithoid,run_frame}` | CLEAN | — |
| README §What "Certified" Means Here, §Representative Vertical Slices, §Limitations, §Evidence Package Lifecycle, crate tables | CLEAN | — |
| CONTRACTs sampled against code: `fs-wasm`, `fs-ivl`, `fs-sos`, `fs-evidence`, `fs-govern`, `fs-checker`, `fs-mesh`, `fs-flowcert-e2e`, `fs-flutter-e2e` | CLEAN | — |
| `fs-neuroshape-e2e/CONTRACT.md` (`safe_radius` invariant, removed fields) | DEFECT | `.39` |
| The other ~125 `crates/*/CONTRACT.md` files | **UNAUDITED** | — |
| README inventory counts (126/127/276 vs 135/136/506) | pre-existing, tracked | `huq.18` |

## Unaudited surface is not clean surface

The largest residual risk in this sweep is coverage, not verdicts. Named
explicitly above and repeated here so a later reader cannot mistake silence for
a clean bill: `fs-geom`'s gated `derived_morphism` (~28k lines), `sheaf_repair`
numerics, `sheaf_merge`, and the router planner; `fs-evidence`'s `identity.rs`
and the bulk of `vv/`; `fs-package`'s JSON emitters and parsers; `fs-adjoint`'s
explain-attribution engines; `fs-solver`'s smoother numerics and block operator
bodies; `fs-flagship-e2e`'s production-scale battery; and ~125 crate CONTRACTs.
A second pass should start there.

## Two blockers found incidentally (not claim-integrity defects, not fixed here)

1. **`xtask` does not compile in the working tree.** Two `E0282` type-inference
   errors in the uncommitted ~2143-line `xtask/src/main.rs` WIP, at pre-existing
   lines 3069 and 3398 whose enclosing functions are byte-identical to HEAD —
   something else in that WIP broke their inference. Reported to its owner; not
   touched.
2. **`fs-wasm` does not compile against the working tree.** It reads
   `report.safe_radius`, which the in-flight `fs-neuroshape-e2e` change removes
   in favour of `safe_step: SafeStepDerivation`. That propagation gap is the
   code half of finding `.39`.
