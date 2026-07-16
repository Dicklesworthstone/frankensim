# fs-opdsl CONTRACT

## Purpose and layer

Layer: **L3 FLUX** (deps: fs-feec L3, fs-material L3, fs-la/fs-sparse
L1, fs-tilelang/fs-qty/fs-math L0/UTIL, fs-rep-mesh L2). The typed
operator IR (plan patch Rev C — "the single most important
implementation addition"): physics operators represented symbolically
over FEEC building blocks and lowered so that PRIMAL apply, JVP,
VJP/discrete adjoint, DWR indicators, preconditioner hints, and MMS
studies come from ONE SOURCE OF TRUTH. Generalizes fs-vskeleton's
`EdgeLaw` seed. "The primal changed but the adjoint didn't" is
structurally impossible: both are derivations of the same tree.

## Public types and semantics

- `Space { degree, n, dims }` — cochain degree (255 = raw vector
  dofs), dof count, `fs_qty::Dims`. Composition requires equality.
- `OperatorDef` — atom/law registry + the unknown's space. Checked
  builders: `apply` (degree AND dims checked → `TypeError`
  `ApplyMismatch`), `add` (`AddMismatch`), `scale` (folds nested
  scales), `pointwise`. An ill-typed operator cannot be represented.
- `Atom` — materialized linear building blocks: `Atom::d` /
  `Atom::d_transposed` (exact ±1 incidence CSR from fs-feec),
  `Atom::mass` (Whitney Galerkin star, symmetric), `Atom::external`
  (the ESCAPE HATCH: hand matrices with `Transpose::Symmetric` /
  `Derived` / `Explicit` — hand atoms pass the same gates and are
  marked `"provenance":"hand"` in every report).
- `PointwiseLaw` — opaque-but-differentiable constitutive nodes
  (value + exact derivative + dims map from one definition);
  `CubicReaction` is the reference law.
- `LoweredOperator` (from `OperatorDef::lower`): `apply` (primal),
  `linearize(u0)` (caches pointwise inputs; slots computed
  STRUCTURALLY in post-order so nested laws index correctly in both
  walks), `jvp` (forward chain rule), `vjp` (reverse chain rule —
  every factor transposed: the discrete adjoint), `materialize()`
  (folds constant linear chains through deterministic
  transpose/spgemm; `None` when nonlinear), `report()` (provenance,
  block hints, tile-kernel metadata, structural fingerprint —
  deterministic).
- `dwr_indicators(r, z)` — algebraic DWR η_i = |r_i·z_i|.
- `mms_poisson_study` — the G1 hook: load from analytic forcing via
  de-Rham + M0, Dirichlet pinning, PCG solve, L2 errors and observed
  orders across a Kuhn ladder, THROUGH the generated operator.
- `fixtures::{poisson, convection_diffusion, reaction_diffusion,
  elasticity, advection_matrix, elasticity_stiffness}` — the
  acceptance operators, defined in the IR. Elasticity's stiffness is
  constitutive-assembled from fs-material's `IsotropicElastic`
  tangent (B-matrix P1, Voigt engineering-shear convention matching
  fs-material).

## Invariants

- d∘d = 0 is honored AT THE IR LEVEL: applying dₖ₊₁ to the image of
  dₖ (also through scale nodes) folds to `Expr::Zero` before any
  float exists.
- Materializing the Poisson chain uses the same spgemm association
  as fs-feec's `stiffness`, so generated == hand BITWISE.
- JVP and VJP are two walks of the same tree with the same
  linearization state: ⟨Jv, w⟩ = ⟨v, Jᵀw⟩ holds mechanically (gated
  per fixture, including the nonsymmetric advection where the gate is
  not vacuous).
- Every scale/add/axpy the executor runs is an fs-tilelang kernel
  with intensity metadata and tier-equivalence twin tests.

## Error model

Construction returns structured `TypeError`s (`ApplyMismatch`,
`AddMismatch`, `UnknownId`) with both sides' spaces. `lower` panics on
ill-typed trees (unreachable through the checked builders). Derivative
calls on nonlinear operators without `linearize` panic with an
actionable message. Shape mismatches in fixture assembly are
structured panics.

## Determinism class

Bit-deterministic cross-ISA: atoms materialize through fs-feec/
fs-sparse deterministic paths; the executor's combinators are
fs-tilelang kernels (bitwise across tiers); plan reports are pure
functions of the registry and tree. Golden FNV-64 over Poisson apply,
convection–diffusion JVP+VJP, linearized reaction–diffusion JVP, and
elasticity apply on kuhn(2): `0x8b28_77cc_cb43_7cbc`, recorded on
Apple M4 Pro (aarch64), verified identical on Threadripper (x86_64).

## Cancellation behavior

Bounded synchronous evaluation; drivers own chunking and Cx poll
points (the workspace L0–L3 discipline). No internal threading.

## Feature flags

None.

## Unsafe boundary

None. `unsafe_code = "deny"`.

## Conformance tests

`tests/opdsl_battery.rs` (11 cases, JSON logging): materialized
Poisson vs hand FEEC bitwise + matrix-free to 1e−13; adjoint identity
⟨Av,w⟩ = ⟨v,Aᵀw⟩ over random vectors for Poisson/convection–diffusion/
elasticity; nonlinear JVP-vs-VJP transpose consistency + generated
chain rule vs fs-ad `Dual64` through the cubic law (both the law's
derivative and the full directional derivative); dd = 0 folded
symbolically (including through scale nodes) and applying as exact
zero; structural rejection of degree, dims, and add mismatches;
generation determinism (regenerated reports byte-identical) with
provenance and kernel-metadata assertions; DWR indicator
localization; MMS Poisson orders ≈ 2 on the n = 4/8/16 ladder THROUGH
the generated operator; elasticity rigid modes (translation +
infinitesimal rotation ≤ 1e−9) and uniform-strain patch (interior
residual ≤ 1e−12); perf comparison hand vs materialized vs matrix-free
MEASURED and logged (sanity bounds asserted; debug-build numbers
documented in the bead close); cross-ISA golden hash.

## No-claim boundaries

- The single-field operator surface remains as documented above; the
  `system` module (bead i94v.1.1.1) adds the multi-field TYPE layer on
  top of it: `SystemDef`/`AdmittedSystem` declare block fields with
  explicit form degree, six-base dims or fs-qty semantic quantity kind,
  basis/frame/orientation references, clock reference, spatial support,
  and state ownership, and admit cross-field structure (sums, atom
  applications, power-conjugate port pairings) BEFORE lowering.
  Ill-typed contractions, mixed frames/clocks, affine-temperature
  misuse (scale/sum/apply on absolute temperature), non-power-conjugate
  pairings, dangling references, non-finite scales, and beyond-cap
  nesting are structured `SystemTypeError` refusals. The admitted
  system mints a `SystemId` (fs-blake3 canonical identity, domain
  `org.frankensim.fs-opdsl.system.v1`) hashing canonical structure
  only: display names are never hash inputs and field/equation tables
  are canonically re-sorted, so renaming/serialization order preserve
  identity while any convention change moves it; byte-identical field
  payloads refuse as ambiguous. NO-CLAIM: the system layer performs no
  lowering, numerical evaluation, or block-preconditioner hinting (the
  solver-stack bead tfz.10 consumes it); pullbacks/clock transfers are
  refusal boundaries, not yet operators; `SemanticType` and `PortKind`
  canonical bytes ride their stable Debug renderings versioned under
  `SYSTEM_IR_VERSION` until fs-qty/fs-couple expose canonical
  encodings; conservation-role bookkeeping stays with fs-couple port
  schemas. All traversals are explicit-stack iterative (depth cap
  refuses work, not recursion). `system::transport` is the versioned
  canonical text transport (magic `fs-opdsl-system-transport-v1`,
  LF/tab records, strict fail-closed parse): any other IR version is a
  `VersionMismatch` refusal pending audited migration; import rebuilds
  a `SystemDef` and RE-RUNS full admission so tampered payloads refuse
  or mint a different identity — the transport carries no authority.
  Round-trip canonicality and the pinned migration golden
  (`sys_014`) hold the identity still; the golden moves only with a
  deliberate `SYSTEM_IR_VERSION` bump and recorded cause. Pullbacks and
  clock transfers are now EXPLICIT OPERATORS, not just refusal
  boundaries: `TransformSignature`/`ClockTransferSignature` register
  named, content-referenced maps between coordinate conventions/clocks,
  and `SystemExpr::Pullback`/`ClockTransfer` are the only admitted
  cross-frame/cross-clock/cross-orientation compositions — wrong
  endpoints refuse with both sides named
  (`TransformEndpointMismatch`), mixed orientations without a
  transform refuse (`OrientationMismatch`), and both tables follow the
  atom identity discipline (content-bearing, canonically sorted and
  remapped, byte-identical payloads refused). NO-CLAIM: v1 transforms
  are space-preserving type-level declarations — the numeric transport
  of values between frames/clocks is the caller's registry's job; the
  IR binds WHICH map was used, never that the map is correct.
- Pointwise laws are scalar diagonal (dof-local). Tensor-valued
  constitutive nodes (hyperelastic energy through `fs_ad::Real`
  generics) and state-dependent laws (plasticity history) are the
  fs-material integration lane.
- MMS hook covers the symmetric-PCG Poisson family; nonsymmetric MMS
  solves and mixed-form studies join tfz.10 (needs nonsymmetric
  Krylov). Convection–diffusion is gated by residual/adjoint checks
  here, not a full MMS solve.
- DWR indicators are the ALGEBRAIC dof form |r_i z_i|;
  element-integrated dual-weighted forms join the higher-order bead
  (tfz.6).
- No common-subexpression elimination beyond scale folding, zero
  elimination, and chain materialization; no symbolic wedge/trace
  nodes yet (they enter with the forms they serve — CutFEM traces
  with fs-cutfem, wedges with the nonlinear FEEC bead).
- Generated kernels: combinators (scale/add/axpy) lower into
  fs-tilelang; SpMV stays in fs-sparse (scatter is a tilelang
  no-claim). Fused matrix-free apply (allocation-free plans) is the
  tilelang-fusion perf lane.
- Perf numbers are debug-build fixture measurements, documented not
  contractual; the ≥80%-of-hand release-build target is the perf-CI
  lane's gate (fz2.4).

## Constitutive-graph adaptation (bead i94v.1.1.2, feature `constitutive-graph`, default-off)

- I01.2 opaque protocol: canonical fs-material law nodes enter ONLY as
  `&dyn LawNode` behind `BoundConstitutiveNode`; the compiler retains its
  OWN `MaterialProvenance` receipt (law/version, state-schema
  version+arity, port dims, compiler-owned differentiability and
  thermodynamic-potential-chart spellings) — proven to survive binding
  exactly with single-field mutation visibility. fs-material never
  imports compiler internals; the dependency arrow points L3-ward only.
- A supplied tangent is EVIDENCE: binding verifies consistent-tangent
  claims against the L3 finite-difference gate before granting the
  `Consistent` lane; `PiecewiseSmooth` routes to a caveated lane;
  `NonSmooth` claiming a tangent refuses typed; unclaimed tangents route
  `DerivativeFree`, where tangent/VJP requests refuse instead of
  differentiating. VJPs are exact transpose contractions of the verified
  tangent. State-owning nodes demand explicit initialization; the state
  codec refuses schema drift with both versions named. The hand-written
  escape hatch binds under identical gates and permanently retains a
  no-generated-consistency marker.
- Batched evaluation (`evaluate_batch`) runs a canonical
  `ConstitutiveGraph` over material-point batches under `fs_exec::Cx`
  with request-drain-finalize: polls at `BATCH_POLL_STRIDE` boundaries,
  drains whole points (never partial), finalizes with a deterministic
  resume cursor whose resumed remainder is bitwise-equivalent to the
  uncancelled run; a defective point refuses the whole batch typed with
  its index named.
- NO-CLAIM boundaries: no SystemExpr grammar splice yet (bound nodes are
  not yet referencable from `AdmittedSystem` equations — that increment
  is identity-sensitive and coordinated with the I01.1 owner); ambition
  [F] promotion (baseline deck, activation threshold, kill criterion,
  reproduction pack) has NOT begun — the feature stays default-off; the
  binding-time tangent gate probes the zero point only (per-application
  re-verification belongs to the generated-program gates).
