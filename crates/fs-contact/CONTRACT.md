# CONTRACT: fs-contact

> Status: ACTIVE (Stage 1, increment 1 — bead tqag). Capability-routed
> body-to-body contact over certified motion.

## Purpose and layer

Blocker B3 (expansion plan, phase E2): body-to-body contact detection
with certificates instead of sampled heuristics. Layer: **L3** (deps:
fs-motion L2, fs-query L2, fs-geom L2, fs-ivl L1, fs-exec L0).
Explicitly NOT a dependency of-or-on fs-solid/fs-mbd solver internals —
those consume adapters; reusable contact protocols live here.

## Public types and semantics

- `normal_patch::NormalLawSource` presents `(U,F,dF)` from an immutable
  `SourcePointContact` or a `FiniteGapResponseCurve`. Point stiffness remains
  derived from the bound bulk/interface model and geometry; finite-gap values
  remain the original curve's integrated force, segment tangent and energy,
  with `Estimate` authority and no extrapolation. No power-law fit or new
  constitutive trait is introduced. Original point/curve errors remain typed.
  Generic storage composition and point lowering live in
  `fs-couple::unilateral_contact::binding`, preserving dependency direction.

- `SpacetimeBody` — a finite body-frame support box bound to a
  `fs_motion::CertifiedMotorTube` (body-to-world). Validation refuses
  non-finite/inverted supports.
- `spacetime_candidates(bodies, window, max_pairs, cx)` →
  `BroadPhaseReport`: the conservative spacetime broad phase. Each
  body's windowed world box is `CertifiedMotorTube::box_action_over`
  over the WHOLE window — an enclosure for every `t`, so a
  non-overlapping pair provably cannot touch inside the window (no
  sampled instants, no tunneling between samples). Deterministic
  sweep-and-prune on world `x` (`total_cmp`, index tie-breaks); output
  pairs sorted; report carries checked/pruned counts and the worst
  motion versor-defect bound, which consumers must carry forward.
- `NarrowRoute` / `narrow_phase(pair, route_a, route_b, iters, cx)` →
  `NarrowVerdict`: capability routing. Stage 1 routes Convex×Convex
  through fs-query's certified `convex_separation` (its semantics pass
  through unchanged: `separation_proven ⇔ lo > 0`, overlap never
  claimed). Any pairing without a compatible declared route refuses
  with `MissingCapability` naming the pair and capability — never a
  guess.
- `ContactError` — typed refusals throughout;
  `CandidateBudgetExhausted` (program risk #2) lists every unresolved
  overlapping pair beyond the budget so the resolved prefix is never
  mistaken for the complete candidate set.

### Finite-patch normal response (bead `frankensim-b8bxd.5.1`)

`normal_patch` is a solver-independent, SI-only finite-patch constitutive
surface. Elastic Hertz sphere/plane and cylinder/plane rungs retain force,
patch dimension, peak pressure, pressure second moment, reversible energy, and
tangent. `PointNormalPatchReceipt` uses `N`, `J`, and `W`; the distinct
`LineNormalPatchReceipt` uses `N/m`, `J/m`, and `W/m`, so zero-valued outputs
cannot erase the dimension distinction. The Hunt--Crossley sphere rung adds
explicitly passive rate-dependent work for a bounded loading/unloading step.
Receipts carry canonical request and
receipt identities, caller input-authority ceiling, uncertainty, and checked
half-space, curvature, small-strain, layer, yield, rate, temperature, and
nonadhesion applicability ratios. `ConstraintBarrierReceipt` is deliberately a
different type with no force, patch, or energy fields, so IPC/barrier output
cannot masquerade as physical compliance.

### Generic normal-patch embedding (bead `frankensim-b8bxd.5.2`)

`normal_patch::embed` maps public constitutive requests into a generic contact
port without accessing law internals or depending on multibody dynamics. A
`SmoothFixed` and converged solver sample alone can publish a typed point or
line action/reaction wrench, zero action/reaction residual, and tangent. The
immutable successor state records a deterministic exactly-once work key; a
checkpoint is the only rollback source. Eventful, stale, future, duplicate, and
nonconverged samples refuse before publishing a port. The embedded transition
retains the law receipt identities, uncertainty, and applicability fields.

### Finite-gap elastic half-space contact

`normal_patch::solve_finite_gap_half_space` is the reference rung for contact
patches whose exact undeformed separation is no longer adequately represented
by a local quadratic tangent. The caller supplies an arbitrary finite relative
gap on a uniform tangent-plane grid and a prescribed rigid approach. Exact
constant-rectangle integrals of the Boussinesq `1/r` surface-compliance kernel
couple the cells; a deterministic active-set solve enforces nonnegative
pressure, nonnegative final separation, and complementarity. The receipt
retains the complete pressure/opening fields, resultant, pressure centroid and
second-moment footprint, reversible elastic energy, active-set count, residual,
and a content identity.

`normal_patch::sample_finite_gap_from_chart` supplies the representation-neutral
geometry bridge. Given any solid `fs-geom::Chart`, a caller-bound geometry
identity, and an orthonormal support frame, it brackets the first
outside-to-inside surface crossing along the inward normal at every grid cell.
Its strict route requires the chart's rigorous continuous trace enclosure,
refuses missing outside/inside brackets or an unresolved zero-straddling
sample, and retains both endpoints of every height bracket rather than silently
collapsing geometry uncertainty into the nominal midpoint. A caller may
explicitly admit a finite nominal `NoClaim` chart for convergence-reference
work; any such sign decision downgrades the complete sampling receipt to
`Estimate`, and its retained intervals cease to be presented as rigorous
bounds. Consequently analytic, SDF, FRep, NURBS-derived, or mesh-derived charts
can feed the same contact solver under an explicit evidence policy; the adapter
contains no Euler-disc profile special case.

This is deliberately a convergence-reference and response-table generator,
not a per-audio-sample dense solve. Positive pressure inside a caller-required
boundary ring is a typed domain-size refusal, so a truncated pressure patch can
never be presented as a physical result. G1 manufactured coverage compares an
independently sampled paraboloid against the analytic Hertz force and footprint.

`normal_patch::build_finite_gap_response_curve` performs that offline-to-hot-path
transition without replacing the constitutive response by a fitted Hertz law.
It runs the dense solver at explicit, strictly increasing approach nodes,
requires nondecreasing force, and retains every underlying contact identity.
Runtime interpolation is piecewise linear in force and analytically integrated
over the same segment, so its tangent and stored energy obey
`dU/d(delta) = F`. The independent dense-solve strain energy remains attached
to every node; construction refuses when the curve-integral discrepancy exceeds
the caller's absolute-plus-relative refinement budget. Evaluation refuses
extrapolation rather than silently extending a table beyond its geometry/load
domain. The curve retains the source gap's authority separately, but every
interpolated response is explicitly `Estimate`: even an enclosed source gap
does not promote a piecewise-linear value between dense nodes without an
interval interpolation remainder.

`normal_patch::FiniteGapResponseFamily` composes curves on one identical
approach grid across a caller-named scalar configuration coordinate. Bilinear
force and stored-energy interpolation preserves `dU/d(delta) = F` at fixed
configuration, while both coordinate and approach extrapolation refuse. This
is the reusable hot-path surface for slowly changing geometry state such as
inclination; it remains `Estimate` and does not hide configuration-grid
refinement error.

### Source-bound JKR adhesion

`normal_patch::adhesive::SourceJkrSphereContact` supplies the short-range JKR
sphere/plane approximation from the original ordered dry interface, two bulk
elastic states and their explicit `adhesion-energy` work-of-adhesion claim.
Construction requires static sources at fixed external conditions, an explicit
attraction range, a Tabor parameter of at least five, and declared small-patch,
half-space, layer, nominal-strain, repulsive-pressure/yield and process-zone
limits. These are model applicability screens, not a certified approximation
error or a finite bound on the singular tensile edge stress. Original sources
remain attached to the owner; responses are model estimates, with explicitly
synthetic inputs retaining their synthetic classification.

The allocation-free evaluator preserves signed approach and force, actual
attached area, the Hertz repulsive contribution, the elastic-plus-surface
potential and its consistent tangent. A bounded safeguarded scalar solve avoids
subtractive cancellation near the displacement spinodal by using an offset from
its radius. `friction_contact_state` feeds the same area and repulsive load
to `fs-tribo` junction shear; a zero or tensile resultant never erases attachment.
`point_port` maps signed action/reaction forces and moments and refuses the
spinodal's unbounded tangent. Existing nonadhesive rungs retain their refusals.

`advance_monotone` follows explicit source-bound history tokens. Attachment at
zero approach and detachment at the displacement-controlled spinodal localize
each jump and report its released energy separately from smooth elastic work.
A prescribed-displacement segment satisfies `work = delta(stored energy) +
jump loss`; accepting the returned token prevents repeated detachment losses.
Load-controlled pull-off is exposed separately and is not used as a displacement
event. Finite machine stiffness, roughness, precontact finite-range attraction,
rate-dependent adhesion, viscoelasticity, mixed-mode peeling and experimental
validation require additional models. The original theory is Johnson, Kendall
and Roberts (1971), <https://doi.org/10.1098/rspa.1971.0141>.

`tests/adhesive_jkr.rs` reconstructs the original physical equations, independently
integrates cycle work, checks potential/force/tangent derivatives, distinguishes
force and displacement pull-off, exercises attachment and repeated separation,
checks zero-load/tensile friction and port signs, and verifies source scaling,
the weak-adhesion Hertz limit and applicability/history refusals.

`interface_binding::adhesive_shear::SourceAdhesiveShear` selects constant
junction shear from explicit `adhesive_static_shear_strength` and
`adhesive_kinetic_shear_strength` claims in Pa at the SAME ordered interface and
query point. Both source receipts are retained. Missing data, incorrect units,
unowned rate/history restrictions, inverted strengths or normal-state identity
drift refuse. No shear strength is inferred from work of adhesion or bulk yield,
and no extra Coulomb term is silently added. `evaluate_jkr` checks the response's
normal-source identity and supplies its real area to the existing passive shear
leaf, including zero and tensile net load. The original input authority remains
unchanged, and frictional power is returned once for the thermal owner.
`tests/adhesive_shear.rs` covers source/refusal cases and the coupled JKR-to-shear
path through zero load, tensile load and detachment. This uncoupled local model
does not evolve shear-induced peeling, partial slip, wear or bond aging.

### Ordered material/interface binding

`interface_binding` retains the complete immutable `InterfaceSystemCard`, its
ordered material and texture roles, and the dependency-light `fs-tribo`
interface identity without promoting authority. Two complete
`IsotropicSolidStatePoint` bundles may bind only when their material roles and
exact physical query points agree with the interface card. The derived Hertz
modulus and limiting yield stress retain both full bulk-property bundles and
their usage receipts. Normal-contact binding additionally requires an explicit
nonnegative, dimensionally correct `adhesion-energy` claim resolved from the
same interface card and exact state point; absence or card/state rebinding
refuses rather than implying zero adhesion.

### Smooth tangential adapter (bead `frankensim-b8bxd.7.2`)

`tangential::smooth` is a solver transaction adapter over the public
`fs_tribo::partial_slip` law. Its regularization and fixed-branch derivative
are local numerical controls, not a smoothness, friction-calibration, or
cross-branch differentiability claim. The creepage deadzone acts on the
tangent-plane slip MAGNITUDE (isotropic), never per frame component, so the
world-frame receipt is covariant under the caller's tangent-frame rotation —
the componentwise form was an executed defect (7e-5 relative world-force
frame dependence) caught when the parked SO(2) conformance test first ran.
The endpoint wrench and its work remain unchanged by energy bookkeeping:
`SmoothResidual` forwards the law's independent work-closure residual and
nonnegative backward-Euler numerical loss. This loss is separate from the
physical microslip heat, including on purely elastic intervals.
The retained constitutive checkpoint includes the accepted tangent frame.
Changing the query gauge re-expresses the preloaded state without changing
its physical wrench. Common material-frame spin passes through unchanged;
the relative-slip regularization must not damp this rigid transport. Prepare,
rollback and restart preserve the same frame/history transaction.

## Invariants

- Broad-phase candidacy is conservative over the query window: a pair
  absent from `pairs` has certifiably disjoint windowed enclosures.
- Output ordering is a pure function of the inputs (deterministic
  sort keys everywhere; no HashMap iteration).
- Refusals leave no partial claim: budget exhaustion returns the
  unresolved remainder, capability gaps name the pair.
- Nonadhesive normal-patch laws refuse non-finite/negative and overflowing inputs, invalid
  curvature/modulus, adhesive, plastic, layered, out-of-rate, or out-of-
  temperature requests rather than silently extrapolating. Solver
  complementarity, residual assembly, and resolved contact state remain
  consumer responsibilities.
- The request declares a local geometry. A sphere law cannot be evaluated as a
  cylinder, and toroidal or highly elliptical contact is a typed refusal rather
  than an effective-radius approximation.
- Embedding does not upgrade a constitutive receipt into a solved contact
  complementarity claim. Its residual is only the explicit action/reaction
  closure of an admitted fixed branch; active-set, global convergence, and
  event resolution stay with the consuming solver.

## Error model

`ContactError` wraps `fs_motion::MotionError` and
`fs_query::QueryError` unchanged (their teaching text passes through)
and adds contact-specific refusals: body-count/support/window
validation, candidate budget exhaustion with the unresolved list,
missing narrow-phase capability, cancellation.

## Determinism class

Bit-deterministic given deterministic inputs: sorted sweeps, fixed
tie-breaks, fs-motion/fs-query deterministic enclosures underneath.

## Cancellation behavior

`Cx` checkpoints per body enclosure and per sweep row; narrow phase
inherits fs-query's cancellation strides.

## Unsafe boundary

None. Workspace lints; no `unsafe` blocks.

## Feature flags

None yet. CCD lanes will gate under features when they land.

## Conformance tests

`tests/contact.rs`, cases ct-001..ct-004: analytic screw-motion broad
phase (approach window overlaps, retreat window disjoint, both against
hand-computed enclosure geometry); determinism replay; budget
exhaustion listing exact unresolved pairs; capability refusal; convex
narrow-phase distance containment at a frozen time against the
analytic value.

`tests/normal_patch_law.rs` supplies G0/G1/G3 independent reconstruction of
Hertz resultants, pressure moments, reversible energy, and numerical tangent;
cylinder approach, energy, tangent, and scaling; Hunt--Crossley
loading/unloading passivity; zero/grazing and adhesion/layer/yield/rate/
temperature refusals; typed point-versus-line receipts; and deterministic
identity/authority replay. Test cards are synthetic and do not admit materials
or a contact configuration.

`tests/normal_patch_embed.rs` supplies manufactured prescribed-approach and
point/line unit checks; action/reaction and moment reconstruction; exactly-once
power/work, checkpoint rollback/retry, fixed-branch tangent difference,
stale/future/duplicate/event/nonconverged refusals, toroidal-refusal
propagation, and embedding identity mutation.

## Certified CCD (bead tqag, increment 2)

`certified_ccd(a, b, window, time_tolerance, max_windows, cx)` proves
clearance or localizes possible contact by conservative window
bisection over `CertifiedMotorTube::box_action_over` enclosures:

- SOUNDNESS (the Sev-0 no-tunneling claim): a subwindow is cleared only
  when the two whole-subwindow image enclosures are disjoint along a
  coordinate axis — no instant inside it can produce contact, with no
  sampling anywhere. Everything not proven clear subdivides to the time
  tolerance and is reported as a possible-contact window; the union of
  reported windows contains every true contact instant. ct-005 drives a
  bullet fully through a thin plate INSIDE the window (both endpoint
  enclosures disjoint — the exact trap endpoint sampling falls into)
  and requires a possible window containing the true crossing,
  localized to under 1% of the window.
- HONESTY: contact is never CLAIMED (box overlap is necessary, not
  sufficient); `ClearWindow` carries a certified lower bound on the
  axis gap; budget exhaustion is a refusal carrying the exact pending +
  unresolved windows in time order, never a truncated verdict (ct-008).
- THE ROOT-GUARD REFUSAL, EXECUTABLE (ct-007): bodies overlapping the
  whole window have no separation sign change for a global-root guard
  to find; certified CCD reports one possible window covering the whole
  domain instead of a false clear. This is why the design bisects
  enclosures rather than guarding roots of `separation(t)`.
- Determinism: LIFO bisection with the earlier half examined first;
  reports replay bit-identically (ct-006).

## Swept-vertex-hull refinement (bead tqag, increment 3)

`refine_possible_windows(a_vertices, a_tube, b_vertices, b_tube,
windows, max_windows, max_iterations, cx)` re-tests each
`PossibleContact` window
for POLYTOPE bodies: every body-frame vertex trajectory is enclosed by
`point_action_over`, and the convex hull of all trajectory-box corners
contains the body's image at every instant (a rigid image of a hull is
the hull of the vertex images). Certified separation of the two swept
hulls (`fs_query::convex_separation`; corner selection is the support
trait's documented exact case) PRUNES the window with a certified gap —
tight exactly where per-instant axis-aligned boxes are structurally
loose. ct-009 pins the trap: two 45°-rotated cubes passing on the
diagonal whose AABBs overlap at EVERY instant (the box route can never
clear them at any tolerance) prune with the analytic edge-to-edge gap.
Soundness: pruning requires `separation_proven` over a SUPERSET of each
swept body, so a window containing true contact can never be pruned
(ct-010 keeps the bullet's crossing window Retained); retention claims
nothing, exactly like the box verdict.

## SDF-obstacle route (bead tqag, SDF increment)

`refine_windows_against_sdf(vertices, tube, obstacle, windows,
time_tolerance, max_windows, cx)` prunes possible-contact windows for a
polytope body against a STATIC exact-distance chart. Soundness: exact
Euclidean signed distance is 1-Lipschitz (the theorem carried by
`TraceStepClaim::ExactDistance` — weaker claims refuse at entry by
capability name), so with the swept-vertex-hull corners enclosed in a
ball of radius `r` around any center `c`, every swept point satisfies
`φ(q) ≥ φ_lo(c) − r`; a positive bound proves the whole subwindow clear
with a certified gap. The center is the overflow-safe component-wise
bounding-box midpoint, and every corner distance is evaluated with
outward-rounded `fs-ivl` operations. Subnormal squared-distance underflow
therefore remains enclosed; arithmetic overflow produces an infinite radius
and disables pruning rather than creating a false-clear certificate. The
center choice affects only tightness. Because certified_ccd MERGES adjacent
possible windows, the route bisects internally with the same
LIFO/tolerance/budget discipline (budget exhaustion refuses with exact
partial state). ct-011 pins the value (a corner-region pass the AABB
route retains forever prunes against the curved surface; a lying chart
with a weaker claim refuses at entry) and ct-012 the soundness arm (a
through-shot's sphere-entry window survives as Retained).

## No-claim boundaries

- Normal-patch receipts are local analytic constitutive evaluations, not
  finite-element contact solutions, material-card admission, friction or heat
  partitioning, lubrication/EHL, plasticity, general adhesion beyond the explicit
  source-bound JKR rung, roughness, calibration,
  experimental validation, Euler-disc target fitting, or any one-millimetre
  optimum/ranking claim. Hunt--Crossley dissipation is a declared point-contact
  model rung, not a measured loss mechanism or coefficient-of-restitution
  shortcut. It augments either the admitted sphere/plane coefficient or the
  actual two-principal-curvature elliptic Hertz coefficient; force, consistent
  tangent, elastic storage, irreversible work, pressure, and patch axes remain
  one constitutive receipt. Toroidal and highly elliptical patches outside the
  elliptic-law numerical envelope have no physical-compliance claim in this
  module; they require a separately validated law. The generic embedding is not a port-schema-v2
  implementation, cancellation proof, global solve, or a replacement for the
  upstream coupled-port contract.
- Material/interface binding proves exact identity, dimensional admission, and
  state-point consistency. Its normal-model bridge selects exactly one
  supported immutable interface-card law (or one explicitly pinned card),
  derives its rate scale, damping coefficient, applicability ratios, and finite
  temperature domain, and rejects missing, ambiguous, out-of-domain,
  dimensionally wrong, or surplus parameter data. This is executable schema
  admission, not physical validation: it does not infer missing
  friction/roughness/damping values, establish a surface as representative of a
  physical specimen, or extend a solid material through a phase transition.

- The finite-gap rung removes the *quadratic-gap* approximation only. It still
  assumes frictionless, nonadhesive, small-strain linear elastic half spaces and
  piecewise-constant pressure cells. Grid/domain refinement, half-space-depth
  and finite-layer applicability, yield/fracture, material uncertainty,
  roughness/load sharing, and comparison with a finite-body FEM/BEM or an
  experiment remain caller obligations. Chart-to-gap bracketing proves one
  first crossing along each requested normal line; it is not a global convexity,
  single-valued-height, support-point, or topology certificate. Its
  second-moment semiaxes summarize the pressure field; they do not assert a
  Hertz pressure shape. A response curve is a passive table approximation to
  its retained dense solves, not a certificate between approach nodes; approach
  refinement and inclination/profile-state interpolation remain caller-visible
  error dimensions.

- Certified CCD verdicts remain ENCLOSURE verdicts: `PossibleContact` /
  `Retained` windows localize in time but never adjudicate contact;
  time-of-impact enclosures tighter than window bounds, refinement for
  non-polytope MOVING bodies (needs swept support maps from fs-motion),
  and MOVING SDF obstacles (the chart is static in this increment) are
  later work in this bead's staging plan. Stage 2 consumes
  simulated-flow tubes through a tube-source-agnostic interface.
- Narrow-phase routes: Stage 1 is Convex×Convex only. SDF-pair local
  gaps (fs-query `ImplicitGapOracle`), nonconvex decomposition,
  interval global optimization, and mixed-route pairings all refuse
  as `MissingCapability` today.
- Penetration depth is never claimed (fs-query's convex overlap
  no-claim passes through); EPA-class certificates arrive with
  fs-query bead hk8f5.
- Rep Router conversion/motion errors do not yet inflate contact
  bounds (fs-query bead fugfk); claims apply to the presented charts,
  not to abstract regions behind conversions.
- The broad phase and CCD prune on certified geometry enclosures, but
  the motion versor defect is REPORTED (`BroadPhaseReport::max_defect`,
  `CcdReport::max_defect`), not folded into the boxes; the fold is
  still open in this bead's staging plan.

## Evolving partial-slip coefficients

`SmoothTangentialAdapter::prepare_with_coefficients` consumes the original partial-slip changing-coefficient continuation without resetting regional history. Actual force, slip, heat, and the successor law use the new coefficients; signed `generalized_work.coefficient_work_j` remains an independent parameter-work port. The coupled caller must balance that port and publish its coefficient driver together with the contact candidate. It is not caloric heat. Existing commit, rollback, restart and work-key deduplication remain active, and ordinary prepare/restore remain coefficient-exact. The smooth-contact regression checks changed reaction, rollback/retry, new-law restart, and duplicate-work refusal through the real adapter.

## Single-contact normal impulse law

`impact::solve_normal_impact` is an allocation-free frictionless Newton law for
one scalar normal contact. Inputs are finite nonnegative gap and activation
tolerance, finite relative normal speed (negative means approach), positive
inverse effective mass `w = J M^-1 J^T`, and restitution `e` in `[0,1]`. Every
negative gap refuses; the activation tolerance never authorizes penetration
repair. Inactive or nonapproaching contacts receive zero impulse.

For an active approaching contact the response is `p = -(1+e)*v/w`, with target
speed `-e*v` and nonnegative loss `(1-e)*(-v)*p/2`. The result includes activation
and impact flags and an independently reconstructed fused impulse/velocity
residual. Arithmetic avoids intermediate overflow when the final response is
representable; nonfinite outputs, unrepresentable positive impulse, and a
restitution-velocity residual above `64*EPSILON` relative to approach speed
refuse. Consistent SI inputs yield N s, m/s, and J. The Newton law and midpoint
impulse work follow [Acary (2014), equations 5 and 21–22](https://arxiv.org/html/1410.2499v1).

The caller owns geometry, contact selection, Jacobian/inertia assembly, actual
state update and global momentum/work checks. This owner neither certifies a
sequence of pairwise impacts as globally admissible nor implements simultaneous
contact, friction, a Moreau–Jean stepper, gap stabilization, or partition of loss
into heat/damage. Six focused unit regressions cover the independent energy
and restitution identities, activation/refusal, scaling and extreme finite
responses. The broad nonsmooth-contact bead remains open.

## Coupled nonsmooth finite-patch tangential contact

`tangential::nonsmooth` solves supplied positive-load finite patches together
using the original unregularized `fs-tribo` regional return maps and the existing
Newton–Krylov owner. The caller supplies B = J M^(-1/2), so mobility B B^T is
passive by construction and may be rank deficient for shared-body contacts.
Slip/spin coordinates have an explicit torsional length and conjugate forces.
Canonical patch ordering makes simultaneous transitions invariant to input
permutation. No branch smoothing is introduced.

Acceptance independently checks endpoint velocity balance, each regional
Coulomb cone and associated-flow complementarity, and mechanical work plus
stored-energy change, physical heat, distinct elastic/mechanical backward-Euler
losses, and signed coefficient work. Coefficient continuation preserves loaded
regional history; parameter work is not thermalized. Prepare/retry leave
accepted history unchanged. Commit consumes all contact work once. Release
requires explicit export of retained elastic storage; recontact starts unloaded.

Seven native physical regression groups cover analytic stick/partial/gross
slide, shared-body recoil and permutation, mixed slip/spin and frame covariance,
checkpoint/release/recontact, cancellation/nonconvergence, timestep refinement,
and changing coefficients. Normal-load/geometry admission, normal complementarity,
impact timing/restitution, Painleve classification, global uniqueness and
classical adjoints remain outside this tangential solver. Regional plastic
increments in `PartialSlipStep` expose the original constitutive flow for the
independent check without changing the return map.

### Low-Tabor adhesion with finite-range attraction

`normal_patch::dmt::SourceDmtSphereContact` supplies an explicitly low-Tabor
sphere/plane model. The contacting DMT-M branch retains Hertz geometry and
subtracts `2*pi*R*W`; an explicitly selected constant-traction Dugdale/Derjaguin
tail supplies finite-gap attraction up to the declared range. Potential and
force are continuous across touch and cutoff. One-sided tangents remain
explicit, and smooth point ports refuse those corners. Negative tail stiffness
remains visible to the enclosing mechanical solver.

Actual contact area excludes the attractive footprint. The source-bound shear
bridge's `evaluate_dmt` therefore produces no junction friction or frictional
heat across a gap, while finite attached area supports shear at zero/tensile
net load. Original source identities and synthetic authority remain intact.
Contact-load inversion includes zero and tensile load. Runtime evaluation is
allocation-free O(1), and conservative endpoint work invents no intrinsic
hysteresis or jump loss.

Tabor <= 0.1 and explicit range/curvature/elastic/stress limits are model-selection
screens, not an accuracy certificate. This is not the intermediate Maugis
crossover, rough contact, a fluid bridge or viscoelastic adhesion. Primary
references: <https://doi.org/10.1016/0021-9797(75)90018-1> and
<https://doi.org/10.1016/0021-9797(92)90285-T>. Seven `adhesive_dmt` regressions
cover original equations, derivative/energy consistency, corners, work cycles,
contact/noncontact friction, scaling and applicability; `adhesive_shear` adds
the original-source DMT-to-shear integration case.

## Simultaneous normal impact admission

`multicontact::solve_normal_impacts` consumes one full contact-space inverse
mass `W = J M^-1 J^T` and solves all geometrically active normal impulses
together through `fs-solver::orthant::active_set`. Active contacts include
initially separating witnesses: an impulse elsewhere may make them approach.
Each declared rebound target is `-e*min(v_before,0)`; the coupled quadratic
imposes nonnegative impulse, endpoint speed above that target and normal
complementarity. Geometrically inactive contacts receive zero impulse while
their cross-contact velocity response is still reported.

Every gap must be finite and nonnegative, including values inside the activation
tolerance. Restitution remains the caller's coefficient in `[0,1]`; positions
are not projected and coefficients are never reduced silently. A solution of
the normal complementarity problem alone is insufficient for acceptance. The
owner reconstructs all endpoint speeds from the original supplied matrix and
checks separation, complementarity and the assembled relative midpoint work
`sum(p_i*(v_before_i+v_after_i)/2)`. Positive relative work beyond the explicit
energy allowance and arithmetic roundoff returns `EnergyIncrease`. Any admitted
positive numerical work is retained as `passivity_residual`; dissipation is not
a hidden thermal-heat assignment. For a physical assembly, relative work equals
kinetic-energy change minus prescribed-support work.

The positive-definite numerical candidate admits at most 64 contact rows.
Redundant/singular or unresolved active contact directions refuse explicitly;
the API does not pick an arbitrary impulse from a nonunique set. Cancellation
returns no accepted partial event. The existing one-contact law is unchanged.
This is a bounded frictionless instantaneous response, not CCD, a Moreau–Jean
time stepper, persistent contact, a Coulomb-cone solve, or a general
Painleve/nonuniqueness classifier. The nonsmooth-contact bead remains open.
