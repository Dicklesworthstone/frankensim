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

### JKR applied load and compliant fixtures

`SourceJkrSphereContact::evaluate_force_controlled` evaluates the stable applied-load
branch from the closed-form JKR radius, including its marginal pull-off point.
It retains the original contact response and rejects loads below pull-off or
outside the original indentation envelope; it does not advance servo/history
state or invent a finite compliance at the force-controlled fold.

`normal_patch::adhesive::fixture::JkrSpringFixture` adds an explicit linear loading
spring between the contact and a prescribed stage. It solves `F = K(D - delta)`
and admits only the branch stable under both contact-radius and fixture
perturbations. The fixture fold satisfies `dF/d(delta) = -K`, interpolating
between force- and displacement-controlled pull-off. Total storage includes
contact plus spring energy; its derivative with respect to stage displacement
is the stage force. The stage port exposes the statically condensed tangent,
with the fold marked nonsmooth.

Immutable source/stiffness-bound history localizes attachment at zero stage
approach and separation at the fixture fold. Each event reports the finite tip
motion and the released total energy. A monotone segment satisfies
`stage work = delta(contact + spring storage) + jump release`; replaying an
accepted detached state cannot release the same energy again. Evaluation is
allocation-free with bounded scalar solves. The original material, geometry
and JKR applicability limits remain in force. The fixture is massless,
quasistatic and linear: frame inertia, nonlinear compliance, servo dynamics
and a dynamic jump solution are not inferred. Unresolvable finite-stiffness
folds refuse rather than silently becoming rigid control. The formulation
follows Deng and Kesari (2019), <https://doi.org/10.1016/j.jmps.2019.07.009>.

`tests/adhesive_fixture.rs` compares the dimensional force/radius equation,
independent full-energy stability Hessian, energy gradients and condensed
tangents, stiffness limits, next-representable values near the two folds,
ports and refusal/history behavior. A complete cycle is compared against an
independent radius-coordinate work integral, including both jump energies.

## Ballistic contact timing and isotropic Coulomb impulses

`ballistic::first_ballistic_contact` finds the earliest closing (or isolated
grazing) zero of a nonpenetrating constant-acceleration scalar gap over an
inclusive finite horizon. It includes an initially departing body's return,
uses a scaled stable quadratic evaluation, and refuses unresolved numeric
ranges. This is an analytic floating-point query, not interval-certified CCD
for arbitrary shapes, motion tubes, or persistent-contact response.

`tangential::impact::solve_isotropic_coulomb_impulse` consumes the existing
`fs-tribo::FrictionLaw`. For equal static/kinetic Coulomb coefficients, a supplied
normal impulse, and positive scalar tangent inverse mass, it returns the exact
two-dimensional maximum-dissipation disk projection. The inverse mass must
include spin, and normal/tangent response must decouple. The normal impulse is
unchanged. Actual endpoint velocity, capacity, nonnegative relative-mode loss,
and independent velocity/impulse-work residuals are checked before acceptance.

Both operations use constant work and storage. G0/G1 tests cover exact impact
and slip responses, spin-aware inertia, grazing/returning trajectories and
missed-endpoint collisions; G3 tests cover tangent rotation and unit scaling.
The mechanics owner supplies geometry, state application, moving-surface work,
and any slip-arrest event splitting. No finite-patch storage, heat partition,
unequal static/kinetic coefficients, general anisotropic or multicontact
friction, calibrated material authority, or certified nonintersection is
claimed by these two primitives.

## Simultaneous non-associated normal/Coulomb point impacts

`frictional_multicontact::solve_frictional_contacts` solves the complete
normal/tangent Delassus coupling in one Newton–Krylov natural-map problem.
Separate nonnegative normal complementarity and tangential disk projection
implement Signorini–Coulomb without the artificial normal dilation of an
associated cone quadratic. Newton restitution and caller-owned nonnegative
normal velocity targets remain explicit, distinct laws. Initial penetration
refuses, and geometrically inactive contacts receive no impulse.

Stable contact identities canonicalize assembly and reduction. Independent
normal primal/dual feasibility, friction-cone, complementarity, sliding and
tangential duality checks follow the nonlinear solve. Actual midpoint contact
work then gates global passivity; a locally valid restitution candidate that
creates assembled energy is refused. All impact loss remains unallocated
mechanical loss, not an automatic thermal or friction-heat source.

This bounded owner admits numerically SPD active rows (at most 21 contacts).
A positive weighted contraction witness for the exact normal/disk map is
required for determined numerical authority. Distinct independently checked
passive roots return `Nonunique`; a checked root without the sufficient
witness returns `UnresolvedDeterminacy`. Rank-deficient manifolds and a complete
Painleve existence classifier remain outside this scope. Configuration
stepping, CCD, contact manifold discovery and a complete Moreau–Jean integrator
are separate owners.

Seven native physical test groups cover analytic oblique coupling, 3D sticking,
explicit target/restitution semantics, shared-body simultaneous response and
permutation/tangent-frame invariance, the existing frictionless oracle,
adversarial global energy injection, independently derived multiple Coulomb
roots, and numerical/geometry/cancellation refusal. The finite-difference inner
forcing floor is 1e-6; the 1e-12 outer target and independent physical gates are
unchanged.

## Straight-line swept-sphere contact timing

`swept_sphere::first_sphere_contact` finds the earliest entering or isolated
grazing contact from actual relative centres and velocities, the sum of the
radii, and a finite inclusive horizon. Common acceleration is admissible
because it cancels from relative motion. Initially overlapping spheres refuse;
touching departing or tangential pairs do not create a new entering event.

The constant-work, allocation-free query uses dimensionless closest-approach
and chord geometry with a rationalized entry root. It avoids squaring physical
distances or speeds, retains a B-to-A normal and relative normal speed, and
refuses unresolvable radius/separation or time ranges. Six G0/G1/G3 tests cover
interior crossings, oblique entry, grazing, misses, initial-contact modes,
near-touch cancellation, rotation, and large/small length scaling.

This is analytic floating-point sphere geometry. It does not certify arbitrary
shapes, accelerating relative trajectories, or general motion tubes; the
mechanics owner still checks the actual event geometry and momentum response.


### Redundant Coulomb contacts with determined velocity

`frictional_multicontact` admits a bounded semidefinite mobility range when a separate Coulomb feedback witness determines its active velocity image. It reuses `fs-la` eigendecomposition, independently checks original-matrix spectral/orthogonality residuals and a separated numerical rank, and never replaces the mobility matrix or adds compliance. All nonlinear, complementarity, cone and energy checks retain the original matrix.

`FrictionalContactDeterminacy::UniqueImpulses` retains the prior independent-row admission. `UniqueVelocity(PsdVelocityWitness)` reports one admissible impulse representative and a unique active velocity image, without claiming unique local pressure, friction heat or independently moving-support work. Friction-weighted normal rows must be observable from the mobility range and the numerical feedback bound must be below one. For an actual mechanical Gram matrix, the unique image implies unique generalized body velocity. General frictionful normal-load indeterminacy remains refused; the witness is a sufficient floating-point admission, not an interval rank/uniqueness certificate or general Painleve classification.

Five new G1/G3/G4 test groups cover finite-height tripod stick/slide with independent force, torque and kinetic-energy balances, an explicit second admissible sticking impulse allocation, exact input permutation, four-point frictionless support, indefinite/ambiguous/unproven-feedback refusals, and cancellation/budget/exact retry. All seven existing core regressions remain passing. The twelve-test native source-selected check does not replace whole-workspace Cargo/DSR verification. Advances `frankensim-ext-contact-nonsmooth-lane-oh0i` without closing its broader Moreau-Jean and classification scope.

### Finite-zone Maugis adhesion

`normal_patch::maugis::SourceMaugisSphereContact` solves the Maugis--Dugdale
cohesive-zone equations (Maugis 1992, DOI 10.1016/0021-9797(92)90285-T;
Leroy and Willert 2021, DOI 10.1002/zamm.201900313, equations 2--8).
Original elastic/interface sources and explicit interaction range select the
constant cohesive traction W/z0. The outer-zone closure is a positive quadratic
at fixed annulus angle; a single bounded scalar inversion selects stable
displacement- or force-controlled contacting equilibria across JKR/DMT limits.

Responses retain signed resultant force, true contact area, the Hertz repulsive
component, analytic tangent, and the finite radial pressure distribution.
The cohesive annulus contributes attraction but no solid junction-shear area.
Source identity and synthetic authority survive unchanged. Geometric smallness,
half-space/layer depth, nominal strain, cohesive stress and compressive pressure
limits are checked on the authored envelope. Smooth ports preserve action/reaction.

A zero-area contacting endpoint can retain attraction; a displacement fold is
also not a detached state. Endpoint smooth tangents refuse. Explicit-budget work
quadrature integrates the original force along the stable contact branch and
returns an error estimate, not an interval certificate. Distinct inputs that
collapse to one representable branch coordinate refuse zero-work inference.
No detached-reference absolute potential, deformed noncontact continuation,
attachment/snap dynamics, roughness, plasticity or universal material calibration
is claimed.

## Source-preserving adhesive held friction

`SourceAdhesiveShear::bind_jkr`, `bind_dmt`, and `bind_maugis` retain the original
ordered normal state, original shear claims and receipts, and actual solid
contact area in `SourceAdhesiveTangential`. Response/source identity mismatch
refuses. The Maugis cohesive annulus and DMT noncontact attraction contribute
no solid junction shear. The explicit source query remains frozen; the held
input must preserve its common contact temperature and repulsive load. A
nominal-pressure source refresh is not inferred. Original synthetic/caller
authority is preserved. The adapter implements the existing tangential trial
trait and neither owns normal-area evolution nor claims mixed-mode peeling.

## Conservative DMT work in native dynamics

`NormalLawSource::DmtSphere` retains the original immutable DMT owner and its
absolute detached-reference potential. Signed negative approach reaches the
finite-range attractive branch. `energy_and_force` remains continuous at
touch and range cutoff, while the smooth-tangent `potential` interface still
refuses those corners. The low-Tabor envelope and source identity remain
authoritative.

`reversible_work_between` integrates the original piecewise force analytically
without subtracting nearly equal endpoint energies. The allocation-free scalar
midpoint numerator uses a centered Hertz series and exact quadratic-hinge
remainders for the same potential. No new damping, heat, snap hysteresis,
Maugis detached reference, or changing-source energy law is inferred.
