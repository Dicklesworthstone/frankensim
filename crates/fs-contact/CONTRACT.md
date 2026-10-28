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
  sweep-and-prune on the axis with the fewest closed-interval overlaps
  (x/y/z ties, then stable body indices); output pairs sorted; report
  carries checked/pruned counts and the worst
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

### Source-backed Hertz--Mindlin tangential contact

`interface_binding::mindlin::SourceMindlinContact` implements the shared
`fs_tribo::TangentialLaw` using an original `SourcePointContact`, both admitted
bulk elastic states and an explicitly selected scalar property from the same
ordered interface card. It obtains the circular patch radius and normal load
from the original Hertz evaluator and derives each shear modulus from its
source E and nu. The coefficient is never inferred from a material name or
copied from the normal model. Original property-use receipts, uncertainty,
authority and the normal operating envelope remain accessible.
`SourceMindlinState` additionally retains the exact original normal-model and
resolved-friction identities. Numerically equal coefficients from a different
ordered interface, specimen history, query or selected claim cannot import
retained shear. Its native constitutive state is exposed read-only.

The caller explicitly selects the classical decoupled, isotropic, constant-load
Mindlin approximation with one local Coulomb coefficient. It does not combine
unequal static/kinetic coefficients. Missing properties, wrong units, rate- or
history-restricted sources, noncircular geometry, nonelastic normal models and
zero-area contact refuse. Stepping retains the original source temperature and
normal load; a nominal-pressure refresh, transverse motion, foreign mechanical
or source checkpoint or excessive step duration refuses before publishing history.
Thermal coefficient changes and changing normal load need a separate physical
history-transfer model; a newly unloaded contact is not such a transfer.

The generic modal tangential consumer can use this source owner directly.
`tests/source_mindlin.rs` covers the independent `8 a G*` stiffness and Hertz
backbone, distinct material/interface substitutions, original receipt retention,
frictionless data, domain refusals and exact retry from retained preload. These
synthetic G1/G3 checks establish numerical/source integration, not calibration
against a measured specimen.

### Source-backed circular Hertz sliding and spin

`interface_binding::mindlin::sliding::SourceHertzSlidingPatch` reuses the
original `SourceMindlinContact` admission, including its explicit ordered
interface coefficient, bulk materials, fixed external conditions and exact
normal response. The native Hertz law integrates local Coulomb traction under
the original pressure distribution for simultaneous two-axis translation and
normal-axis spin. Mixed motion reduces both the resultant sliding force and
spin moment relative to independently saturated limits. The response retains
the original normal-model and resolved-friction identities; complete source
receipts, uncertainty and authority remain available through the binding.

Construction checks a borrowed, source-matching elastic checkpoint with exactly
zero retained shear energy and force. It leaves that checkpoint unchanged and
does not transfer or erase loaded Mindlin history. This selects a distinct
zero-storage gross-sliding model. It does not solve partial slip, static
friction, rolling deformation or a transition between these constitutive laws.
At complete rest its zero wrench is one admissible set-valued selection.

The evaluator requires the original normal load and temperature, rejects a
pressure refresh and enforces the original step envelope. Prescribed slip and
spin are held over the admitted interval: heat is exactly the pressure-integrated
sliding power times that duration. Bodies, changing kinematics and heat
partition remain owned by their mechanical and thermal consumers.
`tests/source_mindlin/sliding.rs` checks actual pure/mixed force, moment and heat,
material/interface/geometry substitutions, complete source retention and
loaded-history, foreign-source and operating-domain refusals.

### Source-backed temperature-dependent Coulomb friction

`interface_binding::thermal_friction::SourceThermalCoulomb` implements the
shared `TangentialLaw` with explicit static and kinetic coefficient requirements
from one original ordered dry interface card. Every candidate evaluates those
same selected claims at its supplied contact temperature and, when explicitly
mapped, nominal pressure. Other admitted external coordinates remain fixed.
Property-use receipts, uncertainty, original selection and authority remain
available; no fitted material-name coefficients or silent source substitution
are introduced. Static pressure needs its quantity semantics, and rate, slip,
area, strain and history restrictions refuse even for constant-valued claims.

The native set-valued Coulomb solver owns traction, sticking, sliding and body
compliance. Temperature and pressure refreshes retain its accumulated slip.
There is no reversible contact potential: coefficient evolution changes the
dissipative law without introducing unreported elastic parameter work. Accepted
states bind the original card, fixed query, selected sources, coefficient roles,
axis mapping and authority. Foreign history, missing/out-of-support conditions,
wrong units, negative coefficients or kinetic friction above static friction
refuse the whole candidate. The mechanical and thermal consumers own normal
load/pressure generation, temperature evolution and heat partition.
`tests/source_thermal_friction.rs` checks the actual temperature/pressure force
and heat response, static breakaway, pinned source support and atomic retry.

### Source-backed Hertz--Mindlin normal loading and recontact

`interface_binding::mindlin::loading::SourceHertzMindlinContact` reuses an
original admitted `SourceMindlinContact` and the shared analytic changing-load
owner. Its initial positive indentation and normal storage are original Hertz
initial conditions; no new normal work is inferred at construction. Shear
storage remains separate from that normal potential in `TangentialLaw`.

Each `normal_step` takes an explicit target indentation and duration. The
original `SourcePointContact` admits the whole monotonic constant-rate segment's
operating box and evaluates its endpoint force, pressure and footprint. Native
history transfer must agree with the original Hertz geometry/energy within
arithmetic roundoff; source-domain bounds are never relaxed. Normal work and
normal energy come from the original source potential. Shear memory/heat come
from the existing analytic Hertz--Mindlin owner, including unloaded recontact.
No private modulus, friction default, coefficient refresh or memory reset is
introduced. Both normal and friction source identities bind every checkpoint.

All selected bulk, normal-interface and friction properties must retain static
source support. Complete query coordinates use the canonical material role
policy, with pressure excluded: a pressure-dependent value, pressure-only
validity bound, or even an extra pressure query cannot silently survive a change
in actual load and footprint. This also refuses unowned force, area, stress,
rate and unknown roles without guessing from axis names. Fixed temperature,
composition, moisture and concentration retain their declared meanings.

The path is explicitly normal-at-fixed-shear followed by held-load shear;
simultaneous body motion needs a separate coupled trajectory owner and temporal
refinement. Temperature stays at the original source condition. Five source
integration groups cover original Hertz work/authority, growth and release,
separation/recontact, typed condition refusal, material/interface substitutions
and exact retry after operating-envelope or foreign-source refusal.

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

`normal_patch::dmt::fixture::DmtSpringFixture` connects that same source-bound
finite-range law to a massless loading spring of stiffness K. With
S=2*pi*R*W/z0, stability requires K+F'(delta)>0. For K>S the attractive-gap
branch is stable and approach/retraction is reversible. For K<S that branch
is unstable: immutable history localizes snap-in at stage D=-z0 and snap-off
at D=-2*pi*R*W/K. The contact's geometrical regime remains separate from the
sign of the total force; an attractive gap never supplies solid shear area.
K=S is neutral with nonunique gap equilibria and refuses explicitly.

The original normal sources, authority and indentation envelope remain in
force. A bounded allocation-free scalar solve supplies contacting equilibria;
gap/separated equilibria and stability thresholds are analytic. Responses retain
spring storage, total potential, force-balance residual and the condensed stage
tangent. Geometric corners have no smooth tangent. Branch jumps report the
finite tip motion and release of contact-plus-spring potential at fixed stage.
The external transient/thermal owner decides where that release goes. Smooth
stage work uses the original exact-piece DMT work and a factored spring-energy
difference; only accepting a successor token commits a release. Total work
closes storage change plus jump release. Machine inertia, tip dynamics,
viscoelasticity and experimental validation are not introduced.

G1/G3 `tests/adhesive_dmt_fixture.rs` checks the original dimensional equilibrium,
the energy curvature of stable and unstable tail solutions, both event
thresholds, independent integrated cycle work, reversible stiff-machine cycles,
potential gradients/tangents, dimensional similarity and source/history refusals.
The fixture is a direct energy derivation for the existing DMT-M/Dugdale model;
the machine-stiffness motivation is Deng and Kesari (2019),
<https://doi.org/10.1016/j.jmps.2019.07.009>.

`normal_patch::dmt::dynamics::DmtInertialContact` resolves finite-mass normal
motion on the same original potential. With explicit mass m, loading spring
K, fixed-reference dashpot c and held external force Q, it advances
`m x'' = Q + K(D-x) - F_DMT(x) - c x'`. K=0 admits load control, c=0 retains
conservative motion, and the prescribed stage D has its own work port.
The source's exact-piece reversible work supplies the discrete contact force;
inertia, spring and damping use work-conjugate midpoint quantities. Snap-in
and pull-off therefore pass through finite-gap motion and kinetic energy,
without equilibrium jumps or automatic conversion of released potential to
heat. Only the explicit dashpot reports dissipation. Endpoint `DmtResponse`
retains the actual evolving solid area, repulsive load, original identities and
authority for tangential consumers. Gap attraction still supplies no shear area.

Each allocation-free proposal has explicit time, iteration, displacement,
velocity, force and source-range limits. `4m/dt^2 + 2c/dt + K > S` is required:
the original force slope is bounded below by -S, so this admits a unique
discrete-gradient root even through the two corners. A safeguarded scalar
solve refuses an out-of-range or unresolved step. Immutable history retains
geometry, velocity, stage, clock and physical model identity; rejected trials
or discarded successors change no accepted state. Separate stage work, held
load work, contact work and damping loss close kinetic/contact/spring storage,
with explicit force and energy residual gates. Mass may be a declared reduced
mass; this scalar owner does not invent the rest of a multibody assembly.

G1/G3 `tests/adhesive_dmt_dynamics.rs` checks detached impulse/work, conservative
attachment and detachment with finite kinetic release, damped settlement to
the original spring equilibrium, separate driven work ports, second-order
convergence to an independent unstable-gap trajectory, dimensional similarity
and transactional refusals. The scalar discrete-gradient method follows
Gonzalez (1996), *J. Nonlinear Science* 6, 449-467,
<https://web.ma.utexas.edu/users/og/PUBLICATIONS/paper_DisHamSys.pdf>.
Energy closure is not a temporal error certificate. Source conditions remain
fixed; adhesion kinetics, elastic-wave propagation, viscoelastic contact and
shear-driven peeling need their own constitutive and mechanical owners.

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

This bounded owner retains a cap of 21 active contacts. Numerically SPD rows
use a positive weighted contraction witness for unique impulses; admitted
semidefinite rows use the separate velocity-image witness described below.
Distinct independently checked passive roots return `Nonunique`; a checked
root without a sufficient witness returns `UnresolvedDeterminacy`. General
unresolved rank-deficient manifolds and a complete Painleve existence
classifier remain outside this scope. Configuration
stepping, CCD, contact manifold discovery and a complete Moreau–Jean integrator
are separate owners.

Seven native physical test groups cover analytic oblique coupling, 3D sticking,
explicit target/restitution semantics, shared-body simultaneous response and
permutation/tangent-frame invariance, the existing frictionless oracle,
adversarial global energy injection, independently derived multiple Coulomb
roots, and numerical/geometry/cancellation refusal. The exact semismooth action
and redundant-contact right preconditioner below supersede directional finite
differences. The 1e-12 outer target and independent physical gates are unchanged.

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

`frictional_multicontact` separates numerical PSD/rank admission from physical Coulomb determinacy. It reuses original `fs-la` eigendecomposition, checks original-matrix spectral and orthogonality residuals and a separated numerical rank, and retains the unmodified mobility for every nonlinear, cone, complementarity and work check. No compliance or dropped physical row is introduced.

`UniqueImpulses` retains the independent-row contraction result. `UniqueVelocity(PsdVelocityWitness)` has an explicit `PsdVelocityCriterion`: `NormalRangeObservable` uses the earlier friction-weighted normal-range feedback bound; `CandidateFrictionDirection` uses subgradients of an audited candidate. For any second Coulomb root, norm monotonicity gives `E + delta_p_n dot C delta_v_t <= 0`, with `E = delta_p^T W delta_p`. When this bilinear feedback annihilates the admitted mobility nullspace and its symmetric energy form is positive on the range, all roots have the same `W*p`, even when individual normal reactions vary. The report retains actual weighted directions, normalized null coupling, energy margin, spectral/orthogonality residuals and traction/subgradient residuals. Near-zero slip uses the declared sticking resolution; sliding directions and null coupling are checked relative to friction magnitude, including tiny coefficients.

A bounded SVD may propose a coherent traction representative, but it must preserve the ORIGINAL mobility image to roundoff and pass the complete original law/work audit before use. Neither representative nor minimum norm assigns a physical pressure distribution, per-contact heat, or independently moving-support work. Distinct checked velocity roots report `Nonunique`; two different impulse gauges alone do not. Unproved determinacy still refuses. These are sufficient floating-point admissions, not interval uniqueness or general Painleve classification.

Sixteen focused integration groups pass: the twelve earlier groups plus physical four-foot cube stick/threshold/oblique sliding with two independently constructed pressure allocations, independent tangent-frame covariance and permutation, an exact heterogeneous-friction example with two passive motions, and tiny heterogeneous-friction refusal. Three original analytic-action/preconditioner tests also pass. Native body motion follows for an actual mechanical Gram matrix; whole-workspace Cargo/DSR verification and broader Moreau-Jean scope remain separate.


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
input must preserve its common contact temperature and current repulsive load. A
nominal-pressure source refresh is not inferred. Original synthetic/caller
authority is preserved. The adapter implements the existing tangential trial
trait and neither owns normal-area evolution nor claims mixed-mode peeling.

`SourceAdhesiveContactUpdate::{from_jkr,from_dmt,from_maugis}` carries the actual
normal response's solid area, repulsive force and original normal identity in a
small allocation-free value. `tangential_for_update` prepares a trial leaf;
`update_contact` publishes only its geometry after checking that exact identity.
The ordered source owners, selected strengths, query and temperature stay fixed,
and their allocations are neither replaced nor copied per update. Source drift
refuses before mutation. This accepts geometry resolved by the external normal
owner; its dynamics, normal work and any junction separation event remain that
owner's responsibility. Source regressions exercise JKR tensile-area evolution,
DMT contact/gap/reattachment, fixed source allocation/temperature, wrong-source
rollback and exclusion of the Maugis cohesive annulus from solid shear.

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


### Exact semismooth action and redundant-contact Newton corrections

The coupled normal/Coulomb natural map now supplies an analytic semismooth Jacobian action through the existing `fs-solver` nonlinear interface. Normal positive-part and circular tangential-projection branches are fixed at each exact Newton iterate. Cached local blocks leave one original mobility-vector product per Krylov action; they do not perturb accepted contact state or turn a mode boundary into a classical adjoint. The unused `difference_step_impulse` configuration field is removed.

For admitted semidefinite mobility only, the owner caches a numerical pseudoinverse of the current Newton Jacobian as a right preconditioner using the existing `fs-la::factor::svd_jacobi`. Its cutoff is proportional to machine precision and matrix dimension. This preconditions the numerical correction only: the original contact mobility, Jacobian action, nonlinear residual, strict outer convergence target and independent physical gates remain unchanged. An unrepresentable right-hand-side component remains visible to the original Krylov residual. No artificial compliance, physical impulse clipping, or acceptance of a stalled solve is introduced.

The evaluation budget counts residual, Jacobian and right-preconditioner calls. Ordinary actions are quadratic; a new semidefinite Newton point may require one bounded cubic factorization, cached across its inner solve, with cancellation before and after. The 63-row cap remains in force.

All twelve physical contact regressions pass against this implementation. Three additional native unit tests independently compare smooth-branch derivatives with central differences, check exact-boundary action linearity and unchanged iterate, and demonstrate that an unfittable right-hand-side component is not hidden by the range preconditioner. Native rigid sustained-support regression and timestep evidence are recorded by the `fs-mbd` consumer. Verification is bounded native source-selected execution, not a full workspace/DSR gate.

## Accelerated relative sphere paths

`swept_sphere::accelerated::sphere_path` admits original relative position, velocity, constant acceleration, radius sum and inclusive horizon. It returns first entering contact and minimum gap on the entire unimpeded quadratic path. The derivative of squared distance is cubic; its quadratic turning points separate every monotone interval. At most three interior extrema and one entering root use bounded bisection (128 iterations each), fixed storage and a nonpenetrating entry bracket. Distances determine gaps without subtracting large squared distances. Exactly zero relative acceleration retains the original straight-line fast path.

Isolated tangent touches are not entering events. A touching stationary pair accelerating inward is reported for the response owner to resolve or refuse; no reaction is invented. Initial overlap, non-finite data and unresolvable numerical ranges refuse. Four G1/G3 groups cover hidden crossings between separated endpoints, accelerated departure/return, three distinct distance extrema, independently manufactured 3-D entry, grazing versus inward support and length/time scaling.

This is floating analytic geometry, not an interval enclosure or arbitrary-shape CCD certificate. A negative minimum after entry describes the unimpeded path; mechanics must stop before entry and independently audit its actual complete segment and witnesses.

## Shared bounded spatial broad phase

`broad_phase` owns one deterministic closed-AABB sweep shared by certified
motion and numerical trajectory consumers. Inputs carry unique stable IDs;
results retain original input indices and sort canonical ID pairs. Touching
boxes overlap. Invalid or nonfinite bounds and duplicate IDs refuse before
enumeration. The input cap is 65,536 boxes.

The sweep counts exact one-dimensional interval overlaps on each coordinate
axis and selects the least populated axis, breaking ties x/y/z. The count
preflights the explicit comparison budget before pair enumeration. Candidate
storage has a separate budget: exhaustion returns a typed error with a bounded
prefix, never a successful incomplete pair set. Cancellation is checked during
validation, counting and enumeration and around bounded sorts. Setup is
O(n log n); enumeration is O(n + k) for the selected axis, followed by sorting
the retained pairs. These are operation-count bounds, not hardware timing
claims. The certified wrapper keeps its complete unresolved-pair error list.

`swept_sphere_aabb` encloses a straight centre segment with an explicit sphere
radius and nonnegative padding, rounding finite bounds outward. It rejects
nonrepresentable enclosures. The caller owns the assertion that its motion and
shape lie inside the supplied boxes; a spatial overlap result alone does not
certify a curved trajectory, collision time, contact manifold or response.
The original certified wrapper still obtains each box and defect from its
`CertifiedMotorTube`; extraction of the shared sweep does not change that
authority boundary.

Seven focused native G0/G3/G4 groups cover 240 independent all-pairs scene
comparisons, permutation/replay, closed touching and signed zero, work/output
budget refusal, outward sphere bounds, malformed inputs/body ceilings, and
cancellation. A 4,096-body vertical fixture performs zero selected-axis pair
comparisons instead of the 8,386,560 required by a fixed x sweep. This is
operation-count evidence. The focused run compiles the production module with
the original geometry owners; it does not establish a full runtime-dependent
certified-wrapper or Cargo/DSR workspace gate.

## Conservative mixed-mode JKR peeling

`normal_patch::adhesive::mixed_mode::SourceJkrMixedMode` selects a circular, ideally brittle, conservative full-stick Griffith extension of the original sourced JKR sphere contact. It retains the original ordered interface identity, query, elastic states, work of adhesion, radius, interaction range and applicability envelope; tangential elastic stiffness uses the shared Mindlin half-space compliance. One potential determines actual contact area, normal and vector tangential efforts, reversible normal/shear/surface storage, energy-release rates and reciprocal coupled stiffness. Force control and displacement control expose their distinct stable branches and marginal peeling limits; bounded solves refuse missing equilibrium, invalid inputs, source-envelope violations and unresolved numerical equilibrium. No slip rule, empirical mode-mixity dissipation, post-instability trajectory or snap heat is inferred. The ideally brittle model can overpredict experimental area reduction and is not a universal adhesion calibration.

`SourceAdhesiveContactUpdate::from_mixed_mode` and `SourceAdhesiveShear::bind_mixed_mode` transfer the actual reduced solid area and original repulsive load into the retained source-bound junction-capacity owner. The caller must recompute the mixed-mode solution consistently with actual tangential loading and separately resolve any slip/peel transition; geometry exchange does not introduce a second elastic shear storage owner. Force control is closed form and displacement control uses two monotone solves bounded to 64 iterations each, without evaluation allocation. Seven original-source physical groups check independent Griffith equations, area reduction, both folds, energy gradients, full stiffness/compliance, vector path work, rotation and scaling; the sourced shear integration and original Mindlin regressions remain passing.

## Fixed-volume planar capillary adhesion

`capillary_bridge::SourcePlanarCapillaryBridge` retains the original liquid and both ordered wetting source states at one frozen physical query. Its declared geometry is a symmetric, translationally invariant liquid bridge between parallel, identically wetting rigid plates, with fixed liquid volume and acute equilibrium angles. Exact circular Young–Laplace profiles supply surface areas, pressure, attractive force, and the derivative of the same Young–Dupré free energy. The explicit topology rule continues a connected profile to zero neck, then retains all liquid in two equal-volume caps; closing caps reform a bridge when they touch. This geometric rule is not a certified physical stability or dynamic pinch-off threshold.

Monotone-gap proposals integrate connected-branch work analytically and return its work-conjugate mean normal force, including event-crossing intervals. Nonnegative topology releases remain separate from mechanical work and are not automatically heat. Source, geometry, quasistatic, energy/work, or cancellation refusal leaves accepted state unchanged. The existing `PointNormalPort` carries equal/opposite endpoint resultants; event work remains with the accepted `BridgeStep`. Evaluations and proposals allocate no memory. Six original-source regression groups independently check geometry, volume, pressure/line force, energy derivatives, event work, scaling, source retention, and refusal. No finite-width end effects, three-dimensional breakup, unequal liquid partition, evaporation, cavitation, dynamic wetting, measured calibration, or full-workspace validation is claimed.

## Capillary mechanical load line

`capillary_bridge::load_line::PlanarBridgeLoadLine` borrows the original sourced bridge and the caller's complete projected held-force compliance. It solves `gap = free_gap + compliance * mean_force` with the original analytic connected-path force, including rupture and coalescence. `BridgeStep` now provides the exact derivative of that mean force with respect to the endpoint gap; it is distinct from the endpoint-force tangent and stays factored at short connected intervals. The solver allocates nothing and uses bounded damped Newton continuation from accepted geometry, with original source, rate, geometry and energy admission on every proposal. A caller-selected positive compatibility-slope margin refuses unresolved attractive branches for timestep refinement; this local solve does not certify global root uniqueness or trajectory accuracy.

No proposal publishes accepted history. Native mechanics must independently reevaluate actual endpoint motion and the complete work account. Four focused source tests extend the six original bridge groups: independent mean-force derivatives on each topology branch; 700 exact held-force two-body steps with conserved momentum and kinetic/surface/release energy through rupture; prescribed and compliant coalescence; and source, range, budget, unstable-slope and cancellation refusal/replay.
