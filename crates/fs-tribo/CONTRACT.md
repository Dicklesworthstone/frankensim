# CONTRACT: fs-tribo

> Status: `frankensim-ext-tribo-dry-baseline-tgbj` baseline. This is not an
> admitted material-data result or an Euler-disc validation result.

## Purpose and layer

`fs-tribo` is an L3 dry-contact constitutive leaf. It provides typed dry
friction, elastic Hertz reference forms, scalar rolling/contour loss laws,
explicit dissipated-work accounting, a bounded flash-temperature candidate,
and caller-owned Archard wear state. It does not depend on contact-solver
internals: the consuming solver owns contact complementarity, tangent
reactions, and residual embedding.

## Public types and semantics

- `InterfaceSystemRef` names ordered surfaces, history, medium, and a nonblank
  caller `source_id`. `InputAuthority` is a caller-declared ceiling
  (`CallerDeclared`, `SyntheticFixture`, or `Estimated`), never an admission or
  verification receipt. Every constitutive response returns this provenance.
- `DryInterfaceSystemCard` is a dependency-independent, immutable
  InterfaceSystemCard-style input. It combines that ordered interface/history
  identity with one `FrictionLaw` and closed SI temperature, nominal-pressure,
  and slip-speed ranges. Its `query` refuses all out-of-domain states and a
  mismatch between the declared slip-speed magnitude and the supplied tangent
  vector. It is not an alias for, or a substitute for, the pending seed-data
  `fs-matdb::InterfaceSystemCard` query/receipt path.
- `ContactFrame` normalizes a finite contact normal. `TangentialSlip` accepts
  only a finite velocity with no material normal component; it refuses rather
  than silently projects a closing/separating velocity into friction.
- `FrictionLaw` provides Coulomb, velocity-dependent, and Stribeck rungs. At
  zero slip it reports static capacity only; it never invents a stick reaction.
- `HertzSpherePlane` and `HertzCylinderPlane` are G1 elastic closed forms,
  expressed in SI. They use caller-supplied radius/modulus and retain the
  caller ceiling; a mathematical closed form does not validate its inputs.
- `ConstantRollingMoment` uses `M = -sign(omega) N a` (N m) and
  `P = N a |omega|` (W). `ConstantContourForce` uses
  `F = -sign(v) F_c` (N) and `P = F_c |v|` (W). They are separate mechanisms;
  callers must not combine them without independent source support.
- `HeatPartition`, `DissipationStep`, `WorkLedger`, and `WearState` have
  private invariant fields and checked getters. `WorkLedger` retains total
  dissipated work and every explicit heat/other channel, and commits all four
  candidate totals only after their closure checks pass. `ArchardLaw` computes
  `dV = k F ds/H` in m3 and commits only a validated candidate state.
- `flash_temperature_candidate` is a typed model-form branch. With complete
  caller thermal conductivity/diffusivity data for every nonzero heat share,
  it reports the uniform-flux semi-infinite candidate
  `2 q'' sqrt(alpha * t / pi) / k`, with `t = traverse_length/slip_speed`.
  Missing properties return a typed `Unknown`, rather than guessed material
  data or a partial temperature result.
- `surface_excitation` admits named, uniformly sampled finite or periodic
  surface-height traces and filters them by the exact one-dimensional marginal
  of an elliptic Hertz pressure footprint. It analytically integrates the
  piecewise-linear trace and its path derivative, then applies only the
  explicitly bounded contact tangent `delta_F = k_n delta_h`.
  `SelfAffinePeriodicProfileSpectrum` can construct a deterministic,
  random-phase, band-limited one-dimensional profile whose spatial PSD follows
  `Phi(q) proportional to q^-(2 H + 1)` and whose Fourier coefficients are
  normalized to an explicit RMS height. Its seed selects spatial phase only;
  it creates no temporal/audio spectrum or apparatus-specific frequency.

### Partial-slip return-map rung

`partial_slip` is a solver-independent, finite-patch Cattaneo--Mindlin-style
*return-map* rung. `NormalPatchView`, `PartialSlipInterface`, and
`GeneralizedWorkOwnership` retain the caller's named normal-patch, ordered
interface/history, and work-owner inputs; their authority remains explicitly
caller-declared. `PartialSlipLaw::advance` has no `fs-contact` dependency and
cannot solve normal contact, evolve patch geometry, or admit a material card.
It yields a reversible tangential/torsional spring core plus a lumped
microslip remainder. Its `PartialSlip` state is not a resolved slipping-area
fraction, pressure field, or traction field. Checkpoints bind all law,
patch, interface, and state data for deterministic replay; equality of those
inputs is not external physical validation. Rolling-deformation loss remains
a separate, zero channel in this rung.

### Rolling and contour-deformation loss rungs

`rolling_loss` provides three immutable, separately sourced candidates:
the `Leine-style Coulomb contour-force rung` (`F_c = -sign(v_c)c_cN`), a
viscous contour alternative (`F_v = -c_vv_c`, with `c_v` in N s/m), and a
finite-patch hysteretic rolling moment (`M_h = -sign(omega_r)N l_h eta_h`).
The last card declares its loss length `l_h` [m] and bounded loss factor
`eta_h`; it is a model-form candidate, not a resolved deformation field.
Every response retains a named patch/load/curvature receipt, ordered
interface/history, temperature/frequency range, source-card authority, and a
dedicated generalized-work key. Principal curvatures are retained directly;
an equivalent-radius substitute requires explicit caller approximation
authority. The rungs have no blend/preference operation, no numerical error
certificate, and no torsional microslip or gross-slide contribution. Pure spin
therefore produces no torque in this slice; a separately sourced torsional
material-hysteresis rung would need a disjoint work channel.

## Invariants

- All public construction and evaluation paths refuse missing identity,
  non-dry medium, negative/non-finite input, and every non-finite derived
  candidate from finite extreme inputs.
- Tangential traction and scalar resistance oppose their declared relative
  rates; output dissipated power/work is finite and non-negative.
- Heat shares are finite, non-negative, and sum to one. A dissipation step's
  channels are finite, non-negative, and close to its total before a ledger can
  mutate. All ledger channels and wear candidate totals are checked before any
  assignment.
- Card state ranges are closed and finite; absolute temperature has a positive
  lower bound. A direct `FrictionLaw` remains a low-level constitutive rung,
  while state-domain enforcement occurs at `DryInterfaceSystemCard::query`.
- A flash result is either complete for both receiving surfaces or an explicit
  `Unknown`; it never silently assigns missing thermal properties.
- Partial-slip refuses nonfinite or nonpositive patch and law inputs, malformed
  tangent frames/kinematics/work ownership, and checkpoint mismatches. It
  computes only its declared lumped return-map state and never creates a
  contact-solver or material-admission receipt.
- Each rolling-loss card refuses non-dry interfaces, malformed patch/source
  receipts, out-of-range temperature/frequency, nonfinite derived work, and a
  wrong work channel. A rolling key conservatively refuses sharing the same
  patch/interval pair with partial-slip because that module does not expose its
  coordinate labels for a finer disjointness proof.
- Rolling/contour candidates are passive: endpoint body power is non-positive,
  heat is non-negative, reversible storage is zero, and the retained work
  identity reconstructs heat over the caller interval.
- No mutable material table, hidden state, or unordered material-pair lookup
  exists. Ordered interface and history identities remain explicit.
- Surface excitation refuses missing texture/source identity, uncovered finite
  footprints, excessive periodic repetition/work, non-finite kinematics, a
  height outside the caller's tangent-linearization fraction, and a
  perturbation that would open the unilateral contact. Both surface sources
  and authority ceilings remain separate in the receipt.
- `PeriodicHarmonicSurface` admits either a measured Fourier decomposition or
  an explicitly declared hypothetical periodic material track. Signed
  sine/cosine height coefficients are keyed by positive integer spatial cycle,
  canonicalized by cycle count, and sampled only with at least eight points per
  shortest period. Duplicate cycles and underresolved grids refuse. An empty
  spectrum is an exact smooth track. No material name, acoustic frequency,
  stochastic realization, seam correction, or authority upgrade is inferred.
- `SelfAffinePeriodicProfileSpectrum` refuses nonpositive/nonfinite RMS height,
  Hurst exponents outside `0 < H < 1`, inverted or zero cycle bands, and
  harmonic bands whose shortest wavelength cannot satisfy the eight-samples-
  per-period rule within the one-million-sample trace ceiling, as well as
  unbounded harmonic counts. Every retained cycle owns a counter-based phase
  stream keyed by the explicit seed and cycle, so extending the band preserves
  shared-cycle phases. Phase changes preserve the declared Fourier power, and
  the continuous periodic profile retains the requested RMS height to
  floating-point roundoff.

## Error model

`TriboError` is a total refusal surface for missing identity, non-dry media,
invalid/overflowing physical input, non-finite vectors, normal slip, malformed
or out-of-domain applicability ranges, invalid partitions, and forged/invalid
dissipation states. A refusal makes no partial state change to `WorkLedger` or
`WearState`. Missing flash thermal data is a typed `Unknown`, not an error or a
candidate.

## Determinism class

Pure scalar operations are deterministic for the same input and ISA. Norms use
the platform `hypot` sequence to avoid avoidable intermediate overflow. No
cross-ISA bit-stability claim is made.

## Cancellation behavior

This bounded scalar leaf has no asynchronous work and no cancellation scope.
Callers scheduling many contact points own cancellation between calls.

## Unsafe boundary

None. The crate forbids unsafe code.

## Feature flags

None.

## Conformance tests

Inline and external G0/G1/G3 tests cover identity/dry/tangent refusals;
ordered-card applicability and slip-speed consistency; analytic rigid-block
equilibrium; independent numerical Hertz values and force-pressure-radius
cross-relations; reversal/scaling/provenance; resistance sign; partition and
ledger-channel closure; typed flash-data insufficiency and power scaling;
forged negative/non-finite work rejection; and rollback/deterministic replay
for work and wear. Surface-excitation tests additionally cover periodic-spectrum
order invariance, exact sample-grid realization, duplicate-cycle refusal, and
the eight-samples-per-shortest-period gate. They also cover deterministic
self-affine replay, exact RMS normalization, the declared PSD slope, seed-only
phase changes, and malformed statistic/band refusal. Test coefficients are
explicitly synthetic fixtures.
The separately owned partial-slip tests cover its scalar return-map admission,
reversal, checkpoint replay, and receipt-mutation refusals.
Rolling-loss tests cover zero load/rate/factor, force and moment reversal,
load scaling and SI dimensions, pure roll/spin/mixed motion, passivity with an
independent work integral, temperature/frequency refusal, work-key overlap,
source-card/patch replay mutation, and a no-target-field source scan.

## No-claim boundaries

This crate does not mint material admission, calibration, the upstream
`fs-matdb::InterfaceSystemCard` query receipt, a flash-temperature error bound
or thermal-port solution, a roughness measurement or statistical roughness
model, adhesion, plasticity, finite-patch partial
slip as a resolved contact field, lubrication/EHL, contact geometry evolution,
wear geometry updates, stop time, Euler-disc ranking, one-millimetre optimum,
sound pressure, loudness, or experimental/video correspondence. The surface
trace leaf consumes a caller-supplied one-dimensional profile and a local
contact tangent; it does not replace a two-dimensional roughness map, nonlinear
re-solution of contact, structural response, or acoustic radiation. Its flash candidate is a declared
semi-infinite uniform-flux model-form estimate, not a temperature measurement
or validation result. A caller may record an estimated or synthetic result,
but must not promote it based on this crate alone.

The self-affine constructor is likewise a geometry realization, not evidence
that a specimen is self-affine or that its RMS height, Hurst exponent,
wavelength cutoffs, or phase realization are correct. Those inputs require
profilometry or an explicitly weak estimate. The one-dimensional profile does
not establish an isotropic two-dimensional areal PSD or resolve asperity-scale
partial contact.

The rolling and contour rungs likewise do not establish deformation-loss
dominance, coupon calibration validity, normal-contact admission, an error
bound, or any apparatus-specific selection between their rival model forms.
