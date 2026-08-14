# CONTRACT: fs-couple

Multiphysics composition through a versioned port-thermodynamic vocabulary: a
lossless interface relation, explicit storage/dissipation/source primitives,
a complete stream-flux bundle, and evidence-bound closed-window physical
accounting without false passivity.

## Purpose and layer

Layer L3 (multiphysics coupling). Dependency-light: `fs-qty` supplies the
canonical six-base `Dims`/power dimension and `fs-iface` supplies neutral field
function-space roles. No domain solver is a dependency.

## Public types and semantics

- `PORT_SCHEMA_VERSION = 2`; `StableId` admits a canonical transport-safe ID.
- `PortValueShape` distinguishes scalar, non-empty vector/tensor, and field
  values; field shapes carry separate neutral `fs-iface::SpaceType` roles for
  effort and flow. `PowerPairing` distinguishes scalar product, Euclidean dot
  product, and field duality with explicit integration-measure dimensions plus
  `FieldMeasureSide` (effort-density or flow-density).
- `PortKind` covers mechanical translation, rotational torque/angular velocity,
  fluid pressure/volume flow, thermal temperature/entropy flow, electrical
  voltage/current, magnetic mmf/flux rate, and chemical electrochemical
  potential/amount flow. Each kind owns canonical generalized effort/flow
  dimensions; semantic kind identity remains distinct even when dimensions
  coincide (for example torque and energy).
- `CoordinateBinding` makes basis, frame, and positive orientation explicit;
  `PortTimestamp` carries a named logical clock and tick;
  `ConservationRole` is canonicalized into a sorted, duplicate-free list.
- `PortSchema::try_new` binds stable ID, legacy/current `PortKind`, six-base
  effort/flow dimensions, value/field shape, coordinates, power pairing,
  timestamp, and conservation roles. V2 admission checks that the contraction
  matches the shape, effort × flow has watt dimensions without exponent
  overflow, the measure-adjusted generalized dimensions match the declared
  `PortKind`, and PR-2 schema-only kinds name their non-energy conserved flow
  where one exists. The original three seed schemas retain their PR-1
  Energy-only role vectors exactly.
- `PortKind::scalar_seed_schema` constructs one canonical scalar coordinate of
  any kind without inventing identity, frame/basis, or clock data; retained
  goldens use it to migrate the existing mechanical, fluid, and thermal seeds.
- `STREAM_PORT_VERSION = 1`; `StreamPort` is not a scalar `PortKind` or a fifth
  thermodynamic relation. It bundles signed mass (`kg/s`), canonically ordered
  species/element amount (`mol/s`), three momentum-rate (`N`), energy (`W`),
  and entropy-rate (`W/K`) values under one `StreamChartBinding`. Its fixed
  roles are Energy, Mass, Amount, LinearMomentum, and Entropy.
- `StreamChartBinding` separates the constituent-basis artifact and its
  explicit, canonically ordered species/element axis from the spatial
  basis/frame/orientation. It also binds the state schema, chemical reference
  state, logical clock/tick, gravity datum, and the closed
  `StreamStressWorkConvention`. Stream admission is owner-outward only until a
  public pullback exists.
- `StreamEnergyChart` structurally selects exactly one of: the canonical
  moving-stream enthalpy chart
  `mdot * (h + |u|^2/2 + g*z) + W_deviatoric`; internal energy plus pressure and
  deviatoric Cauchy work; or one coordinate from an exact mixture
  Euler/Legendre family. Internal-energy, enthalpy, Helmholtz, and Gibbs
  selections all reconstruct canonical enthalpy from the retained conjugate
  terms before transported energy is formed. The caller-declared stream energy
  rate must equal the selected chart bit-for-bit.
- `PressureWorkCrosswalk` recomputes `h = e + p/rho`,
  `volume_flow = mdot/rho`, and `mdot*(h-e) = p*volume_flow` exactly.
  `EulerLegendreCrosswalk` recomputes the mixture Euler identity plus enthalpy,
  Helmholtz, and Gibbs transforms. Both also require an
  `ExactIdentityProofRef` bound to the complete stream context.
- `ChemicalEnergyAccounting` admits either energy embedded in the selected
  state potential or one separately proved species-potential power term. It
  refuses dual ownership, a foreign chemical reference, and a foreign or
  wrong-kind partition proof. Explicit species-potential and Euler/Legendre
  modes require a species-only constituent axis.
- `ConservativeJunction::iterate_added_mass_fixed` /
  `iterate_added_mass_aitken` are schema-bound migration bridges for the
  retained mechanical scalar fixture; they must remain bitwise equal to the
  legacy results and are not general FSI operators.
- `ConservativeJunction` admits two conjugate schemas and evaluates the finite
  scalar seed as `SchemaInterconnection`/`SchemaPort` (shared effort, opposite
  flow). It localizes metadata mismatches and refuses non-scalar/non-finite
  evaluation.
- Four closed primitive variants are distinct in `PortPrimitive`:
  `ConservativeJunction`; `StorageElement` with Hamiltonian/free-energy state
  schema and constitutive-gradient operator; `DissipativeRelation` with a
  constitutive family plus mandatory monotonicity/nonnegative-production
  evidence reference; and `SourceOrReservoir` with an explicit signed
  `AccountingBoundary` carrying the same basis/frame/orientation as its port
  plus included-source/external-reservoir treatment.
- `ACCOUNTING_WINDOW_VERSION = 1`; `AccountingWindowInterval` names an exact
  nonempty interval in one logical clock without treating tick differences as
  seconds. `WindowEvidenceRef` binds a receipt/verifier/statement digest to
  that interval and to one exact `WindowEvidenceRole`: manifest closure,
  manifest row and local port, initial/final inventory, integrated row and
  local port, boundary temperature, element projection, or charge projection.
- `WindowElementSchema` either declares a nonchemical no-claim or retains a
  nonempty canonical element axis plus interval-bound projection evidence.
  Stream and chemical-power manifests require the audited form.
  `WindowChargeSchema` requires direct-coulomb, species-to-coulomb, or
  proven-neutral evidence; neutrality admits only exact zero charge values.
- `WindowManifestEntry` retains contribution, exchange, local port,
  counterparty, topology role, port kind, logical timestamp, optional external
  `AccountingBoundary`, and row-binding evidence. Stream entries additionally
  retain the complete `StreamChartBinding`; every stream in a window must use
  the same state schema, constituent basis/axis, chemical reference, spatial
  coordinates, logical clock, gravity datum, and stress convention. Internal
  exchanges are admitted only as two reciprocal rows of one kind and timestamp.
  Manifest-closure evidence carries the complete canonical row subjects,
  including topology, boundary, and retained stream context, so it cannot be
  rebound to a different or truncated same-interval manifest.
- `WindowBalance` and `WindowInventorySnapshot` carry finite typed energy,
  mass, canonical elemental amounts, charge, and entropy at an exact endpoint.
  `IntegratedWindowTransfer` uses one convention for every signed value:
  positive into the audited system. `BoundaryEntropyBreakdown` makes advected,
  diffusive/chemical, heat, and radiation terms mandatory and retains a
  contribution-bound constant/profile temperature reference whenever heat or
  radiation is nonzero.
- `AccountingWindowSpec` owns the exact manifest, a common stored/stream energy
  reference contract, element/charge/entropy conventions, typed nonnegative
  tolerances, optional reference environment, and closure evidence. An empty
  manifest is legal only as an explicitly evidenced isolated window. Filled
  rows must cover the declared manifest exactly; ordinary power rows may carry
  only balance axes compatible with their `PortKind`, while stream rows carry
  the full multi-balance vector.
- `WindowAuditReport::audit` returns `Err` only for malformed, incomplete,
  foreign, nonfinite, or semantically incompatible input. A physically failed
  audit returns `Ok(report)` with retained `WindowAuditViolation`s. It computes
  `R = (final - initial) - sum(transfer_into)` for energy, mass, each element,
  and charge; computes `S_gen = delta(S) - sum(S_boundary,into)`; applies
  absolute typed tolerances to equalities and a one-sided lower bound to
  entropy generation; and retains canonically ordered per-port and per-exchange
  ledgers. Every internal exchange is independently gated for zero net energy,
  mass, each element, charge, and entropy transfer, so defects cannot cancel
  between pairs.
- `ExergyEnvironment` and `ExergyLedger` optionally report the narrow
  Gouy-Stodola diagnostic `T0 * S_gen` for a finite positive explicit `T0`.
- The raw scalar `Port` remains a backwards-compatible, non-admitting numeric
  container. The legacy `conjugate_to`/`interconnect` migration oracle composes
  only the original mechanical, fluid, and thermal seed kinds and refuses all
  schema-only kinds; raw construction or arithmetic alone carries no
  dimensional, coordinate, clock, identity, or conservation certificate.
- `interconnect(kind_a, kind_b, effort, flow) -> Result<Interconnection,
  CoupleError>` — a Dirac structure (shared effort, opposite flow) whose
  `interface_power` is `0` exactly (power-conserving by construction); refuses
  incompatible ports. `interface_power(&[Port])` = `Σ effort·flow`.
- `EnergyAudit` — `record`, `max_generation`, `is_passive(tol)`: the legacy
  `is_passive` name checks only caller-supplied scalar interface imbalance at
  each recorded exchange. A nonzero balance is a bug alarm, not a proof of
  whole-system or closed-window passivity.
- `AitkenRelaxation::new(omega_init, omega_max)` + `next_omega(residual)` — the
  scalar Δ² dynamic relaxation factor, magnitude-capped.
- `iterate_fixed_relaxation` / `iterate_aitken` — the added-mass interface
  fixed point under fixed vs Aitken relaxation → `FsiResult { converged, steps,
  solution, final_residual }`.

### `vibroacoustic` (bead frankensim-fsim-vibroacoustic-wgkq7)

The crate's first numerical engine: modal structure x rigid-wall-cavity
x exterior-radiation coupling, frequency domain under `e^{-i omega t}`
(matching `fs_bem::helmholtz`), hysteretic stiffness `k (1 - i eta)`.
Sign convention: positive structural deflection points AWAY from the
cavity (composes with the BEM's outward panel velocities).

- `StructuralModes` (mass-normalized, `fs_modal` convention) /
  `CavityModes` / `AcousticMedium` — the data carriers; any cavity
  basis plugs in, `rectangular_cavity_modes` (analytic, with a
  sufficiency-PROVEN enumeration box — an isotropic heuristic silently
  dropped axial modes of elongated cavities, executed regression) and
  `helmholtz_resonator_mode` (lumped lowest mode, `2 (8/3 pi) a` end
  correction; the psi = 1, Lambda = V presentation is EXACT, proven by
  the 3-state reduction) ship.
- `assemble_coupling` — area-integral `C_rq = INT phi_r psi_q dA`,
  pinned against closed-form sin x cos overlap integrals.
- `project_radiation_impedance` — panel-space `Z` (`p = Z v`) reduced
  to modal `Zm` with area weighting on the force side only.
- `VibroacousticModel::frf` / `frf_truncated` — dense complex coupled
  solve; every response carries the complete power breakdown with BOTH
  audit residuals (`input = structural + interface + radiated`,
  `interface = cavity`) — algebraically exact identities, measured at
  1e-13.
- `frf_with_convergence` — truncation error OBSERVED (full vs halved
  bases) and refused (`FS-COUPLE-VIBRO-TRUNCATION-NOT-CONVERGED`)
  above the caller tolerance; near-resonance rows are genuinely
  truncation-sensitive (measured 0.28 on the casebook fixture) and the
  casebook exposes rather than hides that.
- `undamped_natural_frequencies` — exact `x = omega^2` linearization
  (unit-triangular `B`), pinned by the two-oscillator closed form at
  ~1e-16 across a five-decade coupling sweep, a first-order added-mass
  falsifier (quartering 4.00), and a multi-mode independent-determinant
  oracle.
- `VibroError` — stable `FS-COUPLE-VIBRO-*` refusals.

### `modal_acoustic_time`

`ModalAcousticTimeModel` is the generic fixed-rate realization of independent,
mass-normalized structural modes. It consumes generalized force in
`N / sqrt(kg)`, advances viscously damped modal coordinates by the exact
zero-order-hold transition in underdamped, critically damped, and overdamped
regimes, and emits physical observer pressure in pascals. Acoustic observation
uses a caller-supplied complex pressure-per-modal-velocity transfer at each
natural frequency under `exp(-i omega t)`; no material names, digital gains,
or mastering values occur in this layer. Energy/work diagnostics and
transactional state/pressure budgets accompany every sample.
The default initial state is zero displacement and velocity. A caller whose
window begins after a load has already settled may explicitly initialize the
exact static compliance `q = F/omega^2`, `qdot = 0`; using that operation to
erase a real load-on transient is outside the contract.
`observer_pressure_with_transfers_about_static_equilibrium` additionally
removes the exact held-force compliance `F/omega^2` before applying the
complex transfer's displacement-quadrature term. This prevents a static load
from being misreported as DC acoustic pressure while leaving free modal motion
unchanged; the caller must provide the same generalized ZOH force that drove
the structural transition.

The acoustic realization is explicitly narrow-band: one complex transfer value
per mode is exact only at that mode's natural frequency. Broadband radiation,
propagation delay, and feedback impedance require the separately tracked stable
rational-fitting path; this module does not claim them or passivity.

### `broadband_radiation`

Generic offline bridge from solver-neutral sampled scalar-input radiation transfers to
real-tesseral SH filters. Radiation solvers supply complex-SH training rows, direct
held-out fields, source semantics, and diagnostics; `fs-couple` does not depend on them.
It converts `exp(-i omega t)` to `fs-vfit`'s `s=+i omega`, fits prewarped abscissae, and
applies unprewarped Tustin. Direct disjoint held-out fields gate artifacts; stable
proper fits may retain caller-admitted constant feedthrough but always forbid `s*e`.
The runtime owns visible per-filter state, superposes in fixed input order, and refuses
whole samples transactionally. EstimateOnly covers the fixed reference before `1/r`;
moving-source, FW-H, near-field, feedback, room, head, and passivity remain unclaimed.

### `acoustic_realize` / `pcm_wav`

`realize_assembly` turns an `fs_scenario::AcousticAssembly` into observer
pressure by composing nameless primitives: a prestressed Euler–Bernoulli
waveguide (`fs-nlmodal::prestressed_beam_omega` / Kirchhoff–Carrier), a
`bernoulli_aperture`, a `stribeck_friction` contact, a
`traveling_wave_line`, `ModalAcousticTime`, and compact modal radiators.
There is no instrument crate and no instrument algebra. A guitar or
clarinet is one filling of those objects.

- String path, `EA = 0`: triangular-pluck (or bow) modal ICs marched by
  `ModalAcousticTimeModel`, compact transfer
  `H = i ω ρ A_odd / (4 π r √(μ L / 2))` for odd sine modes. Even
  modes have vanishing monopole area and radiate as compact dipoles
  `p ∝ ρ Π̈ / (4 π r c)` with `Π = w ∫ φ (x−L/2) dx`. Fletcher
  inharmonicity `ω_n = n ω_1 √(1+B n²)` when `EI > 0`. Stokes air
  drag from `GasState` plus the authored internal floor. A second
  polarization at `1+detune` is a second member on the same clock
  and shares every plate (it is not an independent unused body).
- String path, `EA > 0`: `fs-nlmodal::kirchhoff_carrier_string`
  (fixed-fixed) or `kirchhoff_carrier_moving_end` (Dirac join)
  + Gonzalez
  pHS step, with the same stiff-string frequencies overwriting the
  linear ω_n. Bridge force `T_eff y'(0)` drives every listed body mode
  (Carcagno-band Sitka pair is a constructor, not a named guitar).
- Bow: MWS regularized friction (steep stiction ramp + falling kinetic
  shoulder). Helmholtz motion is possible with enough modes; it is not
  guaranteed and not a measured rosin curve.
- Reed path: quasistatic or massive Bernoulli valve. Isolated
  cylindrical bores use the `acoustic_chain` ODE with a
  `ViscothermalPin` (massive
  reed = `mass_spring_damper` + jet + face flow). Blow/reed on a
  moving-end string×plate×duct is a leftover `join_port` inlet.
  Moving-end string×duct without a plate is
  `transformer(waveguide, chain, A_inlet)`. Fixed-fixed
  string×duct without a plate shares the ODE clock without
  a force port (`φ(0)=0`).
  Reed lay
  is `fs-dcontact` (not a private Hunt–Crossley).
- `ThinPlate` modes come from `fs-plate` + `fs-modal` (certified
  eigenpairs), not caller Hertz. Mouth pressure × area drives the
  plate; plate volume velocity returns to the waveguide
  (structure–bore).
- String obstacles are `fs-dcontact` power-law potentials on sine
  collocation. On the Kirchhoff–Carrier path the conservative
  potential is `ContactStorage` wrapping the KC Hamiltonian, so
  Gonzalez sees the contact energy. Optional `mu_kinetic` adds
  `fs-tribo` Coulomb traction as a modal port force (not a gradient
  of `H`). The linear modal path still applies both as a force.
  Hunt–Crossley `χ` is a scenario field on the obstacle and a
  dissipative port force (`modal_hunt_crossley_forces`), never a
  term in `H`.
- Reed and structure–bore time port: TMM `R(ω)` sampled on the DFT
  grid and inverse-transformed to an FIR. Isolated linear blow (no
  body) is the same `DelayedFilter` object filled with `IFFT[Z_in]`
  (or the mouth transfer `p/u` when the end is open). A vented
  reflectance FIR does not ring a measurable period; the impedance
  FIR does, so tone-hole shortening is TMM-emergent on one port
  type. There is no one-pole `TravelingWaveLine` fallback on those
  paths. Unflanged-open Nyquist `ka > 1` is refused.
- Plate damping: a single authored viscous ratio is the two-point
  Rayleigh fit (`fs-material::visco`) through that ratio at ω₀ and
  4ω₀ so higher certified modes sit on the stiffness limb. Radiation
  reaction on linear compact radiators is the Rayleigh-integral
  baffled-piston face impedance (same half-space kernel as `fs-bem`,
  written in-tree so couple does not depend on bem) fitted and
  passivity-repaired by `fs-vfit`. `fs-bem` is not a production
  dependency of this crate (cycle through `fs-feec`). Linear plates
  accept in-plane pretension and clamped edges. Von Karman is the
  isotropic SS sine path when `e1 = e2` and edges are free to
  rotate; clamped or orthotropic bending uses DKT-sampled
  displacement with the same sine Airy membrane channel
  (`von_karman_sampled_plate`). The compact far-field observer is
  the baffled on-axis piston `p = ρ A ÿ / (2 π r)`, the same
  half-space as the self-load. A full BEM observer is the
  solver-neutral `DirectionalFarFieldTable` /
  `far_field_observer_pressure` (`p = F e^{ikr}/r`); `fs-bem`
  produces `FarFieldTable` as a **dev-dependency** only (the
  production cycle `couple → bem → solver → feec → couple` is
  refused).
  The finished pressure history is passed through ISO 9613-1
  absorption (`air_path`) with the assembly's explicit humidity;
  Stokes–Kirchhoff is only the fallback outside the ISO window.
  Authored string `ζ` at the fundamental becomes a Prony branch
  (`GeneralizedMaxwell::matching_loss`); higher modes see `η(ω)/2`.
  An optional `HelmholtzCavity` faces the plate monopoles as a
  flow-driven pHS whose damper is compact-mouth `Re Z_rad(ω₀)`.
- Bow roughness: an optional `ContactTexture` is a declared
  self-affine height spectrum. `fs-tribo::surface_excitation` samples
  it; the height perturbs the normal load that Stribeck sees. No
  measured rosin is invented.
- Three-way clock: when `EA > 0` the string member is Kirchhoff–
  Carrier with contact inside `H`. Polarizations share the plate and
  the duct FIR. The certified three-pHS Dirac interconnection is
  `common_flow_dirac(moving_end_waveguide, transformer(plate, load))`
  in `fs-phs` (free-end attachment port, area transformer, Kirchhoff
  star). The load is a Helmholtz cavity or an `acoustic_chain`
  duct — the same object. `PrestressedString.moving_end` selects
  that clock. Bow, obstacle stations, and blow/reed are leftover
  ports of `join_port` (Stribeck and Hunt–Crossley stay port
  forces). Von Karman is the same plate pHS with quartic storage.
  `EA > 0` is `kirchhoff_carrier_moving_end` on the same
  free-fixed port. Fixed-fixed sines keep the one-way bridge
  force; with a plate the duct is still the ODE chain
  (bridge force on the leftover plate port), not the FIR line.
  Without a plate the same cylindrical chain is still the ODE
  clock. A cylindrical bore (uniform or stepped) is the
  `acoustic_chain` LC ladder with an all-regime
  `ViscothermalPin` from `GasState` (`μ`, `γ`, `Pr`;
  series `R` and thermal `G` are three-term Foster networks
  collocated to Bessel Zwikker–Kosten `F(r_v)`, all shear
  numbers). A single linear taper is `spherical_cone`
  (`ψ = x p` plus the Euler near-field shunt); mixed
  cylinder+taper runs stay frustum LC so physical `p`/`U`
  match at interfaces.
  Open tone holes are `AcousticTap` side-branch
  inertances. A quasistatic reed is the Bernoulli port on that
  inlet; a massive reed is `mass_spring_damper` plus the same
  jet, face-velocity flow, and a Hunt–Crossley lay as a
  dissipative port force (not a term in H). With a linear plate
  the mouth is a transformer. Frequency-by-frequency Bessel TMM
  (`LossModel::Bessel`, spherical-wave cones included) is the FIR
  path. 3-D jet broadband remains a no-claim. Plate ×
  Helmholtz in `cavity_phs` is the same transformer, not a
  staggered pair of steps.
- Von Karman geometric nonlinearity is the isotropic simply-supported
  analytic primitive when those hypotheses hold, and the FE-sampled
  Airy construction otherwise. The membrane channel remains the
  sine Airy basis (in-plane movable).
- Reed observer adds the compact far-field dipole of the slit force
  `F = Δp · w · y`. That is not 3-D jet broadband and not the
  fs-aeroac 2-D Curle Hankel kernel. 3-D jet broadband remains a
  no-claim. There is no instrument crate.
- `encode_pcm16_wav` maps pascals through a declared full-scale. It never
  peak-normalizes.

## Invariants

- The Dirac interconnection conserves interface power EXACTLY (to roundoff) —
  the G0 law; incompatible ports are refused.
- Vibroacoustic power identities are exact by derivation and asserted
  at solver roundoff; the input-balance identity is TAUTOLOGICAL for
  any consistently assembled system, so the one-sided interface-normal
  mutation is caught by the `interface == cavity` cross-row residual
  (measured alarm 2.0 vs 1e-16). Repeat solves are bitwise identical.
- Every admitted v2 schema has a shape-compatible power contraction, a checked
  watt-dimensional effort/flow product, and energy. Rotational, electrical,
  and chemical PR-2 kinds additionally require angular momentum, electric
  charge, and amount respectively; legacy seed role vectors remain unchanged.
  Field duality records separate effort/flow function-space roles,
  includes its integration-measure dimensions in the checked product, and
  explicitly assigns the measure to the density side before kind-specific
  dimension comparison.
- Stable relation IDs cannot alias their owned port/boundary IDs. A
  conservative pair requires distinct port IDs, matching physical metadata and
  clock/timestamp, plus outward-from-owner conventions on both sides. PR-1
  refuses common-frame orientations until a public pullback can prove them.
- Storage, dissipation, and source ownership remain separate typed claims. A
  conservative junction cannot stand in for any of them; a dissipative
  relation cannot exist without an evidence reference; a source cannot exist
  without a named accounting boundary whose coordinate/sign binding exactly
  matches the port.
- Every admitted stream has a non-empty duplicate-free constituent axis, one
  shared outward sign convention, finite fixed-dimension rates, fixed bundled
  conservation roles, and exactly one energy chart. Alternate charts cannot
  enter without exact scalar identities and context-bound durable evidence;
  the actual rate axis must equal the proof-bound axis, chemical power is owned
  exactly once, and pressure/deviatoric work uses the single normative
  Cauchy-tension-positive integrated outward-power convention.
- Every accounting report is tied to one exact logical interval, subject-bound
  evidence roles, exact initial/final endpoints, one canonical expected-row
  set, one element axis, one charge policy, and one entropy convention.
  Caller order cannot change arithmetic: rows are sorted by stable contribution
  ID and reduced with a fixed-order compensated sum.
- Global first-law, mass, element, and charge closure use two-sided typed
  tolerances. Entropy generation is never absolute-valued or clamped: only
  `S_gen >= -entropy_tolerance` passes. Individual boundary contributions may
  be negative. Internal pair closure is a separate equality gate on every
  audited axis, and violations retain the exact pair, ports, and rows rather
  than guessing a single cause from a global residual.
- Inventory snapshots and filled contributions retain the exact energy
  reference, element schema, charge policy/basis, and entropy convention under
  which they were admitted; audit refuses same-interval context rebinding.
- A constant-temperature thermal power row must carry only heat entropy and
  satisfy `energy = T * integrated_heat_entropy` bit-for-bit. Radiation on a
  thermal power row requires profile-bound external evidence. Profile-bound
  thermal and stream heat/radiation integrals remain external-evidence claims.
- The energy audit reports an interface-balance failure exactly when some
  caller-supplied exchange has absolute imbalance above `tol` or is non-finite.
- On the added-mass fixture (`μ ≥ 1`): naive staggering (`ω = 1`) diverges while
  Aitken-relaxed coupling converges to `x* = c/(1+μ)`; Aitken never takes more
  steps than a stable fixed under-relaxation.

## Error model

Structured `CoupleError`; no panics. In addition to legacy incompatible kinds,
errors cover invalid IDs/empty shapes, pointwise or field-measure dimension
overflow, non-power products, kind-specific or shape/pairing mismatch, missing
kind-required conservation roles, schema-only kinds sent through the raw
interconnection oracle, localized schema conjugacy mismatch, identity aliasing,
accounting-boundary coordinate mismatch, scalar/non-scalar misuse, wrong
added-mass port kind, and non-finite schema-bound values. Stream errors localize
non-finite fields/components, empty or duplicate constituent axes, non-outward
orientation, proof-bound/actual axis mismatch, species-potential accounting on
an element axis, chart/proof binding or identity mismatches, non-positive
density, chemical double counting, and declared/chart energy-rate disagreement.
Window errors additionally localize clock/interval/endpoint mismatches,
evidence-role/subject mismatches, malformed or one-sided internal exchanges,
duplicate rows/ports/elements, missing or unexpected manifest rows, foreign
stream context, missing stream/chemical element schema, unproved charge policy,
negative inventory/tolerance, nonpositive temperatures, incompatible
power-port transfer axes, absent temperature evidence, nonfinite integrated
values, and element/charge axis violations. Physical conservation or
second-law failures are report violations, not construction errors.

## Determinism class

Fully deterministic: schema admission canonicalizes roles by enum order,
stream admission canonicalizes constituents by typed ID, and window admission
canonicalizes manifest, per-exchange, port, and element ledgers by stable typed
identity. Exact crosswalks use a pinned operation order and bit equality;
window reductions use a fixed-order compensated sum. Interconnection, relation
construction, audit, and iterations are pure functions of their inputs.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`src/broadband_radiation.rs`: G0 complex/real-SH round trip/reconstruction against
`fs-bem`; admission, neutral replay, properness/direct-feedthrough retention, and
held-out mutation; causal transactional runtime superposition and bitwise replay.

`src/vibroacoustic.rs` unit tests: closed-form overlaps, two-oscillator
split + dropped-coupling mutation, added-mass first-order falsifier,
exact energy balance + one-sided normal-flip mutation, reciprocity,
truncation refusal, cavity/resonator pins (incl. the elongated-cavity
regression), multi-mode pencil vs independent determinant, refusals,
determinism. `tests/vibroacoustic_casebook.rs`: the fs-plate ->
fs-couple -> fs-bem composition — box-with-flexible-top vs independent
perturbation with the sealed-cavity stiffening direction, coupled FRFs
with observed per-frequency truncation deltas, and the radiation solve
whose per-period energies audit GREEN through a real
`WindowAuditReport` (with the `entropy_on_nonthermal_power_port`
refusal retained as doctrine: dissipation leaves through THERMAL
ports); plus the guitar T(1,1) comparison against Carcagno et al. 2018
with the exact two-DOF product invariant and the material-vs-measured
damping inequality.

`tests/couple.rs` (38 cases): v2 scalar-seed migration goldens for all three
legacy kinds and bitwise whole-result migration of the added-mass fixture;
rotational/electrical/magnetic/chemical watt-dimensional goldens, required-role
admission, raw-oracle refusal, and semantic kind mismatch refusal; schema
fail-closed metadata; localized junction mismatch; non-scalar refusal by the
scalar evaluator; field-duality measure dimensions, density-side assignment,
measure-application overflow, and distinct effort/flow spaces; all four
primitive descriptors and identity-alias refusals;
complete moving-enthalpy stream admission and canonical constituent ordering;
one proved explicit-chemical contribution;
bit-exact enthalpy/internal-energy pressure-work equivalence; one-ULP and
foreign-context crosswalk/stress-evidence refusal; wrong-identity and
non-positive-density refusal; all four exact Euler/Legendre coordinates
reconstructing the same canonical enthalpy and chemical double-count refusal;
proof-bound axis mismatch and species-mode/element-axis refusal;
empty/duplicate/non-finite/wrong-orientation stream refusals and
declared-energy mismatch;
green first-law/mass/element/charge/entropy closure with explicit signed
advected/diffusive/heat/radiation rows and `T0*S_gen`; retained red first-law
and negative-entropy reports; exact missing-row and foreign-evidence-role
refusal; exact manifest/element/charge evidence-subject refusal; same-interval
energy/accounting-context rebinding refusal; stream/chemical element-schema
and typed power-axis refusal, including cancelling forbidden entropy slots,
inexact constant-temperature energy, and constant-temperature radiation;
boundary-temperature and nonfinite fail-closed admission; all six
permutations of a cancellation-sensitive compensated sum; an explicitly closed
isolated entropy window; one-sided internal-pair refusal; and a hidden internal
entropy source localized to the exact two rows, ports, and exchange while an
external row masks its pair-local energy defect from the global first law;
legacy power-conjugate ports; exact interface power conservation and
incompatible-port refusal; energy-audit imbalance and non-finite alarms; the
Aitken Δ² factor; naive staggering diverges where Aitken stays stable; Aitken
accelerates over stable fixed relaxation; light added mass converges naively;
determinism.

## No-claim boundaries

- `acoustic_realize` is a composed description→waveform, not a named
  instrument product and not a 3D jet. Bow is regularized friction, not
  measured rosin or Helmholtz-guaranteed. A linear `ThinPlate` is DKT +
  certified modes + compact monopoles. `geometric_nonlinearity` selects
  the von Karman modal pHS (SS sine or FE-sampled DKT `w`), not a
  full shell or an in-process BEM body.
  String+duct+plate share one sample clock and, when `EA > 0`, a
  Gonzalez string whose contact lives in `H`. The three-pHS Dirac
  join lives in `fs-phs`; this realize path still uses the FIR
  duct and one-way bridge force. Isolated blow is an impedance
  FIR, not a reflectance FIR. 3D broadband jet noise remains a
  no-claim. WAV encoding is a physical-scale dump, not mastering.
- `vibroacoustic` is frequency-domain only: time-domain realization is
  the vector-fitting bead's scope; non-rectangular cavity bases beyond
  the lumped Helmholtz mode arrive through the same `CavityModes`
  carrier (numeric Laplacian producers are not shipped here); no
  modal-density/SEA regime; the casebook's box fixture is validated by
  independent perturbation and exact power audits; the guitar fixture
  compares the coupled T(1,1) pair against Carcagno et al. 2018
  (CC-BY, Table I) through a DECLARED shape-surrogate top (deviations
  +10-17% recorded, product invariant exact) — measured mobility
  MAGNITUDES exist only as figures in the license-compatible
  literature and are not compared; `condition`-free dense solves inherit `fs_la::eigen_complex`
  boundaries.

- The FSI fixture is the classic LINEARIZED added-mass interface map
  (`H(x) = −μx + c`) — enough to reproduce the instability and its fix; a full
  nonlinear FSI solve over real fluid/structure subsystems is the consumer.
- `AitkenRelaxation` is the scalar Δ² relaxer; the vector INTERFACE
  QUASI-NEWTON (IQN-ILS) accelerator and MULTIRATE co-simulation are staged.
- PR-2 admits the rotational/electrical/magnetic/chemical vocabulary and
  kind-specific generalized dimensions. It does not by itself prove causality,
  DAE index, source closure, a constitutive law, or a physical port adapter.
- `ChemicalPotentialAmountFlow` describes species electrochemical potential
  and species amount flow. Reaction affinity/extent rate requires an explicit
  stoichiometric crosswalk and is not represented by that kind.
- The scalar evaluator does not execute vector/tensor/field or general
  multi-port Stokes–Dirac operators. `StorageElement` and
  `DissipativeRelation` carry durable public operator/evidence references; this
  crate does not execute or validate the referenced domain law in PR-1.
- Field-duality admission checks dimension arithmetic including the declared
  measure; it does not prove quadrature, trace pullback, orientation, or the
  numerical duality operator.
- `AlongFrame`/`AgainstFrame` can describe schemas and audit boundaries, but
  PR-1's conservative junction refuses them rather than assuming an unproved
  effort/flow pullback. A later neutral transfer API must make that transform
  explicit before admission.
- PR-3 admits `StreamPort`, one selected energy-accounting chart, exact numeric
  crosswalk identities, and context-bound durable proof references. It does
  not execute the referenced verifier, equation of state, stress operator,
  constituent-map artifact, or chemical partition proof.
- PR-3 freezes a stationary or co-moving accounting surface. It does not model
  boundary velocity, Reynolds-transport terms for a moving/deforming control
  boundary, or cross-frame kinetic-energy transforms.
- PR-4 audits direct mass, exact caller-projected element amounts, direct or
  externally projected coulombs, first-law energy, boundary entropy, internal
  pair closure, and optional `T0*S_gen`. It does not audit momentum, species
  production/reaction progress, mass-versus-molar consistency, or derive
  element/charge values from stoichiometry, valence, molar mass, or Faraday's
  constant.
- `fs-couple` checks evidence interval, semantic role, and subject identity but
  does not execute the referenced verifier or discover omitted sources. A green
  report means the supplied, closure-evidenced scalar algebra closes; it does
  not independently prove quadrature, physical-time integration, source
  discovery, EOS validity, species-to-element/charge projection, or the truth
  of the retained evidence statement.
- Logical clock ticks are ordering coordinates, not seconds. Every window row
  must arrive pre-integrated in physical time with an integrated-transfer
  receipt; endpoint sampling is not promoted into an integral.
- Constant-temperature thermal power rows check `Q = T*DeltaS_heat` exactly
  and therefore refuse radiation; radiative thermal rows require the profile
  evidence lane. Profile temperature, stream heat/radiation, and
  nonequilibrium entropy decompositions are retained external-verifier claims;
  this crate does not execute their quadrature or infer an average temperature.
- PR-4 currently admits a conservative homogeneous stream-context window:
  external stream rows must share state schema, constituent basis/axis,
  chemical reference, spatial coordinate binding, clock, gravity datum, and
  stress convention. Heterogeneous inlet/outlet contexts require an explicit
  upstream projection/crosswalk into one common audit context; this PR does not
  define or execute that crosswalk.
- Included-source entropy means entropy transferred into the audited system.
  Internally generated entropy must remain in `S_gen`; supplying it as an input
  would subtract it and is outside the admitted convention.
- The exergy ledger is only Gouy-Stodola destruction relative to the named
  environment, not a complete open-system availability-flow balance.
- Neither `WindowAuditReport::is_green`, `StreamPort`, nor an
  `AccountingBoundary` is a subsystem, discretization, co-simulation, or
  whole-machine passivity certificate.
- The two `fs-iface::SpaceType` entries record neutral effort/flow field roles;
  they do not supply a mortar/Nitsche/harmonic transfer operator, certify that
  the declared pair is dual, or prove inf-sup compatibility for a particular
  domain adapter.
- The energy audit's balances are supplied by the caller each exchange; wiring
  them onto the ledger is the coupling driver's integration.
- Dirac interface losslessness does not establish passivity of component
  storage/dissipation/source laws, spatial or temporal discretization,
  interface transfer, nonlinear iteration, multirate windows, or the coupled
  system. Those obligations require a signed, closed-window energy audit.
