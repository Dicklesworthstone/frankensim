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

### `render`

Block render API (music bead `frankensim-music-v8-root-3ez8g.2.1`): the
callback-shaped hosting layer for performance images. `RenderContext`
owns admitted voices (`ReedBoreVoice` — the `realize_reed_bore` physics
restructured for blockwise stepping, with the one-shot realizer now a
thin wrapper over it so the two paths share one loop body;
`ModalStringVoice` — the exact-ZOH modal runtime hosted verbatim) and
advances them one block at a time into a pre-sized scratch buffer.
INVARIANTS: a block boundary performs no arithmetic, so block size is a
pure loop-partition choice — the one-shot path and every block partition
are bitwise-identical (tested at 64/571/full); typed `ControlDelta`s
apply only BETWEEN blocks, transactionally, and are logged as
`ControlRecord`s carrying their D17 lift description (empty = pure input
move); cancellation (`render_under_gate` over `fs_exec::CancelGate`) is
polled only at block boundaries with drain semantics — an in-flight
block completes, the rendered prefix is a whole number of blocks, and a
resumed context continues bitwise-identically (tested). ALLOCATION:
the massless-reed voice with an empty plate bank is allocation-free per
block (counting-allocator gate); the massive-reed lay path
(`dissipative_modal_forces` returns a `Vec` per sample) and the modal
voice (`step` builds a per-sample energy frame) are DISCLOSED allocating
voices — fusion candidates (bead 3ez8g.15), not silently admitted.
Determinism: one-host bitwise, inherited from the hosted kernels.
No-claims: no device/OS audio host (Franken-only policy; live output
would be a quarantined adapter decision); no image hopping yet (D17
lifts beyond input parameters land with the articulation/track beads);
refusal mid-block poisons the context (mid-sample state is not rewound)
and is documented rather than hidden. Boundary facts a consumer must
know, both discovered by this bead's battery: the characteristic-line
realization refuses `UnflangedOpen`/`FlangedOpen` terminations whose
Nyquist `ka` exceeds 1 (at 48 kHz that is every bore wider than
~2.27 mm radius — the measured-load `Termination::Tabulated` lane,
bead zolja, is the lift), and the Zwikker–Kosten Bessel path's former
`r_v ≈ 160` continued-fraction ceiling was an iteration-budget bug fixed
in `fs-phs` (budget now scales with `|z|`).

### `music_render` (binary)

The music lane's render CLI (bead ib15w): `music_render <reed|string>
<out.wav> [--seconds S] [--block N] [--full-scale-pa P]`. Pinned fixture
compositions render through the block API and encode through the ONE
pascals→PCM owner, `pcm_wav::encode_pcm16_wav` — the recorded seam
decision (beads ib15w + h7xu5.7.8): the cinematic stereo/receipt-hashed
encoder stays cinematic, no third RIFF writer exists in the music lane.
Laws, all e2e-tested through the real binary: same arguments →
bit-identical WAV + provenance sidecar; existing output/sidecar paths
REFUSE (evidence is never overwritten); full-scale is a MAPPING, never
normalization — the same physics renders identical `peak_pa` at any
full-scale, PCM peaks scale inversely, and an undersized full-scale
CLIPS with the count reported (the reed fixture peaks ≈10.7 kPa —
mouthpiece pressures are kPa-scale, so full-scale choices must be too);
sample rate pinned at 48 kHz (ecosystem coherence with fs-psycho). The
provenance sidecar (`frankensim-music-render-provenance-v1`) carries
fixture, rates, clip count, peak/rms pascals, and the domain-hashed
blake3 of the WAV bytes — no wall-clock, no commit stamp (git history
of committed artifacts carries those). No-claims: fixtures are pinned
compositions, not a project/assembly loader (that is fs-cli product
territory); mono only; no mastering, room, or loudness processing of
any kind.

### `bakeoff`

Bake-off receipt schema and protocol for the music-program claims registry
(bead `frankensim-music-v8-root-3ez8g.1.4`). `BakeoffReceipt` records two
contender images run on one shared-card fixture against a caller-supplied
reference QoI map, with logical budgets (states, steps, solver iterations —
never wall-clock, so bytes are host-stable), observed failure modes, and a
reviewed outcome. `BakeoffOutcome` has exactly three variants —
keep-both, keep-for-subset, refuse-newcomer — and structurally cannot
express deleting an admitted image (doctrine D21: menus, not winners).
Receipts serialize to a canonical line-oriented byte encoding
(`frankensim-bakeoff-receipt-v1`, strict round-trip decoder) and are
content-addressed via `fs-blake3` under
`org.frankensim.fs-couple.bakeoff-receipt.v1`. The harness MEASURES;
the outcome is a reviewed judgment passed in by the caller — residual
tables cannot adjudicate claim scopes alone. The executed
`tests/bakeoff_string.rs` fixture (modal-ZOH vs the Gonzalez-stepped
pHS modal bank, analytic damped-oscillator reference) golden-pins its
receipt at `tests/receipts/string-modal-zoh-vs-phs-bank.bakeoff`; the
pHS contender's phase error is gated against the implicit-midpoint
dispersion model `(2/dt)·atan(omega·dt/2)` rather than a blanket state
band, because pointwise residuals near zero-crossings are unboundedly
phase-sensitive. No-claim: a receipt proves the comparison ran and what
it measured, never that a claim scope is correct — scope promotion is
the per-track gates beads' review job through the instrument-claims
registry.

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
  paths. Unflanged-open Nyquist `ka > 1` is refused; a
  flanged mouth uses the Rayleigh piston and is not.
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
  (`ψ = x p` plus the Euler near-field shunt). Any
  multi-section chain that contains a taper stitches
  ψ-lines onto the LC ladder with transformer `x` at
  each interface.
  Open tone holes are `AcousticTap` side branches on
  the ODE path: a compact neck is one inductor; a long
  neck (`kℓ > 0.2` or `ℓ > 4b`) is the same 2-cell LC
  line plus flanged mouth the TMM chimney already is.
  The same `ViscothermalPin` sits on a lumped neck or
  on each chimney cell of a line.
  A scenario `LocallyReactingWall` is the ODE `WallPin`
  (per-cell shunt LC, `A_w = 2π a dx` with slant on a
  taper) and the TMM `Y' = 2π a slant / Z'` addend on the
  FIR `Z(ω)` / `R(ω)` path, including tone-hole chimneys.
  Same numbers both clocks. The one-pole `TravelingWaveLine`
  is not that path.
  A closed pad is the TMM cavity compliance on the station
  cell; `foster_branches > 0` adds the bore's thermal
  Foster on that remaining `C`. A scenario `ToneHole.open_fraction`
  σ is the same admittance mix
  `Y = σ Y_open + (1−σ) Y_closed` the TMM `HoleState::Vent`
  already is.
  A quasistatic reed is the Bernoulli port on that
  inlet; a massive reed is `mass_spring_damper` plus the same
  jet, face-velocity flow, and a Hunt–Crossley lay as a
  dissipative port force (not a term in H). With a linear plate
  the mouth is a transformer. An unflanged ODE mouth is the
  same compact `(R, X)` as the TMM: `Re Z_rad` plus mass
  `Δℓ = 0.6133 a` on the last flux. A flanged mouth is
  `0.8216 a` (Rayleigh piston above `ka = 1` on the TMM
  path; the unflanged Nyquist ceiling does not apply).
  `foster_branches > 0` replaces the pin `Re Z_rad` with
  a Foster match of that load across a band (unflanged
  samples stay under `ka = 1`; flanged may use the
  piston). Tabulated `R(ω)` remains the FIR path.
  An open chimney with
  `foster_branches > 0` carries the bore's Foster series
  on a lumped neck, or series plus thermal Foster on
  each cell of a long-neck line (extra states). Frequency-by-frequency
  Bessel TMM (`LossModel::Bessel`, spherical-wave cones
  with local-radius lossy substations) is the FIR path. 3-D jet broadband remains a
  no-claim. Plate ×
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

- `mm_line` module (bead 3ez8g.4.2): the multimodal characteristic-line
  runtime for brass. `MmLineBank` realizes one exact-FIR plane-mode
  reflectance line PER VALVE COMBINATION from the fs-duct multimodal
  authority (eight round trips of impulse response — an executed
  lesson: a four-trip tail left the lowest peak 92 cents flat from
  wrap-around aliasing), terminated into the tabulated bell load
  inside its support with a DISCLOSED analytic fallback outside it
  (splice bins counted in the per-combo receipt). Passivity (`|R| <=
  1` on the realization grid) is enforced ON THE STORED TAPS — a
  valve switch rebuilds from them, so enforcing only the live line
  would resurrect an unenforced reflectance mid-performance. Valve
  switches at block boundaries use the carry-outgoing-history lift
  (in-flight waves persist; the new taps govern future reflections)
  plus an optional linear CROSSFADE from the old line — an
  instantaneous reflectance swap is genuinely discontinuous physics,
  and the crossfade is the no-click lane; every switch is logged with
  its lift name and carried sample count. `MatrixFirLine` is the
  full N^2-FIR bake-off arm with per-mode incoming-energy logging.
  `cup_junction` + `CupState` is the D18 mouthpiece-cup lumped
  compliance, TRAPEZOIDAL so the runtime matches the analytic shunt
  `Z/(1 + i omega C Z)` to `O((omega dt)^2)`. Realization doctrine
  (the recorded bake-off): at a brass input plane the higher modes
  are deeply evanescent, so the MM-derived `R_00` — which already
  contains every interior mode-conversion round trip — is what a
  plane-driven source sees; the receipt (ml-002) measured the matrix
  arm within 1 cent of the dominant line with a plane-mode energy
  fraction above 0.999 at 9x the cost.

- `brass_loop` module (bead 3ez8g.4.3): the composed brass playing
  loop — no `Trumpet` type. `BrassVoice` = a 1-DOF
  `fs_phs::mass_spring_damper` lip island (outward-striking; parameters
  at reduce-lab card scale with provenance labels on every input) x the
  per-valve multimodal characteristic lines x the D18 cup junction.
  PITCH IS NEVER ASSIGNED: the control surface is blowing pressure,
  lip tension, valve combination — structurally enforced by the
  battery's source grep. THE NUMERICS LESSON (executed): an explicit
  Bernoulli flow computed from the previous sample's bore pressure
  carries a one-sample delay whose phase flips the aperture's
  small-signal conductance ANTI-dissipative above fs/4 — a crook combo
  locked onto a 7.3 kHz parasitic scream; the fix is the closed-form
  IMPLICIT flow/junction solve (the linear cup relation substituted
  into the Bernoulli law gives a quadratic in the flow), leaving no
  delay in the flow path while the lip gap stays explicit (lip mass
  filters it). Passivity on the line taps is additionally enforced on
  a 4x-oversampled DTFT grid (inter-bin Gibbs overshoot makes a
  grid-only |R| <= 1 claim insufficient). Per-block diagnostics carry
  the lock estimate, gap statistics, blow/bore work, and the lip
  island's supply defect.

- `wind_line` module (bead 3ez8g.6.3): the wind articulation runtime —
  `WindLineBank` holds per-fingering exact-FIR characteristic lines
  over a typed `FingeringTable` with the carry-history + crossfade
  switch (the MM bank's D17 lift, plane-wave edition), and the
  char -> VFIT hop: a settled note may leave the exact-FIR phrase
  image for the cheaper rational hold image via the REPLAY-PRIME lift
  (the incoming IIR line is fed the recent outgoing history, outputs
  discarded, before the crossfade), hopping back on gesture
  resumption. Hop-readiness is GATED: entering an image whose
  registry gate is not green refuses structurally
  (`ImageNotGated`). "When to hop" is MEASURED DATA
  (`data/claims/wind-hop-policy.tsv`: settle-detector parameters +
  click-vs-hop-timing table, minted by the fixture, gated on every
  run), not doctrine. Fundamental estimation in the fixtures is
  AUTOCORRELATION — three estimator traps were executed en route:
  a global FFT peak returns the dominant harmonic; even-factor HPS is
  structurally blind on odd-harmonic (closed-pipe) spectra (the
  |X(2f)| factor lands on the absent even harmonic); odd-factor
  products cannot disambiguate f0 from 3f0 (the odd series is
  self-similar under 3x). The period domain has neither ambiguity.

- `piano_vertical` module (bead 3ez8g.5.1): the piano vertical
  composition — still no `Piano` type. Three unison strings (exact-ZOH
  modal images whose mode series carries the stiff-string
  inharmonicity law `f_n = n f0 sqrt(1 + B n^2)`; tension is state) +
  a duplex segment + a small modal board coupled through the
  POWER-CONJUGATE bridge (`F = sum c_k v_k` in, `g_k = -c_k v_b`
  back) + the 87zbd felt hammer island + pedals as coupling states
  (dampers = viscous drag through the strings' own force ports; una
  corda = struck-string count). The Weinreich aftersound EMERGES:
  in-phase unison motion pumps the lossy board (fast early decay);
  the surviving quasi-antisymmetric configuration barely couples
  (slow aftersound). MEASUREMENT LAWS learned executing the gates
  (each an executed misread first): the board-velocity spectrum
  carries the BOARD's own lines — the string law must be read from
  the string state; the board envelope is structurally BLIND to the
  antisymmetric survivors — the two-stage fit runs on the string
  ensemble energy; strong detune SUPPRESSES energy exchange (real
  physics) — the beat gate reads the AUDIO envelope, whose chorus
  modulation persists; and envelope smoothers must be scaled between
  the coupling slosh (~33 Hz, measured) and the detune beat (2-4 Hz).

- vibration gates (bead 3ez8g.7.1): the string/plate/nonlinear
  registry review. The TRUNCATION GATE makes modal honesty executable:
  a row's retained series must have converged (top-retained-mode energy
  share below the authored gate), and the deliberately under-truncated
  falsifier (N=3 at 220 Hz, whose top mode sits mid-band where the felt
  hammer still delivers ~1e-1 of the energy) FAILS the same gate the
  disclosed N=12 fixture passes at 6e-20. The linear<->nonlinear
  SELECTOR THRESHOLD is committed data with fixture provenance
  (`data/claims/vibration-selector-thresholds.tsv`), re-derived every
  run from the KC glide law the fs-nlmodal battery pins at 1e-12 — a
  hand-edited threshold fails. THE GONZALEZ CAVEAT travels with the
  KC/von-Karman promotion: conservation under discrete-gradient
  stepping is structurally blind to gradient errors, so nonlinear
  evidence cites FD-gradient/trajectory oracles, never energy
  conservation.

- Dispersive-waveguide bake-off (bead 3ez8g.7.2): the [F] string
  alternative measured against modal ZOH on ONE shared music-wire card
  with the analytic stiff law `f_n = n f0 sqrt(1 + B n^2)`,
  `B = pi^2 E I/(T L^2)` derived from the card, as the oracle. The
  executed receipt decided KeepBoth: modal keeps every certified claim
  (D21 asserted against the live registry in the gate test); the
  waveguide earns the budget-constrained hero-string scope
  (B_hat/B = 1.014, worst partial 9.8 cents with the linear
  fractional-interpolation phase as the named residual owner, ~10x
  fewer flops per sample, state cost flat in partial count — while
  modal holds FEWER states on this card, recorded honestly). The
  DISPERSIONLESS control on the same loop misses the law by 41 cents:
  the allpass stage, not the loop, carries the stiffness.

- Xylophone bar (bead 3ez8g.12.1): the cheapest FULL filling — a pure
  composition of live parts (analytic free-free bar modes + exact-ZOH
  runtime + a Hertz mallet island; no new physics, only binding and
  gates). THE BOUNDARY CONDITION IS LOAD-BEARING: the free-free
  `beta L` values are the `cosh cos = 1` roots computed in-fixture by
  deterministic Newton (self-verified at machine zero, never
  transcription-trust), giving the non-harmonic ladder
  f2/f1 = 2.7565, f3/f1 = 5.4039; the executed falsifier proves the
  gate discriminates (the PINNED n^2 family measures f2/f1 = 4.000 and
  fails). Material: the Indian-rosewood matdb pack (E_L = 13.53 GPa;
  density = SG*1000 as an authored Estimate, disclosed). The named
  successor is the ARCHED bar (ratio tuning toward 3-4, a
  varying-section solver) — v1 claims the uniform prismatic bar
  exactly and only that.

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

`src/mm_line.rs` `mm_line_tests`, cases ml-001..ml-006 — held-fingering
line fidelity vs the MM authority (0.39 cents worst over three valve
combinations; the explicit scattering LOOP additionally carries its
one-sample junction skew, measured 38 cents at a ~127-sample period and
gated at a disclosed 45), the realization bake-off receipt, the valve
no-click gate (carry+crossfade within 3x the control's slew; the
hard-cut falsifier is at least 2x worse), per-line passivity on the
stored taps, the cup shunt vs its analytic transform, and
causality/splice/refusal arms (the causality gate is the executable
witness of the e^{-i omega t} conjugation discipline).

`src/brass_loop.rs` `brass_gate_tests` + `brass_observer_tests`, cases
bl-001..bl-008 — the EMERGENT gates: lock 238 cents above the
cup-loaded column peak (outward-striking band) and 281 cents from the
lip alone at 10.4 kPa; tension walks the lock 208.9 -> 445.1 Hz across
distinct column peaks; valves drop the whole series (-125/-200 cents,
monotone, no pitch anywhere in the inputs); slotting (+78 cents of
continuous lipping inside a slot, 1310-cent jumps between); the cold
horn plays -50.4 cents flat vs the -48.6 sound-speed prediction with
temperature entering ONLY through the gas card; the no-click valve
change (envelope held through the switch; the fresh-voice hard-cut
falsifier collapses to 40% RMS); structural honesty + refusal arms;
and the radiated-observer demonstration from the SAME bake's SH
directivity table at the lock frequency (compact-mouth near-omni is
the correct physics at ka 0.09).

`src/wind_line.rs` `wind_articulation_tests` + `wind_hop_tests`,
cases wa-001..wa-004 + the committed-policy gate — the sigma slur
(+265 cents with the envelope held through the change), the register
vent (MEASURED, not the textbook sketch: with an INERTIVE 1.8 mm /
6 mm-chimney vent the twelfth-class flip happens at 0.70L, x3.36,
18 cents from the TMM authority's own peak, while the 0.30L vent
leaves the lock in the low regime — position matters AND the lock
follows the impedance peaks), the structural gate refusal, and the
hop policy (settled hop click 0.504 vs the early-hop falsifier 0.864;
the note survives on the hold image; the committed policy artifact
re-gated every run with the measured nuance that the settle detector
is conservatively late — the right side to err).

`src/piano_vertical.rs` `piano_vertical_tests`, cases pv-001..pv-005 —
THE TILT CONTRAST on the composed vertical (felt trend 159.7 vs the
matched linear spring's 1.000: the hysteresis is audible, not
asserted); the inharmonicity law measured in the audio (worst 4.23
cents over partials {1..5,7} with the B=0 falsifier missing by 15.04;
partial 6 disclosed as the duplex collision line, partial 8 as
weakly-excited at the strike point); the two-stage decay (string
ensemble 3.28/s early vs 1.15/s late) with audio-envelope beats
scaling with detune (0.370 s -> 0.239 s at doubled detune); pedal
topology (dampers collapse the note; una corda energy ratio 0.67 with
the third string still fed through the bridge); and the dissipation
ledger + bitwise replay (worst per-window growth 0.0).

`tests/vibration_gates.rs`, cases vg-001..vg-003 — the truncation gate
with its executed under-truncation falsifier (6.2e-20 vs 1.0e-1 across
a 2e-3 gate), the selector-threshold re-derivation pin (relative 1e-8
against the committed artifact), and the gate-summary enumeration
(every vibration-filling registry row must appear in
`data/claims/vibration-gate-summary.tsv` with a status MATCHING the
registry; the plate 1-port's honestly-ungated row is asserted by
name).

`tests/bakeoff_dispersion.rs` — the committed-receipt gate (bytes
decode, hash, KeepBoth outcome, D21 modal-row-untouched against the
live registry, the newcomer's own row present) and the
passivity/control gate (loop bounded and decaying under g < 1 with the
unit-magnitude stage; constructor refusals; the dispersionless
control's partial 8 sits within 20% of harmonic while the law demands
+18 cents).

`tests/xylophone_bar.rs`, cases xb-001..xb-004 — root self-verification
(residual at machine zero, classical table at 1e-8, asymptote check) +
named refusals; the rendered ratio gate within 1% with the executed
pinned-family falsifier; strike-position physics emergent (mode-1 node
found in-fixture at 0.2242 = the cord-mount point; striking it kills
the fundamental ~8 orders, the center strike kills partial 2, and the
non-noded partials keep speaking); and the listening chain gate
(WAV -> sidecar -> receipt, Unadjudicated).

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
- The multimodal line inherits the fs-duct modal image's disclosed
  boundaries (matched-mouth higher modes, no tone holes, no mean flow)
  and adds its own: the explicit scattering junction carries a
  one-sample loop skew shared with the reed voice (disclosed, measured,
  not a line error); a hard (fade = 0) valve swap is genuinely
  discontinuous and is the caller's choice; the crossfade lift is an
  engineering lane, not a claim about real valve port geometry; and
  the realization cost is O(n_fft) MACs per sample per line (the
  matrix arm is N^2 of that) with no real-time claim until a budget
  row exists.
- The brass loop is Estimate end to end (authored card-scale lip; the
  claim decision is the brass-gates bead). Its own v1 boundaries: the
  lip pressure-collection face area is authored ABOVE the
  orifice-strip area (mouth-side lip face; disclosed in provenance);
  no WallPin lane through the multimodal image yet; the radiated
  observer is demonstrated per-direction at the lock frequency, not as
  a broadband convolution (the broadband lane is
  `broadband_radiation`); valve changes are discrete combos with the
  crossfade lift — a trombone's continuous slide needs the wind
  epic's fractional-delay bead, a dependency of claim, not of code.
- The wind articulation runtime inherits the char image's plane-wave
  scope; the hold-image realization here is a fixture-grade 12-pole
  fit (the certified hold lane is the clarinet casebook's); the
  register-vent finding is geometry-specific (an inertive chimney
  vent inverts the classic pressure-release node rule — recorded as a
  measured hypothesis, not a general law).
- The piano vertical's board is an authored spruce-scale modal set (the
  full fs-plate orthotropic chart with rib stiffeners is the
  piano-gates upgrade lane); bridge geometry is authored v1 (the
  ingest lane upgrades later); the explicit board coupling is
  semi-implicit at audio rate (the ledger gate bounds it); cabinet,
  room, and JCA linings stay out by the bead's own boundary. Claim
  promotion happens in piano-gates (.5.3), not here.
- The vibration gate review promotes claims on EXISTING evidence (the
  .1.4 bake-off receipt, the fs-plate Olson-Hazell battery, the
  fs-nlmodal analytic pins); it mints no new physics. The plate
  vfit-driving-point row stays ungated with the missing item named (no
  executed plate-1-port fit artifact); bowed and vibroacoustic rows
  carry evidence but await their own mechanism/corpus gates (.7.3-.7.5,
  .7.4).
