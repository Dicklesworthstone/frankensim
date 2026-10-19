# CONTRACT: fs-couple

M4.1-W1 extends `unilateral_contact` with explicit physical span measures [m],
per-length stiffness, static/kinetic coefficients and a declared friction ramp
speed (mandatory for nonzero kinetic friction). No 1/n normalization or hidden
1e-3 m/s scale remains in the adapter; constitutive refusals retain owner/field
in a typed DContactError. Loss-port Results propagate. `binding::bound_point_obstacle`
lowers SourcePointContact verbatim and is also used by MaterialContact.
`BoundNormalStorage` composes the existing point or finite-gap potential with
the pHS Storage contract, preserving original U/F/dF and source errors in its
fallible API; the infallible Storage view signals invalid trials as nonfinite
values for the numerical owner's refusal. Whole finite-gap footprints already
integrate their areas and are not multiplied by another area weight. Sources
must have an unloaded zero-approach origin. Point Hunt--Crossley loss is a
separate port evaluated at actual discrete effort and explicit dt; its envelope
is checked before projection. No time integration, cancellation framework,
material data, thermal evolution or interpolant authority is invented.

Multiphysics composition through a versioned port-thermodynamic vocabulary: a
lossless interface relation, explicit storage/dissipation/source primitives,
a complete stream-flux bundle, and evidence-bound closed-window physical
accounting without false passivity.

## Purpose and layer

Layer L3 (multiphysics coupling). Domain assembly depends on physical
solvers. The neutral vocabulary is owned by `fs-ports` and reexported here
without conversion or duplicate types. That lower boundary depends only on
`fs-qty` and `fs-iface`; FEEC imports it without depending on domain assembly.

## Public types and semantics

### Typed gesture consumer bindings

`gesture_path::bind_paths` binds each description-validated parameter path to
an explicit existing pressure voice or actual material-string bow member.
Every source path requires a destination; incompatible/unconsumed paths return
`GestureError::NoConsumer { path }`. `BoundPaths` retains the original resolved
v2 schedule, alias-migration identity, description and exact destination map for
render provenance. V2 schedules remain typed throughout lowering and call the
original pressure/friction sampling loops without changing time rounding, ramps, station
checks, integration, cancellation, force transfer or stored physical state.
Numeric u8 alias decoding is isolated in fs-scenario's compatibility module;
typed destination members use usize without a legacy u8 conversion.
The existing single-string, material-bank and simultaneous-bow entry points have
path-bound counterparts returning that provenance alongside their renderer.
Legacy file entry points remain v1 compatibility adapters; a typed schedule
requires description admission and explicit path bindings. No new consumer or
silent fallback is invented for topology events, frets, valves, pitch, prestress
or terminations. G0 `gesture_schedule::unconsumed_targets_reach_a_typed_consumer_or_refusal`
checks every alias and bitwise pressure-control equivalence; original render
equivalence tests remain unchanged. Addressing alone establishes no physical
validation or conservative topology-change claim.

### Declared per-note control resolution

`description::per_note` resolves `MinimumThresholdPressure` using a model-only
monotone onset predicate and explicit SI domain, pressure tolerance and maximum
evaluations. Its API takes no audio/reference input. The non-oscillating lower
and oscillating upper bracket, every probe, model identity, objective, target,
note, domain and budget enter a deterministic content-addressed receipt.
`record_resolution` binds the receipt to the corresponding Performance control;
`fs_scenario::performance::lower` retains it in render provenance. Other declared
objectives remain representable but refuse this pressure-specific resolver.
Model identities must bind initial conditions and the predeclared observable;
evaluation cancellation is the model's responsibility. The search is bounded
and synchronous and publishes only a completed result. Monotonicity is an
explicit model obligation; a sampled bracket is not a proof of global
monotonicity or asymptotic bifurcation onset. G1 `per_note_rule::pn_001` uses
the existing reed fixture's 0.25 s sustained-pressure observable with 18.75 Pa
bisection width; G0 `pn_002`–`pn_004` cover API shape, provenance and refusals.
The authored four-note `music_lane.sh --fixture midi-reed` exercises MIDI,
lowering, model-only resolution, schedule rendering, replay and named refusal.
This fixture is plumbing evidence and makes no reed realism claim.

### Port-schema adoption contract

`PortedModule` exposes stable local schemas, coordinate charts, internal N-way
junctions and the last accepted physical-time transfers for CoupledModalSystem,
SharedCavity, GasPiston, DoubleActingGasPiston, DynamicAperture, ModalAperture and
the ReedBoreVoice characteristic inlet. All schemas pass PortSchema::try_new.
Graph embedding must namespace local identities. The chart timestamp names the
clock origin; each audit binds its own exact accepted interval. Held effort work
uses the original conjugate displacement increment, never logical ticks as time.

Adapters observe the existing solvers; they do not feed back into their arithmetic
or acceptance. Typed transfers and window evidence are materialized lazily on
request. Exact-network connection work is minus the spring-energy change and
the independently integrated damper loss along the same full trajectory, never
endpoint reaction times displacement. Multiple damped columns require bounded
per-column integration scratch; this is not an allocation-free observer.
Rejected trials preserve the preceding account. An independent observation
refusal is retained as `ModuleObservationFailed`, and observational overflow as
a typed transfer error; neither publishes an invented zero record or rewinds
an already accepted step.

SharedCavity retains its mass-normalized modal chart and signed volume projection.
The canonical schema is the explicitly transformed reference-coordinate pair:
F = f_q sqrt(M_ref), v = q_dot/sqrt(M_ref), with M_ref=1 kg, not a specimen modal
mass or a point-force measurement. Raw owner units remain N/sqrt(kg) and
m sqrt(kg)/s. Its compression junction uses the original signed volume pullbacks
to pressure/volume-flow coordinates, including neck and compression storage.
Each CoupledModalSystem column similarly names its B_i^T q_i work charts and
one connection-storage endpoint. DoubleActingGasPiston's two gas faces share a
single common-flow wall junction; internal face work is not an exterior supply.

Available energy windows use the owner's actual storage change, prescribed work,
loss and existing joule tolerance. No balance is reconstructed as a leftover.
Piston heat ports use the source-backed discrete gas-midpoint temperature and
Q/T_midpoint entropy quadrature; gas entropy change remains independently
measured. This energy projection claims neither gas species/element balances,
exact isentropy nor external-reservoir thermodynamic validation. Negative
numerical entropy generation remains a visible WindowAuditReport violation.

MagentaBadger's ruling in Agent Mail 2636 explicitly leaves two audits unavailable:
DynamicAperture has no owner-declared energy tolerance, and ReedBoreVoice /
DelayedFilter has no physical line storage/loss owner or joule tolerance. Their
real transfers remain available, but ModuleWindowAudit::Unavailable names the
missing inputs. Those gaps keep IRP .162 open; unavailable is never a passing
window. No stepping or music-lane golden is intentionally changed.

G0/G3 tests in tests/port_schema_adoption.rs cover admitted schemas for all seven
owners, available/unavailable accounts, sign-flip refusal, two-port parity,
three-component/shared-cavity N-way descriptors, the common-wall descriptor,
and accepted-record preservation on refusal. Each module emits a native-run
JSONL row with residual/tolerance or explicit null/unavailable fields. Native
verification and reed/string baseline hashes remain central-run requirements.

### Integrator doctrine and current audit boundary (IRP M5.2-W1)

Render-path linear parts must use exact held-input propagation or explicitly
frequency-matched energy-consistent propagation. Nonlinear storage retains its
discrete-gradient owner; generalized-alpha is not admitted on render paths.
Held effort exchanges use the same displacement work at both ports, rate
changes require antialias filtering, and phase/contact/solver budgets must be
admitted explicitly. A midpoint path that remains is admitted only inside the
fs-phs one-cent band (`admit_midpoint_band`, pole-rate bound times dt <= 0.04);
the table below says which paths are exact and which are banded.

`tests/integrator_warping_audit.rs` contains the current path inventory and G1
one-step-map measurements at 1/5/10 kHz and 48/96/192 kHz. The common fs-phs
linear owner is measured for its callers; that is not a complete nonlinear or
connected render. Public ImpactSystem and Dynamic/Paired/ModalAperture seams
are exercised directly, including the uniform scalar condensation. The exact
modal owner is the control. Midpoint phase is compared to
`2*atan(omega*dt/2)/dt` within 1e-6 cents. The closure test
`every_render_path_is_listed_with_integrator_and_measured_warping` accepts
a row only if it is exact, admits no more than the one-cent band, is the
first-order network row with its measured order ladder, or is a PianoAudio
factor whose companion measurement stays within 1 cent. A path outside its
band refuses with a typed `TimeResolution` error; nothing is retuned.

| Current path | Linear owner | Current band / audit boundary |
| --- | --- | --- |
| KC string (`realize_coupled_kc`) and nonlinear plate (`thin_plate::VkBody`) through fs-phs `step_render` | Exact ZOH when the storage is exactly quadratic, otherwise Gonzalez midpoint | Midpoint admits only the one-cent band on the split's pole-rate bound and refuses with `TimeResolution` |
| Dirac joins (`realize_dirac_join` / `step_descriptor`) | Gonzalez midpoint | One-cent band on the descriptor pole-rate bound; typed `TimeResolution` refusal |
| Glottis fold island, PlateBank cavity held-flow island, coupled ODE chains, blown and reed ODEs, brass lip held-force island, plate-cavity transformer | Exact ZOH (fs-phs `linear_exact`) | No warping; nonlinear forces are held over the step |
| TsDriver speaker coil/cone/panel | Full coupled exact ZOH | Held voltage, independently integrated current/work/loss; dense reference step, not an allocation-free callback |
| ImpactSystem | Gonzalez midpoint | One-cent band on the complete retained block (modal stiffness, volume-spring dyads and resistance), typed `ImpactError::TimeResolution`; the per-mode 0.9*pi guard remains a separate validity check |
| ModalAperture, uniform condensation, and Dynamic/PairedAperture with `BladePropagator::Midpoint` | Midpoint response | One-cent band (pole-rate bound times dt <= 0.04, DEF-04) refuses at construction and at every step; the material-modal adapter additionally caps its caller's bound at 1 |
| PairedAperture with `BladePropagator::ExactZoh` | Exact ZOH per elastic blade (scaled exponential; Van Loan viscous loss) | No warping at any omega*dt (audit rows ~1e-11 cents up to omega*dt 1.31). Contact/Bernoulli forces are held over the step and solved at the midpoint opening; blade material memory refuses |
| DynamicAperture with `BladePropagator::ExactZoh`, authored or plate-bound (closure, Maxwell arms), and every plate-valve performance file | Exact ZOH of the opening plus its Maxwell arms (scaled exponential; Van Loan viscous and arm loss), through `reed_bore::contact::step_with_opening` | No warping (audit rows <= 3e-10 cents up to omega*dt 1.31); arms checked against an independent RK4 integration. Contact/Bernoulli forces held over the step and solved at the midpoint opening |
| CoupledModalSystem | Exact ZOH of the full retained bilateral pencil (component modes, connection springs and dashpots), split into connected K/C islands (cd7f65bc7); nonlinear contacts keep held discrete-gradient reactions around the same affine response | Coupled algebraic Nyquist bound; measured against the exact two-body eigenproblem; the two-body sticking traction matches an independent exact-flow integration within 1e-10 N |
| ModalAcousticTimeModel | Exact ZOH | Caller-declared Nyquist guard; no frequency-warping claim for its observer transfer |
| Source-backed GasPiston | Midpoint gas/momentum | One-cent band on the gas-spring rate sqrt(rho c^2 A^2/(m V)) at the start, midpoint and end states; typed `PistonError::TimeResolution` |
| Circuit/thermal/network piston | First order in general | Observed order 1 +- 0.1 on the four-rung companion ladder (`FS_IRP_NETWORK_AUDIT_LOG`); no exact-linear claim |
| Live PianoAudio at 1/2/4x | Exact local blocks plus finite-step coupled solves | Companion-measured coupled-board frequency within 1 cent, window discrepancy included, at each factor (`FS_IRP_PIANO_AUDIT_LOG`) |

The coupled-mode test uses the exact symmetric stiffness eigenproblem and a
three-rate phase-error ladder. The coupled pencil steps exactly, so the scalar
midpoint atan oracle does not apply to it. Piston tests reuse the actual
source-backed thermochemical and thermal/electrical fixtures. The speaker
migration has a held-voltage time-partition regression and a public-state pole
audit. The audit and its companions pass natively (M5.2-W1 closed on ts1 b46 at
a4b6c4c22); consumer fixtures that request clocks outside a midpoint band
refuse with TimeResolution and must move to a resolved clock or an exact block.

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
`held_force_for_port_velocity` solves the exact-ZOH scalar endpoint constraint
using free modal motion plus unit-port-force mobility, with conjugate force
and velocity projections. It refuses invalid/nonpositive mobility and
nonfinite inputs/results without changing state. The ordinary step still
admits state and energy budgets. The bowed runtime's set-valued contact uses
the same free-motion plus unit-force endpoint response, cached per bow station
(see the Bow bullet under the string paths).
G1 compares an independent oscillator solution
and demonstrates the old pinning error; G3 checks mixed-sign multimode ports
in all damping regimes and full/half/subsample durations. This enforces only
endpoint velocity: continuous no-slip motion, exact within-step event times,
and a passive contact discretization are not claimed.
Each frame additionally retains mode-local viscous work-minus-energy differences
and their scale-aware roundoff allowances. Negative loss beyond a mode's allowance
refuses the entire step, even when another mode could mask it in the aggregate.
The aggregate work/energy calculation is unchanged; summing the mode-local
differences can differ by floating-point roundoff. Neither an undamped mode's
energy drift nor a negative value within its allowance is physical heating.
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
apply only BETWEEN blocks, transactionally.
`SetModalForce` changes one mass-normalized held force [N/sqrt(kg)] in an
admitted modal voice, preserving displacement, velocity and its physical clock.
Zero force releases the input without silencing or resetting the resonator.
All voice/mode indices and force values in a batch are checked before any
input or control log changes. This exposes generalized excitation, not a
spatial pluck/contact law or an automatic score scheduler. Controls are logged as
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
refusal mid-block permanently poisons the context (mid-sample state is not
rewound). The failed block's output is partial and must be discarded.
Subsequent block or control calls return `RenderError::Poisoned` before
touching output, voices, or control logs; reconstruction is required.
Request sizing and control-admission refusals do not poison a healthy context.
`ScheduledRenderer` admits immutable `ScheduledControl` schedules.
`ScheduledRenderer::from_pressure_gestures` connects complete pressure tracks
directly to hosted voices, requiring their actual sample periods to match the
declared audio rate. It validates destination bindings even for an empty horizon,
requires context sample zero, and preserves supplied physical state. Compilation
and all admission finish before controls or simulation state advance.
The renderer
splits caller blocks exactly before their integer sample indices. Equal-index
events retain input order; end-boundary events wait for the next nonempty block.
All controls are validated and log capacity reserved before stepping. The
constructor refuses invalid schedules or budgets; no resampling or
score-to-excitation compilation is claimed. Its cancellation gate polls before each
caller block and consumes no events on cancellation. G3 tests compare scheduled
force changes, release and re-excitation against direct modal stepping across
multiple block partitions, with cancellation/resume at callback boundaries.
`pressure_gesture_controls` binds a selected blowing-pressure gesture track to
that scheduler. It samples existing gesture ramps at the declared control rate,
holds each value, and places ticks at the first audio sample at or after their
time using integer clock conversion. Admission bounds evaluated ticks and storage;
the caller supplies the destination voice's audio rate and a half-open horizon.
G3 checks a nonintegral clock ratio and pressure ramp against per-sample control
application across callback partitions. Other gesture targets and continuous-time
interpolation between control ticks are not implemented by this binding.
Boundary facts a consumer must
know, both discovered by this bead's battery: the characteristic-line
realization refuses `UnflangedOpen`/`FlangedOpen` terminations whose
Nyquist `ka` exceeds 1 (at 48 kHz that is every bore wider than
~2.27 mm radius — the measured-load `Termination::Tabulated` lane,
bead zolja, is the lift), and the Zwikker–Kosten Bessel path's former
`r_v ≈ 160` continued-fraction ceiling was an iteration-budget bug fixed
in `fs-phs` (budget now scales with `|z|`).

### Aperture-junction solver modes (bead frankensim-2s4i5)

IRP M5.5 adds `render::ladder::RealizedResolution` and
`ReedBoreVoice::{new_with_fft_length, realized_resolution,
set_junction_step_tolerance, junction_residual_max, initial_ir_tail_energy}`.
The constructor shares the existing `with_line` admission/initialization; it
only passes the explicit request to `characteristic_line_dense`. Actual FIR
history length exposes the 8192 clamp. Strict mode has **no** pressure-step
tolerance: it still terminates on an exact zero or adjacent binary64 endpoints.
The massless FastNewton seam admits an explicit relative tolerance in (0,1],
keeps 1e-9 as the existing default, and preserves every corner/fallback rule.
Massive reeds refuse this tolerance seam rather than reporting a tolerance the
coupled structural solver did not use. The observed volume-flow defect is
recorded in m3/s, and the initialized periodic FIR's last-quarter tap norm is a
dimensionless diagnostic. Neither bounds pitch, onset, or partial-level error.

`tests/convergence_ladder.rs` binds the existing M-CL0 fixture to fs-ladder,
with at least three measured sample-rate, FFT, and Newton-tolerance rungs and
the joint group at fixed physical IR duration. The legacy observer's backward
jet-force difference supplies the declared smooth-signal first-order hypothesis;
this is not a global-order assertion for the oscillator. The fixture measures
raw pitch/partials and a **finite-horizon onset** (0.25s, 1e-3 Pa RMS, nondecaying
last-quarter RMS, 8 bisections in the declared 0..4800 Pa bracket). An unbracketed
onset refuses. This is not an asymptotic bifurcation/extinction measurement.
G1/G3/G5 report construction is exercised, but Q2 accuracy authority remains
unavailable where no conditioning map supplies a bound in the acoustic QoI's
units. Spreads and signed finest-reference joint interactions remain data.
`scripts/e2e/music_ladder.sh` runs a prebuilt native test executable and retains
JSONL, a raw-scaled PCM render/provenance sidecar and its fs-blake3-bound report
in a fresh directory. It performs no build and never overwrites evidence.

Both modes share the closed-slit characteristic constraint
`(1-r0) p_plus - p_minus_history = Zc Ubody`. The reflection denominator
is not clamped, and closure is evaluated at the body-loaded bore pressure;
body motion can reopen the slit. Singular or nonfinite candidates do not
qualify for this shortcut. G1 junction tests cover signed body flow,
reflected incoming waves, and reopening.

The reed aperture junction has two solver modes on
`ReedBoreVoice::set_solver_mode`. `Strict` (default) is the
deterministic path: a local sign-changing bracket, up to 32 symmetric
bracket expansions if needed, then up to 2100 bisections (covering the
binary64 exponent range, including subnormals). The former 96-step bound
could refuse near-zero roots solely because their pressure scale was far
below the initial bracket width. It stops at an
exact residual zero or adjacent floating-point pressure endpoints,
returning the endpoint with smaller flow residual. Nonfinite residuals,
failure to bracket, or exhausted refinement return errors; a grid guess
is never returned as a solution. This removes the former dimensionally
invalid pressure-scaled flow tolerance. Root selection is deterministic,
but multiple quasistatic roots and the Bernoulli dead zone remain model
limitations; this is not a uniqueness or exact-arithmetic certificate.
G1 tests independently reconstruct volume-flow balance across body flow,
reflection, and large incoming pressure, and test refusal without a root.
`FastNewton` is a DECLARED FAST MODE (the
fs-rand ziggurat precedent): DAMPED island Newton on the analytic
Jacobian of the shared Bernoulli-plus-wave-junction residual — each
step is halved up to four times until the residual magnitude strictly
decreases, and a sample that never improves, hits a sqrt kink or
non-finite state, exhausts its iteration budget, or shows a vanishing
slope is handed untouched to the strict path and counted in
`FastSolveStats`. Convergence requires an UNDOUBLED full step inside
`NEWTON_STEP_TOL`, so a still-descending line search never reads as
a converged root. It is NOT bitwise-equal to the strict path and must
never become the default until it earns the same proof. Receipts (bead
2s4i5): the junction-level battery lives inline (`fast_mode_tests`) —
a 64-point open/interior/closing stimulus sweep checks both modes against
a volume-flow scale `w H sqrt(2 Pc/rho)` and checks pressure agreement
within `1e-8 * (1 Pa + |p_strict|)`. The integration lane
(`tests/reed_newton_fusion.rs`) owns the fallback-hit-rate receipt on
a full nominal render (bounded structurally at 0.5; measured 12-22%
warm-seeded) and before/after budget rows stamped with build profile —
debug rows are diagnostics forever, per the budget-lane doctrine.
Recorded finding: end-to-end rendered waveforms are NOT comparable
between modes beyond onset — the self-oscillating loop amplifies
per-sample root differences into phase-divergent (O(kPa)) waveforms,
so waveform equality is never a claim either mode makes.

### Material-driven viscous slit flow

`bernoulli_aperture::viscous::ViscousSlit` owns the shared reduced channel law
`dp = 12 mu L U/(w h^3) + rho U|U|/(2 w^2 h^2)`. `from_gas` retains one
coherent `GasState`, derives viscosity from it, and refuses a different density
at evaluation. `PairedAperture` and `ModalAperture` can select that law from
their retained material gas plus the physical channel length. `DynamicAperture`
accepts the same source-bound law; a plate-valve performance's optional
`slit_length` record uses its actual regional inlet gas, not the exterior gas.

All three mechanical owners use this flow inside their same-sample pressure,
contact and swept-volume solve and account for the channel loss as `dp * U`.
Independent blades and signed modes retain their own conjugate pressure areas;
the shared slit width is counted once. Spatial single/symmetric/independent
plate closures integrate the same law analytically along each positive P1 edge
segment, with an analytic pressure tangent. Closed portions contribute zero.
An area-equivalent height is never substituted into the cubic viscous term.
Spatial pressure condensation uses safeguarded Newton with a passive bracket,
96 evaluations and a 32-epsilon relative pressure-balance target; failure
refuses the trial before any mechanical, material-memory or wave state advances.
Law attachment is allowed only before the first accepted step.

G1/G3 tests in `viscous_aperture`, `modal_aperture::viscosity`, the paired spatial
module and the constitutive owner cover flow/work closure, material-gas swaps,
the cubic-gap laminar limit, independent point quadrature, partial closure,
edge subdivision, paired/modal equivalence and full modal damping. This remains
a fixed-state, parallel-channel, quasi-steady model with one common streamwise
length. It does not resolve transverse/entrance flow, channel inertance,
pressure recovery, gas evolution, compressible flow or measured instrument
accuracy; it is not the plan's resolved-flow acceptance rung.

### `music_render` (binary)

The music lane's render CLI (bead ib15w): `music_render <reed|string>
<out.wav> [--seconds S] [--block N] [--full-scale-pa P]`. Pinned fixture
compositions render through the block API and encode through the ONE
pascals→PCM owner, `pcm_wav::encode_pcm16_wav` — the recorded seam
decision (beads ib15w + h7xu5.7.8): the cinematic stereo/receipt-hashed
encoder stays cinematic, no third RIFF writer exists in the music lane.
Reed runs accept `--temperature-k K` and `--ambient-pressure-pa PA`; these
feed the shared `GasState` into the existing duct/reed calculation.
Optional `--relative-humidity RH` supplies a dimensionless fraction in `[0, 1]`
to the shared `GasState::try_new_moist_air` mixture owner. Default zero follows
the exact dry-air path. Nonzero humidity changes the gas properties consumed
by the duct and reed calculation; it is not an output pitch/gain adjustment.
For nonzero humidity the owner refuses temperatures outside its saturation-fit
window and excessive vapor fractions. This does not model condensation or
moisture uptake in solids. The ambient values, relative humidity and gas
specification are recorded in the sidecar. Other fixtures refuse these
overrides. This is a stationary gas-state change, not
solid thermal expansion, heat transfer or phase change. Reed runs also accept
`--duct-length-m M` and `--duct-radius-m M` for the cylindrical duct, feeding
the declared geometry to the existing acoustic solver and recording it in the
sidecar. These dimensions must be finite and positive; other fixtures refuse
them. Optional `--duct-outlet-radius-m M` selects a linear cone when different
from `--duct-radius-m` (the inlet); omission or equal radii retains the cylinder.
Both cone endpoint radii are recorded for replay. Other fixtures refuse
them. Existing acoustic admission also applies: a positive radius can still
exceed the selected radiation model's `ka` limit and must refuse.
Geometry remains fixed during the performance. No automatic wall-material
response is implied. Gas-model validity
limits do not certify the fixed solid specimen at extreme temperatures.

The file-driven path also accepts `frankensim-flue-performance-v1`:
`music_render wind crates/fs-couple/examples/flue.performance OUT.wav --block 512`.
The fixed-order SI records explicitly supply mechanical/output clocks, the output
sample window and pascal scale/ceiling, bore and exterior gas, cylinder dimensions
and prescribed terminal reflection, mouth dimensions, jet card/gain/onset seed,
jet limits, wave memory, fractional phase band/budget, optional viscothermal cells
and band, receiver position/resolution/band, and control compilation limits. The
suffix is the existing canonical, typed blowing-pressure gesture schedule.
`jet_card tonal-interim` explicitly selects the retained stage-I source card;
`jet_card json {CANONICAL_JSON}` supplies a complete one-line card. No filename or
instrument name selects physical dimensions, material properties or calibration.
`gas` and `exterior` accept the same legacy, named-mixture and source-backed NASA-9
records as other wind sources; `exterior same` shares the supplied bore gas. A
separate exterior must have the same static pressure. Gas changes enter jet speed,
Reynolds number, propagation, selected viscothermal loss and physical radiation.
The optional `windway none`, `windway slit LENGTH_M` or
`windway unsteady LENGTH_M TRANSVERSE_CELLS CONVECTION_CELLS` record after
`flue` selects the shared ideal, quasi-steady or unsteady source-gas windway
before performance setup. Actual channel length, model, viscosity and resolved
cell counts remain in the output sidecar. Omitted/`none` retains the ideal jet.
A finite channel is explicit physical geometry; no lip-channel length is
invented for a transverse flute.

`bernoulli_aperture::slit_channel` resolves fixed, fully developed planar
Newtonian momentum from the retained `GasState` density and viscosity. Its
2..=4096 cell-centred finite volumes share a no-slip wall diffusion operator;
backward Euler evolves the transverse profile from rest under signed pressure.
An optional uniform kinetic-head exit is solved with the mean flow in the same
implicit update. Kinetic storage, pressure work, physical viscous heat, exit
work and numerical dissipation remain separate. Standalone refusal leaves the
accepted channel state intact. No fitted relaxation poles or alternate gas
formulas enter this owner.

`FlueVoice::with_unsteady_windway` connects that accepted flow to jet speed,
source volume, Reynolds coverage and the original card's growth law. Pressure
controls accelerate the fluid and release retains its momentum. A separate
2..=4096 cell semi-Lagrangian convection field transports the existing feedback
across the actual mouth cut-up at the evolving speed, including stop and
startup, without a minimum-speed history clamp. The voice's original failure
poisoning applies to the combined sample. Channel and acoustic mouth work
ledgers retain their distinct ports; numerical damping is not material heat.
Single-tube and regional-gas branched file sources use this same builder and
the same physical receiver. Independent steady Poiseuille and transient
sine-series regressions exercise space/time refinement, signed release,
reversal, gas interventions and callback partitioning; native execution is
required before treating them as passing evidence.

These are explicit numerical axes for a fixed, incompressible, fully developed
planar channel and the original estimated jet card. Side-wall/entrance/axial
profile development, compressibility, thermal transport, moving walls and
turbulence remain outside the channel. The exit uses uniform bulk kinetic
head despite the resolved internal profile. This does not establish a resolved
jet, moving-reed channel or experimentally validated instrument claim.

This loader constructs the existing `FlueVoice::from_fractional_network` and
`FlueOutletPerformance::with_multirate_radiation`. Every mechanical outlet-flow
sample reaches the shared differentiating Rayleigh receiver before the sole PCM
owner. The file declares 48 kHz CLI output and an integer mechanical/output ratio;
`--decimate` refuses a second observer. The sidecar retains both clocks, original
input/card identities, source gas coefficients and the actual receiver delay.
This is an **EstimateOnly** composition of a tonal jet and stationary rigid
cylindrical pipe with a prescribed reflection and one-way infinite-baffle outlet
radiation. The selected card keeps its narrower laboratory scope. It supplies no
resolved turbulence, matched radiation feedback, measured flute geometry, wall
material response or calibrated SPL claim; `--acceptance` refuses this source.

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

`broadband_radiation::modal` and `modal_acoustic_time::broadband` add an opt-in
one-way pressure receiver over complete mass-normalized modal velocity fields.
Computed training/independent held-out samples use exp(-i omega t), with explicit
producer identity, SI units, pure delay and producer discretization estimate.
Continuous vector fits are conjugated for the fs-vfit owner; proper stable
state-space filters join each original mode in an augmented held-force ZOH
exponential, reusing fs-phs's checked deterministic exponential. There is no
Tustin warp or velocity difference. The original modal owner still advances
mechanical states and work/loss records, preserving those bits exactly. Its
legacy pressure coefficients are zeroed only in the observer's private clone,
so the admitted receiver owns the pressure ceiling. A failed step publishes no
filter/delay/mechanical state. Block render failures poison their stream wrapper.

The pure delay uses exact integer taps or an explicit centered Blackman/sinc
fractional kernel; all taps must be causal within the physical delay. History
starts at zero. This approximates a bandlimited sampled pressure history, not
continuous retiming of all held-force spectral images. Held-out error includes
the actual kernel and retains EstimateOnly authority; it is not a continuum
certificate. Producer error remains separate. Constant samples with zero delay
use the exact legacy equilibrium-relative quadrature operation order.

Caps: all 4096 modal inputs are retained in <=256-input cold batches, each grid
has <=512 frequencies, <=4194304 transfer cells total, per-mode order<=128,
iterations<=128, delay-kernel half-length<=4096, delay history<=1048576 samples,
and total augmented matrices<=4194304 entries. No mode/pole truncation or implicit
budget enlargement is allowed. Host cancellation is at bounded pressure-block
boundaries; no hard-real-time, continuum passivity or load-model claim follows.

The v1 modal-file optional `broadband` block reuses its original parser and
ProjectedForceProgram; v2-v5 coupled/nonlinear trajectories refuse this
independent-mode receiver. CLI sidecars retain every producer/fit/error/delay
field, and the legacy path remains the default. The native G0/G1 battery is
`tests/broadband_modal_transfer.rs`: independent analytic held-input sphere
endpoint response at five off-resonance drives (1e-6 amplitude/radian bands),
constant compatibility, bitwise mechanical energy/work/loss, all-mode cap,
held-out refusals, integer delay, and three-rung delay/pole/grid ladders. The
pole/grid axes use the independent nonrational Bessel/Struve piston owner:
the delay-stripped sphere is already first-order rational and cannot establish
pole-count liveness. Numerical/native qualification is pending. The audible
replay/refusal sub-lane is `scripts/e2e/music_modal_broadband.sh`, requiring the
native qualification log and a supplied binary; it makes no instrument-realism
or general source-mortar/loading claim. M-BAR-Q2/M-GT-Q2/M-VNP-Q2 must supply their
own computed mode-basis samples explicitly before adopting this observer.

Generic offline bridge from solver-neutral sampled scalar-input radiation transfers to
real-tesseral SH filters. Radiation solvers supply complex-SH training rows, direct
held-out fields, source semantics, and diagnostics; `fs-couple` does not depend on them.
It converts `exp(-i omega t)` to `fs-vfit`'s `s=+i omega`, fits prewarped abscissae, and
applies unprewarped Tustin. Direct disjoint held-out fields gate artifacts; stable
proper fits may retain caller-admitted constant feedthrough but always forbid `s*e`.
The runtime owns visible per-filter state, superposes in fixed input order, and refuses
whole samples transactionally. EstimateOnly covers the fixed reference before `1/r`;
moving-source, FW-H, near-field, feedback, room, head, and passivity remain unclaimed.

### Material string bank rotational bridges

`bowed_string::runtime::reciprocal::bank` admits explicit small-angle
`StringBodyRotationConnection` laws between a sourced torsional string mode
and a body rotation. The same polar inertia and GJ/shear-viscosity sources that
drive bow surface motion supply the string angle map. A body rotation can come
from two retained bridge footprints in the original material plate basis and
their physical separation; this is a secant rotation, exact for rigid rocking,
with a spatial approximation for a deforming plate. It is not an independently
tuned body mode. Raw angular modal maps retain their producer's authority.

The existing coupled-modal owner receives the exact work-conjugate conversion
`x = r theta`, `k_x = k_theta/r²`, `c_x = c_theta/r²`, using the sourced outer
string radius. Reciprocal torque, rotational spring storage and dashpot loss
therefore enter the same solve as bending, all other strings and the shared
body/cavity. Its force limit bounds `|torque|/r`. Body acceleration supplies the
existing pressure observer after bow release; twist receives no direct audio
gain. The new binding must precede playback and string motion. Original source,
body, gas, pressure and interface bindings survive multiple-bow scheduling.
G2 checks an independent two-oscillator angular-spring solution; G3 exercises
torque-only release, shear and body-density substitutions, original spatial
footprints and joint bowing; G0 covers rotational preload energy and refusals.
This does not establish nonlinear bridge contact, moving endpoints, local
plate-slope convergence, large-angle torsion, exterior radiation loading, or
measured violin/cello/bass fidelity.

### `string_specimen`

`acoustic_realize::LinearMaterialStringRuntime` incrementally advances an explicitly
linearized pinned string from a `ResolvedStringSpecimen`, using the existing
exact-ZOH modes, loss admission and compact strip acceleration/jerk observer.
Material/operator state, modal displacement/velocity, observer history and sample
position publish at one epoch after all admission checks and a final `Cx` poll.
`step` returns a `MaterialStringFrame` pairing that epoch and the acoustic frame
with `StringDissipation`: material Kelvin-Voigt heat, oscillating-cylinder air
loss, and authored damping without an identified thermal destination. Coefficients
are retained when the modes are constructed and rebuilt with each rebind.
For constant additive viscous coefficients during a sample, each mechanism's
energy is its coefficient fraction times the mode-local viscous loss: all share
the same integral of squared velocity. This follows the viscous work law in
[MIT 8.01, section 23.5](https://ocw.mit.edu/courses/8-01sc-classical-mechanics-fall-2016/mit8_01scs22_chapter23.pdf).
The existing Rayleigh option still replaces the entire damping model; its loss
is authored, with no silently added air or material contribution. Tiny negative
mode-local differences and undamped drift remain in a signed roundoff residual,
not in heat. That residual must fit the combined modal/aggregate/reduction
roundoff allowance before the runtime commits. Diagnostic microphone pressure
adds no second radiation-energy debit. The returned partition does not itself
debit reservoirs or establish a closed thermal feedback/entropy model.
Refused steps and rebindings preserve the complete accepted runtime. Cloning an
accepted runtime retains the acceleration history needed for identical continuation.

`project_point_forces` maps prescribed transverse loads in newtons at material
stations to the current mass-normalized sine basis. `point_motion` uses the same
participation to return accepted displacement [m] and velocity [m/s]. The transpose
pair preserves power and held-load work; endpoint participation is exactly zero
for pinned supports. The resulting forces feed the existing mechanical and thermal
steps, including after refinement and material updates. Projections must be used
before the next epoch change. This does not supply contact laws, opposite-body
reactions, moving-load work or compliant-support dynamics. G1 checks an independent
physical oscillator and work; G3 checks multipoint power/work through thermal
updates; G4 checks endpoints, invalid inputs, overflow, cancellation and retry.

`initialize_point_load_equilibrium` connects the same physical force projection
to the modal owner's static solve `q_k = f_k / omega_k²`. It initializes only an
unexcited runtime before its first sample, establishes zero prior acceleration,
and advances the epoch without changing material state or the sample clock.
Holding those forces retains equilibrium; releasing them excites the same live
mechanics and pressure observer. Stored initial energy assumes a settled loading
history before the simulation window; it is not a simulated actuator-work or
contact receipt. Existing motion cannot be overwritten. G1 checks static compliance
and release response; G4 checks repeat, invalid-force and cancellation refusals.

`rebind_ambient` resolves a new prescribed gas through the existing gas owner
and rebuilds air drag and compact observer weights atomically. Material receipts,
modal coordinates, previous acceleration and accepted sample count are retained;
the epoch advances with zero parameter work because this model has no gas-added
mass or gas-dependent elastic storage. The thermal string exposes the same update
with expected-epoch checking and unchanged material enthalpy. Boundary heat must
still enter through thermal transport. This is a piecewise prescribed ambient
change, not an atmosphere transient, fluid-loading work or propagation-delay solve.
G1/G3 checks use cylinder drag and compact pressure formulas including retained
jerk history; G4 checks invalid gas, stale epochs and cancellation preserve retry.

Rebinding requires the exact same immutable card, length, linear density and
sine-mode count. This makes the mass-normalized basis identical, with no sorted-mode
matching, displacement reset or sound-level rescaling. Tension, EI, damping and
observer weights rebuild together; the source state remains inspectable by other
consumers. Signed parameter work is the change in modal elastic energy at fixed
displacement/momentum, checked against an explicit caller budget. It is returned
for the thermal/external parameter owner's opposite transfer, not booked as damping.
The axial prestress reservoir, heat transport and closed thermomechanical balance
are not supplied here. Coefficients are externally scheduled and piecewise constant.
The acoustic observer retains the preceding accepted acceleration across updates;
abrupt coefficient changes can produce a transient, without a crossfade.

Explicit `refine_modes` increases the retained sine basis at fixed geometry without
advancing time. Existing coordinates and observer acceleration history are retained;
added modes start at rest with zero previous acceleration. Shared admission enforces
mode, Nyquist and material-frequency limits before publication. The thermal wrapper
checks the expected epoch and preserves enthalpy. G1/G3 checks continuation energy,
pressure and independent forced response; G4 checks atomic refusals and cancellation.
This does not reconstruct omitted history or estimate truncation error, and does not
support coarsening or transfer between different geometries or eigenbases.

The rung omits transverse-motion-induced axial stretching, moving supports, extra
polarizations, Prony memory, basis/geometry transfer, phase change and radiation
backreaction. It uses prescribed piecewise ambient gas and fixed listener distance, with no atmospheric
path filter or common propagation delay. G1 uses a synthetic thermal-expansion card
and an independent damped-oscillator/pressure calculation. G3 cycles a same-card
density/temperature update at prescribed mass through changed cross-section, bending
stiffness and pressure without changing the retained sine basis. G4 exercises material-work
and late observer refusals, cancellation and exact suffix replay. This is focused
implementation evidence, not measured string fidelity or full MR05/MR34 completion.
The loss-partition G1 test independently integrates squared analytic modal
velocity and feeds only material heat into the existing equilibrium enthalpy
owner. Its synthetic constant-coefficient solid absorber checks the combined
mechanical/thermal/air balance and rejects duplicate heat, air credited to the
solid, and a missing mechanical term. It stays wholly on the solid branch;
thermoelastic feedback and phase dynamics are not tested. G3 changes loss laws
without resetting vibration and checks that the returned destinations change.
Existing G4 refusal/retry tests compare the complete frame, including loss fields.

`acoustic_realize::ThermalMaterialStringRuntime` owns that incremental string,
its immutable material card and a uniform equilibrium enthalpy state. It closes
the dissipative temperature/viscosity feedback for the prescribed-radius or
prescribed-mass, fixed-tension, constant-density, constant-Young's-modulus solid rung. Those
restrictions keep elastic and axial prestress storage independent of temperature;
there is no omitted stiffness-change or expansion work. Density must agree with
the solid branch of the thermal chart. The initial absolute temperature and card
must agree across both owners. Original query-axis schemas and selected claims
are retained; subsequent queries pin those claims, with no automatic substitution.
Constant properties and linearly interpolated temperature-only laws are admitted;
exact-only temperature tables do not authorize an evolving path between samples.
Viscosity must be nonnegative throughout its supplied source curve, including
interior knots whose values cannot be inferred from positive endpoints.
For `FixedMass`, the prescribed mass and `m/L` remain exact through thermal
updates; constant density and loaded length keep the resolved radius and basis
unchanged. Thermal area and volume use `ResolvedStringSpecimen::radius_m()` and
its cross-section, rather than reconstructing mass through a rounded radius.
This does not admit changing density, expansion or geometric state transfer.

Each sample advances mechanics at its initial viscosity, credits material loss
and explicit external heat to total enthalpy, then resolves the new material
temperature and rebuilds damping. Final pressure and retained observer history
use the updated acceleration. This is a first-order partitioned coupling; exact
mechanical substeps do not make the coupled solution exact. A positive caller
temperature-increment budget constrains every accepted step, and temporal
refinement must assess coupling error. Air loss and mechanical roundoff remain
outside material heat. The returned first-law residual includes work, both storage
changes, outgoing losses and numerical loss, with a scale-aware arithmetic
allowance including the enthalpy reference. It is not a physical-error certificate.
All candidates publish at one epoch after the final cancellation poll. Reusing
an accepted input epoch refuses before consuming another heat input or sample;
thermal/material refusal after the mechanical substep leaves both owners intact.
Cloning preserves exact continuation, including thermal state and source pins.
The direct `step` takes external heat in joules. `step_with_current` instead
resolves the existing `fs-material::conductor` from the specimen's pinned
electrical-resistivity claim, current temperature and actual length/area. An
ideal source deposits `I² R(T_initial) dt` once in the external heat account;
material vibration loss remains separate. The returned `ElectrothermalStringFrame`
retains initial/final sourced resistance and the same accepted thermal/mechanical
epoch. Final electrical resolution precedes publication; any refusal preserves
the accepted state and retry position. This is first-order partitioned uniform
Ohmic heating with insulated boundaries, not circuit dynamics, skin/contact
resistance, a charge-field solve, Lorentz forces or motional emf. Synthetic
G1/G3/G4 tests check analytical temperature refinement, current-sign symmetry,
resistance/damping/pressure response and atomic retry, not experimental fidelity.
`step_with_current_and_thermal_transport` composes the same current source
with the existing transport callback. The thermal owner receives vibration
plus Joule heat as internal deposition and returns actual boundary heat.
`ElectrothermalStringFrame<R>` retains Joule heat, boundary heat and the
concrete transport report separately; its coupled external-energy account is
their sum, so electrical work is neither lost nor counted twice. Both paths
share current admission, and final electrical resolution still precedes
publication. Zero current reduces to the ordinary thermal transport trajectory.
The independently prescribed fluid/enclosure reservoirs retain the existing
uniform-body, Biot and infinite-reservoir restrictions.
`step_with_voltage` and `step_with_voltage_and_thermal_transport` share these
same transactions under `OhmicDrive::Voltage`. Each sample freezes the current
`V/R(T_initial)` and deposits its Joule energy; the next sample uses the updated
resistance. The frame retains the authored drive, initial sample current and
current demanded at the final accepted resistance. Final electrical evaluation
must succeed before publication. Positive temperature coefficients therefore
reduce fixed-voltage heating while increasing fixed-current heating. The
partition remains first order; no within-sample circuit solve is claimed.
`step_with_thermal_transport`
accepts a thermal proposal through the shared `fs-material::phase` input/result
carriers. Passing `|input| environment.advance(cx, input)` with
`fs-conduction::lumped::LumpedThermalEnvironment` uses the existing
`solve_lumped_enthalpy` for one
backward-Euler step over the same audio sample. It derives the exposed cylinder
area `2 pi r (L+r)`, volume and Biot length `V/A` from the resolved specimen,
with no independently authored thermal mass or geometry. Convection and surface
radiation act on this same complete boundary; supports are thermally insulating.
The environment is a prescribed infinite, uniform reservoir and enclosing
radiation field with independent fluid and enclosure temperatures. Opposing
convection and radiation fluxes contribute their signed heat once. It neither changes nor debits the separately
declared acoustic gas. Transport may be declared or card-backed; the conduction
owner enforces its existing card identity and temperature-coverage rules.
The callback must not mutate or debit an external source before acceptance.
`step_with_thermal_state` also accepts caller-owned state by exclusive borrow.
Its callback reads that accepted state and returns an owned candidate alongside
the normal thermal proposal. Only after all coupled checks succeed does the
method replace the caller's state. A finite support therefore advances with
the string without a separate caller commit. G3/G4 finite-support integration
checks heat debit, temperature/damping/pressure feedback, late foreign-chart
refusal and exact retry of both owners. The transport owner validates its
opaque state; callbacks must not mutate shared interior state or external
resources. This is synchronous publication through exclusive borrows, not a
cross-thread or durable transaction.
`step_with_drive_and_thermal_state` provides the same publication boundary
for either `OhmicDrive::Current` or `OhmicDrive::Voltage`, reusing the existing
electrothermal path. External thermal state remains a candidate until final
resistance/current resolution also succeeds. The callback receives Joule plus
material heat; its boundary transfer excludes both. Finite-support tests check
this separation, current/voltage callback parity, cancellation and retry.
The coupling owner verifies that the proposal uses its original phase chart
and finite heat/residual budgets, then checks the coupled energy balance. This
keeps fs-conduction out of fs-couple's production dependencies: FEEC's optional
topology feature already depends on fs-couple. The cross-crate test has a
development-only conduction dependency.

Material loss enters once as mean internal power for that sample. Boundary heat
uses the actual endpoint convection/radiation powers, not a subtraction from
the computed enthalpy change that could hide thermal solve error. The returned
transport trajectory retains those powers, its Biot value and energy residual.
The caller's thermal energy-residual budget is checked separately and appears
as `thermal_solve_tolerance_j`, zero on the direct-heat path. Coupled first-law
acceptance uses that budget plus the reported arithmetic allowance; neither
certifies time-discretization error. A late thermal, phase or material refusal
preserves mechanics, temperature, source receipts, epoch and observer history.

The fixed corpus Biot ceiling is applied, including the conservative radiative
coefficient over the chart. This is a necessary uniform-body check, not proof
of initially uniform temperature or homogenized internal heating. Rapid local
sources and nonuniform support transfer require spatial transport regardless
of step size. The analytical convection reference and lumped assumptions are
described in [MIT's transient heat-transfer notes](https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node129.html).
The thermal chart retains its existing source authority. Finite-support contact
uses the conduction adapter's fixed resistance and uniform-body limits; spatial
conduction, entropy production certification, thermal
expansion, changing stiffness/density, and solid-liquid transfer are not claimed.
G1 compares the coupled trajectory with an independently integrated continuous
ODE and checks temporal convergence. Synthetic heat capacity and viscosity make
feedback observable in a bounded test; they are not measured material data.
G1 additionally checks analytical heating/cooling and first-order refinement
for two radii, plus separate and combined convection/radiation inputs changing
the material damping and diagnostic pressure. G4 tests Biot, thermal residual,
source-domain, phase and cancellation refusals, duplicate epochs, retry and
cloned continuation. These exercise the existing solver directly, not a mock.

`with_uniform_circular_material_state` lowers one homogeneous isotropic circular
specimen to the existing `PrestressedString` realizer. It consumes the shared
`ResolvedMaterialStatePoint` API, requiring only positive SI density and Young's
modulus. Area, area second moment, linear density, mass, EA, EI, and diameter
follow from the same radius and length. All source receipts and the complete
material state remain available in `ResolvedStringSpecimen`; its specimen hash
binds geometry and material resolution, excluding mechanical loading and solver
options. It is not the complete scenario identity.

The original binding holds the template's geometry and tension fixed.
`with_uniform_circular_material_and_prestress` additionally takes an explicit
`StringPrestress`: fixed tension, fixed extension from a stress-free length, or
target primary undamped fundamental frequency. Fixed extension computes `T=EA*strain`
and enforces the caller's declared linear strain limit; that assumption is not a
measured yield boundary. Target frequency inverts the existing pinned-pinned
Euler–Bernoulli dispersion relation, including EI. Slack/compressive solutions and
extension/target prescriptions on the moving-end basis refuse. The prescription
remains separate from the material receipts, available through `prestress()`.

`with_uniform_circular_material_and_constraints` also accepts a
`StringGeometryConstraint`: fixed radius or fixed total mass at the template's
current loaded length. Fixed mass derives `r=sqrt(m/(rho L pi))`, retains the
prescribed mass, and sets `mu=m/L` directly. It never reconstructs either mass
quantity through the rounded radius. Thus changing density at fixed mass/length
keeps exactly the same numerical modal basis while changing EA, EI and observer
width. The radius-derived area and `rho A` agree with these mass constraints to
floating-point roundoff. Geometry and prestress constraints are independent and both remain
available on the result, separate from material authority. The specimen hash
binds resolved geometry, mass and linear density, not the authored comparison
policy. The v2 identity distinguishes prescribed masses even when rounding maps
them to the same radius.

`with_uniform_circular_thermal_extension` consumes the existing
`IntegratedIsotropicThermalExpansion` and current elastic resolution from the
same exact card, query coordinates and axis conventions. It computes the
small-strain eigenstrain relation `T = EA * ((L-L_ref)/L_ref - integral(alpha dT))`,
using one reference configuration, as in
[MIT 16.20, section 3.7](https://web.mit.edu/16.20/homepage/3_Constitutive/Constitutive_files/module_3_no_solutions.pdf).
The current geometry/mass constraint remains explicit; both elastic and expansion
receipts survive in the result. `FixedThermalExtension` records the mechanical
prescription, while `thermal_expansion()` distinguishes the sourced path from
raw authored eigenstrain. A subsequent ordinary binding requires the thermal
path to be supplied again; it cannot silently retain old thermal evidence.
Total, thermal and positive elastic strains must stay within the authored limit.
Slack, compression and domain exit refuse; no buckling solver is implied.
Fixed tension remains a separate force-controlled prescription. G1/G3 check a
piecewise-linear expansion curve across temperature knots, heating/cooling force
changes, mass and receipt retention, and pitch from the real pressure solver.
G0 checks current-state/source mismatch, invalid reference length, unsupported
moving supports and strain-domain refusals. This initializes independent uniform
specimens; it does not evolve temperature, exchange heat, transport strain or
solve lateral contraction during playback.

All policies preserve the selected loss law, supports, and modal budget. These compare
specimens at resolved states; thermal evolution, transverse contraction,
yielding, anisotropy, or melting are not solved. Target pitch describes an undamped
linear eigenfrequency, not a damped or finite-amplitude pressure frequency. Diameter
supplies the existing compact strip-radiator width, without a new claim of accurate
circular-wire radiation. The density and modulus requirements use exact
dimension-only quantity schemas; equal-dimension semantic kinds or waveform
forms cannot be erased when lowering to the solver. Admission refuses missing
or quantity-schema-mismatched properties,
invalid mechanical inputs, and nonfinite/zero derived products before publication.
The finite scalar calculation is deterministic; no inner-loop database access is
added. `tests/string_specimen.rs` checks G0 admission, G1 mass,
EA/EI and beam modes, and G3 geometric/material scaling plus pressure-measured pitch
through the actual acoustic realizer. G1/G3 cases also cover extension-derived
tension and material changes under explicit extension/target-frequency constraints,
plus fixed-mass diameter/stiffness changes and their bending-sensitive pressure pitch.
Synthetic data check numerical behavior, not physical validation of a measured material.

`ResolvedStringSpecimen::with_kelvin_voigt_bending_loss` selects a constant
Kelvin–Voigt bending law using the retained `kelvin_voigt_bending_viscosity`
property [Pa s]. Density, Young's modulus and viscosity must be dimension-only
validity-wide scalar constants or temperature-only curves at fixed positive
absolute T. Frequency samples are not promoted to a broadband law. The
viscosity must declare a frequency domain, intersected with any
frequency restrictions on density and modulus. Other state coordinates stay at
their admitted query point. The descriptor carries `eta I`, the intersected band,
and the material bundle identity. Rebinding an already viscous specimen recomputes
all three from the new geometry and material. No source uncertainty is improved.

All three material loss binders (Kelvin–Voigt string/plate and Prony string)
use `ResolvedScalarProperty::angular_frequency_band_rad_s`. One static typed
frequency coordinate can have any name and be expressed in Hz or rad/s;
legacy untyped `omega` remains rad/s. The normalized band intersects every
participating coefficient's applicability before the existing modal checks.
Multiple frequency coordinates refuse because they may be independent physical
variables. No source receipts or interpolation decisions are rewritten. G3
tests compare actual pressure bit-for-bit across the two unit conventions and
legacy inputs; narrowing only density support rejects an out-of-band string.

The shared linear/Kirchhoff–Carrier/moving-end modal path projects
`mu y_tt - T y_xx + EI y_xxxx + eta I y_xxxxt = 0`, giving
`zeta_k = eta I k^4/(2 mu omega_k)`. Thus tension contributes stiffness but no
viscous loss. This is the constant-coefficient Kelvin–Voigt term in
[Sakthivel et al. (2023)](https://arxiv.org/abs/2301.07931), with tension added
and projected onto the existing sine/cosine bases. The explicit law replaces
the authored constant modal ratio; simultaneous Rayleigh
or nonzero authored damping refuses. Air resistance is still added separately.
Every retained reference frequency, including the secondary polarization, must
lie in the descriptor band. G1 checks `eta I` and rebinding; G3 checks actual
pressure decay against an independent diameter-form expression in all three
solver paths, including zero-viscosity air-only decay. G0 covers unavailable
data, wrong quantities, sampled frequency curves, duplicate losses and out-of-band modes.

This Kelvin–Voigt choice models linear bending strain-rate loss only. It does not supply axial
viscosity, material memory, amplitude-dependent loss, time-evolving material
state, thermoelastic/joint/radiation dissipation, or measured material validation.
The modal applicability check concerns reference frequencies; nonlinear
frequency shifts and generated harmonics are not an in-band certificate. The
existing moving-end cosine approximation and compact observer remain unchanged.

`with_prony_bending_loss` binds a causal generalized Maxwell bending spectrum
from explicit modulus [Pa] / relaxation-time [s] property pairs and
`equilibrium_young_modulus` [Pa]. `with_standard_linear_solid_bending_loss`
selects the single `relaxing_bending_modulus` / `bending_relaxation_time` pair
through that same binding. The explicit equilibrium modulus must equal the
Young's modulus used to derive EA/EI; a sampled dynamic modulus cannot silently
be treated as an equilibrium value. All selected coefficients, including density,
must be dimension-only validity-wide constants or temperature-only curves at
fixed positive absolute T. Every relaxation-time claim must
declare a frequency band, intersected with all other coefficient bands. An empty
spectrum is explicitly elastic and requires an equilibrium-modulus band. Geometry
or material rebinding resolves every retained source pair again, recomputing
each delta EI and replacing the source identity and time constants. Missing or
mismatched selectors refuse rebinding; duplicate source pairs refuse binding.
These are resolved inputs, not measured-material validation.

Each mass-normalized modal coordinate receives generic `fs-phs` relaxation arms
with `a_j = delta EI_j k^4/mu`: `q'' + c q' + omega0^2 q + sum_j a_j(q-v_j)=force`,
`v_j'=(q-v_j)/tau_j`. This is the internal-strain generalized Maxwell formulation
described in [Bleyer's viscoelasticity tutorial](https://bleyerj.github.io/comet-fenicsx/tours/nonlinear_problems/linear_viscoelasticity_jax/linear_viscoelasticity_jax.html),
projected onto the existing beam modes. Nonlinear axial storage, contact ports,
body joins and pressure observation use the same augmented pHS state. A held
pluck initializes every `v_j=q`, representing full relaxation before release.
Positive-strength branches add memory states in branch-major, mode-major order;
zero-strength branches retain source selection but add no state or active pole.
The owned spectrum is borrowed throughout the audio loop. These
material arms produce storage dispersion and dissipation together; no PCM decay
envelope or substitute modal damping ratio implements it.

Admission checks both relaxed and instantaneous reference frequencies against
the material band and Nyquist guard, using the sum of all branch stiffnesses
and including fixed-end detuned polarizations.
Moving-end secondary polarization refuses on this path. Kelvin–Voigt, Rayleigh
and nonzero authored internal loss cannot be combined with the Prony descriptor;
air drag remains separate at the relaxed reference frequency. G0/G1 tests cover
sample-clock admission (`dt/tau_j <= 0.1` for every nonzero arm and
`dt*omega_instant <= 0.1`, otherwise increase sample rate),
source/quantity/band refusals, geometry rebinding and actual pressure pitch/decay
against independent characteristic-polynomial roots in fixed linear, nonlinear
and moving-end realizations. Multi-arm G1 pressure is compared with independently
integrated hereditary-strain dynamics and checked for second-order refinement;
G3 checks equivalent branch splitting/permutation, zero-arm omission, detuned
polarizations, and body/contact composition in fixed and moving-end paths.
This remains linear bending memory at a fixed state. No measured spectrum fit,
axial creep, evolving temperature, phase change, arbitrary pre-release loading,
or in-band certificate for nonlinear harmonics is supplied.

### `driving_point::regional::field`

`RegionalDuct::from_gas_field` discretizes the shared `fs_material::gas_field::AxialGasField`
at explicitly refinable midpoint cells (at most 64 total), retaining the
existing cylinder/cone TMM and pressure/volume-flow interconnection. It adds no
propagation or loss law. The exact x=0 state feeds source density and the inlet
characteristic impedance, and x=L feeds the existing terminal radiation load;
an independent exterior can still be supplied explicitly. Uniform fields keep
the original complete geometry and arithmetic. Nonuniform tone holes and smooth
radius profiles refuse pending the corresponding continuous-field duct support.
`ReedBoreVoice::new_with_gas_field` feeds this same regional line into the existing
reed runtime and compact-jet observer. No instrument constants, tuning law,
evolving breath or calibrated exterior observation are added.

The focused field tests compare uniform responses bitwise, exercise cone
radiation and endpoint states, and refine a supplied 1 K axial temperature
gradient over 4/8/16/32/64 cells. Fundamental resonance convergence is second
order; the third-to-first register interval is checked against independent
energy-weighted perturbation. A rendered-pressure test checks uniform legacy
identity, exact source density and the gradient's nonzero audio effect. These
are numerical tests of prescribed fields, not measured instrument qualification;
native results are recorded separately.

### `acoustic_realize` / `pcm_wav`

`realize_assembly_with_gas_and_absorption` accepts a caller-resolved gas and
an explicit `GasAbsorptionModel`. Its observer always uses the selected law;
invalid ISO-air conditions refuse without a classical fallback. The original
entry points preserve legacy absorption arithmetic. `RayleighMedium::from`
copies density/sound-speed bits without reevaluating an atmosphere. G3 tests
cover explicit ISO/classical selection, invalid ISO humidity, exact medium
views and legacy history equality. GasState itself carries no source receipts.

File-driven bowed, reed and plate-valve performances expose
`from_bytes_with_gas_sources` and retain their authored `GasRecord`s. Ambient,
regional and exterior declarations use the same L3 reader. The CLI loads
bounded source artifacts itself: repeated `--gas-component SPECIES NASA9
TRANSPORT-CLAIMS`, explicit `--gas-water none|ID` and `--gas-phase
frozen|liquid-water`. A primitive reed render additionally supplies
`--gas-record FILE`; wind/ensemble gas records live in their performance files.
No thermochemical route falls back to a preset. Source-aware sidecars retain
composition, temperature, pressure and all three pack identities per species.
Omitting the new options preserves legacy numerical and sidecar paths.

`realize_assembly` turns an `fs_scenario::AcousticAssembly` into observer
pressure by composing nameless primitives: a prestressed Euler–Bernoulli
waveguide (`fs-nlmodal::prestressed_beam_omega` / Kirchhoff–Carrier), a
`bernoulli_aperture`, a `stribeck_friction` contact, a
`traveling_wave_line`, `ModalAcousticTime`, and compact modal radiators.
There is no instrument crate and no instrument algebra. A guitar or
clarinet is one filling of those objects.

- Ambient state comes from `GasState::try_new_moist_air(T, p, RH)` once,
  before realization. Excitation, duct impedance/resonance, radiation loads,
  and observer absorption share that state. RH=0 uses the original dry-air
  constructor bit for bit; RH>0 requires 253.15–323.15 K and water mole
  fraction at most 0.15, with model-domain failures returned as `Ambient`
  errors. Invalid RH refuses admission. This is an estimated ideal mixture:
  viscosity and conductivity retain the dry-air fits, the saturation
  enhancement factor is omitted, and no exact humid-medium prediction or
  experimental validation is claimed. G0/G3 assembly tests check the dry
  boundary, refusals, independent mixture density/speed/impedance, deterministic
  pressure, and waveform-derived duct pitch. The fixed eight-cell Cauer duct
  retains its first-order spatial boundary error; the pitch check uses its
  analytic LC frequency scale and does not establish 5% continuum accuracy.
  ISO absorption remains the sole
  bulk absorption inside its window; no extra classical term is added.
- String path, `EA = 0`: triangular-pluck (or bow) modal ICs marched by
  `ModalAcousticTimeModel`. The shared fixed-end string observer reads actual
  acceleration: odd modes give `p = rho A_odd q_ddot / (4 pi r)` in physical
  displacement coordinates. Even modes have zero monopole area and signed
  moment `M_even = -w L²/(n pi)`, giving `rho M_even q_dddot/(4 pi r c)`.
  This follows the spatial expansion of distributed retarded monopoles
  ([Wood, eq. 6](https://euphonics.org/11-7-1-compact-sound-sources-monopoles-dipoles-quadrupoles/)).
  Linear acceleration includes held forcing and viscous loss; nonlinear
  acceleration uses the admitted pHS storage and `(J-R) grad H + G u`, so
  stretching, conservative contact, and applied forces reach emitted pressure.
  Jerk is a backward difference of accepted accelerations: first-order at the
  reported endpoint, second-order at the interval midpoint for smooth motion.
  The initial acceleration is evaluated from the initial state and first held
  force; no zero-history startup spike or ideal release impulse is synthesized.
  The observer retains the linear path's 1e7 Pa pressure ceiling.
  This is an uncoupled compact strip approximation with the even-mode observer
  on the positive string axis. Common propagation delay, angular listener
  geometry, cylindrical-wire radiation, resolved impact impulses, and radiation
  back-reaction are not supplied by this observer. G1 tests check signed SI
  pressure against a damped oscillator, even-mode jerk and refinement, nonlinear
  stretching-force pressure, and silence under a held static load. Fletcher
  inharmonicity `ω_n = n ω_1 √(1+B n²)` when `EI > 0`. Stokes air
  drag from `GasState` plus the authored constant modal ratio. The reusable
  `air_path::oscillating_cylinder_air_resistance_per_length` computes the
  Desvages (2018) section 3.2.2 cylinder resistance approximation
  `R=2 pi mu (1+r sqrt(2 rho omega/mu))` [N s/m²]; both string paths use
  dimensionless `zeta=R/(2 mu_linear omega)`. Width is interpreted as circular
  diameter. The approximation retains a constant low-frequency continuation;
  added mass, end effects, confinement, rarefaction, turbulence and nonlinear
  fluid drag remain unclaimed. Explicit Rayleigh damping still replaces the
  complete loss sum when no material bending law is selected. The explicit
  Kelvin–Voigt path above instead combines only bending viscosity and air drag.
  Without a selected constitutive law, `damping_ratio` is a finite nonnegative
  constant per-mode reduced-model input. Zero explicitly assumes no internal
  loss; elastic EI does not create dissipation. No hidden `2e-7*omega` addend
  or relaxation spectrum inferred from a single ratio remains. The authored
  ratio is not a material measurement; physical loss requires admitted data
  and a selected Kelvin–Voigt or causal Prony law. G1 checks the third mode's
  actual SI pressure against the independently damped oscillator, and G3 checks
  equivalence to explicitly zero bending viscosity across linear, KC and Dirac
  paths. Invalid authored ratios refuse instead of being clamped.
  G1/G3 tests check the resistance formula and
  pressure ring-down under changed gas pressure. A second
  polarization at `1+detune` is a second member on the same clock
  and shares every plate (it is not an independent unused body).
- String path, `EA > 0`: `fs-nlmodal::kirchhoff_carrier_string`
  (fixed-fixed) or `kirchhoff_carrier_moving_end` (Dirac join)
  + Gonzalez
  pHS step, with the same stiff-string frequencies overwriting the
  linear ω_n. Bridge force `T_eff y'(0)` drives every listed body mode
  (Carcagno-band Sitka pair is a constructor, not a named guitar).
- Bow: a set-valued Stribeck contact (`stribeck_friction::StribeckLoadLine`).
  Each substep solves the exact-ZOH endpoint load line `v = v_free - m F`
  against the unregularized kinetic law of `fs-tribo::FrictionLaw::Stribeck`
  and its static cone `mu_s N`. The residual has at most one local maximum and
  one local minimum, so its stable roots are bracketed exactly. Where a stick
  reaction and sliding roots coexist, the contact keeps its branch
  (McIntyre–Schumacher–Woodhouse hysteresis) and changes only at a fold, a
  slip reversal or a lift. No capture window or ramp scale enters the bow
  path. G1 checks root selection against a dense-scan oracle, the hysteresis
  loop against an independent fold scan, and the Coulomb closed form. The
  bowed gates measure a 4..64 substep ladder at the playable point (playable
  on every rung, steady pitch within 12 cents). The former capture window
  (`max(stiction, 0.02)` and a 0.15 m/s flip basin) is gone; at 3.9 N, its old
  gate point, the regime was not converged in substeps under either contact.
  Helmholtz motion is possible with enough modes; it is not guaranteed and
  not a measured rosin curve. For explicit ODE embeddings elsewhere,
  `StribeckFriction` delegates traction to
  `fs-tribo::FrictionLaw::regularized_traction_1d`, reversing the port's
  driver-minus-body velocity to the owner's body-minus-driver convention. The
  shared owner supplies deterministic exponential evaluation and refuses
  malformed coefficients, negative/non-finite normal loads, and
  unrepresentable force. G1/G3 checks compare the force with an independent
  Stribeck formula and the work of equal-and-opposite contact forces. This
  wiring admits no sourced interface data or authored modal losses.
  `BowedStringCard::from_material_specimen` lowers a resolved circular specimen
  using the shared pinned-beam dispersion, Kelvin–Voigt loss and cylinder-air
  drag owners. It requires a matching loss/material identity and frequency
  band, rejects memory, moving supports and polarization splitting, and applies
  runtime bandwidth admission. The caller retains the source receipts; this
  fixed-state, small-amplitude card grants no new material authority.
  Elastic and viscous support shear are evaluated at the bridge, `x = 0`:
  `T y'(0) - EI y'''(0) - eta I y'''_t(0)`. Explicit cards must include the
  viscous bending contribution in their modal damping ratios. Bow velocity
  and bridge-force histories observe the same completed sample endpoint.
  G1 tests compare independent held-load motion, energy and support shear;
  G3 checks shared pluck/bow losses and sensitivity to material viscosity and
  ambient pressure. This does not evolve material state, route damping heat,
  model nonlinear stretch, or establish bow-contact passivity; substep
  convergence is measured at the playable point only.
  `FrictionIsland::InterfaceStribeck` consumes an ordered dry-interface card
  through `fs-tribo`, retaining its identity, history and caller authority.
  The complete sliding law, including velocity strengthening, drives force;
  other friction families refuse. Fixed temperature, load/nominal-area pressure
  and actual slip must remain inside the card domain at admission and both
  endpoints of each substep, including sticking. These sampled checks do not
  bound interior extrema. G1/G3 compare signed, load-dependent motion to an
  independent oscillator and check runtime domain departure. This path supplies
  neither measured-data admission nor heating, wear or evolving contact area.
  Bowed plate outputs report signed volume velocity in m³/s, not point velocity.
  Plate observers admit positive finite listener distance before stepping and
  refuse non-finite volume flow or pressure. G0/G3 checks cover malformed
  distances and inverse-distance pressure scaling with identical motion.
  This compact observer adds no propagation delay or general near-field claim.
  The material string bank can retain both transverse displacement planes
  through `MaterialStringSource::with_transverse_polarizations(TransverseDirections::try_new(bow, bridge)?)`.
  Both planes use the same original circular section, installed tension, density,
  bending stiffness, material loss and gas drag. Their physical directions project
  displacement and force with the same work-conjugate map; shared body compliance
  determines polarization mixing without an authored detuning. Coordinate order is
  first transverse plane, second plane, then optional material torsion. Surface
  bow/pluck ports include the original radius-times-angle contribution; centreline
  loads and unilateral stops carry no invented radial torque. Stops follow the
  source's declared centreline/bow axis. Existing explicit torque/rotation bridge
  connections compose with both planes and retain their own source-derived angular
  maps and passive laws. Material body/cavity pressure and scheduling remain in the
  original joint solver.

  `frankensim-material-strings-v3` extends the v2 section/support grammar with one
  required `polarizations single` or `polarizations two BOW_X BOW_Y BRIDGE_X BRIDGE_Y`
  record after each string section. Both directions must be finite unit vectors.
  Existing v1/v2 files and optional `rotational-bridges` records retain their
  semantics. `polarized-pluck.material-strings` exercises the data-to-PCM route.
  This circular, frozen-state linear model does not supply arbitrary
  six-degree-of-freedom bridge motion, independent stop-normal geometry, shared
  nonlinear two-plane tension, helical winding/slip, or a measured-instrument timbre
  guarantee. Nonlinear material source conversion uses the existing shared
  two-plane Kirchhoff–Carrier stretching law.

  Seven focused physical tests cover material-mode identity, rotational
  covariance, directional sound changes, force work conjugacy, constructor
  ordering/refusals, combined rotational-bridge response and two-plane unilateral
  contact. Three example tests cover source input, material/angle swaps,
  PCM/callback identity and admission. Native verification is recorded separately.
  `ReciprocalStringBank::with_cylinder_gas_loading` selects a causal gas
  operator for the same circular material strings. The shared
  `air_path::cylinder` owner derives the unsteady incompressible Stokes impedance
  `Z(s)=pi*rho*a^2*s+4*pi*mu*z*K1(z)/K0(z)`, `z=a*sqrt(rho*s/mu)`, from
  actual outer radius and the retained complete gas. The continuum follows
  [NASA/CR-2007-215050, equation (7)](https://ntrs.nasa.gov/api/citations/20080006726/downloads/20080006726.pdf).
  Its positive spectrum per log relaxation rate,
  `8*mu/[pi*(J0^2+Y0^2)]`, supplies physical free inertias joined by passive
  dashpots in the original exact reciprocal system. Both transverse planes
  carry instantaneous added gas mass and viscous memory; original solid
  density, bending, installed tension, material viscosity and torsional inertia
  remain distinct. The active transverse card replaces sampled gas resistance
  with intrinsic material loss plus the explicitly collapsed slow tail, so
  fluid damping is not counted twice.

  All material-string performance versions accept the optional
  `cylinder-gas-loading X_MIN X_MAX INTERVALS` row before an optional
  `prepared-linear` row and then `cavity`, after any section, polarization and
  rotational-bridge records. The same
  original geometry/gas supplies the operator; no extra material coefficients
  are authored. This binding precedes additional bows, stop profiles, static
  preload and the physical pressure/PCM schedule. File construction counts all
  mechanical/fluid coordinates and string/body/torque/cavity/tensor-loss ports
  against the explicit mode, connection and complete setup-work allowances.
  The source fixtures and regressions compose material/mixture substitutions,
  regional body loss, both transverse planes, torsion, cavity, stops and release.

  These fluid components append after physical bodies. Full tensor body loss,
  cavity connections, pressure observation, dry interface and torque maps
  survive binding, either cold torque/fluid order, and multiple-bow scheduling.
  Pressure and force maps follow the changed transverse mass normalization;
  fluid coordinates receive no fictitious direct acoustic pressure weights.
  Added kinetic storage, reciprocal reactions and viscous heat enter the same
  complete work balance. Static initialization admits exactly zero-stiffness
  positions only when their applied static force is exactly zero. Their
  already-zero positions remain gauges while their velocities remain physical
  free inertias; no artificial spring or general singular preload solve enters.

  `CylinderFlowResolution` declares independent slow/fast dimensionless
  wavenumber cutoffs and log-grid intervals: `0<X_MIN<=0.5`,
  `X_MIN<X_MAX<=3000`, and `1..=4096` intervals. The exact integrated slow
  tail is collapsed to resistance with dynamic error bounded by
  `R_low*lambda_min/omega` at positive frequency; the fast tail is truncated.
  Refine all three axes and the actual coupled response. Every branch adds one
  coordinate and one connection per transverse mode. The shared mechanical
  connection ceiling is 4096, with caller mode/connection/setup budgets still
  enforced before allocation; that port capacity alone does not raise the
  ordinary exact connected-island limit. An explicit prepared linear budget
  admits the complete larger system as described below. Repeated/live bindings and loaded elastic
  modes outside their original material bands refuse.

  Operator regressions compare complex force with an independent modified
  Bessel reference and exercise every cutoff/grid axis. Material-bank tests
  exercise pressure-dependent inertia, viscosity-only changes in actual body
  pressure, both transverse planes, unchanged torsion, fluid kinetic energy,
  torque binding order and preload/release. Source review is not native test
  execution. This local long-cylinder, small-motion Newtonian incompressible
  model supplies no finite-amplitude convection/wakes, axial/end/confinement
  or neighbouring-string interaction, rotary drag, rarefaction, compressible
  radiation, evolving atmosphere or measured violin/cello/bass fidelity claim.

- Reed path: quasistatic or massive Bernoulli valve. Reeds with explicit
  positive stiffness derive their effective pressure area from the
  declared closing relation `A = k H / Pc`. Characteristic-line motion and
  ODE pressure force/swept-face flow share that area; a fixed 25 mm face
  length no longer overrides explicit stiffness and closing pressure.
  This is an effective modal area, not measured geometry. The zero-stiffness
  legacy reduction still derives stiffness using a 25 mm face length, and
  damping is now an explicit finite nonnegative ratio, consumed as
  `c = 2 zeta sqrt(k m)`. Reduced valve cards pass `zeta = eta/2`, preserving
  their existing `c = eta sqrt(k m)` model. This removes the runtime's fixed
  0.35 substitution; authored fixture ratios remain explicit assumptions,
  not measurements. G1 covers actual velocity decay; G4 rejects invalid ratios.
  G1 tests exercise the force balance and actual moving-reed acceleration.
  The massive characteristic-line junction now uses the current opening and
  swept-face flow `-A v_mid`, plus body volume flow, in the pressure solve. Positive
  velocity widens the opening: the ODE and characteristic paths share this
  sign, conjugate to `F = -A (Pmouth - Pbore)`. It corrects the ODE's former
  positive swept-flow sign; see Silva et al., arXiv:0810.2870, Eq. 27.
  Its fixed-opening residual is monotone, bracketed between zero jet flow and
  zero pressure drop, and uses the same signed Bernoulli primitive as the ODE.
  The massless solver and its Newton Jacobian are not used for massive reeds;
  fast-mode massive samples count as strict fallbacks. G1 tests cover closed
  apertures, reverse flow, and instantaneous pressure/face work cancellation;
  G3 covers mass-dependent bore pressure at equal static compliance.
  Characteristic-line linear mechanics now uses midpoint position and velocity,
  with that same average velocity in swept flow. Eliminating the linear reed
  response produces a positive effective load for the existing aperture solve.
  G1 checks the discrete spring/kinetic energy change against pressure work,
  Bernoulli loss and viscous loss without contact, plus energy preservation and
  second-order phase convergence of the isolated oscillator.
  Aperture geometry and elastic lay force are held at the start of the step.
  Single-slit Hunt–Crossley loss uses the same midpoint velocity as mechanics:
  the active reaction is `elastic - elastic*chi*v_mid`; a second zero-reaction
  solve handles the nonadhesive unloading clamp. The coefficients come from
  the admitted fs-dcontact obstacle, restricted to one unit opening coordinate.
  G1 exercises loading, unloading and within-step reversal; a lagged-velocity
  negative control injects energy where the accepted loss removes it.
  This proves the sign of loss work against the held elastic reaction, not
  that the latter matches the change in nonlinear contact potential.
  The isolated ODE reed retains its own pHS stepping scheme. This is not a
  nonlinear contact or whole-bore energy certificate: contact work, coupled
  geometry/time refinement and onset validation remain MR68 obligations.
  Isolated
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
- Card-driven acoustic compilation: `material_assembly::compile_material_assembly`
  resolves exact `MaterialCard`/`QueryPoint` inputs with explicit claim selection,
  then calls the existing circular-string, uniform-plate and regional-chart
  specimen binders. Only consumer-needed scalar properties are resolved: rho/E
  for elastic strings, the selected in-plane plate constants, and additional
  viscosity or equilibrium/modulus/time fields for explicitly selected sourced
  string/plate bending loss and explicitly selected isotropic plate thermal loss.
  Inactive property families are not queried. Loss selectors
  are retained; old numerical loss coefficients and citations are rebuilt.
  `StringPrestressBinding::ThermalExtension` selects the existing alpha(T)
  integrator and sourced thermal string binder. Its reference point and alpha
  claim selection are explicit; the immutable card and current point come from
  the same string source as density, elasticity and selected loss properties.
  Mechanical/loss pins and expansion pins cover their separate requirements.
  Piecewise-linear expansion curves are integrated across all crossed knots;
  current alpha times temperature change is not substituted for that integral.
  Fixed supports use the reference stress-free length to derive elastic strain
  and tension. Fixed mass or current radius remains an independent constraint.
  The specimen retains both mechanical and integrated-expansion receipts, and
  the ordinary acoustic realizer uses the resulting tension and rebuilt loss.
  `Prescribed` keeps explicit force/extension/tuning constraints and raw authored
  eigenstrain; it does not query or claim sourced expansion. Invalid paths,
  missing/ambiguous data and slack/compressive states refuse through the existing
  owners. No thermal transport, finite deformation or melting is added.
  `CompiledMaterialAssembly` keeps the numeric assembly and source-bound specimens
  immutable, and its `realize` method calls the existing ordinary/chart dispatcher.
  Unbound components remain authored. Material conditions remain separate from
  ambient-gas conditions; this compiler does not infer thermal equilibrium or
  upgrade the authority of the whole assembly. Geometry, constraints, forcing,
  supports, damping choices and solver controls retain their existing meanings.
  Resolution/binding errors preserve component context and the owning diagnostic;
  failure leaves the input and prior compiled specimens intact. A thermal plate
  selects the existing six-property thermoelastic resolver; optional bending
  viscosity is requested in that same atomic bundle, including exact claim pins.
  The resulting state rebuilds the fixed-mass/thickness specimen and both loss
  mechanisms. Stale numerical coefficients and citations are replaced together.
  An orthotropic material selection with this isotropic thermal law refuses.
  Specimen temperature must be an explicit positive absolute T coordinate; the
  surrounding gas does not supply it implicitly. Compilation is
  synchronous with no cancellation guarantee; ordinary runtime admission runs at
  realization. This is a Rust description API, without a new wire format, global
  scenario identity, optical binding or conservative material evolution.
  G0/G1/G3 exercise explicit-binding pressure parity, consumer-minimal exact-source
  retention, simultaneous string/plate card substitution at fixed mass, regional
  chart substitution through moving-end string/cavity pressure, active Kelvin–Voigt
  and Prony coefficient rebuilding, and missing-member/property/domain/law refusals.
  Thermal compiler tests vary specimen temperature through source curves and
  check changed pressure in linear, sine von Karman and sampled-FE plate paths,
  with and without separately selected viscosity. A fixed-mass card swap checks
  rebuilt thickness and both damping addends against independent Zener and
  Kelvin–Voigt arithmetic; conflicting heat-capacity claims require the selected
  pin to survive compilation. These fixtures validate implementation, not an
  experiment or a transient heating simulation.
  String thermal compiler tests independently calculate integrated strain,
  fixed-mass area/inertia, tension and stiff-string frequency, then measure
  emitted-pressure pitch at two temperatures with and without sourced viscosity.
  The viscosity fixtures supply scalar coefficients at explicitly supported
  temperature states, preserving the existing frequency-wide constant-law
  requirement; sampled coefficients are not promoted to broadband loss laws.
  They also check exact alpha selection independent of mechanical pins, force
  control without alpha resolution, and unchanged inputs after path/slack refusals.
- Uniform plate material binding: `with_uniform_plate_material_state` consumes
  the same `ResolvedMaterialStatePoint` carrier as the circular-string and
  disc adapters. Isotropic elasticity requires rho/E/nu; principal-axis
  orthotropic plane stress requires rho/E1/E2/nu12/G12. It retains every original
  property receipt and checks exact scalar quantity schemas. Card axes map to
  section material axes, with reciprocal nu21 on exchange. The explicit
  `ThinPlate::material_angle_rad` rotates that uniform section counterclockwise
  from rectangle x, including D16/D26 bending/twisting coupling. Geometry,
  supports and forcing remain fixed. The angle participates in the v2 specimen
  identity. No out-of-plane constants are invented in either
  linear or FE-sampled acoustic plate construction. A fixed-thickness or
  fixed-mass prescription derives thickness and mass coherently; geometry and
  material/model identity are retained separately from supports, pretension,
  authored damping and solver controls. These are independent frozen specimens,
  not a conservative thermal evolution. The purely elastic binder refuses stale
  thermoelastic addends. `with_uniform_isotropic_thermoelastic_material_state`
  instead rebuilds the mechanics and loss inputs from the existing resolved
  thermal bundle, retaining its complete receipts and optional viscosity data.
  The source-bound orthotropic adapter also refuses geometric
  nonlinearity: the existing sampled nonlinear plate's isotropic membrane
  approximation is not a resolved orthotropic membrane law. Anisotropic thermal
  loss and phase changes remain unavailable here. The raw
  nonlinear path checks E, nu and G before selecting isotropic sine modes.
  Its rotated FE bending still uses the documented isotropic membrane
  approximation. The modal search uses a sine trial Rayleigh quotient;
  it is not an exact off-axis frequency or a completeness certificate.
  G0/G1/G3 cover receipts, reciprocity, fixed-mass geometry,
  semantic aliases, model boundaries and real pressure changes from material
  and axis swaps, with aligned eigenfrequencies checked against Navier theory.
  An arbitrary angle changes actual assembly pressure; a quarter-turn recovers
  the pressure and eigenfrequency of reciprocal axis exchange.
- Uniform material-bound plate viscosity: `ResolvedPlateSpecimen::with_kelvin_voigt_bending_loss`
  selects proportional isotropic Kelvin–Voigt bending. The dimension-only
  `kelvin_voigt_bending_viscosity` [Pa s] claim must be nonnegative and carry an
  explicit frequency band. Rho, E, nu and eta must be validity-wide scalar
  constants or temperature-only curves at fixed positive absolute T; their
  frequency bands intersect. Temperature, moisture and other
  coordinates remain at the admitted state point, with original uncertainty
  and observation status retained. The assumed viscous plane-stress tensor has
  the same Poisson ratio as the elastic tensor: the moment law is
  `M = D (curvature + (eta/E) curvature_rate)`, giving `C_b = (eta/E) K_b`.
  This is a proportional constitutive specialization of
  [Kelvin–Voigt strain-rate damping](https://doc.comsol.com/6.3/doc/com.comsol.help.sme/sme_ug_theory.06.029.html#3441080),
  not a measured general bulk/shear viscosity tensor.
  The ordinary linear, sine von Karman and sampled-FE plate paths consume it.
  Each retained reference-mode frequency must lie in the admitted band, even
  for a zero-strength viscosity claim. Only bending energy receives this loss;
  geometric prestress does not. For FE modes the diagonal coefficient is
  `(eta/E) phi^T K_b phi / (phi^T M phi)`. Off-diagonal damping coupling is
  omitted, a modal approximation when prestress makes damping nonproportional.
  The nonlinear path retains bending-rate damping but supplies no nonlinear
  membrane viscosity or guarantee that finite-amplitude frequency shifts and
  generated harmonics remain in band. Authored scalar damping alongside this
  law refuses. Explicit thermal loss remains separately additive; its inputs
  must exclude this viscosity's contribution to avoid double counting.
  `compile_material_assembly` resolves the selected viscosity with the plate's
  elasticity and density. A material rebind recomputes its coefficients and
  source identity rather than carrying the old material's viscosity forward;
  missing source data refuses. Thermal rebinding of an already viscous raw
  plate refuses because it would otherwise retain stale viscosity inputs.
  G1 checks the actual sine-mode damping operator and ringdown against the
  analytic Kelvin–Voigt plate equation with timestep refinement. Independent
  geometric-energy integration checks prestressed DKT and sampled-FE damping,
  including additive thermal loss. G3 changes only viscosity through the card
  compiler, retaining elastic eigenfrequencies while changing actual assembly
  pressure in all three plate paths. Fixed-temperature curve tests retain the
  source interpolation receipts and compare plate frequencies/damping, every
  Prony string branch and actual pressure against independently evaluated
  coefficients. The thermal string compiler also exercises interior T values
  with simultaneous integrated expansion and bending viscosity. These are
  frozen-state calculations, not constitutive evolution during heating.
  These synthetic implementation fixtures
  are not experimental material/ringdown validation.
- Regional linear plates: `ResolvedPlateSpecimen::section` exposes the rotated
  numerical section for assignment to `PlateChart` elements. The caller retains
  the resolved specimens for source receipts; a numerical chart does not mint
  material authority or a regional provenance identity. The full-rectangle
  specimen mass is not reused as regional mass: fs-plate integrates local rho*h
  and rotary inertia from each actual section.
  `compile_plate_material_chart` additionally binds disjoint named regions to
  exact resolved material states, constitutive/axis choices and fixed-thickness
  or fixed-region-mass prescriptions. The immutable `ResolvedPlateChart` retains
  these inputs alongside its numerical chart. It integrates actual triangle
  areas to derive h = m/(rho A_region); `with_region_material` rebuilds sections,
  mass and identity using the retained constraints. Failure leaves the original
  binding unchanged. Empty/duplicate names, overlapping or missing assignments
  and out-of-range triangles return `PlateAssignment`; geometry and section
  refusals preserve the owning `PlateError`. No fallback material is inserted.
  The specimen identity includes mesh geometry/connectivity, boundary-node set
  and per-element source/model/thickness/angle in mesh order. Names, assignment
  order and repartitioning with identical resolved per-element inputs do not
  alter it. It is not a complete acoustic-scenario or authored-constraint identity.
  Compilation is synchronous, with no cancellation guarantee or material-state
  evolution; cloned numerical charts do not automatically retain source authority.
  `certified_chart_radiators` consumes that chart with `PlateChartRadiation`
  (explicit eigenvalue search interval, mode budget, support/prestress, authored
  damping and nodal unit-force footprint). All force weights must be finite
  and sum to one; supported forces react into the support without redistribution.
  `realize_assembly_with_plate_chart` feeds the resulting modes into the same
  assembly dispatcher, body ports, piston loads, string/cavity/duct coupling
  and pressure observer as uniform plates. A simultaneous `assembly.plate`
  refuses. Extra authored lumped bodies still compose. The modal interval
  need not include the fundamental; certification applies only in that window.
  It does not certify stability or exclude buckling modes outside that window.
  No spatial thermal-loss law, nonlinear membrane, remeshing, phase evolution
  or full-wave radiation is claimed. Interfaces are perfectly bonded, with
  thickness centered on one flat mid-surface as specified by fs-plate.
  G3 checks uniform/element partition pressure equivalence, legacy uniform
  path parity including moving-end/cavity ports, and pressure/mode changes
  from local density, stiffness, thickness and grain on an unstructured chart.
  G0/G1/G3 additionally check assignment refusals, exact source retention,
  integrated regional mass, region-order/partition pressure equivalence and
  actual pressure changes after named material replacement under both mass
  and thickness constraints. These use synthetic materials, not experimentally
  validated material data or measured sound.
  G0 checks force footprints, modal controls and conflicting descriptions.
- Linear plate projection: the eigenvector, its modal mass `phi^T M phi`,
  signed radiating area `integral(phi dA)`, and signed bridge-force projection
  remain in one basis. Surface integrals use triangle-area P1 nodal quadrature;
  the declared drive footprint is uniform traction over the generated
  rectangle's first cell strip, normalized to one unit of total force.
  There is no mass inferred from drive amplitude, absolute-value radiation,
  or artificial area floor for cancelling modes. Generalized modal force is
  `b_drive F + A_mode p_face`. The pHS ports divide both coefficients by
  `sqrt(modal_mass)`, and its observers convert acceleration back by the same
  factor. Thus eigenvector sign/scale cannot change physical pressure or work.
  This follows modal virtual work and the compact limit of the
  [Rayleigh integral](https://euphonics.org/4-3-2-the-rayleigh-integral-and-the-baffled-piston/).
  The optional diagonal piston self-load uses the physical panel area to set
  its equivalent circular radius; input velocity is `A_mode q_dot/A_panel`
  and generalized reaction is `-A_mode p_face`. Modal normalization cannot
  change that geometry. The approximation omits mutual radiation impedances,
  higher multipoles and higher-order DKT trace reconstruction; eigenvalue
  certification does not certify radiation accuracy. G1/G3 check manufactured
  quadrature/mass, zero-monopole cancellation, mode sign/scale invariance with
  and without the actual fitted self-load, conjugate pHS ports, and density
  scaling through the actual assembly pressure output.
- Plate damping: a single authored viscous ratio is the two-point
  Rayleigh fit (`fs-material::visco`) through that ratio at ω₀ and
  4ω₀ so higher certified modes sit on the stiffness limb. Both linear and
  von Karman banks add Zener thermoelastic loss only from explicit specimen
  inputs. `with_isotropic_thermoelastic_state` binds a resolved six-property
  bundle to the plate, moving elasticity, density, thermal properties, and
  temperature together and retaining its state identity. There is no
  density-threshold material guess or fixed-temperature substitute. Missing
  thermal inputs leave that loss unavailable/unclaimed; requesting the scalar
  law on an orthotropic plate refuses (including a mismatched shear modulus
  when E1=E2). The law estimates the first through-thickness thermal relaxation
  mode of a homogeneous isotropic thin beam/plate; it does not solve a thermal
  field or establish an anisotropic loss law. Callers must exclude this addend
  when authored damping or resolved thermal transport already includes it.
  With prestress, the thermal addend is weighted by the mode's bending/total
  strain-energy ratio: `phi^T K_b phi / phi^T (K_b + K_G) phi` for DKT modes,
  and `D k^2 / (D k^2 + T)` for the isotropic sine modes. Authored loss is
  unchanged. The same projection reaches the sampled-FE von Karman bank.
  This is weak-loss, diagonal modal bending damping about the prestressed
  equilibrium; it does not model finite-amplitude membrane thermal loss or
  off-diagonal thermal coupling. The energy separation follows
  [Fedorov et al., Eq. (1) and Sec. II](https://arxiv.org/abs/1807.07086).
  G1 tests independently integrate geometric energy and compare actual
  unforced sine-mode ringdown to the analytic damped oscillator under refinement.
  G0/G1 tests inspect both assembled damping operators against independent
  Zener arithmetic, vary equal-density specimens and temperature, and cross
  the former density threshold without switching laws. Radiation
  reaction on linear compact radiators uses the sole Bessel–Struve piston
  owner in `fs-phs`, converting acoustic p/U to face p/v by the physical
  area, then fitting and passivity-repairing with `fs-vfit`. Evaluation and
  fitting refusals propagate; there is no silent low-ka fallback or unloaded
  assembly success. `fs-bem` is not a production dependency of this crate.
  G0 `thin_plate_refuses_instead_of_low_ka_fallback` exercises this refusal.
  The existing compact one-branch load keeps ka<=0.5 and dt*fmax<=0.1;
  its analytic leading resistance and inertance coefficients are tested at
  1e-14, with a 512-frequency grid targeting <2.2% complex and <4.2%
  resistance error versus the closed form (native measurement pending).
  These are sampled discrepancies, not continuum enclosures. Linear plates
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
  `CompactBody::drive` reuses the shared exact-ZOH modal transition under
  held force, converting through mass-normalized coordinates. Acceleration,
  volume velocity and pressure describe the completed endpoint. G1/G3 tests
  compare analytical damped/undamped motion and held-force subdivision,
  including a step beyond the former Euler stability limit. Radiation-filter
  reaction remains sampled and held, not an exact continuous coupled solution.
  If motion refuses after the filter advances, the caller must discard the run;
  the filter memory is not rolled back.
  The finished pressure history is passed through ISO 9613-1
  absorption (`air_path`) with the assembly's explicit humidity;
  Stokes–Kirchhoff is only the fallback outside the ISO window.
  Authored string `ζ` stays a constant modal ratio; a relaxation spectrum
  is present only when the caller explicitly selects `PronyBending`.
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
  numbers). A single widening taper is `spherical_cone`
  (`ψ = x p` plus the one inlet Euler near-field shunt).
  Any other taper, alone or in a multi-section chain, is
  frustum LC cells (exact frustum inertance and volume),
  frankensim-d1j6s. A `Closed` termination is the chain's
  rigid outlet (`ChainOutlet::Rigid`, `U = 0`), matching
  the characteristic line's `Termination::Closed`; before
  frankensim-gs17m it was the pressure-release end, an
  octave away.
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
  Its modal damping inputs are complete: it never selects aluminum/steel from
  density or takes the solid's temperature from the gas. `PlateCavitySpec::bind_plate`
  uses the ordinary `certified_radiators` owner to rebuild frequency, material
  damping, force participation and monopole area together. Force and flow ports
  are normalized by the same square root of modal mass. Source-bound callers
  pass their `ResolvedPlateSpecimen::plate()` and retain that specimen's receipts;
  this numeric adapter does not mint material authority. Thermal and viscous loss,
  including prestress dilution, are included once by the plate owner. Nonlinear
  plates refuse at this linear adapter. Failed admission preserves previous modes.
  G1/G3 tests carry compiled source curves through this pHS pressure path at two
  independent solid temperatures on both sides of the former density threshold,
  compare damping with Zener plus Kelvin–Voigt arithmetic, and exercise repeat
  determinism and failed rebinding. They are synthetic implementation checks;
  this adds neither thermal evolution nor exterior plate-radiation back-reaction.
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

DEF-05 modal apertures now admit relative tolerances strictly below 1/128.
Work and flow closure use the actual declared relative ceiling, with no 128x
post-solve allowance. `ModalApertureFrame` records the work and flow
`SolveCertificate`s (joules and m^3/s), including the observed residual, scale,
requested/effective budgets, gamma-n arithmetic floor and the actual scalar
residual-evaluation count. These explicit budgets never auto-raise below their
floor; a refusal retains `fs_phs::PhsError` through
`AcousticRealizeError::AlgebraicSolve`. No state or sample clock is published.
Work's floor uses `gamma_(4*modes+12)` times the absolute energy/work operands;
flow uses `gamma_(2*modes+8)` times its absolute flow operands. These are
conditional arithmetic diagnostics, not bounds on material, flow-model or
time-discretization error. G0/G4 `modal::refusal_tests` tests the 1/128 boundary
and a 1% unaccounted external-force work defect through the production gate.

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
  authority. Its length comes from the shared `LineLadder`
  convergence ladder (DEF-02, IRP .82), not from a fixed round-trip
  multiple; an executed lesson is that a four-trip tail left the
  lowest peak 92 cents flat from wrap-around aliasing. The ladder
  starts at four round trips and doubles until every combo's
  fourth-quarter tail energy is below `tail_tol` and its first three
  peaks move less than `peak_tol_cents`. It refuses with
  `LadderUnconverged` at the ceiling. The line is terminated into the
  tabulated bell load. A spliced analytic fallback outside the
  table's support is a spectral jump, so the ladder refuses it by
  name. Passivity (`|R| <= 1 + 1e-12` on the 4x-oversampled grid) is
  checked on the stored taps without clipping and refuses otherwise.
  A valve switch rebuilds from those taps. Valve
  switches at block boundaries use the carry-outgoing-history lift
  (in-flight waves persist; the new taps govern future reflections)
  plus an optional linear CROSSFADE from the old line — an
  instantaneous reflectance swap is genuinely discontinuous physics,
  and the crossfade is the no-click lane; every switch is logged with
  its lift name and carried sample count. `MatrixFirLine` is the
  full N^2-FIR bake-off arm with per-mode incoming-energy logging.
  It is sized by the same ladder:
  - DC is the matrix's own low-frequency limit;
  - the convergence tail is the Frobenius sum over all N^2 entries;
  - peak shift and passivity are checked on the plane entry;
  - its `MatrixLineReceipt` records each entry's tail.
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

- `flue_loop` module (bead 3ez8g.10.3): the composed flue playing loop
  — no `Flute` type. `FlueVoice` = the jet-drive island x the wind
  line. NO CARD -> REFUSE: the island is parameterized only from a jet
  card minted by fs-aeroac's lab (`JetIsland::from_card`), and only a
  card whose claim class can drive a tone is admitted (the tonal
  interim card; the broadband-refusal-boundary card refuses with
  `CardClassCannotDrive`). Card-backed numbers: convection ratio
  `c_p/U = St (h/delta) / (stage + 1/4)` from the stage-I lock, the
  profile ratio `theta/b`, the Reynolds validity band, the claim class;
  authored and labeled: the spatial deflection gain, the onset seed,
  the mouth radiation resistance. PITCH IS NEVER ASSIGNED: blowing
  pressure and geometry gestures (cut-up, labium offset, jet angle)
  are the controls. Model: Bernoulli jet speed; receptivity = the
  leaky INTEGRAL of the delayed mouth velocity (its 1/omega cancels the
  drive derivative's omega, so the pipe's resonance alone selects the
  mode); the sinuous mode's amplification BAND from the card's theta
  (two low-pass sections at Michalke's tanh neutral `omega_c = U/(4
  theta)` and a high-pass at `omega_c/2` because growth vanishes at
  low Strouhal); labium cut `Q_j = Q_max (1 + tanh((eta - y0)/b))/2`;
  series drive source `rho (Delta_m/S_m) dQ_j/dt` in an implicit mouth
  loop with the mouth inertance, radiation resistance, and the line;
  an EXACT per-block energy ledger (source = stored + numerical +
  radiation + pipe; closure ~1e-15). THREE EXECUTED TRAPS on the way
  to a note: (1) a frequency-flat deflection with the omega of the
  drive derivative locks the highest stage the delay allows (9.4 kHz);
  (2) a continuous onset seed at 1e-3 U, integrated, flaps the jet by
  its own half-width and the labium slope averages to nothing (no
  lock; the seed-independence falsifier in the fixture is the guard:
  a limit cycle's RMS does not move when the seed drops tenfold);
  (3) an amplification band with no low edge lets the integral gain
  win at 6 Hz where the open pipe has no impedance (a relaxation cycle,
  not a note). MEASURED on the flute-class fixture (6 mm bore, 12 kHz
  line for the radiation law's `ka < 1`, 10 x 1 mm flue, 4 mm cut-up,
  tonal interim card): 10-50 Pa lock the fundamental (0.90-1.01 f1),
  110-240 Pa the octave (1.85-2.02 f1), 360 Pa the twelfth (2.84 f1);
  a 5 mm cut-up locks one mode below a 3 mm one at equal pressure; four
  rank voices are bitwise their solo selves in any stepping order;
  replay is bitwise and keyed by Philox voice identity. The card's
  Reynolds band (144-264) does not cover playing pressures (Re ~1500),
  so `claim_check` refuses every claim from this fixture: that boundary
  is inherited by the flue-gates bead, not papered over.

  Optional `with_viscous_windway(length_m)` binds the same material-gas
  `ViscousSlit` used by reeds to the actual flue exit width/height. Its solved
  bulk speed reaches volume injection, Reynolds number, convective delay,
  amplification and the physical outlet receiver. Schedule admission and
  accepted pressure controls check that actual speed against history/Mach
  bounds; cold attachment also rechecks a previously held pressure. Omission
  preserves the ideal jet. `flue_windway` tests cover viscosity-only material
  swaps at unchanged density/sound speed, actual microphone pressure, control
  refusals and callback partition invariance with simultaneous mouth loss.
  The supplied rectangular channel is a quasi-steady parallel-plate reduction,
  retaining the card profile; exit kinetic head is not acoustic mouth loss.
  No entrance/profile development, windway inertia, turbulent jet, transverse
  flute lip geometry or measured instrument fidelity follows from this option.

  `flue_loop::network::performance::file` also admits explicit branched tube
  graphs, section-local material gas records, addressed viscothermal sections,
  passive endpoint/interior R-L-C loads, and multiple coplanar radiating
  openings. The jet and windway bind the actual inlet section's complete gas,
  independent of source order. The same physical graph, fractional/loss owner
  and multirate receiver retain the original node addresses through loss
  subdivision. CLI provenance describes every supplied node, section gas and
  selected opening. The single-tube grammar retains its original construction.
  `music_render_wind::flue::network` supplies G0/G3 source/direct-construction,
  gas/branch/receiver sensitivity, callback/PCM and refusal regressions. This
  adds data access to existing physics; it does not add finger-pad mechanics,
  junction end corrections, mutual radiation feedback, gas evolution or
  measured-instrument authority. See `examples/FLUE_NETWORK_PERFORMANCES.md`.

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

- Cane reed card (bead 3ez8g.6.2, completing T-Wind): "shaped like a
  reed, made of cane" — a tapered vamp mesh + the CC-BY cane packs
  through the fs-solid reduce lab (f1 = 2709 Hz, the measured
  clarinet-reed class) driving the EXISTING realize_reed_bore loop;
  the authored-parameter path remains the Estimate debug image. The
  LAY is a measured-geometry chart (tip opening + facing length ->
  circular-arc radius) consumed as the rest gap, rolling-stiffness
  hardening the named successor; contact stays the fs-dcontact
  Obstacle pattern. Provenance chain in the card: mesh digest,
  material resolution identity (heel E_L dry; the blank-zone density
  mismatch and the 10x E_T anisotropy DISCLOSED; eta = delta/pi from
  the damping pack), reduction options, retained-compliance
  disclosure. Measured sensitivities: more-open lay +46% threshold;
  harder cane brighter at matched breath (+1.9%) with the ordering
  inverting near the soft threshold (both regimes recorded); the
  stiffness threshold shift is below the ladder rung (logged, never
  overclaimed).

- `glottis` module (bead 3ez8g.8.2): the glottis as the SAME object
  as a lip reed and relief valve — two islands on the menu (D21): the
  1-DOF image with the TITZE SURFACE-WAVE mechanism (the flow sees the
  gap one mucosal propagation time late; a bare quasi-steady one-mass
  valve measured NO phonation to 6 kPa — the delay IS the mechanism;
  wave speed authored Estimate) and the two-mass image (4-state
  coupled pHS with the coupling INSIDE the Hamiltonian, corrected
  Steinecke-Herzel driving, asymmetric seed escaping the unstable
  symmetric rest, dcontact collision). FIVE energy/mechanism defects
  were EXECUTED en route and their causes are recorded in-code where
  each was fixed: the short low-loss tract's FIR wraps non-passively
  (4x-oversampled enforcement, the MM-bank lesson); the junction
  Picard iteration does not contract when Zc U rivals the drive
  (closed-form quadratic solve, the brass lesson); a stale one-sample
  line push adds junction phase slip; the S-H recovery factor's SIGN
  is load-bearing (inverted, the cycle locks antiphase onto the tract
  resonance); and the symmetric rest state is a driving-pressure
  deadlock (S-H's own displacement seed). Fold parameters:
  masses/stiffnesses/length RECEIPTED from the licensed two-mass pack;
  damping, rest gap, collision, mucosal speed authored Estimate — the
  honesty ladder's permitted rung, card-driven either way.

- `tract` module (bead 3ez8g.8.1): the tract as A(x,t) with LICENSED
  charts (Assaneo 2011 PLOS CC-BY Table 2 — the authors' own 10-tube
  /a/ and /u/; Spanish, acoustically-inverted, formant-valid register
  disclosed; Story 1996 is unlicensable and the Dresden CC0 meshes are
  the named MRI-grade upgrade vein). D22 as a WHOLE-CURVE gate: char
  impulse-spectrum envelope vs the TMM |Z_in| envelope of the SAME
  walled chart. STRUCTURAL FINDS the lane earned: characteristic_line
  consumers must pass VOLUME-normalized zc (rho c/S — an area-free
  rho*c flattens R to +1 and silently erases the tract shape; the
  glottis module carried this defect and its gates were re-based under
  the corrected load); a rigid-wall short tract realizes ACTIVELY (its
  ringing outlasts the FIR window and time-aliasing makes |R| > 1 —
  the tissue walls are required for the REALIZATION, not just the
  physics); a uniform passivity clip that fires flattens the shape
  (backstop only); dense realization via characteristic_line_dense
  (n_fft override). The morph carry needs p-minus CONTINUITY or every
  rebuild micro-clicks.

- Vowel gates (bead 3ez8g.8.3, completing T-Voice): the first voice
  claims — the vowel from the licensed area function, the source from
  the self-oscillating folds, their INTERACTION through the load at
  the glottis. The corpus loop closes (rows
  acoustic-assaneo-2011-tract-{a,u}): the charts' TMM F1 lands in the
  published Spanish vowel CLASS bands (827 / 340 Hz; a class check,
  never a per-subject match) and the length-shifted falsifier fails
  the same gate at 1104 Hz. THE FLAGSHIP EMERGENT GATE: sweeping the
  fold frequency toward F1(/u/), the phonated f0 rides ~35% above the
  free-fold frequency (load stiffening) and REGISTER-BREAKS at the F1
  ceiling (deviation -147 Hz at the top of the sweep) — the physics a
  source-filter vocoder structurally cannot have, measured. The
  1-DOF-vs-two-mass bake-off receipt is EXECUTED (KeepBoth:
  cheap-ensemble vs solo-realism; slope -6.6 vs -11.5 dB/oct on the
  identical fixture); the ooh/aah ladder is minted for the owner's
  ear.

- `pickup` module (bead 3ez8g.9.2): the Faraday pickup as a
  PURE-READER PORT — D20 as structure (pose + per-mode window gains
  only; `emf()` reads the string owner's modal velocities; no step, no
  state, no re-discretization anywhere). The raised-cosine aperture
  window is unit-normalized so the point-pickup limit is EXACT; the
  window is what makes bridge-vs-neck voicing physical. v1 provenance:
  linearized flux with an authored dipole-like height profile
  (Estimate, disclosed); the MAGNETOSTATIC BAKE (offline B-map lab,
  the lips/jet authority pattern) is the recorded follow-up that
  upgrades provenance without changing the runtime shape.

- `speaker` module (bead 3ez8g.9.4): the Thiele-Small driver + sealed
  cabinet as ONE port-Hamiltonian — the Bl gyrator is an antisymmetric
  J pair (power-exact by construction), the sealed box is
  `Sd^2 rho c^2 / V` extra cone stiffness (exactly why the classic
  resonance-shift algebra holds), and one fs-plate-class panel mode
  couples through the shared cavity term (the cabinet sound as
  physics; X-Consist with ONE retained mode, disclosed). Excursion
  past the card's x_max REFUSES — the suspension nonlinearity is not
  modeled and silent linear forte is the lie the refusal prevents.
  Card = authored datasheet-class Estimate; low-ka piston load
  disclosed; cabinet IR packs are FORBIDDEN and absent.

- Electric chain (bead 3ez8g.9.5): `tests/electric_chain.rs`
  composes string -> pickup -> tone RC -> triode island -> TS driver
  with every stage a named, receipted object. The composition laws:
  ONE string owner (a bare string stepped alongside stays BITWISE
  identical — the chain reads, never owns); voicing survives the
  chain (per-harmonic bridge/neck output ratios = pickup window-gain
  ratios on string-carried lines); the harmonic ladder is
  DEVICE-made (single-mode string, so 2f0/3f0 at the speaker is
  constitutive curvature); the box composes transparently
  (chain-measured excursion-per-volt = driver-alone). The
  split-vs-full-DAE bake-off receipt
  (`tests/receipts/electric-split-vs-full-dae.bakeoff`, KeepBoth)
  pins the split image to the analytic gain and to the monolithic
  1e-14 authority. MEASUREMENT LAW recorded: harmonic line
  projections need INTEGER-CYCLE windows — at a non-integer cycle
  count the fundamental leaks ~1.7% into every harmonic line,
  drive-independent, and read as a 27% violation of a ratio law
  that holds.

- Composed guitar body (bead 3ez8g.7.4): `tests/guitar_body_gate.rs`
  gates the braced-top + ladder-braced-back + Helmholtz model's
  BRIDGE-ADMITTANCE triad against the Carcagno 2018 BR-guitar corpus
  rows (F1 +5.9%, F2 +9.9% in 15% envelopes; F3 -32.6% in an
  authored 40% envelope with the mechanism DISCLOSED — the flat
  rectangle surrogate has no X-bracing/arching, so its third
  breathing family sits low: the named weakest claim). Executed
  falsifiers: the midline-brace order flip (frequency-order breathing
  selection picks a 6e-15-volume mode; select by max |INT phi dA|)
  and the MEASURED refutation of nodal-line folklore (an ECCENTRIC
  brace on the nodal line shifts the mode 11.6% via torsion +
  eccentric-membrane coupling). The wolf note is a pHS pair: tune
  EFFECTIVE frequencies (k = m w^2 - k_c; the coupling spring detunes
  a naive pair by 10% and no hybridization happens), then the string
  on the breather beats at DETRENDED depth 0.046 (the estimator
  lesson: divide out decay before measuring modulation) vs 0.990
  clean off-resonance. Truncation sensitivity is PER-FREQUENCY
  (spread 23.6x) and its direction is fixture-specific: here
  between-peaks suffers more than on-peak because low-index resonant
  modes survive basis halving.

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

`tests/cane_reed_card.rs`, cases cr-001..cr-003 — the provenance-chain
card mint (identity reproducibility, retained-compliance floor, the
reed-class f1 band, the lay's facing-scale radius); the sensitivity
sweeps with a SUSTAINED-oscillation speak detector (tail at 2% of
blowing pressure and non-decaying — the absolute-threshold detector
fired on ringing at every rung, twice, before this); and refusals by
name (a card without a modulus cannot resolve; lay-chart nonsense;
an empty orifice face).

`src/glottis.rs` `glottis_tests`, cases gl-001..gl-005 — both islands
self-oscillate against the live tract with onset curves in the
physiological range (1-DOF 188-276 Pa, two-mass 404-865 Pa, rising
with stiffness, none on the ladder floor — the detector demands
DETRENDED sustained output after a settling ramp read as phonation on
the first run); the collision/damping ledger (undriven kick decays,
per-step dissipation non-negative); the two-mass vertical phase
(upper lags lower 0.306 rad at the fundamental) with faster-closing
flow; the bake-off FIXTURE QoIs on the identical tract + card (f0
152/173 Hz, OQ 0.78/0.82, slope -6.6/-10.0 dB/oct — the receipt is
vowel-gates' to mint); refusals by name.

`src/tract.rs` `tract_tests`, cases vt-001..vt-004 — the whole-curve
D22 gate (envelope correlation 0.87//a/ and 0.90//u/ with dominant
resonances within 50 Hz; identity /a/ >> /u/); the yielding-wall F1
raise (340 -> 401 Hz, nearest-peak comparison past the wall's own
~100 Hz mode); the u->a morph (100 Hz control rate, carried, clicks
18x less than the five-cut waypoint control on the REFLECTED wave,
F1 track 353 -> 806 Hz logged); licensing-first refusals (unlicensed
data refuses by name). The glottis module's gates were RE-MEASURED
under the corrected load in the same change (two-mass onsets
404-715 Pa physiological and monotone; the 1-DOF's ~5.8 kPa
non-monotone onsets recorded as the fixed-delay mechanism's real
behavior; flow skew logged, not asserted).

`tests/vowel_gates.rs`, cases vw-001..vw-004 — the formant class gate
with its executed length-shift falsifier (a uniform AREA scale would
not move formants; the falsifier must shift LENGTH); the sung
interaction sweep (ten stiffness points crossing F1 with per-point
f_free/f0/deviation logging); the committed bake-off receipt gate
(KeepBoth scopes + both rows alive in the registry); and the
vowel-ladder listening chain (Unadjudicated).

`src/pickup.rs` `pickup_tests`, cases pk-001..pk-005 — station
voicing with every per-mode EMF line matching an INDEPENDENT
Simpson-integrated analytic gain within 5% (bridge centroid 452 Hz vs
neck 328 Hz); two mixed pickups combing per the gain-sum prediction
with an executed notch; the exact point-pickup limit (5e-4) plus the
pure-reader law (identical repeated reads; EMF recomputable from the
public gains and the same state slice); refusals by name; and the EMF
driving the .9.1 circuit DAE end-to-end (RC-load rolloff matching the
first-order analytic response at 0.25x and 4x cutoff with supply
defects at 1e-21).

`src/speaker.rs` `speaker_tests`, cases ts-001..ts-005 — the sealed
box shifts resonance per fc = fs sqrt(1+Vas/Vb) (measured 64.7 vs
65.1 Hz analytic); the electrical impedance matches the analytic TS
motional curve to the printed digit at four frequencies; the panel
colors the response near its resonance and NOT far away (executed
contrast against the rigid box); a 60 V resonance drive trips the
excursion refusal while a gentle drive never does; refusals by name
plus the supply-defect ledger at 3e-19.

`tests/electric_chain.rs` (7 + 2 mints): ec-001 compose + logging +
the one-owner witness; ec-002 voicing ratios within 6% (gain-floored
to string-carried lines); ec-003 H2 0.07% -> 0.86% across grid drive
0.03 -> 0.40 V RMS with H2 >> H3; ec-004 composed-vs-direct box
transfer within 0.04% at 65/260 Hz; ec-006 committed bake-off
receipt re-asserted (split gain 59.950 = authority 59.950, analytic
59.923; pluck L2 4.6e-12); ec-008 listening digest chain
(electric-clean-vs-driven, Unadjudicated); ec-009 the electric
gate-summary artifact enumerates every registry row.

`tests/guitar_body_gate.rs` (5 + 1 mint): gt-001 triad + product
invariant (defect 0.0) + damping honesty (material-only Q 39/54/78
above measured 34/18/36 as it must be); gt-002 the executed
midline-brace negative control; gt-003 wolf note; gt-004 truncation
spread; gt-006 guitar-strum listening digest chain (Unadjudicated).

## No-claim boundaries

- TOMBSTONE-3ez8g.13.5: `traveling_wave_line` is dead code with zero
  in-tree consumers, superseded by the exact-FIR characteristic line;
  every public item is deprecated + doc(hidden) and the module header
  denies new consumers. No music image may cite it (standing law).
  Deletion stays available on explicit written owner permission.


- The electric chain claims NOTHING about the string (one owner; the
  vibration epic's rows carry it) and nothing about perceived tone
  until the listening receipt is adjudicated. The preamp input gain
  is a measured-calibration volume knob (disclosed, not physical);
  budget rows wait on the budget lane covering this chain (D25);
  live-default stays no.

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
- The glottal islands claim PHONATION MECHANISMS, not vowels: no
  formant filters, no vocoder, no glottal-waveform playback (the
  waveform emerges from the valve + tract loop or not at all). The
  uniform 17.5 cm tract is a fixture, not an articulation model
  (.8.1's lane); the bake-off receipt between the two islands is
  minted in vowel-gates (.8.3), not here.
- The tract module claims chart realization and morph mechanics, not
  vowels-with-a-glottis (.8.3 composes them); no fricative jet, no
  nasal branch (named refusals); mouth radiation is IdealOpen v1 (the
  compact-mouth Z bake is the recorded refinement — the radiation-fit
  ka ceiling forbids UnflangedOpen at full audio band).

### IRP radiated observation admission (M2.6-W0)

`pcm_wav::observation::ObservationPolicy` distinguishes `Acceptance` from
`Diagnostic`. `ReedPerformance::from_bytes_with_policy` refuses `listener_m`
under Acceptance; its physical outlet still requires the unchanged flanged
feedback load and uses one exterior gas for feedback and reception. The raw
voice, positional reed CLI fixture and legacy baseline remain diagnostics.
`music_render wind --acceptance` selects the flanged-development path explicitly.

`TwoSidedPressure::with_observation_policy` and
`AperturePerformance::new_with_observation_policy` require the actual installed
reaction binding under Acceptance. They check the load's identity, radius,
clock, band and exterior density/sound speed (relative tolerance 1e-12).
`from_baffled_environment` uses that checked path. Observation adds no second
reaction load, changes no bore gas and never changes a termination to admit it.
Existing one-way builders retain Diagnostic scope.

G1/G2 analytic tests in `baffled_pressure_qualification.rs` check on-axis piston
pressure at z/a=2,10,100 and signed 2J1(x)/x directivity at ka=0.5,2,5 within
1e-3, with radial/angular refinement. The on-axis fixture uses a=25 mm to
respect the existing z>=50 mm domain; the far-field fixture uses a=7 mm.
The band parameter is admission metadata: its ladder checks invariant kernels
and live spatial/temporal refusals, not waveform convergence as the band changes.
G0 checks the finite-difference 0.1-rate and differentiating 0.45-rate caps at
cap and cap+1 Hz. G3/G5 adapter tests compare the same accepted terminal flow
against the shared owner and repeat pressure bitwise across callback sizes.
These are authored infinite-baffle fixtures, not real-instrument fidelity or a
full-band reaction match. Native qualification remains required; rustfmt alone
does not establish the numerical claims.

## Frequency-domain branched chart and network port (IRP M2.1-W8)

`tract::network::TractBranch` and `TractChart::network_with_branches` turn a
chart and explicit branch geometry/source identities into `fs_duct::network`
data. The resulting ports are glottis, oral and the supplied branch labels.
All propagation, junction equations and graph elimination stay in fs-duct;
there is no nasal-specific solver. The former no-nasal-branch documentation is
lifted for this frequency-domain API. TractVoice's time-domain series image is
unchanged; dynamic branch-line realization remains M2.2-W3c.

`driving_point::network::NetworkDrivingPoint` is the explicit FD port seam for
player-airway/bore consumers. Construction requires the graph, active label,
every other port's load, and a caller-supplied scalar `FluidPressureFlux`
PortSchema, with owner-outward orientation. Gas, omega and loss remain explicit
at each evaluation. `harmonic_port` performs the actual complex transfer
`p=Z U_in`, `U_out=-U_in` with that schema. It does not claim a time-window
energy audit, glottal dynamics, or a completed M3.1-W7 player-airway realization.

`tract::network::tests::network_branched_tract_nasal` uses Assaneo's licensed
/a/ chart plus a DECLARED SYNTHETIC 0.13 m x 4 mm-radius branch after section 2.
The mouths are explicitly ideal pressure release; this is a pole-zero diagnostic,
not a measured nose or an acceptance voice render. G3 checks the lower first
formant and added oral-transfer antiformant, with 8/4/2 Hz frequency grids
(last-pair feature positions within 4 Hz); all sampled N-port matrices must be
reciprocal and have PSD real Hermitian part to 1e-12 normalized tolerance.
Qualitative side-cavity context: Dang, Honda & Suzuki, JASA 96,2088 (1994),
<https://pubmed.ncbi.nlm.nih.gov/7963023/>. No human-specimen frequency oracle is
inferred from that citation. Tests print formants/antiformants and provenance.
Additional cross-crate tests exercise the explicit typed transfer, missing-load
and wrong-schema refusals, and unbranched graph/chain equivalence (1e-10).
All new numerical assertions remain pending central native verification.

## Assembly model-domain admission (IRP M6.1-W1b)

`description::validate(&AssemblyDescription, &ClaimDomain)` is the mandatory
admission seam for description-render consumers (M5.4-W2). Its private
`Admission` token borrows the exact immutable description and claim, exposes the
unchanged requested band, and records each selected registry row, module,
configuration, evaluated axis, refinement route and resolved material-state
identity. It does not render or narrow a band. Existing raw component render
APIs are unchanged; the description renderer must accept Admission directly.

`ClaimDomain` declares gas, dynamics, an existing `selector::RegistrySnapshot`
and explicit `ModuleDomain` selections by EntityId. This reuses
`fs_evidence::ValidityDomain`; every required or extra constrained axis must be
provided. No missing-axis-as-unconstrained rule applies here. A refused or absent
registry row cannot admit a model. Ungated rows remain visibly ungated in the
limitations; domain admission is not a gate promotion or Q3/Q4 acceptance.

Primitive parametric geometry vocabulary (all dimensioned SI; a
`test-fixture:` prefix marks authored fixtures without changing semantics):

- `cylinder.v1`: radius, length, evanescent_relative_error.
- `tone-hole.v1`: hole_radius, bore_radius, chimney_height, edge_radius.
- `circular-opening.v1`: radius.
- `lumped-reed.v1`: rest_opening, width, closing_pressure, attack, mass,
  stiffness, damping_ratio.
- `circular-string.v1`: radius.

Construction and every named ConfigurationState are checked. Assignments must
name actual parameters with matching dimensions. Unsupported imported geometry,
other part/model combinations, incomplete coverage, duplicate part/module
selections and mismatched material card/query bindings refuse. Required valve
and string elasticity/density resolve through MaterialSource and the existing
state-point owner; explicit Unknown refuses before card resolution. No property
is defaulted. Geometry-source citations and authored bounds are declarations,
not a proof of the bound or of real-specimen equivalence.

Domain policies explicitly supply frequency_hz for every module, ka for compact
radiation, edge_displacement_ratio for compact holes,
length_over_shock_distance and evanescent_relative_error for ducts,
mouth_pressure_pa for reeds, bow_force_n/pluck_force_n for strings. These
policies intersect owner limits. The §1.2 rounded audit cutoffs are
1.841*c/(2*pi*a) and 3.832*c/(2*pi*a), with no symmetric-only escape. Segment
geometry, Bessel/wide-tube scope, compact radiation, reed primitives and gas
coherence reuse their existing owner admission functions. A strict hole/bore
boundary is recorded with upper_inclusive=false.

Edge displacement uses the declared lowest excited frequency and conservative
pressure envelope. Inverse shock distance uses beta=(gamma+1)/2; linear string
drag requires KC=2*pi*displacement/diameter <= 1. A zero pressure gives zero
inverse shock distance without an infinity sentinel. The numerical policies for
“much less than” and evanescent error must come from the caller's declared
evidence; this validator does not invent them. These algebraic evaluations are
deterministic for identical inputs and perform no stochastic work or hot kernel.

G0/G1 tests `description_validate::dv_001` through `dv_010` cover analytic
cutoff agreement at 1e-12 relative, exact/next-value boundaries, ka=1 and declared
ka=0.5 policies, strict hole radius, unknown materials, slit refusal, mandatory
axes, retained limitations, dynamics-driven edge/shock/KC and configuration
coverage. Native execution is pending central verification. No audio golden or
wire schema changes; no specimen-fidelity, resolved anatomy/flow or model-form
equivalence claim follows from this token.

## Material-derived circuit storage

`MaterialCircuit::new_with_storage` binds every R/C/L branch to its material
owner, lowering dielectric/core coefficients through `fs-material::field_storage`
into the existing Kirchhoff descriptor. Original branch ordering, conserved
charge/flux linkage, source receipts and per-resistor Joule heat are retained.
Changing dielectric or core data changes the actual electrical evolution.

`rebind_storage_candidate` proposes a zero-time material/geometry change at the
exact same charge and flux linkage. It returns each branch's field-energy change
and the total signed parameter work. The coupled physical owner must supply or
receive that work when accepting the candidate; it is not counted as resistor
heat or source-port work. Incompatible ideal-source constraints, incomplete or
duplicate bindings, unrepresentable energy and cancellation refuse atomically.
C/L remain fixed within an ordinary circuit step. This linear lossless update
does not supply continuous thermal/mechanical feedback or switching impulses.

## Reciprocal material capacitor motion

`material_circuit::motion::MaterialCapacitorMotion` composes the existing
source-bound RLC circuit, uniform dielectric capacitor and regional rigid-body
inertia. The active capacitor's original quadratic term is replaced exactly once
by q²d/(2 epsilon A); all other C/L storage and original J/R/G topology are retained.
The same descriptor solve advances charge, flux, gap and momentum. Voltage and
attractive force are derivatives of the same field energy; mass comes only from
`ResolvedBody`. Changing supported dielectric or electrode material data changes
the coupled motion and electrical history without changing a physics law.

Original electrical ports keep their order and an external force port appends.
Bulk conductor and electrical contact losses reuse the existing native-effort
accounting and remain separately sourced heat destinations. Explicit refreshed
contacts must retain their IDs/branch assignments. Their losses are counted once.
The combined field/kinetic storage, original source work, force times actual
displacement, and heat must close before any joint candidate is returned.

The endpoint retains exact proposed charge/flux and rechecks the original circuit
constraints at the new gap. Incompatible ideal clamps, unsupported source data,
open gap bounds, unrepresentable values, balance failures and cancellation refuse
without partial publication or charge projection. Native domain-aware finite
differences shrink unsupported probes rather than expanding material validity.
A validity-wide constant dielectric coefficient, fixed area, ideal electrodes
and a rigid translational guide are required. Thermal dielectric curves require
field entropy before electrocaloric coupling; fringing, breakdown, impact,
dielectric loss and electrode deformation are not inferred.

Focused regressions cover data-only density/permittivity swaps, both charge
signs and zero-charge motion, nontrivial driven RLC/contact topology, a 100 µm
gap at 100 V, and atomic constraint/domain refusal. Native execution is pending;
source review and independent numerical oracles are not a Rust test pass.

### Material two-plane nonlinear string composition

`MaterialStringSource::into_nonlinear` admits both transverse planes to one existing fs-nlmodal Kirchhoff–Carrier stress channel, with one source-derived EA, elongation, live tension and cross-plane stretching work. Original EI, prestress, mass, material/gas loss, torsion and directional centreline/surface/torque maps retain their ordering. `NonlinearMaterialString::prepare` retains that source while composing explicitly supplied receiver/contact/felt/volume participants in the original `PreparedImpactSystem`; `MaterialStringRun` implements the existing `ImpactSource` for the existing causal pressure renderer. The continuous vector-slope bound is invariant to transverse basis rotation. This is a frozen, moderate-slope reference model. Scalar translating supports still refuse two planes; full rod/longitudinal dynamics, transfer from a nonlinear bow-friction bank, acoustic reaction, full instrument validation and exact-linear/SAV bandwidth work remain open.

`ReedPerformance::from_bytes_with_gas_field` accepts the same source-retaining `AxialGasField` for authored geometry and gestures. Its exact inlet state drives reed density and impedance; the resolved midpoint cells supply one existing regional TMM for both feedback and terminal pressure. The exact outlet field state supplies the default exterior, and an explicit exterior record replaces that single binding for both flanged loading and baffled observation. The field remains frozen during rendering. Common pressure, full length, geometry/refinement and conflicting per-segment source declarations refuse through the existing owners. Standalone and ensemble sidecars can inspect the original primitive stations and source-pack identities.

The flue performance file accepts optional `walls COUNT FMIN_HZ FMAX_HZ` plus addressed `wall NODE isotropic PACK T_K OMEGA_RAD_S THICKNESS_M stress-free|strain-fixed elastic|kelvin-voigt` records after the original duct geometry, followed by `radiation prescribed|compact-baffled`. The CLI loads caller-relative NormalizedMaterialCardPack artifacts (8 MiB each, 64 MiB aggregate), while the library resolves original cards at independent solid T/omega states with SingleClaimOnly property selection and retained source/usage receipts. Shared artifact paths load once; geometry, branches, regional gases, existing losses and outlet addresses remain with their original owners. Each material wall occupies the existing degree-two node's physical dual strip; unselected and endpoint half-strips stay rigid, and actual area/coverage is reported. Compact radiation binds every declared observed outlet through FlueVoice::from_baffled_environment, then validates each installed binding against the lowered graph and the same exterior/receiver declaration. Material substitutions reach wall storage, jet-returning bore pressure and the common physical microphone; exterior substitutions affect both reaction and propagation when selected. This remains the existing frozen small-motion isotropic local-ring image plus independent low-ka baffled loads: no full shell coupling/exterior wall radiation, full-band unbaffled or mutual reaction, resolved turbulence or measured Boehm claim. The composition remains EstimateOnly and these additions do not close full M0.5 or M2.5-W3c/W4b qualification beads.

### Fixed supply-passage momentum in material reed instruments

The shared `bernoulli_aperture::slit_channel` exposes its retained no-slip
profile as the exact next-step discrete affine boundary `dp = Z*Q - history`
for a fixed `StaticPressure` passage. Its source gas, dimensions and clock
supply the impedance; history is fluid momentum, not a second authored pressure
source. `TwoSidedAperture::with_supply_channel` adds this impedance and history
to the existing single, independent-pair or modal reed load. The same accepted
through-flow drives the channel and both acoustic ducts; physical blade-face
pressure subtracts the actual channel pressure drop. All participant trials
and the complete flow/energy account pass before channel, reed and paired
waveguide states are published.

The supply profile retains initial kinetic energy, release and signed reverse
flow. `TwoSidedFrame::supply` separately exposes storage, pressure work, wall
heat and backward-Euler numerical dissipation; both losses enter the numerical
account once. The material-wall and shared baffled-environment constructors
compose the same supply with original regional wall sources, Maxwell histories,
fractional/viscothermal propagation, radiation reaction and the physical
microphone. Existing schedules consume the resulting source directly. The
`two_sided_reed` example selects explicit passage geometry and mesh with
`--supply-slit LENGTH_M WIDTH_M HEIGHT_M TRANSVERSE_CELLS`.

The passage is additional fixed, planar, fully developed incompressible
Newtonian geometry, distinct from an existing duct segment or the moving reed
clearance. Closing-wall remapping, squeeze work, compressibility, axial flow
and evolving atmosphere remain outside it. A kinetic-head terminal refuses
affine condensation. Source/clock/load mismatches and complete-step failures
leave accepted state unchanged. Added regression sources cover pressure/flow
condensation, material single/paired reeds, initial fluid momentum, viscosity
interventions, physical microphone output and transactional refusals; native
execution and measured clarinet/oboe validation remain separate requirements.

`SteadyAxialGas::solve` in fs-material can feed its source-retaining `AxialGasField` directly to the existing regional duct, reed-bore and scheduled ReedPerformance field constructors. Wall cooling and explicit gas injections therefore change the same EOS/transport medium used for bore propagation, reed inlet loading and physical termination. The shared thermal owner supplies stationary enthalpy/species balances; the emitted field is frozen for acoustics. Its grid, caloric, conservation and low-Mach limits remain explicit and do not imply resolved transient breath evolution or measured clarinet/oboe equivalence.

### Explicit prepared propagation for complete modal networks

`PreparedModalBudget` selects cold preparation of the complete, fixed-period
linear mechanical operator. `CoupledModalSystem::new_prepared` and
`new_with_columns_prepared` admit initially larger networks; consuming
`with_prepared_linear_propagation` on the coupled, friction or reciprocal
string-bank owner preserves its original state and sources. The budget has
positive `max_modes`, `max_dense_bytes` and `max_multiply_terms` fields, with
a shared ceiling of 512 mechanical modes. Checked conservative additional
dense-storage and cold multiply-work counts precede prepared allocations.
Original setup-work, physical, material-band and transaction checks remain.

The prepared state retains every displacement, velocity and held force.
Original signed spring/damper columns, rest loads and material damping form
the same mass-normalized stiffness and damping matrices. A scaled full-state
exponential gives the held-force transition; separate component and individual
connection loss Gramians retain physical resistance ownership. Short-interval
Van Loan integration and Gramian doubling avoid a stiff inverse exponential
over the full sample period. Independent transition refinement, symmetry and
the complete quadratic energy identity gate cold acceptance; original sample
energy, per-port work and force/pressure ceilings gate publication. Preparation
does not truncate modes, replace the nonlinear friction solve or change clocks.

The bank retains the selected budget across cold torque and gas reconstruction,
and its complete owner continues through additional bows, stops, preload and
pressure/PCM rendering. All material-string performance versions accept
`prepared-linear MAX_MODES MAX_DENSE_BYTES MAX_MULTIPLY_TERMS` after the optional
cylinder record and before `cavity`. The file retains its 128-coordinate cap
and selects preparation before attaching torque/fluid. Its initial bank must
still fit the ordinary constructor; use the direct prepared constructors for
an initially larger generic network. Omission retains the ordinary path and
its existing 42-mode connected-island limit.

The supplied two-plane layered-string fixture contains 54 total mechanical
coordinates, including torsion, body, cavity and 48 free gas-memory inertias.
Its explicit 64 MiB / 200-billion-multiply-term budget authorizes costly cold
preparation, not a runtime or real-time guarantee. Regression sources cover
ordinary/prepared equivalence, noncommuting matrices and rest loads, individual
loss/work integrals, timestep partition, an analytic 43-mode connected system,
stiff free inertia, complete resource refusal and actual prepared material-string
preload/pressure. Native execution, numerical convergence of the chosen physical
basis and measured violin/cello/bass fidelity remain separate requirements.

`ReedPerformance` accepts an optional authored `gas_field` block after segments, selecting the existing continuous Magnus or midpoint field owner with explicit resolution, station positions, temperatures, canonical source IDs and mole fractions. Numeric resolution remains midpoint shorthand. Only caller-loaded species/NASA-9/transport packs supply the material model; declared components may be exactly zero at the inlet or throughout the field. A positive dyadic composition is used solely to admit the complete source catalogue and is never evaluated physically. Original knots and pack identities remain retained, with active knot receipts on the existing gas-record/provenance surface. The same field supplies acoustic propagation, exact inlet reed coupling and default exterior loading/pressure observation. Existing exact-length, common-pressure, source-domain, phase, geometry and aggregate-work refusals remain, and explicit segment or injected-field coexistence is refused. Wind and ensemble sidecars distinguish continuous steps per smooth span from midpoint cells and authored segment samples. Regression coverage includes source/order/zero/seed invariance through actual midpoint and continuous pressure, callback parity and the unchanged absent-block path. This remains a prescribed frozen field through the existing compact reed/duct/radiation model, not a transport/reaction solution or measured instrument qualification. The checked-in N2/CO2 example includes the existing source-pack construction and rendering commands.

### Material fluid films and finite heat

`render::plate::impact::squeeze::material::MaterialSqueezeFilm` resolves the same `FluidTransport` snapshot used by flow and heat, retains the concrete source and its receipts, and compiles viscosity into the original `fs-tribo::ResistiveFilm`. Its pressure and analytic tangent remain inside the existing joint mechanical equation. `HeatedImpactSystem::new_with_films` routes selected contact and incompressible-film work to the same finite caloric stores. Film heat is evaluated by the original pressure owner at the actual discrete-gradient velocity and midpoint geometry, then deposited once before mechanical acceptance. Existing thermal-network subcycling, conduction, external sources and rollback apply unchanged. Network reports distinguish `deposited_heat_j` from true combined-system `external_heat_j`.

The film is frozen, incompressible, Newtonian, no-slip and quasistatic. Caller-declared thin-gap, pressure and applicable inertia/cavitation/rarefaction restrictions remain; no continuous viscosity-temperature feedback, added fluid inertia, phase-changing geometry or implicit wall/fluid heat partition is inferred. Compressible-film free-energy loss and unrelated mechanical/acoustic losses are not relabelled finite heat. The source binding remains inspectable while mechanics owns its immutable compiled pressure operator. Regression sources cover data-only liquid viscosity/caloric swaps, analytic pressure/work, unchanged mechanical stepping, finite temperature and network conduction, and atomic refusals. Native execution was unavailable in this session; the independent numerical reduction is not a Rust test pass.

### Source-driven thermoelectric transport

`circuit_enthalpy::thermoelectric::step_candidate` composes the existing prepared material network, Kirchhoff descriptor and enthalpy residual. Each constant-absolute-Seebeck solid leg binds an original voltage-source branch in exclusive oriented series with its material resistor; bulk distances and areas must match its two source-derived conductive links. One pinned source claim covers both junctions and the actual Joule-heated bulk. The native midpoint source current drives Kelvin heat `[-S*T_a*I, S*T_b*I]` into the original junction cells while the existing resistor owner deposits Joule heat once in its bulk cell. Shared junctions sum the signed coefficients across materials, yielding thermocouple voltage and Peltier heating without material-name presets. Internal Seebeck work and Peltier heat are removed from external-supply accounting before the total electrical/thermal balance gate. All source, geometry, numerical and cancellation refusals preserve both histories. This bounded zero-field, fixed-geometry, single-phase leg model rejects variable-S/Thomson, tensor/AC and electrical-contact extensions; other thermal cells retain their original phase and boundary behavior. Endpoint coefficient coupling and lumped Joule deposition require time/spatial refinement. Physical convention reference: https://doc.comsol.com/6.4/doc/com.comsol.help.heat/heat_ug_theory.07.087.html .
