# CONTRACT: fs-cli

The stable command-line membrane for the Cooling 0.1 product workflow (bead
`frankensim-extreal-program-f85xj.6.2`). The binary is named `frankensim`;
the package is `fs-cli` so the workspace retains its flat `fs-*` crate
convention.

## Purpose and layer

Layer L6 (HELM). `fs-cli` turns command-line arguments and project bytes into
deterministic result records and structured diagnostics. It owns presentation
and exit semantics, not project-schema, solver, report, or package authority.
Those remain with `fs-project`, `fs-session` and the cooling pipeline,
`fs-report`, and `fs-package` respectively.

It also owns the product-level geometry-import orchestration checkpoint for
bead `frankensim-extreal-program-f85xj.6.3`: exact caller-supplied raw bytes
flow through `fs-io` quarantine/promotion, `fs-project` persistent assignment,
and one atomic `fs-ledger` operation. Lower layers remain free of filesystem
and L6 policy.

As of bead `frankensim-extreal-program-f85xj.6.5` (slice 1) it additionally
owns the deterministic solve-orchestration driver: content-derived run
identity, session budgets derived from the project, staged ledgered
execution, cancellation at stage boundaries and fixed solve-owned evidence
checkpoints, durable stage checkpoints
through the fs-exec legacy v1 snapshot envelope, and honest
budget-exceeded/stage-gap terminals.

## Public types and semantics

The v0 grammar is intentionally small:

```text
frankensim [--json] validate <project.fsim|project.json>
frankensim [--json] import <project> <source>... <ledger.db> --unit <unit> --max-hole-edges <n>
frankensim [--json] import <project> <source>... <ledger.db> --unit <unit> --step-root <id> --target-h <spacing>
frankensim [--json] solve <project.fsim|project.json> <ledger.db>
frankensim [--json] solve --resume <run-id> <ledger.db>
frankensim [--json] report <run-id>
frankensim [--json] package <run-id>
frankensim [--json] compare <left-run> <right-run> [<ledger.db>]
```

`report` (and `run`) also export `<run>.field.vtu` (driver v33+): the
published tetrahedral mesh with nodal `temperature` and per-cell
`region_label`. The bytes are the conduction stage's retained
`solve-conduction-field-vtu` artifact, which the conduction receipt names as
`field_artifact`, and they are re-hashed on export.
`import` takes one source per project geometry row and binds them in
declaration order; a count mismatch refuses `cli-import-source-count`, and a
swapped pair refuses on the pinned source hash. A faceted-STEP import names
one root entity and takes exactly one source.
`--json` may appear once at any position. Unknown flags, duplicate/missing
operands, mixed mesh/STEP policies, non-integer repair/root values, and
non-finite or non-positive STEP spacing are refused. Project inputs are capped
at 16 MiB before parsing.
`.fsim` selects the canonical s-expression spelling and `.json` the canonical
JSON spelling; unknown extensions are refused rather than guessed.

`validate` invokes the strict `fs-project` reader and all of its recognition
and semantic checks. A successful result reports the canonical project hash,
schema version, zero findings, and the exact authority class
`structural-project-admission`.

`report` and `package` execute against a run whose seven stages are all
sealed (an incomplete run refuses `cli-report-run-incomplete` and names the
next stage). Before either verb reads a byte of the retained report or
package, it proves the run by **sealed evidence**: the driver-state chain is
re-attested exactly as resume does (row identity, driver version, exact
lineage edges, intact edge seals, recovered card packs reproducing the run
identity), every retained receipt is re-hashed against the hash the sealed
driver state recorded, and the conduction operation's retained outputs are
taken from the sealed operation itself with their kinds checked. Exports never
replay physics; only `solve --resume` re-executes retained stages, because it
continues computing from the verified state. The report result discloses this
as `"verification":"sealed-evidence"` next to
`"authority":"projection-of-retained-receipts"`. MEASURED 2026-09-02: with the
replay, `report`/`package` on the heatsink example cost twice the solve
(3.4 s of a 4.2 s validation was re-solving conduction).

`compare <left-run> <right-run> <ledger.db>` diffs two completed runs of the
same project by their retained receipts and, like the other export verbs,
never replays physics. Both runs go through the resume loader and the report
receipt binding above; the rows are read from the retained QoI receipt (values,
colours, requirement outcomes, the eight budget terms), the material-resolve
receipt (pack-set root and per-kind card identities), and the seven stage
receipts. A stage row is `unchanged (same receipt)` when the hashes agree,
`unchanged (same inputs; differs only by binding keys)` when the two receipts
differ in no top-level key but `run`, `project_hash`, or `import_op` (which
run, project, and import op the receipt belongs to, not what the stage
computed), and `changed` otherwise; the differing keys are listed so the
reader can see why. A design change (another card pack, limit, fidelity, or
geometry import) is a different canonical project, so the two project hashes
are reported side by side with `same_project` and a changed hash is one of the
rows. The result carries `"authority":"projection-of-retained-receipts"` and a
no-claim: a changed value is a change between two Estimated candidates, not a
verified effect. Refusals: no ledger operand (`cli-export-ledger-required`), a
missing ledger (`cli-export-ledger-missing`), an unknown or incomplete run (the
solve loader's own codes), two runs with no common requirement to diff
against because their QoI names or requirement identities differ
(`cli-compare-project-mismatch`, naming both project hashes), and a receipt
whose shape cannot be read (`cli-compare-receipt-shape`). Comparing a run with
itself reports `changed:false` and no row marked changed. MEASURED 2026-09-02:
the hostile foam twin of the reference project (rebound to a 0.04 W/mK card
under the same chemistry key) shows exactly one material card change, a hotter
maximum with a smaller nominal margin, both verdicts `indeterminate`, and five
of seven stages changed (`import-verify` and `assign` differ only by binding
keys); a merely renamed project shows the hash change and a bit-identical QoI.

### Cooling-network total-enthalpy transient

The experimental `cooling-network <request.json>` binary command admits an
opt-in `transient.enthalpy` storage model in its existing
`frankensim.cooling-network.v1` request. It calls the production
`HeterogeneousEnthalpyBackwardEuler` owner (also for one uniform chart), including
its ambient or closed-enclosure radiation endpoint when requested, and the
existing quasi-steady mixed-air coupling. This command
remains a nominal estimate outside the native `.fsim`/ledger workflow.

The enthalpy object requires explicit Newton/Krylov/backtracking budgets and
one uniform or nodal specific-enthalpy initial state. The original uniform
chart fields remain valid. Alternatively, `materials` supplies 1–256 named
records, each with a positive fixed reference density, 2–4096 equilibrium
chart knots and caller-declared material identity/source; `element_materials`
assigns exactly one name per tetrahedron. These storage assignments are
independent of conductivity and contact material labels. Uniform fields and
the assigned representation cannot be combined. Each enthalpy vertex belongs
to exactly one named record: different charts require separate interface
vertices and an explicit thermal contact, even if their numbers coincide.
No chart blending or temperature-to-enthalpy interface conversion is inferred.

Each chart's optional `phase` is `solid-liquid` by default, preserving the
existing mixed-phase ingress. Explicit `solid` and `liquid` call the material
owner's single-phase constructor: every knot has the corresponding liquid
fraction 0 or 1, and temperature rises strictly with enthalpy. An ordinary
sensible solid therefore needs no fictitious melting endpoint.
It excludes the temperature/capacity initializers and `transient.nonlinear`.
The material owner validates the chart; temperature and phase are evaluated
from accepted specific enthalpy, including unsmoothed latent plateaus. The
declared identity does not mint a source or material-resolution receipt.
Storage is `sum_i m_ref_i h_i`, with tetrahedral reference mass lumped at the
vertices. Phase change never changes that reference mass or adds latent heat
again. Existing conductivity assignments and supported contacts remain in the
spatial operator.

Each uniform or named enthalpy material may also declare `phase_conductivity`
with `law: "linear-liquid-mass-fraction"`, finite positive `solid_multiplier`
and `liquid_multiplier`, and a nonempty `source`. At each nonlinear trial,
each element uses `K_eff = [s_s + (s_l-s_s) f_bar] K_base(T_bar)`, where
`f_bar` and `T_bar` are arithmetic means of its nodal liquid **mass** fractions
and temperatures. The multiplier follows enthalpy material assignments;
the base tensor still follows the independent conductivity assignments,
including their temperature dependence and anisotropy. This caller-declared
phenomenological law does not infer a volume-fraction mixture or material
calibration. An omitted law has unit multiplier. Robin, radiation and contact
transfer laws are not multiplied.

Newton and endpoint adjoints include both `dT/dh` and the direct `df/dh`
conductivity term. Latent enthalpy can therefore affect neighboring
temperatures even when its local `dT/dh` is zero. Fixed/adaptive forward
trajectories, repetition, sizing and supported fixed-grid adjoints use the
same prepared law. The multipliers remain fixed control data; no derivative
with respect to their values is reported. A fraction-slope corner used by a
variable multiplier refuses a classical endpoint derivative, as a
temperature-slope corner already does. The supplied law and source are
retained beside the corresponding output material/chart fields.

Every air-reference and radiative iteration holds the same accepted `h_old`
fixed. Physical time advances only after the full implicit residual, solid
energy, coupling heat and complete endpoint energy gates pass. Radiation uses
the actual endpoint patch mean and remains separate from the air convection
load. Accepted input, air gain, radiation loss and storage enter the window
once. Exhausted work, chart-domain errors or cancellation publish no partial
trajectory. The Newton worst-case inner-column product must fit the existing
linear budget, and `linear_restart <= min(solid vertex count, 256)` also bounds
Krylov workspace; fixed step-count and wall budgets cover the complete invocation.

The final JSON includes `solid_specific_enthalpies_j_kg`,
`solid_liquid_mass_fractions`, existing final temperatures and the scheme
`backward-euler-total-enthalpy`. `transient.enthalpy` retains initial/final
total reference enthalpy and the chart/source/density/control/work declarations.
Assigned charts retain their named `materials`, `element_materials`, and
derived `vertex_materials`; each history observation resolves its own chart.
Each accepted history row also carries specific-enthalpy extrema and a
reference-mass-weighted mean liquid fraction. A new schedule may restart from
the returned nodal enthalpy; no temperature-to-enthalpy inverse is inferred on
a plateau. The forward JSON does not store full nodal history.

On fixed grids, `transient.adjoint` admits `final` or `sampled-peak` objectives
and an explicit `max_checkpoint_bytes`. The dedicated enthalpy tape carries
the physical h cotangent through every earlier step; it never divides by
`dT/dh` on a latent plateau. It binds checked endpoint implicit derivatives
with complete mixed-air and declared ambient or enclosure radiation feedback. Retained
h/temperature/reference checkpoints and control accumulators are charged
before allocation. Reconstruction requires exact accepted h and temperature
bits and rechecks physical residual/energy gates. The cap is not a bound on
all solver workspace. Reverse solves share the original invocation deadline;
failed replay, derivative gates or cancellation publish no partial gradient.
Fixed-count `transient.repeat` retains one chronological tape across every
cycle. Its initial-h derivatives refer only to the original state; each
interval workload, inlet and fan control affects every occurrence. Checkpoint
admission accounts for all retained endpoints, and reverse carry crosses
cycle boundaries without a temperature-to-enthalpy inverse or a cold reset.

The report's initial-state keys are
`dtemperature_dinitial_specific_enthalpies_k_kg_j` and
`dtemperature_duniform_initial_specific_enthalpy_k_kg_j`. Interval
`dtemperature_dpower_multiplier_k` and graph-node-ordered
`dtemperature_dinlet_temperatures` retain the ordinary trajectory conventions.
For an admitted single affinity-scaled fan bank, interval
`dtemperature_dlog_fan_speed_ratio_k` includes uniform air-capacity scaling,
complete solid/air/radiation feedback and supported Reynolds dependence of
convection. Declared scalar HTC has zero speed dependence. A convection card
without an admitted smooth Reynolds derivative leaves the fan derivative
unavailable (`null`), while the complete forward candidate remains usable.
Optional `component_power: true` and `contact_resistance: true` use the same
accepted trajectory and coupled adjoint. Component reports retain the original
consistent-P1 PowerMap footprints, including overlaps and zero-watt components;
interval rows differentiate applied watts at every occurrence of that base
interval. Base-watt rows affect only `power_scale` intervals, leaving explicit
`component_powers_w` overrides fixed. Contact reports differentiate the persistent
log resistance through every earlier endpoint. The source-density pullback and
watt-valued nodal-load multiplier already include the storage discretization;
these contractions add no timestep factor or per-control solve. Their retained
accumulators share `max_checkpoint_bytes`. Missing component footprints or
contact interfaces refuse before trajectory publication. The report keys and
units match the ordinary fixed-capacity trajectory controls.

Charts, reference densities, conductivity, geometry, fan curves, quadratic
loss coefficients, convection law data and radiation controls remain fixed;
contact resistance is fixed unless its sensitivity is explicitly requested.
Material slope corners and validity endpoints refuse classical endpoint derivatives;
plateau interiors are supported. Ties select the existing earliest sampled
maximum/active-vertex branch without claiming a unique derivative.

`transient.power_design` and `transient.fan_speed_design` use the existing
workload- or fan-multiplier sizing searches and are mutually exclusive.
With a sampled-peak adjoint it can use its checked directional sensitivity;
omitting `transient.adjoint` selects the existing derivative-free search.
Fan candidates bind actual scaled speeds into an immutable schedule shared by
forward solves, reports and adjoint replay. Flow and correlated HTC are
recomputed at each candidate; none is frozen at the original speed. Every
candidate restarts from the same initial enthalpy and runs the complete
physical trajectory with residual and energy acceptance before sampled
feasibility is classified. The reported feasible workload or fan schedule applies
to sampled endpoints; no continuous-time peak bound or global optimum is
claimed. `examples/cooling-network/enthalpy-power-sizing.json` is a runnable
synthetic phase-changing request with coupled air and ambient radiation;
`examples/cooling-network/enthalpy-fan-sizing.json` exercises the fan control.

`transient.repeat` admits fixed `cycles` or bounded `until_periodic`, with the
existing `max_total_steps` cap. Actual accepted nodal h and T cross each cycle
boundary. The last-cycle `transient` report remains in local time; cumulative
`repeated_cycles` contains every cycle's energy, accepted work and sampled
peak, including warm-up. Sizing restarts each candidate from the original h
and uses the complete repeated experiment. Periodic stopping requires both
`temperature_tolerance_k` and explicit positive
`specific_enthalpy_tolerance_j_kg`, applied to maximum absolute nodal
start/end differences for the declared consecutive cycles. Temperature
agreement alone cannot establish a settled latent state. Temperature-storage
requests reject the enthalpy tolerance. A failed periodic gate, total-step
budget or cancellation publishes no partial result. Periodic stopping does
not admit adjoints. Optional `repeat.fan_controller` uses the existing hysteretic
controller, sampled from accepted temperature at each cycle start and held
through that cycle. Its multiplier composes once with the prescribed fan
schedule and any outer sizing multiplier; every reachable speed is admitted by
the actual fan bank. Periodic passes require both field residuals and an
unchanged next controller state. The result includes sampled sensor values,
applied multipliers and switch counts. Controller state advances only at cycle
boundaries, never on nonlinear or rejected adaptive trials. Adjoint requests
with hysteretic control explicitly refuse. The synthetic
`enthalpy-controlled-fan.json` example exercises phase storage with closed-loop
cycle-rate fan control.
The cycle-map residual is not distance to an infinite-cycle limit or a bound
on future peaks. The synthetic `enthalpy-repeated-pulse.json` example carries
both latent and sensible storage through three cycles with a global adjoint.

Total-enthalpy schedules also admit `transient.adaptive` with positive
`absolute_tolerance_k` and an explicit positive
`absolute_specific_enthalpy_tolerance_j_kg`, plus the existing
`relative_tolerance`, `minimum_trial_step_s` and `max_trials`.
A complete coupled full step is compared with two coupled half steps. The
maximum normalized discrepancy across every nodal temperature and every nodal
specific enthalpy must pass; relative scales use changes from the old field,
so a reference offset cannot relax acceptance. The second half step starts
from the first half step's h. Only the accepted pair contributes history,
phase summaries, sampled peaks and energy; coarse and rejected trial solves
remain charged to solid work. There is no extrapolation.

The trial ceiling respects both `max_step_s` and each interval's requested
spacing `duration_s / steps`; trials stop at workload/fan events. The minimum
trial step bounds refinement, while shorter interval tails or explicitly
requested spacing remain possible. Both accepted endpoints must fit
`max_steps` and the remaining repeated `max_total_steps` budget before a
trial starts. Trial statistics reset each cycle; `transient.adaptive` describes
the final cycle, while `repeated_cycles` owns cumulative accepted work and
energy. Adaptive forward repetition, periodic stopping and derivative-free
workload/fan sizing use this same h history. Adaptive adjoints explicitly
refuse. The reported `estimated_local_error_ratio` appears only at the second
endpoint of an accepted pair; it is a local endpoint discrepancy estimate,
not a global, midpoint or continuous-time error bound.
`examples/cooling-network/enthalpy-adaptive-pulse.json` supplies a synthetic
mixed sensible/latent field. `tests/cooling_enthalpy/adaptive.rs` adds analytic
plateau/repeated-energy checks, accepted-grid replay and bounded refusal cases.

Prescribed top-level `recirculation` is supported throughout fixed/adaptive
enthalpy trajectories, repeated cycles and workload/fan sizing. Each solid/air
endpoint solves its mixed intake from the current exhaust with the declared
return fractions held fixed. Accepted history reports mixed supply temperatures,
true fresh-air/exhaust heat gain and mixing residuals. The final top-level report
retains the declared topology and full supply/return streams. Separate endpoint,
whole-window and all-cycle gates close stored enthalpy against source energy,
actual fresh makeup/remaining exhaust and net radiative heat. Only
accepted endpoints enter `fresh_exhaust_energy_gain_j` and
`fresh_exhaust_energy_residual_j`; repeated-cycle summaries expose both per-cycle
and cumulative values. Coarse/rejected adaptive trials contribute work only.
Fixed-grid inlet sensitivities refer to fresh makeup and include implicit return
feedback. Return fractions are fixed controls; their transient derivatives or
design remain unsupported. This is instantaneous adiabatic thermal mixing with
an imposed pressure reset, without return-duct hydraulics or air residence time.
`examples/cooling-network/enthalpy-recirculation.json` combines mixed-phase
storage, a prescribed return, component watts and a final-temperature adjoint.

Top-level `radiation.enclosure` admits the same explicit gray emissivities and
reciprocal, closed view-factor matrix as the other cooling paths. Each patch
emits at its P1 area-mean temperature and applies a uniform radiative flux over
its exactly matching convection region. All radiation and mixed-air trials
retain the previous accepted h field. Acceptance reassembles the complete
nonlinear h residual, checks fresh radiosity against applied heat, and closes
stored enthalpy against the original convection and signed enclosure residual.
Only convection enters air transport. Opposite patch exchanges change local h
but cancel from whole-system heat; `radiative_out_w` and accumulated
`radiative_energy_loss_j` are the signed enclosure closure residual, not the sum
of outgoing patch heat. Final patch rows retain both `outward_heat_w` and
`applied_heat_w`. This path composes with adaptive steps, fixed/repeated schedules,
prescribed returns and workload/fan sizing. Fixed-grid history adjoints include
the complete reflected-radiation h Jacobian, including latent zero temperature
slopes and phase-dependent conductivity. Emissivities, view factors and geometry
remain fixed controls. `examples/cooling-network/enthalpy-enclosure.json` and
`tests/cooling_enthalpy/enclosure.rs` supply a synthetic two-patch analytic heat
reference, mixed-phase trajectory finite differences and accepted restart checks.

This CLI consumer admits fixed or adaptive timesteps and fixed reference
densities with uniform or explicitly assigned equilibrium charts. Steady design
and controls outside these two transient searches, time/mesh studies,
and return-fraction design explicitly refuse.
No moving geometry, phase advection,
fluid storage, phase kinetics or inter-step peak/error certificate is claimed.
`tests/cooling_enthalpy.rs` exercises the real binary with independent analytic
nodal/whole-window latent storage, hot/cold reservoir and unchanged air-load
references, endpoint summaries, deterministic replay, restart parity and
unsupported/data/budget/cancellation refusals. Its heterogeneous matching-contact
case checks every nodal reference mass, two distinct latent temperatures,
internal-transfer cancellation and both hot/cold radiative reservoirs against
independent algebra. Single-phase sensible diffusion is compared with an
independent four-node linear solve, and a latent insert/ordinary-solid case
checks joint transport. The fixtures are synthetic:
`examples/cooling-network/enthalpy-phase-pulse.json` and
`examples/cooling-network/enthalpy-contact-materials.json`; requests and numerical
references are documented in `examples/cooling-network/TRANSIENT_COOLING.md`.

### Material discovery

`frankensim [--json] discover <request.json> <pack>...` reads the five existing
canonical material-data pack families and queries the real in-memory
FrankenSQLite `MaterialStore`. It writes no input file or persistent database.
Duplicate exact packs collapse; conflicting pack IDs refuse transactionally.
Candidate order follows the store's canonical order, independent of file paths.

Request schema `frankensim.discovery.v1` requires these fields:

- `target`: `materials`, explicitly unbound `properties`, or an object with
  ordered `surface_a` and `surface_b` material-state objects. Each state has
  nonempty `chemistry`, `phase`, `process`, and unsigned `revision`.
- `properties`: scalar `{name, unit, kind}` requirements. Units use the
  existing fs-qty grammar; output values use coherent SI. `dimensional` is
  explicit and is never a wildcard for semantic kinds. Supported real-scalar
  kind names include `absolute-temperature`, `temperature-difference`,
  `energy`, `torque`, `pressure`, `stress`, `heat-capacity`, `entropy`, and
  fs-qty's material-convention names. Unsupported kinds, tensor/test-context
  fields, and all other unknown fields refuse. Full tensor/hardness-test
  requirements remain available through the typed store API.
- `models`: `{law, version, pin?}` requirements; version is a positive exact
  integer, and an optional pin is a member-card content hash. No implementation
  or material association is inferred from a name.
- `domain`: `{mode, axes}`. `local-state` axes have `{name, kind, value}`;
  `envelope` axes have `{name, kind, lower, upper}`. Coordinates are quantity
  strings such as `293 K` or `26.85 degC`, normalized by fs-qty. Equal endpoint
  dimensions and ordered values are required; duplicate axes refuse. An
  explicit axis kind of `legacy` requests the source's untyped coordinate
  profile; this cannot prove source-axis units or satisfy a typed coordinate.
  `dimensional` and semantic kinds preserve their exact descriptors.
- `selection`: `single-claim-only` or `prefer-observation-backed`, delegating
  unchanged to the owner. Duplicate property/model requirements refuse, and
  at least one property or model is required.

`discover_paths` accepts an actual `LawRegistry` for embedders. The binary
currently has no registered constitutive graph factories, so associated model
requirements retain `UnknownLaw` implementation gaps; merely loading model
data never enables a law. Model membership, source selection, card and node
domains, and factory admission use the existing store implementation. A
successful report does not initialize state or bind a connected solver graph.

JSON output (`frankensim.discovery-result.v1`) carries complete/partial/
unavailable candidates, exact pack hashes and named identities, the requested
state profile, property descriptors, selected claims and SI endpoint values,
all property/model gaps, and globally unknown requirements. Gap details retain
the owner's diagnostic and witness; their Debug rendering is explanatory text,
not a stable serialized Rust error type. Text output reports the same candidate
and requirement outcomes. Exit0 means discovery completed, including a report
with gaps or zero candidates. It does not mean the requested simulation is
admitted. Read/encoding/size failures exit3; request, canonical-pack, and store
refusals exit4; malformed CLI grammar exits2.

Admission caps:64KiB request,32 packs,16MiB each and64MiB aggregate,64 entries
per requirement/axis array. Reads enforce actual byte counts as well as file
metadata. Processing is synchronous and bounded by these caps; it has no OS
signal/cancellation integration or peak-memory guarantee. No source TSV parser,
source acquisition, model evaluator, or new catalog schema lives in the CLI.

The actual lead example uses the existing source compiler:

```sh
cargo run -p xtask -- matdb-pack --manifest data/matdb/seed-v1/lead-pure-nbs-c447/manifest.tsv --out /tmp/lead.fsmatpk
cargo run -p fs-cli -- --json discover examples/material-discovery/lead-heating.json /tmp/lead.fsmatpk
```

That request explicitly searches unbound properties: the source pack combines
point-specific physical values with distinct processed hardness conditions,
so naming it lead does not define a coherent specimen. The requested293–650K
envelope exceeds density's declared293K point; conductivity, specific heat and
latent heat are absent under the requested names. A local density request at
293K returns the original11340kg/m³ with its selected claim. These checks
verify source transport and declared coverage, not measurements or an operating
heating/melting solver. `tests/discover.rs` exercises the actual binary's
positive named-material result plus domain/semantic/input/model gaps;
`xtask/tests/matdb_pack_cli.rs::g0_cli_lead_heating_discovery_explains_actual_missing_inputs`
uses the real source compiler and the same example request through CLI dispatch.

### Euler cinematic static admission

The `cinematic` command is the stable user membrane for bead
`frankensim-h7xu5.8.1`:

```text
frankensim [--json] cinematic inspect <config.fscine> (--trajectory <artifact>|--run-reduced) [--dry-run]
frankensim [--json] cinematic storyboard <config.fscine> (--trajectory <artifact>|--run-reduced) [--dry-run] <host-resources>
frankensim [--json] cinematic daily <config.fscine> (--trajectory <artifact>|--run-reduced) [--dry-run] <host-resources>
frankensim [--json] cinematic representative-4k-frame <config.fscine> (--trajectory <artifact>|--run-reduced) [--dry-run] <host-resources>
frankensim [--json] cinematic final <config.fscine> (--trajectory <artifact>|--run-reduced) [--dry-run] <host-resources>
frankensim [--json] cinematic resume <config.fscine> (--trajectory <artifact>|--run-reduced) [--dry-run] <host-resources>
frankensim [--json] cinematic verify <config.fscine> --trajectory <artifact> [--dry-run]
frankensim [--json] cinematic mux <config.fscine> --trajectory <artifact> [--dry-run]
```

`<host-resources>` is the complete tuple `--memory-bytes`,
`--free-storage-bytes`, `--wall-time-s`, `--workers`, and
`--paths-per-second`; a partial tuple refuses. Launch modes require it because
their named profile estimate is admitted against caller-declared host
availability. `inspect`, `verify`, and `mux` may omit it and then report
`resource_admission=not-requested`; that is not a host-capability claim.

The command decodes the complete 1 MiB-capped `.fscine` document, binds both
budget references to the exact named profile identity, resolves regular-file
assets relative to the canonical config directory, rejects path escape and
stale bytes, and hashes admitted assets under per-file and aggregate byte
ceilings with cancellation checkpoints. `--trajectory` additionally decodes
the canonical bounded Euler trajectory artifact, checks the configured
trajectory schema version and exact artifact identity, and leaves the artifact
reader to validate its distinct wire-codec version. `--run-reduced` records an
unverified request in the write-free plan; it does not schedule, persist, or
pretend that the campaign has run.

Every successful inspect/plan admission emits deterministic
composition/partition identities, the exact profile and estimate (if
admitted), trajectory verification facts, planned stages, authority, and
no-claim text without physical locator paths. Static admission never creates
the configured artifact root. Until their producer Beads land, non-dry
launch/resume, whole-bundle verify, and mux calls fail closed with the owning
dependency; storyboard/daily execution belongs to `frankensim-h7xu5.8.3`,
representative/final/resume to `.8.2`, verification to `.8.4`, and quarantined
mux to `.8.5`.

### Geometry import

The library surface exposes `RawGeometryLibrary` and
`import_project_geometry`. A raw-source library binds one caller path/label,
byte payload, length unit, repair or STEP root/sampling policy, and optional
named face groups to the strong identity of an exact canonical project geometry
row. Physical labels are provenance only: `.fsim` continues to identify
imported receipt/content rows, never machine-local paths. STL, OBJ, and PLY use
the quarantine/promote route. Strict triangular faceted STEP uses the
caller-selected positive `FACETED_BREP` root, validated unit ID, and finite
positive target spacing, then assigns on the exact repaired soup returned by
the lower topology/SDF handoff. Its retained wrapper preserves the separate
native decoder and tessellation-import receipts without merging their
authority.

The `import` command executes this path for the product reference shape: exactly
one canonical project geometry row, one raw source file, and one SQLite ledger
destination. Mesh callers must explicitly select a hole-repair cap; STEP
callers must explicitly select root and spacing; both supply the source
coordinate unit. Multi-source and named-group adapter construction remain on
the library surface until a bounded source-manifest grammar is ratified.
Project validation and bounded source-file reading finish before the ledger is
opened; format admission and import refusals are then retained in that ledger.
The command derives deterministic execution seed from the project and caps raw
source bytes by both the import default and declared project memory budget.

On success the library atomically retains:

- exact hostile source bytes as input lineage;
- the exact `fs-io` promoted receipt in both an artifact and the Imports
  extension table;
- a deterministic lossless PLY spelling of the promoted finite mesh;
- each `fs-project` assignment report;
- a complete orchestration receipt and terminal successful operation.

Parse, source-hash, parser-version, promotion, and post-promotion assignment
refusals retain all evidence available at the refusal stage and finish one
terminal error operation. Project-admission, resource-envelope, and
pre-cancellation failures occur before ledger side effects.

### FrankenScript executor v0 (`run <program.fs|program.fs.json> <ledger>`)

`run` routes an operand named `*.fs` (s-expression) or `*.fs.json`
(canonical JSON) to `frankenscript::run_program_path` (bead
rc-root-q61wp.38). The executor parses the program, admits it through
`fs_ir::admission::admit` (the only authority for units, budgets, capability
grants and the explicit pillars), re-recognizes the LOWERED canonical program
as a `Study`, and binds every clause statically before anything runs:

- `(let <name> (cooling.project "p.fsim" [:hash "<project hash>"]))` binds a
  project; the hash pin, when present, must equal the project's canonical
  hash, the study seed must equal the project's `seeds.root`, and a declared
  `(budget (wall …) (mem …))` must cover the project's solve-time and memory
  budgets.
- `(cooling.import <project> :sources ("…") :unit "…" :max-hole-edges <n>)`
  (or `:step-root <id> :target-h <length>`), `(cooling.solve <project>
  [:materials (…)] [:interfaces (…)])` and `(cooling.run <project> …)` call
  the same stage drivers as the `import`, `solve` and `run` verbs, in program
  order; the first non-success step stops the program.
- `(study.run "study.fsim" [:budget <n>] [:hash "<blake3 of the file>"])`
  calls the canonical `study` driver (Journey B: the 2-D marquee). The study
  file's `seeds.root` must equal the program seed; its own budgets are
  enforced by the study driver. Study run ids chain through receipts that
  embed wall time, so parity with `frankensim study` is every receipt field
  except `consumed_wall_s`, `predecessor` and the run id.
- Any other verb or `let` binding refuses as `frankenscript-not-executable`
  ("not executable: verb X has no stage binding"); argument errors refuse as
  `frankenscript-argument`. Admission findings render as `admission-<check>`
  diagnostics (warnings on stderr with `"severity":"warn"`). Every refusal
  happens before a ledger is opened.

Paths resolve against the program's directory; card packs come from the
program (a `--materials`/`--interfaces` flag on a program `run` is a usage
error). The result record carries the study name, the program hash
(`hash_domain("frankensim.frankenscript.program.v0", lowered canonical)`),
each bound project's path and hash, and every executed step's own result
record. `examples/heatsink-fan/heatsink-fan.fs` (Journey A) and
`examples/marquee/bracket-2d.fs` (Journey B) are the tracked programs.

### Solve orchestration

The library surface exposes `run_solve`, `resume_solve`, `SolveRunId`,
`SolveStage`, `SolveDriverState`, `SolveOutcome`, and `SolveRefusal`. The
`solve` verb wraps them with bounded project I/O and a real monotonic clock.

Run identity is content-derived before any side effect:
`hash_domain("org.frankensim.fs-cli.solve-run.v1", project canonical hash ||
constellation || workspace || root seed || driver version)`. Budgets travel
inside the project hash, so raising a budget starts a new run whose completed
artifacts still deduplicate by content. Driver semantics version 7 binds the
current executable prefix: version 3 added retained card packs and material
resolution, version 4 added flow-network execution, version 5 conditionally
added conduction when the project carries its explicit conduction setup, and
version 6 added exact matching-P1 finite-contact lowering. Version 7
length-frames every variable-width interface-resolution identity field so
distinct declared names cannot collide by shifting delimiters. Older run ids do
not collide with or silently resume under version 7. Every solve
operation carries the 32-byte run identity as its ledger `session` value; the
run's own operations are its index. Resume and import discovery walk globally descending
`visible_op_ids_page_controlled` pages under a first-page high-water mark,
then reconstruct each bounded candidate through
`read_op_fields_controlled`. Retained field hashes are compared with the
fixed operation-content sidecar only after cheap row/IR/Five-Explicits
prefilters; a malformed sidecar on an unrelated row therefore cannot block
discovery. No ledger schema extension is required.

The session capability token is derived from project budgets: wall seconds
from `:solve-time`, memory bytes from `:memory-bytes`, cores pinned to 1
(the slice-1 driver is single-threaded), core-seconds = wall × cores, and one
verb per stage. The driver measures its own wall time (the governor meters
what it is told) and charges the `fs-session` governor after every stage.
Warnings fire at declared consumption fractions (0.5 and 0.9 of the wall
grant); any `Throttled`/`Paused` enforcement stops the run after the current
stage with the honest `budget-exceeded` terminal (exit 6) carrying the
completed stages and their durable artifacts. The partial is never presented
as complete.

Stages execute in pinned order: `import-verify`, `assign`,
`material-resolve`, `flow-network`, `conduction`, `qoi`. Each completed stage
is one atomic ledger operation: frozen stage IR, the Five Explicits, the
stage receipt artifact, lineage links, and the sealed driver state. The first
stage additionally retains the exact canonical project source and every
admitted card pack as input lineage. `import-verify` streams every retained
import artifact through the ledger's content-hash verifier (row presence is
not authority), checks the frozen import byte/count envelope against retained
raw, canonical PLY, and canonical assignment-report bytes; `assign` binds
declared targets to that verified evidence; and `material-resolve` builds the
complete card library from the admitted packs, resolves every declared
material and interface binding through `fs_project::resolve_bindings` under
`BindingRequirements::thermal_steady_v1()`, and retains each selected claim's
replayable usage receipt as a `solve-material-usage-receipt` artifact.
`flow-network` lowers the declared fan system and enclosure network to its
interval-certified operating point. `conduction` executes only when
`cooling.conduction` explicitly supplies one interior seed per region and a
partition of physical boundary laws. It replays the retained promoted meshes
and assignments, constructs and audits a labeled tetrahedral volume, binds
region labels to material-card conductivity, distributes declared regional
power, resolves declared matching-P1 contact faces against their interface
cards, solves the heterogeneous steady problem, and retains the canonical
stage receipt, a `solve-conduction-solution` temperature-field artifact, and
interface-resolution evidence when contact is present.
Projects without that optional declaration refuse with
`cli-solve-conduction-undeclared` (exit 4), and a project with no
`temperature-max` requirement refuses at `qoi` with `cli-solve-qoi-undeclared`
(exit 4): project defects named at the stage, not stage gaps. `cli-solve-stage-gap`
(exit 5) is reserved for a stage with no producer at all (none today); such
refusals name their producer bead and are retained as
terminal error operations.
The conduction receipt (schema v6) discloses the PUBLISHED mesh's quality
census (tets, vertices, min dihedral, max radius-edge, sliver and flat-tet
counts), the base mesh's flat-tet repair (found / repaired / unrepaired /
rounds — fs-mesh removes zero-volume tets by edge removal before the audit) and
the recovery budgets and statistics that produced it (memory budget, Steiner
cap and depth, segment/facet recovery rows). Radius-edge is disclosed, not
enforced: constrained refinement is bridge plan B2c. The stage REFUSES
(`cli-solve-conduction-mesh-quality`, exit 4) when a zero-volume tet survives
the repair or the smallest dihedral of any solved rung is below 1° (measured:
the reference heatsink sits at 7.3° after repair, 6.4° on its second uniform
rung); slivers between 1° and 5° solve and are disclosed. The Steiner cap
derives from the declared memory budget with the fixture default as a floor,
so identical inputs at the fixture budget mesh identically.

**The uniform h-ladder** (receipt block `ladder`, driver version 16). When the
project declares `solver.fidelity = "ladder"`, the stage solves the audited base
and then up to two further rungs, each one uniform 1→8 refinement of the
labeled complex (fs-mesh CONTRACT item 16: walls split in place with their
parent facet, labels replicated, volume preserved), taken while the next rung
still fits the declared memory budget (`memory_bytes / 256` tets). Declared
matching-P1 finite contacts are resolved again on every rung: the two
region-owned traces must have exact matching subtriangles, each contained in
its retained source face and covered by the same declared interface/card.
An unbound, ambiguous or nonmatching descendant refuses; no contact is dropped
or replaced by perfect contact. The receipt carries each rung's interface-pair
count and the final field's regenerated interface evidence and heat/jump report.
Other fidelities except `adaptive` solve the base once. Adaptive fidelity has
its own study block, described below. Other single-rung paths
report `"stop":"fidelity-single-rung"`. The block carries one row per rung (tets,
vertices, `h_m = (volume/tets)^(1/3)`, min dihedral, the QoI stage's functional
`t_max_k` — the nodal maximum over the ThermalLimit region — nonlinear and
Krylov iterations, final residual) and a `richardson` estimate over the last
three rungs by the Eça–Hoekstra procedure: `observed-order` (monotone
differences, observed order in [0.5, 2]: GCI half-width `1.25|f3−f2|/(2^p−1)`
and the extrapolated value), `data-range` (anything else three rungs can show:
`3 × (max − min)` of the QoI over the rungs, a bound on the observed variation
rather than an asymptotic estimate), `converged-exactly`, or `single-rung` /
`two-rungs` (no estimate). The published field, temperature range, energy
balance and quality census are the FINEST rung's. Cost is disclosed, not
hidden: three rungs are about seventy base solves (MEASURED 2026-09-02, debug
build, heatsink example: 685 → 5,480 → 43,840 tets in 0.5 s, 7.8 s and 130 s;
T_max 301.99578 → 301.99615 → 301.99610 K, i.e. `data-range`, half-width
1.1 mK), which is why the ladder is a declared study and not the default. The
linear solve's true-residual gate is two decades below the declared
`tolerance-rel` (floor 1e-13); the crate default of 1e-12 refused the second
rung at 1.37e-12 after 1,256 Krylov iterations.

**Adaptive thermal goals** (`solver.fidelity = "adaptive"`, receipt block
`adaptive`). Each candidate mesh is solved, uniformly enriched and solved
again. The enriched goal comparison uses the production residual and actual
transposed Jacobian, including nonlinear `k(T)` terms and fixed finite contact.
Declared independent airflow branches participate through their analytic
wall-to-reference feedback at the flow-network operating point. Each mesh
rebuilds its own wetted areas and card-derived coefficients, closes the primal
exchange, and binds the exact converged Robin references to its adjoint. Both
the enriched field and prolonged coarse field are marched through the same
fine-mesh air model before forming residual differences. The full coupled
transpose residual is checked after a bounded IQN-ILS interface solve; the
primal iteration trace is never differentiated. Branch reference and watt
gates remain independent. The history records branch count, interface work,
residual and tolerance; hydrodynamic feedback and mixing are outside this model.
The reference primal and dual both pass independent residual checks. Adaptive
primal stopping is tightened to satisfy the comparison's residual gate; step
stagnation alone cannot provide an accepted reference. The explicitly observed
linearization remainder is retained, including remaining dual/arithmetic error.
The goal is temperature at the enriched maximum's vertex. The exact owner-change remainder
`I(T_coarse)(v_max) - max(T_coarse)` is retained, so the reported goal remains
the declared region's actual nodal maximum. Acceptance requires both the
absolute dual change plus both absolute remainders and the independently measured
maximum change to fit `accuracy_rel * abs(T_max_coarse)`. This is observed
two-space agreement, not a continuum error bound; it adds no discretization
half-width or Verified authority to the QoI receipt.

When unresolved, incidence-distributed absolute dual residual contributions
mark cells accounting for half their sum. Conforming longest-edge-star
bisection preserves regional labels and source-face ancestry; every resulting
solve still passes the refined-mode 5-degree dihedral floor, and a bisected
candidate whose worst radius-edge ratio exceeds both 2.0 and its uniform
enrichment's own ratio is treated as a quality failure (bisection may not
degrade shape; admission of the Ruppert-refined base itself is dihedral-only
because radius-edge is disclosed, not a P1 element-quality measure). A zero marking
signal requests a global refinement, never success. The receipt retains each
comparison's mesh sizes, contact-pair counts, maximum, quality, primal/dual
residuals, linearization and owner remainders, nonlinear-Jacobian selection,
marked count, actual solved-mesh count, and stopping reason. Memory-derived tet
limits and a 17-primal-solve cap retain the last independently probed candidate
when further work cannot fit. These output counts are not an allocator peak
memory guarantee. Coupled comparisons have at most 32 interface sweeps, each
with the declared inner Krylov budget, and use a 1e-11 K primal reference gate.
Temperature-dependent contact resistance remains refused by material lowering.
An active material kink or invalid temperature span refuses the unique tangent;
there is no frozen-conductivity substitute. Marked and uniform meshes regenerate
the same declared matching contact and boundary/source ownership before solving.
If a marked candidate fails the quality floor, the already solved admissible
global enrichment is reused only when its own next probe fits; the receipt
counts these uniform quality fallbacks. Otherwise the last fully probed field
is retained with `mesh-quality-refinement-budget`. No low-quality candidate is
solved, and unrelated solver/material failures still refuse.

The remaining declared wall budget is checked between numerical operations.
Expiry records `cli-solve-conduction-adaptive-wall-budget` through the ordinary
stage-refusal path; a timing-dependent partial physics receipt is not published.
Prior completed pipeline stages remain recoverable, while unfinished adaptive
conduction restarts. An individual primal or dual solve does not preempt at a
wall deadline. Successful physics receipts replay without a timing cutoff.

The QoI receipt (schema v2) carries the ladder's estimate as the
**Discretization term** — the first measured term of the eight — through
fs-airflow's discretization-receipt seam: `"state":"interval"` with the
half-width in kelvin, the term's provenance role
`thermal-qoi-discretization-receipt`, and a `derivation` naming the method
(`richardson-gci`, `eca-hoekstra-data-range` or `bitwise-agreement`), the
ladder status, order, rung count, refinement ratio, safety factor and the
conduction receipt it was read from.

**Declared-input propagation** (conduction receipt block `propagation`, driver
version 24). When the project declares an envelope and exactly one
temperature-maximum requirement, the conduction stage re-runs flow-network and
conduction at base fidelity for each vertex of the declared interval inputs:
inlet/reference temperature at the envelope bounds crossed with fan-curve
pressure x (1 -/+ declared tolerance) for the Boundary-conditions term; every
card-derived coefficient x (1 -/+ the card's discrepancy allowance), inside
the conjugate fixed point, for the Model-form term; one joint worst corner
whose excess over the summed half-widths is added to the boundary term. Each input half-width is
the largest deviation of the region maximum from the base-fidelity nominal;
monotone response per input is assumed and disclosed. A refused vertex, or a
declared (non-card) coefficient, keeps its term NO-DATA under a specific
reason. Measurement is always negligible (the pipeline ingests no observation
data).

**Published linear solver error** (`solver_algebraic`, driver version 24).
For fixed linear conductivity, Robin references and matching contact, the
Solver-algebraic term now uses the actual published mesh and temperature field,
including an adaptive or ladder endpoint. `fs-conduction` computes outward
residuals and a verified inverse bound. Every competing region vertex is
enclosed using the whole-field bound, so the maximum may relocate without
inheriting an invalid selected-node claim. No second primal or selected-node
dual solve is needed. The receipt retains the maximum interval, residual and
inverse bounds, stability-proposal diagnostics and work counts. Missing inverse
evidence or a finite-range/structural refusal leaves this term NO-DATA. The
residual's evaluation roundoff is already included, not a second independent
uncertainty contribution. This encloses the stored floating-point linear
system only; the engineering term remains Estimated and establishes no
assembly, continuum, material or physical-validation claim. Nonlinear
conductivity and ambient radiation retain the explicitly named
100x tighter base-solve comparison as an estimate, without a linear enclosure
claim. Old driver checkpoints cannot resume into the changed analysis.

**Maximum-goal corrections** (`solver_control`, driver version 25). Each
supported fixed-reference linear mesh solve allocates 10% of `budgets.accuracy_rel` times its
pre-correction temperature rise above the coolest declared reference to
algebraic error. With no declared reference, the basis is explicitly the
pre-correction absolute temperature. This is a frozen requested-accuracy
allocation, not a measured discretization term. The maximum-controlled driver
performs bounded corrections only when the field misses that Kelvin target.
The correction budget is the existing linear iteration cap, with at most two
defect retries and eight iterations between checks; all are retained alongside
actual work, initial/final bounds, the allocation scale and the stop reason.
Every adopted field independently passes the original physical residual gate;
energy, Robin and contact fluxes are recomputed from it. Adaptive comparisons,
ladder rows, roundoff and QoI all consume this accepted field. A rejected
physical candidate preserves the original field and records the failed gate.
Missing inverse evidence, exhausted work, zero/unrepresentable accuracy scales
and unsupported nonlinear models never report goal completion. The
original physical solve still runs; this correction path makes no speedup claim.

**Coupled solid/air maximum error** (`solver_control`, driver version 26).
For temperature-independent conductivity with fixed flow, geometry and heat
transfer coefficients, the actual ordered air paths lower to their full affine
reference law. The published field is checked against `(A - B C) T = b + B d`,
including the contribution of prescribed temperatures and matching contact.
The Solver-algebraic term uses `outward-coupled-linear-maximum-enclosure`; its bound
and inverse include both the complete feedback residual and response-solve
errors. The verified coupled inverse may use whole-state contraction or the
port Schur check. Missing inverse evidence stays NO-DATA with its specific
reason; it never selects a frozen-solid enclosure or a tolerance
comparison. Port, transfer, response-storage, lowering and verification work
are capped, and all response columns share one linear-iteration allowance.
This path analyzes the field already accepted by the physical conjugate solve.
Its `assessment-only` receipt separates response preparation from zero
primal correction work and discloses whether the field meets the requested
accuracy allocation; coupled goal corrections remain unsupported. Preparation
refusals retain unknown response work as `null`, while evaluation refusals keep
work already spent. Missing-bound reasons fit the QoI uncertainty contract,
with full numeric diagnostics retained in `solver_control`. The bound
concerns the stored affine coefficients, excluding their assembly rounding,
exact exponential coefficients, nonlinear physics, hydraulic uncertainty,
discretization and physical validation. Driver identity prevents old receipts
from being resumed under these changed solver-error semantics.

**Verified inverse proposals** (driver version 28). Fixed-linear maximum
analysis can verify a numerical inverse when positive diagonal scaling cannot
establish a bound. The outward `I - A R` check includes every proposed column's
residual and rounding; storage, verification work and the single shared
stability-iteration budget remain bounded. This enables consistent contact
and enriched Robin systems without assuming the proposed inverse is exact.
The additional route admits at most 256 free unknowns and retains NO-DATA
when its structural or verification conditions cannot be established.

**Adaptive algebraic balance** (`solver_control`, driver version 29).
Each supported fixed-linear coarse/enriched comparison now allocates one tenth
of its current prospective Discretization term to the outward sum of both
mesh solutions' maximum-error bounds. The term is still the existing Estimated
`2 * max(abs(estimated_change), abs(measured_change))`, with its assumed
order-one scope. The sum rounds upward and the allowance rounds downward.
If necessary, either rung is corrected toward half the pair allowance.
Every changed field is independently checked against its original physical
residual gate, then prolonged and compared again before marking or stopping.
Each pair permits at most four correction rounds; each rung shares its original
primal-iteration cap across all corrections and permits at most four retarget
calls, with cumulative stability preparation capped by five original allowances.
History records the final bounds, pair allowance, correction rounds, physically
accepted field updates and additional correction iterations. Small fixed-linear
systems can establish a missing inverse using independently verified numerical
inverse columns, with at most 256 free unknowns and explicit storage/traversal
caps; all column solves share the existing stability-preparation allowance. The final `solver_control` uses `measured-adaptive-discretization` and
keeps the last internal correction target separately from the final allowance.
A covered pair with missing bounds, exhausted work or an unattainable allowance
stops as `solver-algebraic-budget` and supplies no Discretization term. It retains
the accepted physical field and its matching solver-error evidence. Unsupported
coupled, nonlinear and radiating models retain their prior adaptive comparison
with an explicit unsupported balance status; this change claims no continuum
bound or new correction support for those models.

The Roundoff term is fs-conduction's componentwise bound
`γ_k Σ|λ_i|(|b_i| + (|A||T|)_i)` on the published solve. It uses the exact
operator behind the field (including a radiating solve's combined partition),
with λ the adjoint at the hottest vertex of the requirement region. It
excludes air-network arithmetic and tracks only assembled magnitudes, not
element-level cancellation, so it is Estimated. Propagation re-solves skip it.
The Parameters term moves every material binding's declared relative
conductivity tolerance (fsim v6) to its lower and then its upper bound
together, re-solving at base fidelity with the card's sampled knots scaled
over the same temperature span. Its half-width is the larger deviation. With
no declared tolerance it stays NO-DATA, because the bound cards state no
conductivity uncertainty. The Geometry term moves every exterior boundary
vertex of regions assigned from an artifact that declares a surface offset
(fsim v7) by +/- the offset along its area-weighted outward normal, keeping
the same topology. The nominal lowering therefore stays valid, while the
solve, the source (watts conserved over the perturbed volume) and the
airflow wetted areas use the perturbed geometry. A tet collapsing below a
tenth of its volume refuses, and contact-coupled bodies refuse. While any
term is NO-DATA the verdict stays Estimated / indeterminate and the `no_claim`
text counts the NO-DATA terms. The
report projects the interval term's magnitude into the uncertainty table and
the QoI's `discretization_error`, and renders the ladder rows and estimate as
its convergence section (fs-ladder's vocabulary via fs-report); a single-rung
run shows no study.
A run id whose ops exist in the ledger but none of which is admissible under the
running driver (another driver version, a stage-receipt schema this driver does
not read, a non-deterministic or off-branch op) refuses with
`cli-solve-run-incompatible` (exit 4), listing how many ops were seen and why
each class was rejected; `cli-solve-unknown-run` is reserved for a run id with
no ops at all.

Card packs reach the run through the repeatable `--materials <pack>` and
`--interfaces <pack>` flags. Admission decodes each pack through its own
canonical `fs-matdb` envelope, orders the set by pack root, collapses
byte-identical duplicates, and refuses when two distinct packs reconstruct
one card — the card library is keyed by card identity, so there is no
last-one-wins path. Caller order and caller paths are not semantic: they
never reach the run identity or any receipt. The canonical set root **is**
bound into the run identity, so a different pack set is a different run
rather than a second answer for one run. Because the packs are retained
against the run's first operation, resume recovers them from the ledger
whichever stage was interrupted, and re-deriving the run identity from the
recovered set is what makes that recovery an attestation rather than a
restatement.

Stage checkpoints use the fs-exec **legacy v1** snapshot envelope
deliberately: the v2 envelope's expectation token is in-process-only and its
post-restart authorized resume path is explicitly unfinished (fs-exec
CONTRACT), while the v1 `LegacySnapshotExpectationV1` is reconstructible from
durable fields alone. The sealed driver state (run, project hash, consumption
totals, completed stages) is retained as a `solve-stage-state` artifact; its
ledger content hash and `open_expected` admission prove only bounded codec
integrity. Public `SolveDriverState`/`CompletedStage` construction and a valid
legacy envelope grant no resume authority.

Before opening a governor or publishing progress, resume collects every
retained driver state of the run (bounded checkpoint decode plus shape
validation), then independently re-attests candidates longest prefix first:
the first fully verified prefix wins, an equally long valid competitor refuses
as ambiguous, and strictly shorter candidates are re-attested only when every
longer one failed, because they cannot otherwise affect longest-prefix
selection. `solve --resume` re-attests by replay (each retained stage is
re-executed and its receipt must be bitwise the retained one); the export
verbs re-attest the same chain by sealed evidence (re-hashed receipts, sealed
outputs) and say so. Resume discovery requires the exact supported
solve-stage schema, field order, run, driver version, stage/ordinal pairing,
diagnostic state, and logical clocks before reading that operation's bounded
edge set; unrelated successful same-session operations therefore cannot
obscure a valid checkpoint by carrying a wide output set. A stage-shaped
candidate must then have positive, distinct, increasing operation ids; finite
nonnegative one-core totals (`core_s == wall_s`); exact canonical stage rows,
Five Explicits, deterministic main-branch execution context, and typed
operation-content identity sidecars; complete nontruncated edge sets with the
writer's exact edge-count seal; exact canonical stage receipts; one checkpoint
per stage; and a direct predecessor-checkpoint chain. Stage zero must consume
the exact strictly decoded canonical project and one exact import summary. Its
import operation must carry the exact project, Five Explicits, canonical
versioned IR field order and complete limits/policy shapes (including
named-group arrays), deterministic main-branch context, typed content sidecar,
exact typed edges, and strict whole-input-consumed receipt grammar. Unrelated
successful import operations are rejected by the cheap
row/IR/Five-Explicits prefilter before any bounded edge scan. Import evidence
is capped at 255 sources so the complete worst-case `4*N+1` edge set fits the
ledger's 1024-edge scan. Raw artifacts are streamed under their frozen
per-source/aggregate caps; opaque promotion receipts are streamed under a
4 MiB per-artifact cap. PLY and assignment-report evidence that solve must
parse has a 64 MiB per-artifact cap. Before payload verification, solve
requires the exact nine-line ASCII header emitted by `fs_io::ply::write_ply`
within a 256-byte header window, canonical unsigned vertex/face counts within
the frozen import limits, and a matching number of newline-terminated body
records. Solve then makes two bounded passes over that exact writer grammar:
the first decodes finite vertex coordinates and in-range `u32` triangle
indices; the second requires every numeric token to equal Rust's canonical
writer spelling with exact single spaces and newlines. No generic
`fs_io::ply::read_ply` or full-buffer `write_ply` call remains in solve
verification. Before reading any source payload, solve preflights the checked
total of the summary and every raw,
promotion, PLY, and report artifact against the smaller of the project memory
budget value and a 512 MiB hard admitted-input/work envelope. This is a
deterministic per-candidate resource-admission bound, not a claim that peak
allocations fit within that value. One invocation examines at most 8,192
visible operation IDs in at most 128 64-row pages and charges a single 1 GiB
cumulative work ledger. This is a deterministic conservative byte-equivalent
admission/work measure, not exact CPU accounting: controlled
operation/artifact input bytes are charged once, accepted UTF-8 copies and
direct byte-comparison passes add their processed bytes, derived render/receipt
outputs add their exact output bytes, and fixed page rows, sidecars, edge
items, descriptors, and project/entity/assignment items add documented proxy
charges. If captured history continues beyond either cap, or cumulative work
would exceed 1 GiB, solve returns
`cli-solve-work-envelope`; it does not report a false not-found result or
publish a refusal row that would worsen the same history bound. A genuine
candidate outside the separate per-candidate solve envelope produces
`cli-solve-import-envelope`. Competing independently valid longest checkpoints
are ambiguous and refuse rather than winning by row order. Only after this
complete attestation does resume re-charge recorded consumption, so the budget
continues instead of resetting.

## Output and exit contract

### Adaptive 3-D topology study

With the `sdf3-study` feature, `study <study.fsim> <ledger.db>` admits the
separate `fsim-sdf3-study :version 1` producer demonstrated by
`examples/marquee/bracket-3d-adaptive.fsim`. It executes the raw implicit-domain
3-D elasticity/SIMP solver with independent SI body loads, a continuation
schedule, goal-guided octree refinement and raw-density transfer. Each new
stage checks compliance and physical-volume directional derivatives and
restores the projected-volume constraint before optimization. Geometry and
Krylov budgets are cumulative; tree, per-solve and cooperative wall budgets
are explicit. A stopped run retains its last accepted solved model.

Sealed receipts use the existing study/report/package routing. Reports export
the retained cell densities, connectivity, physical nodal displacements for
each load, stage histories, gradient discrepancies and two-grid corrections.
Export does not run the optimizer. The package is a structural provenance
envelope; it confers no mechanics certificate. The admitted implicit boundary
is fixed, and comparisons are discrete numerical observations. There is no
continuum bound, optimum, manufacturing or peak-RSS guarantee. Durability is
at invocation completion; optimizer resume and the generic `--budget` override
are explicitly refused. A binary without the feature refuses this producer
with an actionable feature message.

### Fixed-background 3-D minimum-volume stress study

With `sdf3-study`, the same `fsim-sdf3-study :version 1` surface accepts
`(objective :type volume-fraction :sense minimize :unit "1")` with the explicit
`stress-limited-simp` optimizer in `examples/marquee/bracket-3d-stress.fsim`.
It uses `fs_topopt::sdf3::design::StressDesignStudy3` and its exact-discrete
stress adjoints, graph filter, projection and projected augmented Lagrangian.
The initial octree remains fixed; `maximum-level` must equal `initial-level`.
The same independently integrated body, reference-pressure and traction laws,
SI domain and clamp declarations as the compliance study apply.

The limit is the normalized volume-and-load-weighted qp von Mises aggregate,
not a sampled or continuum maximum-stress bound. The report retains both
sampled maxima separately. The density floor must stay above the ersatz-model
stress turnover. Two heterogeneous directional finite differences check the
initial stress and volume derivatives before optimization. These probes and
all optimizer trials spend the declared evaluation and Krylov allowances;
`max-stress-points` bounds the retained point family in one evaluation.

Every accepted update is sealed before further optimization. Accepted AL
iterates may be infeasible, so exports include both the last accepted state
and the least-volume accepted feasible incumbent. Only a converged, feasible
current iterate reports `completed`; iteration/evaluation/linear limits report
`budget-exhausted` with the actual constraint violation and finite retained
fields. KKT residuals apply to the last accepted iterate. `--budget N` caps
updates in this invocation. `report` and `package` export sealed artifacts
without solving. `study --resume study-<receipt-hash> <ledger.db> --budget N`
continues from the accepted optimizer state under the original source and the
same executable. It preserves multipliers, penalty, spectral step, inner
tolerance, previous outer violation, all accepted history and the distinct
least-volume feasible incumbent. Rebuilt geometry and one accepted-endpoint
stress re-solve (two when the incumbent differs) spend the original geometry,
Krylov, wall and evaluation allowances; earlier optimization updates and the
initial gradient gate are not replayed. The physical endpoint responses must
match the retained source before continuation. Reports count restoration
evaluations explicitly. The invocation cap adds updates to the retained count
without extending the original target or any resource allowance. Completed
or target-exhausted receipts return unchanged. Historical receipts without
optimizer state remain exportable but cannot resume. Failed restoration and
work lost after a crash cannot be durably charged; previous checkpoints remain
intact.
The structural evidence package grants no continuum, manufacturing, global
optimality or physical-validation authority.

### Fixed-grid 3-D robust load-space minimum-material study

With `sdf3-study`, `fsim-sdf3-study` also admits objective
`dilated-volume-fraction` and optimizer `load-space-limited-simp`. This uses
`fs_topopt::sdf3::robust::envelope::load_space::design::limits::LoadSpaceLimitedDesign3`:
the complete compliance Gram matrix and its checked spectrum cover the authored
zero-centered force ellipsoid `f=Fq`, `||q||_2<=1`. Every basis weight must be
one; force amplitudes and correlations belong to the physical columns.

The minimized objective is dilated material fraction; equilibrium and the
hard compliance cap use eroded material. The smooth spectral constraint does
not relax the physical cap. Initial spectral/material gradients pass the
existing finite-difference gate. Only accepted cap-feasible iterates may become
the least-material incumbent. KKT diagnostics continue to describe the last
accepted AL iterate, which can differ from the exported incumbent.

Geometry, supports, material, physical regions and projection thresholds stay
fixed. The original wall, Krylov and quadrature budgets cover all solves;
`max-optimizer-evaluations` excludes the separately reported initial gradient
probes. Accepted prefixes are committed before further updates and standard
report/package exports contain their actual fields without a re-solve.
`study --resume study-<receipt-hash> <ledger.db> --budget N` restores the shared
projected AL state, preserving its multipliers, spectral step, penalty schedule,
accepted history and least-material feasible incumbent. The original source
and executable are pinned. Geometry is rebuilt and accepted/incumbent responses
are independently re-solved; one or two endpoint evaluations are charged to the
original optimizer allowance, and all geometry/Krylov/wall work is cumulative.
Earlier search updates and the initial gradient experiment are not replayed.
The invocation cap adds updates without extending the original target. Historical
receipts without optimizer checkpoints remain exportable but cannot resume.
Failed restoration and uncheckpointed work lost after a crash cannot be durably
charged; previous results remain intact. Authority remains Estimated, with no
stress, buckling, continuum, mesh-independent optimum, geometric manufacturing
tolerance or experimental-validation claim.

### Scientific field export for 3-D studies

`report <study-run> <ledger.db>` also exports `<study-run>.vtu` for the adaptive
SIMP, stress-limited and robust load-space 3-D producers. The generic exporter
uses the sealed selected-design artifact and needs neither original source
files nor a new solve. It retains physical node coordinates, VTK voxel corner
ordering, raw and projected cell densities, cut-cell volumes and the available
displacement fields. Stress exports retain sampled relaxed stress; load-space
exports retain all robust density realizations, basis displacements and the
worst-combination displacement. The source artifact hash and explicit export
status accompany the path in text and JSON output. No field is invented when
there is no solved baseline or no feasible selected design.

These cells are the active octree background boxes, not a reconstructed cut
surface or manufactured solid. Export preserves the selected incumbent's
Estimated authority. Bounded allocation and validation reject malformed field
dimensions, nonfinite data and inconsistent voxel geometry.

### Native cooling probability study

`study <study.fsim> <ledger.db> [--budget N]` recognizes
`fsim-uncertainty-study :version 1` and executes the
`native-cooling-uncertainty-v1` producer. The example is
`examples/native-uncertainty/study.fsim`. It binds 1..=32 explicitly independent
uniform inputs to an admitted native cooling project, with 2..=256 fixed-count
Monte Carlo samples and an explicit sampling seed and total wall allowance.
Every accepted observation comes from an ordinary geometry import, material
resolution, staged solve and retained `temperature-max` QoI. Probability support
does not widen the base's operating envelope or card domains, change its solver
policy, or convert an engineering interval into a probability distribution.

Each observation retains its ordinal, parameter values, child run, project hash,
QoI receipt and kelvin value. Only a completed fixed-count study publishes its
mean, sample standard deviation, descriptive standard error, empirical quantiles
and pass fraction against the requirement's limit minus its explicit margin.
Partial, cancelled and refused runs retain accepted observations with
`statistics: null`. A refused physical sample terminates the study at that
ordinal; it is never clipped, counted as a pass, replaced or skipped.

`--budget N` caps additional samples in this invocation; zero admits and retains
the model without evaluating it. `study --resume study-<receipt-hash> <ledger.db>`
loads the bounded statistical checkpoint, verifies accepted child QoIs, and uses
retained project, geometry and card bytes without reopening the source paths.
The original total sample and cooperative wall allowances continue across
resume. Completed and refused terminal receipts return unchanged. Study outputs
are committed together at invocation completion; completed child solve evidence
has its own normal solve-stage durability. Input bytes obey an explicit fraction
of the base memory budget; this is an admission cap, not measured peak RSS.

The existing `report` and `package` verbs export retained study artifacts.
Replay on the same numerical implementation reproduces the completed report;
invocation wall use and checkpoint ancestry remain separate receipt data.
All version-1 statistics are Estimated, advisory descriptions of the declared numerical
input model. They do not supply confidence sequences, optional-stopping
decisions, physical validation, continuum-error certificates or compliance
signoff, and they do not replace the child's engineering uncertainty budget.

Version 2 requires a predeclared `bernoulli-mixture` compliance policy with
required probability, alpha and minimum sample count. The model's raw
`temperature-max <= limit - margin` indicators enter fs-eproc's fixed
Beta(1/2,1/2) likelihood mixture. At or above the declared minimum, a lower
confidence bound at least the target means `meets-probability-target`; an upper
bound strictly below it means `below-probability-target`. Either stops with
`decision-reached`, termination `probability-target`, and exit 0. An unresolved
lifetime sample cap returns `budget-exhausted`, termination
`lifetime-sample-budget`, and exit 6. Resume cannot change the retained policy
or add to the declared lifetime cap. Version-1 bytes and semantics are unchanged.

The version-2 `compliance` report contains both probability-confidence endpoints,
the empirical frequency, policy, count and decision; `statistics` stays null.
All-pass or all-fail samples retain nonzero probability uncertainty. Mathematical
time-uniform coverage requires the declared Bernoulli model; deterministic
floating-point inversion is Estimated and is not an outward-rounded proof.
No mean control variate, QMC-point confidence claim, replacement of failed
samples, physical-error bound or safety-signoff authority is introduced.

Version-3 frozen mean controls also compose with `:method latin-hypercube`
and its explicit replicated layout. Coordinate-secant or nominal-adjoint
calibration finishes before probability observations; its completed probes
and frozen coefficients survive the existing native resume path. Only complete
independent hypercubes enter the controlled mean and between-replicate error.
The raw `lhs` report, retained child temperatures and numerical compliance
remain unchanged. `mean_control.estimate.statistical_unit` explicitly names
`complete-independent-latin-hypercubes`. Neither adjusted samples nor their
errors acquire physical-QoI, confidence-sequence or compliance authority.
See `docs/NATIVE_LATIN_HYPERCUBE.md` for the executable study declarations.

For a native transient project, version-3 `nominal-adjoint` mean controls select
the existing final-time accepted-trajectory adjoint. Calibration binds the sealed
run's storage law, final temperature, horizon and fine-grid state count before
using its accumulated derivatives. Power coordinates retain the authored duty;
a workload schedule that replaces static power gives that static coordinate a
zero derivative. Prescribed convection coefficients convert the sampled
`dT/dlog(h)` to `dT/dh`. Natural reservoirs, air branches, fan banks and radiating
surfaces retain their distinct physical identities. Constant-capacity storage
also supplies prescribed Robin temperature/coefficient controls; the enthalpy
producer does not, and varying unsupported coordinates refuse with a secant
alternative. Enthalpy storage supplies independent uniform prescribed-bath
temperature controls held through the full selected trajectory. Calibration
checks the exact nominal kelvin value, boundary ownership, available status and
fine-grid coverage; shared prescribed vertices cannot supply independent
coordinates. Constant-capacity fixed-temperature controls remain unsupported.
The native bath regression checks a fresh physical finite difference, unchanged
raw sample statistics and source-independent resume of the frozen calibration.
Sampled-peak and steady adjoint requests are never silently reinterpreted as the
final-time objective.

Raw probability trajectories, compliance and retained temperatures are unchanged.
The frozen mean control and calibration receipt use the existing budget and
source-independent recovery path. The derivative remains a selected numerical
branch and does not certify smoothness at a tied maximum or physical accuracy.
`data/reference-project/cooling-transient-uncertainty.fsim` is a runnable example
with declared power and convection laws over the existing transient model.

### Normalized thermal study

`study <study.fsim|study.json> <ledger.db> [--budget N]` executes the existing
scalar Poisson/CutFEM radius optimizer. The canonical example is
`examples/marquee/thermal-2d.fsim`: unit source, zero temperature on all outer
and hole boundaries, unit square, dimensionless compliance, fixed hole
centers and area equality. The executable reader rejects fields that do not
survive canonical recognition, including ignored geometry or unit declarations.
The separate free-boundary elasticity producer executes the real level-set
optimizer; q61wp.16 remains partial against its broader error-control goals.

For the plain or projected-volume single-load `elasticity-2d` producer, append
`(assessment :type elasticity-dwr :max-solves-per-attempt 2)` after the optimizer
section to assess the exact final accepted design. This explicit opt-in requires
at least 256 MiB admitted memory for the coarse/enriched state pair and the
existing bounded mesh-level envelope 2..=5. It authorizes two assessment solves
per attempt, each capped at 60,000 CG iterations with a 1e-12 recomputed
Euclidean relative-residual target. Memory is an admission allowance, not
measured peak RSS. Single-load projected-stress mode still refuses this
assessment request. An explicit weighted-sum load family is supported below.

The final optimizer update is durably retained before assessment starts.
Cancellation and lifetime-wall-budget checks bracket the indivisible DWR phase.
An interrupted attempt preserves the accepted geometry, records elapsed wall
time and leaves assessment pending; a numerical refusal retains its reason and
elapsed charge with `numerical-failure`. Resume can retry assessment on the
already complete optimizer checkpoint without another geometry update. A
completed terminal assessment is reused without solving again. Disabled
assessment leaves existing output unchanged.

The retained report, receipt continuation and HTML expose signed estimated
goal error, absolute cell-indicator sum, coarse/enriched compliance on the
same bilinear geometry, residual contribution breakdown, actual solve
residuals/iterations and geometry identity. Absolute indicator sum is marking
mass, not an interval radius. No certified continuum-error bound, adaptive
optimization, physical validation or optimum is claimed. An unmet material
constraint remains `constraint-unmet` even if error assessment succeeds.

Projected-volume final assessment binds the actual feasible current design,
including the prepared baseline when a bounded search accepts zero updates.
The existing constraint evidence records an optional `optimizer_terminal`
before the assessment starts. This distinguishes a finished
`no-feasible-descent` search from a paused optimizer: cancellation or numerical
refusal can retry only the assessment, without repeating projection, optimizer
recovery solves, candidates or accepted updates. An assessed stalled design
remains `no-feasible-descent`; estimated numerical error does not establish
optimization convergence. Prescribed regions and refinement-origin geometry
are retained unchanged, and the assessment uses the declared final grid.
Disabled assessment omits this terminal marker and preserves existing behavior.

The `final_dwr_assesses_the_retained_design_and_resume_reuses_completed_work`
binary regression exercises numerical reporting and final-checkpoint recovery.
Projected-volume regressions additionally cover a genuinely accepted feasible
endpoint, a zero-update feasible stall, stops before and after the real DWR
phase, retry without search/recovery callbacks, and completed-result reuse.

The existing projected-stress independent-load study also accepts this same
assessment when its aggregate is `weighted-sum`. Its explicit
`:max-solves-per-attempt` must be exactly twice the declared number of load
cases (including the primary case and every zero-weight case), so two cases
require `4`; the native limit of 16 cases caps an attempt at 32 solves.
Cases are assessed sequentially on the unchanged final level set with the
same per-case traction support/vector, material and solver as the optimizer.
The 256 MiB minimum and mesh-level limit 5 apply; only indicator evidence is
retained between cases, not previous displacement systems. Coarse compliance
must reproduce every retained case and the weighted sum bitwise before the
assessment is published. Results retain per-case and weighted coarse/enriched
compliance, signed estimates, absolute indicator mass and residual breakdowns,
plus each actual solver residual and iteration count. Weights are not
probabilities; forces are never summed before equilibrium. Worst-weighted-case
objectives refuse this linear-sum assessment rather than hiding active-case
switches. Compliance DWR does not estimate or certify stress error.

Only an area- and sampled-stress-feasible `completed` or `no-feasible-descent`
endpoint is assessed. A zero-update feasible baseline is eligible. Paused
optimization, exhausted work during stress restoration and infeasible
restoration stalls perform no assessment solves; the latter retain an explicit
assessment refusal with the existing optimizer status. The load-family history
retains its `optimizer_terminal` before assessment, and cancellation is checked
before each two-solve case and before final publication. Each case remains
indivisible. Resumption retries only assessment of the retained endpoint, even
with zero recovery-solve allowance; it neither restarts candidate search nor
re-solves optimizer recovery. Interrupted attempt wall time remains charged,
and completed assessments are reused byte for byte. Numerical assessment
refusals retain the accepted geometry and permit assessment-only retry.

Every accepted transition retains the source, iteration log, exact parametric
geometry, SVG/HTML report, JSON summary and format-9 package in one ledger
transaction. A `study-<receipt-hash>` run ID names that immutable checkpoint.
`report` and `package` accept this ID plus its ledger and export retained bytes
after checking their hashes, kinds and producing-operation lineage. Package
export re-runs the structural checker. A partial study report remains visibly
partial and Estimated; neither a checker pass nor iteration exhaustion proves
an optimum or a guaranteed PDE error bound.

`study --resume <study-run-id> <ledger.db> [--budget N]` uses the returned
`run_id` (`study-<receipt-hash>`), not the stable `receipt.study_id`, and replays the retained prefix
and requires the same trace before advancing. It recharges prior wall use and
also meters replay; it does not reset the original budgets. `--budget N` caps
new iterations in this invocation, while the declared total iteration limit
continues across resume. Replay checks remaining wall time between solves and
refuses with the retained checkpoint if it cannot reconstruct the prefix within
budget. Successful resume charges replay time; an aborted replay publishes no
new checkpoint, so its failed-attempt cost is not durably accumulated.
Wall/core time uses the session governor. The current
memory envelope is admission-only (at least 128 MiB, mesh level at most 5,
at most 32 holes and 256 steps), with no measured RSS enforcement. Cancellation
is polled between bounded transitions through the internal gate seam; the
binary has no OS signal handler. Pre-cancellation publishes no operation.
Budget exhaustion exits 6, cancellation exits 130, and ledger/export failures
cannot report success. A geometry artifact is the exact circle/plate
parameterization, not a sampled SDF grid. Full .20 acceptance remains open.

- stdout carries final result records only;
- stderr carries diagnostics and, for solve in JSON mode, one
  `frankensim.cli.solve-progress.v1` JSON-line per completed stage and per
  budget warning;
- JSON mode emits one complete object per line in deterministic field order;
- text mode emits stable `key=value` result rows and `ERROR`/`FIX` diagnostic
  pairs;
- exit `0` is success, `2` usage, `3` input I/O/encoding/size, `4` project or
  run refusal, `5` unavailable product stage or solve stage gap, `6`
  budget-exceeded honest partial, and `130` explicit cancellation observed
  before publishing a cinematic result. The older solve driver retains its
  existing refusal-class cancellation surface; this bead does not silently
  redefine that separate compatibility boundary.

Diagnostic codes and fix text are machine-facing compatibility surface.
Human prose may improve without changing a code or exit class.

Each machine code maps to exactly ONE exit class regardless of which layer
catches it. In particular the card-pack resource caps
(`cli-solve-card-pack-count`, `cli-solve-card-pack-size`) always exit `3`
(input/resource cap): such invocations match the documented grammar — the
card-pack flags are unbounded repetition — so `2` (usage) would misreport
correct syntax, and the cap is enforced at parse, read, and admission layers
without changing the class (bead p63op).

## Error model

Command-line failures are `Diagnostic` records with stable codes on stderr
plus a result record naming the status. Import failures are
`GeometryImportRefusal` values; solve failures are `SolveRefusal` values
(code, optional stage, what, fix, optional owning-bead dependency, optional
run identity, optional recorded ledger operation). Both refusal families are
retained in the ledger when a ledger is open and the refusal happened inside
a stage; recording failure degrades loudly into the refusal text rather than
silently. Budget exhaustion is not an error: it is a successful
`SolveOutcome` with the `budget-exceeded` status, because the work performed
is real and durable.

Solve refusal codes: `cli-solve-project-invalid`, `cli-solve-budget`,
`cli-solve-session`, `cli-solve-ledger`, `cli-solve-ledger-transaction`,
`cli-solve-import-evidence`, `cli-solve-import-envelope`,
`cli-solve-work-envelope`, `cli-solve-assignment`, `cli-solve-capability`,
`cli-solve-stage-gap`, `cli-solve-conduction-undeclared`, `cli-solve-qoi-undeclared`,
`cli-solve-run-incompatible`, `cli-solve-conduction-mesh-quality`,
`cli-solve-cancelled`, `cli-solve-run-id`,
`cli-solve-unknown-run`, `cli-solve-resume-identity`,
`cli-solve-resume-complete`, `cli-solve-resume-budget`,
`cli-solve-ledger-path`, `cli-solve-ledger-open`, `cli-solve-budget-exceeded`,
`cli-solve-usage`, `cli-solve-card-pack-read`, `cli-solve-card-pack-size`,
`cli-solve-card-pack-count`, `cli-solve-card-pack-source`,
`cli-solve-card-pack-decode`, `cli-solve-card-pack-identity`,
`cli-solve-card-pack-conflict`, `cli-solve-material-nonfinite`,
`cli-solve-material-receipt-envelope`.

The `material-resolve` stage additionally propagates `fs_project::bind`'s own
violation codes verbatim (`project-material-card-unknown`,
`project-material-state-mismatch`, `project-material-unbound-region`,
`project-binding-domain-uncovered`, and the rest of that vocabulary). This is
deliberate: flattening every binding cause into one CLI code would destroy
the diagnostic identity the binding layer already defines.

## Invariants

- Argument order never changes semantic output except for the documented
  position-independent `--json` flag.
- A successful validation has exactly zero `DecodedProject::findings()` and
  no lenient default or canonicalization receipt, because the CLI uses strict
  readers.
- Every refusal has a non-empty code, message, and suggested fix.
- User-controlled strings are escaped before JSON emission; every JSON record
  is one line.
- No unavailable `report`, `package`, or `compare` invocation writes a run, report, package,
  checkpoint, or ledger artifact. A solve run writes only through its staged
  operations, and a stage the driver cannot execute refuses with retained
  evidence instead of substituting a skeleton stage.
- Solve run identity is derived before any side effect, printed in every
  result record, and bound as the `session` of every operation the run
  writes; the same project in the same workspace always derives the same run.
- Each completed solve stage commits atomically: stage receipt, sealed driver
  state, exact sealed lineage links, and terminal outcome land together or not at all
  (fs-ledger's crash battery is the durability authority).
- The driver state sealed after stage N lists exactly stages 0..=N with their
  positive, distinct, increasing operation ids and receipt hashes. Resume uses
  the unique longest fully re-attested retained prefix; a malformed candidate
  cannot lend authority to its decoded fields, and equally long valid
  competitors refuse as ambiguous.
- The slice-1 single-core driver records finite nonnegative consumption with
  `consumed_core_s == consumed_wall_s` in every checkpoint. A non-finite or
  backward fresh clock refuses before a checkpoint or lineage seal is
  published.
- A pre-cancelled solve publishes nothing; cancellation between stages leaves
  the completed prefix durable and resumable with bit-identical stage
  evidence relative to an uninterrupted run. Cancellation during fresh or
  resume evidence verification returns `cli-solve-cancelled` directly with no
  refusal operation, stage receipt, checkpoint, seal, progress line, or
  terminal-success publication. A zero-stage prefix recommends a fresh retry;
  `--resume` guidance is emitted only after at least one durable stage exists.
- Budget continuity survives resume: recorded consumption is re-charged
  before any new stage runs, and a resume that would already be exhausted
  refuses rather than resetting the meter.
- Import accepts exactly one source for every exact project geometry row and
  no extras; declaration order, not insertion order, determines retained rows.
- Raw bytes must reproduce the project row's FNV source hook and exact
  `fs-io` parser version before promotion.
- Import refuses a caller-owned ledger transaction so its artifacts,
  extension rows, lineage, and terminal outcome commit or roll back together;
  solve refuses a caller-owned transaction for the same reason.
- Every import and solve operation freezes project-derived units, seed,
  budgets, versions, and capabilities in the ledger Five Explicits. Frozen
  import IR also binds every import/assignment resource limit and, in project
  declaration order, exact source-row identity, source unit, repair cap or
  STEP root and target-spacing bits, and ordered named-group mappings. Caller
  path labels do not enter semantic identity. Frozen solve stage IR binds the
  stage name, ordinal, run identity, project hash, and driver version.
- Solve admission checks frozen raw-source byte caps against streamed actual
  lengths (including checked aggregate addition), assignment request and
  predicate-work caps against the exact project and promoted PLY face counts,
  mesh vertex/face and named-group face-range caps against canonical retained
  PLY, and selected-face totals against canonical retained assignment-report
  rows.
- A STEP success retains both lower receipts and writes/assigns the exact
  repaired soup whose counts and fingerprint appear in the import receipt.
- Caller-supplied named-group face ordinals are never laundered across
  face-removing repair. Duplicate/degenerate removal with non-empty groups
  refuses until an adapter supplies an explicit remap or callers use geometric
  selectors. Orientation-only repair, vertex compaction, and appended hole
  faces preserve existing face ordinals.
- Successful import never truncates assignment results: resolver/report count
  must equal the prepared geometry count.

## Determinism class

Argument parsing, validation formatting, unavailable-stage refusals, and
geometry import identities are pure functions of arguments and input bytes
except for the explicit file/ledger boundaries. They read no clock, RNG,
network, or machine state.

Solve stage receipts, driver-state payloads, run identities, and operation IR
contain no wall-clock values; operations use logical times, so deterministic
stages reproduce identical artifact content hashes across independent runs
and ledgers (conformance-tested). Wall time enters only budget accounting
(the meter and the `wall_s` reporting fields) and is explicitly not identity.
Whether a run ends `completed` or `budget-exceeded` is wall-clock dependent
by design — a budget is a wall phenomenon — but every artifact the run
retains on either path is deterministic.

## Cancellation behavior

Validation is bounded by the 16 MiB CLI input cap but has no asynchronous
cancellation surface.

Geometry import has explicit source-count, per-source-byte, aggregate-byte,
and assignment-work caps. It polls the supplied `fs-exec::Cx` before source
work, per source, before and after promotion, and before ledger publication.
A pre-cancelled attempt publishes nothing. Once the atomic SQLite transaction
begins, the bounded ledger calls finish or roll back; cancellation does not
leave a partial successful operation.

The solve driver owns a caller-supplied `fs-exec::CancelGate`. It polls at
every stage boundary and throughout solve-owned evidence materialization,
stream verification, incremental UTF-8 validation/copy, JSON cursor
advancement, PLY body scanning, canonical payload decoding, writer-token
verification, named-group face-reference and duplicate validation, retained
project comparison, import-IR canonical comparison, and reconstructed
stage-receipt comparison. Solve-owned artifact readers poll at phase entry
before ledger dispatch, before and after every shared controlled tile, and
again after the controlled reader returns before classifying a missing,
stopped, corrupt, or complete result. Byte scans and comparisons use an
input-byte checkpoint cadence no wider than 64 KiB and repeat the current
checkpoint after inspected work before returning a mismatch or parse error.
It also polls after every successful stage body immediately before clock,
charge, or ledger publication. Cancellation returns `cli-solve-cancelled`
directly; it never passes through the retained-refusal writer, so the already
durable prefix is unchanged.

The 64 KiB value is an input-byte checkpoint cadence, not a wall-clock latency
claim. The controlled fs-ledger reader hashes and delivers one shared
at-most-64-KiB tile before invoking solve; a callback break prevents later
tile hashing/delivery. Its SQL engine may already have materialized the
containing storage row of up to 4 MiB. JSON cursors poll at phase entry, at
each crossed fixed boundary, and at the final current position before their
wrapper classifies success or syntax failure. UTF-8 conversion and output
reservation results, byte-comparison mismatches, named-face range failures,
and PLY line/header/payload failures are likewise classified only after the
post-work checkpoint. Import named-group duplicate checks preserve the
writer's face and group order: they reserve each complete duplicate set once,
scan face references in tiles representing at most 64 KiB of `u32` values,
and bound every solve evidence label at 4 KiB. The enclosing operation IR
inherits fs-ledger's guarded `MAX_OP_IR_BYTES` 1 MiB ceiling.

Allocation and the named duplicate-set insertion itself remain bounded opaque
intervals. Retained project parsing through `fs_project::parse_sexpr` is one
opaque call over at most `MAX_PROJECT_BYTES` (16 MiB), bracketed by phase-entry
and completion checkpoints; solve makes no intra-call latency claim for that
lower-layer parser. Canonical import-project rendering and canonical stage
receipt construction are also captured as bounded opaque calls bracketed by
entry/completion checkpoints, followed by tiled byte comparison. State-envelope
decoding is bounded by 4 MiB and bracketed before its result is classified.
The PLY preflight, payload walk, and writer-spelling wrappers each repeat a
bounded post-result checkpoint, including on short malformed input. Canonical
JSON integer and finite-float token widths are capped before standard-library
numeric conversion; the finite writer bound is 327 bytes, including the
negative smallest subnormal. The binary does not yet install an OS signal
handler; gate wiring exists at the library seam and is conformance-tested
there.

## Unsafe boundary

No unsafe code.

## Feature flags

- `sdf3-study` enables the numerical 3-D compliance and stress-constrained
  study producers described above.
- `thermal-verification` enables an explicitly requested
  `temperature-volume-mean-bound` report on native conduction solves. Declare
  exactly one report without an output `:region`; the selected volume is the
  existing `temperature-max` requirement region. The separate mean output
  does not replace or certify that maximum-temperature requirement.

The mean producer bounds the nominal linear PDE on the final published
polyhedral mesh domain using the existing equilibrated primal and dual flux
verifier. It retains the original field, source, boundary data, material cards
and exact pinned property receipts. A global constant representation is
admitted only from an explicitly unconstrained scalar constant source claim;
finite material validity spans and flat sampled curves remain refused.
Matching finite-resistance contacts preserve the original paired faces,
independent temperature traces and source query receipts. Their
resistance-weighted temperature jumps enter the existing broken-H1 continuum
majorant; the selected region is bounded in the complete contacted assembly.
Contact resistance must likewise come from an explicitly unconstrained scalar
constant claim. Finite validity spans, sampled/state laws, caller-supplied maps
and nonmatching contacts require their own source-validity admission.
Transient, nonlinear, radiation and natural/forced-convection feedback models
require their own bounds and are refused by this producer.
Cell and iteration allowances bound work, and cancellation is cooperative;
no intra-kernel wall-time guarantee is claimed.

The HTML/JSON report projects the producer's Verified nominal numerical
interval as a distinct QoI, while all unknown physical uncertainty stays
unknown. The standalone package retains the exact interval and source
receipt address as an Estimated statement: its existing structural checker
does not independently authenticate or reverify the continuum producer.
Builds without the feature explicitly refuse this output request.

Native steady mixed convection is available without a feature flag. Separate
`natural-convection` and `airflow-convection` targets exchange heat through the
same solid on auto, ladder and adaptive fidelities. Every natural
coefficient update encloses a complete existing multi-branch IQN-ILS air
exchange, including its final solid response and decomposition check. The
existing conductivity, matching-contact and ambient-radiation solvers remain
active within each response. This is separate natural reservoirs and forced
air on disjoint boundary patches, not an empirical combined coefficient for
one surface or a buoyancy modification of the flow network.

The outer loop allows at most 80 evaluations; each inner air exchange retains
its 100-iteration cap plus one final solid response. Acceptance requires every
natural h to change by at most 1e-10 relatively and every natural patch's
actual `h(T_wall) A (T_wall-T_ambient)` to match its applied Robin heat within
`max(1e-12 W, 1e-6 * max(|Q_actual|, |Q_applied|))`. Each air branch keeps
its own independent heat gate. Natural reservoir heat and radiation never
enter air enthalpy transport. Heated and cooled Churchill-Chu walls retain
their signed flux; equilibrium and card-domain refusals remain explicit.
Cancellation and the remaining invocation wall budget are checked between
bounded responses; failure publishes no partial conduction field.

The `natural_convection.coupling` receipt records outer/inner work, the admitted
coefficient tolerance, and final coefficient and natural heat discrepancies.
The returned derived boundary is the one actually applied to the published
field. Inner conductivity/radiation work retains its separate existing budgets.
The linear inverse-error certificate is unavailable and the whole-model
tolerance estimate remains. The result retains Estimated authority.

Driver 79 admits adaptive mixed-convection goals, including optional ambient
radiation, on at most 64 uniform Robin traces. Each enriched comparison
re-evaluates the original Churchill-Chu law at both the solved fine field and
prolonged coarse field, and marches both fields through every physical air
path. The existing core FGMRES transpose includes the signed natural
coefficient slope and the full streamwise air-reference Jacobian together;
radiation also differentiates its combined coefficient and weighted reservoir
reference. Natural and forced-air targets must remain distinct. The shared
primal exchange uses the adaptive air-reference gate and tightens the natural
coefficient tolerance below the goal residual tolerance. The independently
reconstructed primal and dual residuals must pass before a comparison can mark
cells or supply the observed discretization estimate.

Nonradiating mixed history rows name `natural-air-full-transpose`; radiating
rows retain `area-mean-radiation-full-transpose`. Both disclose
`air_reference_feedback` and `natural_coefficient_feedback`. The existing
nonlinear and maximum-owner remainders, refinement quality gates and
last-probed-field budget stopping remain active. `cli-solve-conduction-boundary-adaptive`
refuses inconsistent retained boundary laws or insufficient feedback storage.
This is an Estimated two-space comparison, not a continuum bound or a
linear inverse-error certificate. Natural reservoir and radiative heat never
enter the air enthalpy balance.

Driver 63 admits steady mixed nominal adjoints through the existing complete
material/radiation/air transpose. Signed natural wall and ambient coefficient
partials share that operator with stream-wise air-reference feedback; no
coefficient or air reference is frozen for differentiation. The unchanged
accepted field must pass the complete constitutive residual gate. Power,
material, fan-speed, air-inlet and radiation controls use the resulting dual;
the boundary-extended request also includes prescribed-temperature lifts.
These are local Estimated derivatives without inverse or gradient-error
certificates. Two physical re-solve regressions in
`tests/solve/mixed_convection/adjoint.rs` cover both natural heat signs,
smooth conductivity, radiation and the prescribed-temperature lift. `tests/solve/mixed_convection.rs`
supplies signed-reservoir, independent correlation/energy, radiation
composition, adaptive feedback, deterministic replay and unresolved-budget
regressions.

Native transient conduction is available without a feature flag through the
schema-v13 `conduction.transient` declaration. The runnable input and commands
are in `data/reference-project/cooling-transient.fsim` and `TRANSIENT.md`.
Both backward-Euler grids start from the declared initial state; their
combined step cap includes the coarse and nested fine solves. The ordinary
field, thermal QoIs and requirement margin describe the fine final-time
state. Temporal comparison is separate from the unknown spatial error, and
no steady algebraic/roundoff certificate or steady input-propagation solve
is attached. See the retained step/storage receipts for energy and work.
Bound temperature-dependent conductivity selects the shared endpoint
Newton/FGMRES solve with its full k′(T) tangent. Each accepted endpoint retains
nonlinear work and the joule residual/threshold; material-domain and iteration
failures publish no partial field. The default storage model uses declared
constant capacities. The native driven-face test checks an
independently derived transient balance and a frozen-conductivity twin.

An explicit `transient.enthalpy` selects the shared heterogeneous P1
total-enthalpy solver instead. Each conduction region supplies a bounded
equilibrium chart, fixed reference density, attributed card identity, source
text and initial specific enthalpy. Sensible and latent energy use the same
specific-enthalpy state; no separate latent account or apparent heat capacity
is added. A temperature plateau preserves the declared liquid mass fraction
through its enthalpy coordinate. Pure solid and pure liquid charts are also
admitted. Distinct regional charts require separate interface vertices and a
declared contact law; shared-vertex mixtures are explicitly refused.

Schema v15 permits `:initial-temperature from-regional-enthalpy` for this
storage model. Every body starts from its own declared specific enthalpy,
resolved through its own chart; separate bodies may start at different
temperatures. The historical kelvin quantity remains a uniform-temperature
assertion and must agree with every regional enthalpy. Constant-capacity
storage still requires that quantity. No initial temperature or latent state
is inferred from an ambient reservoir or contact neighbor.

Initial regional temperatures are checked against the sourced base
conductivity before time integration. The enthalpy receipt names
`initial_state`, records the actual initial temperature on every material
row, and uses `initial_temperature_k: null` at the top level for regional
mode. Both grids restart from the same physical enthalpy vector. Matching
contact transfers heat between the bodies without adding an external energy
source. Driver version 62 distinguishes this initialization and receipt
semantics from older retained solves.

Every coarse and fine grid starts from the original regional enthalpies.
Each accepted endpoint uses actual endpoint conductivity and the original
physical old enthalpy throughout Newton, FGMRES and line-search trials. The
shared nonlinear residual and independent joule-energy gate must both pass
before history advances. Prescribed temperature, Neumann/Robin laws, matching contact and
event-aligned regional power schedules remain active. The declared Newton,
backtracking and existing linear allowances bound each response; checked
products and memory estimates include both trajectories and retained evidence.
Cancellation or any failed domain, work or energy gate publishes no partial
conduction field.

`fs-cli-transient-enthalpy-v1` retains the actual material and chart identities,
sources, reference densities, accepted step energy, true nonlinear/Krylov
residuals and work, mass-weighted phase fractions, and final nodal enthalpy,
liquid mass fraction and reference mass. The normal field and report use the
fine-grid final temperature. Coarse/fine enthalpy and phase differences remain
visible even when both temperatures agree on a plateau; equal temperatures
do not become a zero-error temporal claim. This is Estimated fixed-domain
transport, with source attribution rather than authentication of supplied
phase data. Prescribed temperatures require a unique enthalpy on every incident
phase chart; latent-plateau ambiguity and out-of-domain baths refuse. The
original initial enthalpy remains physical history. Fixed endpoint constraints
therefore charge their complete startup storage to the bath once, including
when the whole body is prescribed. Step rows retain both the applied
`dirichlet_in_w` and fresh `physical_dirichlet_in_w`; radiation, natural and
airflow audits use the latter with physical reservoir fluxes and actual air
heat gain. The final solver report counts only free degrees of freedom.
Solve driver 80 binds these bath and constrained-row semantics into both solve
and native uncertainty identities; earlier retained evidence cannot resume
under the changed physical and derivative model.
Native trajectory adjoints are admitted under the scope below. Expansion,
motion and phase kinetics are not inferred.
`data/reference-project/cooling-enthalpy.fsim` and `TRANSIENT.md` give a complete
synthetic melting example through validate/import/solve/report. Driver version
57 separates these semantics from older retained solves.


Native total enthalpy also supports the existing quasi-steady airflow
convection path, with optional sourced ambient radiation. The branch-major
IQN-ILS owner retains the fixed flow, inlet and transport properties. Every
solid, radiation and air trial receives the same accepted old specific
enthalpy; only an accepted time step advances storage.

The last counted solid response retains its applied Robin references. A
separate shared endpoint evaluation uses the air march's updated references
and actual nonlinear radiation. Acceptance requires the complete nodal
residual, physical energy residual and a separate actual-air energy balance.
Air enthalpy gain, off-path convection and radiation each enter storage once.
The standard branch decomposition describes the applied convection; separate
physical reference and flux rows describe the endpoint audit. No additional
solid solve is hidden after acceptance.

The residual target is frozen from the initial enthalpy under explicit
branch-inlet seed references, labeled
`initial-enthalpy-with-branch-inlet-seed-boundary`. This seed is a boundary
guess, not a claim that the old air state is converged. Each air callback
audits that old state at its applied references and tightens both positive
Newton tolerances. Air-reference convergence reserves an independent share
of the nodal residual and joule tolerances. These gates do not replace the
existing independent branch heat checks.

Per-step and aggregate `transient.conjugate` evidence records air iterations,
all solid responses, Newton/Krylov work, radiation trials, direct endpoint
audits, physical energy and actual-air energy. The 100-air-iteration policy
multiplies the radiation and Newton caps in admission; direct audits are
bounded by 102 per step without natural convection. Memory admission includes
air staging, retained per-segment receipts and final inner Newton histories
per accepted endpoint; discarded coupling trials reuse bounded workspace.
Enclosure radiation still refuses. Driver 70 adds enthalpy trajectory
adjoints with the complete natural and air boundary feedback described below.
Driver version 65 separates this physical coupling from older retained solves. `examples/heatsink-fan/heatsink-fan-enthalpy.fsim` is the
native fan-cooled synthetic phase-storage example.

Driver 66 closes separate natural reservoirs inside each native enthalpy/air
callback, including ambient radiation, regional initial states and phase-aware
conductivity. All trials use the same old specific enthalpy. The final audit
uses both actual natural coefficients and updated air references at the last
counted field; the original residual and energy gates remain in force.
Natural-reservoir and radiation heat remain outside branch air transport.
Applied boundary rows stay attached to the returned field and are distinct
from the independently evaluated physical boundary.

Natural work and endpoint audits include all air trials under checked products
of the 100 air, 80 natural and declared radiation/Newton/Krylov limits. Per-step
natural/radiation rows describe the final air trial; trajectory counters include
all trials. The mixed audit limit is 16102 per step, including the seed, every
old-state audit, all natural audits and the final physical audit. The synthetic
`data/reference-project/cooling-enthalpy-airflow.fsim` demonstrates a warm
natural reservoir, a cold radiation reservoir and fan-driven solidification.
The mixed G1/G3 regressions in `tests/solve/transient_enthalpy_conjugate.rs`
check signed heat, integrated latent storage and phase fraction against
independent balances and exercise that input through the native pipeline.
Air storage, fan dynamics and moving phases are not inferred. Native enthalpy
trajectory adjoints are admitted under the scope below.

Driver version 67 admits the existing `temperature-final-adjoint` and
`temperature-sampled-peak-adjoint` report outputs also support native
enthalpy trajectories with prescribed exterior laws and fixed matching
contact. They use the temperature-maximum requirement region and the
completed fine grid; the peak selector includes the initial state and
retains its earliest maximizing branch. Tied maxima are marked explicitly.

Only requested adjoints retain full accepted fine endpoints. Storage
admission includes those checkpoints, source bases, reverse staging and
derivative rows. The heterogeneous shared `linearize_accepted` owner
recomputes physical residual, energy, chart and material smoothness before
each reverse step. Revalidated h, temperature and phase must match their
retained fields. No nonlinear endpoint is re-solved or replaced, and no
primal field or energy total comes from the reverse calculation.

The `fs-cli-native-enthalpy-adjoint-v1` object appears under
`transient.adjoint`. Its controls are initial nodal and regional specific
enthalpy, absolute regional delivered watts per fine interval, interval
source multipliers, and regional reference density plus its uniform
multiplier. The shared source projector is the same one used for forward
power, including zero-baseline regions. Density differentiation changes
fixed reference storage mass while holding the phase chart, conductivity,
contact and geometry fixed; it models no density-driven expansion.
Initial enthalpy gradients are not labeled initial-temperature gradients.

`transient.adjoint.prescribed_temperature` reports independent uniform bath
temperature controls held through the selected trajectory. Each declared
boundary retains its target, nominal kelvin value and
`dtemperature_dprescribed_temperature_k_per_k`. The shared transpose retains
the fixed-node objective and free-node coupling; source, history, density,
conductivity, contact and convection contractions exclude fixed residual
rows while retaining source basis columns that load neighboring free nodes.
Each bath trace vertex contributes once. Declarations sharing prescribed
vertices report an explicit unsupported reason without a numeric derivative.
An initial sampled peak has zero future bath influence. Chart endpoints and
slope-changing branches continue to refuse classical derivatives.
The native prescribed-bath tests compare fully fixed storage and radiation
reactions, a one-free-node P1 recurrence and freshly perturbed trajectories.

The reverse carry remains in h coordinates through latent plateaus;
temperature seeds are pulled back only at the selected objective state.
A selected initial peak uses a checked local chart slope and has zero
future source or density influence. The shared transpose checks its true
residual within the admitted per-solve linear work cap. Reports retain
revalidated endpoint counts, total dual iterations, worst residual and
bounded reverse work. Chart/material derivative corners refuse rather than
inventing a classical gradient. Those history/source/density coordinates hold
prescribed boundary parameters, constitutive laws, phase multipliers, contact,
geometry and the time grid fixed. Additional material and contact coordinates
are supplied by driver 74 under the scope below.
Driver 70 includes full natural-convection and airflow feedback in this
native adjoint path. Driver 68 includes complete ambient-radiation
state feedback. Driver 69 also returns per-surface emissivity and reservoir
temperature controls, with query temperature and source support held fixed. This is a local Estimated
derivative, without a gradient-error, continuous-peak or validation claim.
`data/reference-project/cooling-enthalpy-adjoint.fsim` supplies a synthetic
latent-to-sensible heating example with independent analytic expectations.

For ambient radiation, the heterogeneous shared
`linearize_accepted_storage_with_ambient_radiation` owner revalidates the
original nonlinear residual and energy before retaining the complete
area-mean radiation tangent. The original convection boundary is supplied
once; no converged radiation coefficient is frozen or reapplied. History,
source and regional reference-density controls use the same physical
transpose. Admission additionally bounds `2 * vertices * patches` feedback
entries, at most 64 patches, and one active radiation workspace. The optional
`transient.adjoint.radiation` row records this method, feedback cap and
additional workspace. Capacity-storage radiation adjoints are admitted by driver 75 below.
`tests/solve/transient_enthalpy_radiation_adjoint.rs` covers hot and cold
reservoirs against independent scalar roots and fresh perturbed native
solves, a frozen-secant negative control, latent-to-sensible history and an
initial sampled peak. `cooling-enthalpy-radiation-adjoint.fsim` supplies the
native radiative sensitivity example.

Driver 75 admits ambient-radiation history adjoints for native heat-capacity
storage. Each forward radiation endpoint retains its original physical
residual target alongside the accepted temperature. The reverse pass audits
the actual gray-radiation residual and energy, then binds the same unchanged
field to `C/dt + J(T)` without another primal solve. Complete mean-dependent
coefficient and weighted-reference factors enter the true transpose; the
consistent Robin face trace is preserved. All initial-temperature, regional
absolute-power, common-capacity, original-convection, emissivity and reservoir
controls use the same complete multiplier in watt residual coordinates.
`revalidated_solid_endpoints` counts consumed physical audits; reconstructed
endpoints and reconstruction Krylov iterations are zero. Surface source and
parameter identities use the same `radiation.surfaces` fields described below.
Admission additionally bounds four nodal factors per Robin port, one active
endpoint workspace, aggregate control rows and a scalar target per fine step.
The existing `tests/solve/transient_radiation.rs` checks all seven controls
against independent scalar implicit derivatives and fresh native perturbations
for both radiation heat signs, exact forward-field preservation and zero
future influence at an initial sampled peak. The synthetic
`cooling-radiative-adjoint.fsim` example requests the new native report.

Driver 76 extends native heat-capacity adjoints to the in-domain Churchill-Chu
vertical-plate natural law, alone or with ambient radiation. The actual
coefficient is rebound at the saved Robin wall mean, including signed
wall/ambient and film-density dependence. Combined radiation receives both
coefficient slopes and both weighted-reference terms. Reverse admission uses
the outer natural endpoint's original physical target, never a nested frozen
radiation target, and the same physical energy tolerance. The complete
watt-coordinate multiplier feeds every history, power, capacity and cooling
control without another timestep factor or primal solve.
`transient.adjoint.convection.natural_reservoirs` reports independent target
and ambient identities with `dtemperature_dambient_temperature_k_per_k`;
the parameter changes both the explicit reference and natural coefficient.
Natural ports are omitted from ordinary fixed-Robin control rows. Correlation
data, pressure and geometry remain fixed, and unavailable derivatives at
equilibrium or domain endpoints retain explicit refusal. An initial sampled
peak keeps zero future controls. Storage admission charges one bounded active
natural-law workspace before either grid executes. Existing native natural
tests compare hot/cold histories, with and without radiation, to independent
implicit recurrences and fresh physical-control perturbations. Their unequal
natural/radiative reservoirs exercise both weighted-reference contributions.
`cooling-natural-radiative-adjoint.fsim` supplies a synthetic native example.

Native heat-capacity adjoints also admit quasi-steady airflow, including
multiple segments per branch and composition with distinct natural traces and
shared ambient-radiation traces. The accepted endpoint binds the complete
`C/dt + J(T)` transpose with every upstream air-reference derivative and the
weighted-reference derivative on a radiating air port. The same watt-coordinate
multiplier drives history, regional power, capacity, natural/radiative
parameters, branch inlet temperatures and independent absolute fan-bank speeds;
no extra timestep factor or primal solve enters these contractions.

Air controls use the existing `convection.air_inlets` and `fan_banks` identities
and derivative fields described below. Inlet edits hold the hydraulic density
envelope, pressure and transport properties fixed. Available fan coordinates
include hydraulic flow, heat-capacity rate, smooth in-domain Reynolds/card
dependence and downstream air response. An unavailable bank retains its
identity and explicit unsupported reason throughout the selected prefix;
neither partial derivatives nor substitute zeros are reported. Derived air
ports are omitted from fixed-Robin control rows. Initial sampled peaks have
zero future controls and require no transpose.

Every consumed air endpoint rechecks its actual marched references, complete
physical residual and separate actual-air energy balance against the original
accepted target and joule tolerance. The outer natural target remains
authoritative when it encloses air/radiation. Revalidation preserves the
counted field and applied-operator evidence without replacing either.
Admission bounds 64 Robin ports and 64 fan banks, four nodal feedback factors
per port, the bounded reference/control workspace and the fine target tape
before either integration grid runs. No dense nodal Jacobian is allocated.
Existing `tests/solve/transient_conjugate.rs` adds independent uniform
backward-Euler/exponential-air and Gnielinski/absolute-speed derivatives, fresh
native physical-control probes, unchanged-field and initial-peak checks, and
a nonuniform two-segment air/natural/radiation control direction. Radiation on
the downstream air wall exercises the weighted off-diagonal reference terms.

Driver 69 retains the radiative parameter owner alongside its bound reference
storage. The shared combined pullback solves one complete physical transpose
and contracts h, source, density, emissivity and reservoir controls with that
same multiplier. The native report sums the two surface controls over the
selected prefix; future steps and an initial sampled peak contribute zero.
`transient.adjoint.radiation.surfaces` preserves declaration order and gives
each surface's name, target, nominal emissivity and reservoir temperature,
card and property-usage receipt identities, fixed query temperature,
`dtemperature_demissivity_k` and
`dtemperature_dreservoir_temperature_k_per_k`. Surfaces sharing a card retain
independent coefficients. These controls perturb the absolute model inputs
throughout the trajectory; they do not differentiate the immutable source
record or propagate its uncertainty. Emissivity remains in (0,1], so a nominal
value of one reports only nonpositive feasible increments. Reservoirs require
positive kelvin and may lie outside the wall-property support interval.
The time grid remains fixed; no horizon or time-step control is reported.
Memory admission adds two aggregate arrays, their headers and bounded surface
rows, independent of fine-step count. The existing radiative trajectory tests
compare both controls to independent analytic derivatives and fresh native
solves with new immutable emissivity packs, including latent-to-sensible
history and an initial sampled peak.

Driver 70 completes the native total-enthalpy reverse path for in-domain
Churchill–Chu natural convection and independently supplied quasi-steady air
branches, including their composition with ambient radiation and phase-aware
conductivity. Each consumed checkpoint retains the original forward physical
residual target. The reverse reconstructs actual card coefficients and air
references at the unchanged accepted temperatures, checks the full physical
residual against that retained target, and independently checks both physical
energy and actual-air enthalpy gain with natural/radiative heat kept separate.
No nonlinear solve, history update or convergence-target relaxation occurs.

The shared typed Robin attachment adds the signed natural wall derivative,
including film-density dependence, and the complete branch-major upstream air
reference matrix to the already complete enthalpy/radiation tangent. The same
multiplier contracts initial h, regional watts, reference density and existing
radiative parameters. Capacity-storage natural and airflow adjoints are
admitted under the complete endpoint transpose described above; other natural
cards require their own differentiated constitutive law.
Driver 72 adds physical reservoir, branch-inlet and independent fan-speed
controls under the scope below. Geometry and transport properties are fixed.
`transient.adjoint.convection` names the method, port counts, original-target
rechecks and bounded workspace. Additional low-rank factors are capped by
`4 * vertices * convection_ports`, plus existing radiative factors; at most
64 physical convection ports are admitted. Tests in
`tests/solve/transient_enthalpy_convection_adjoint.rs` compare heated/cooled
natural histories with independent scalar derivatives, forced-air controls
with fresh native finite differences, and mixed natural/two-segment-air/
radiation histories with a coupled finite difference.

Driver 72 accumulates constant-through-trajectory convection parameters with
that same completed joule-residual multiplier, scaled by each endpoint's
duration exactly once. `transient.adjoint.convection.natural_reservoirs`
retains each target, nominal ambient kelvin and
`dtemperature_dambient_temperature_k_per_k`; it includes the explicit
reference and full signed ambient/film-density coefficient partial.
`air_inlets` retains each independently supplied branch, nominal inlet kelvin
and `dtemperature_dinlet_temperature_k_per_k`. This inlet control holds the
hydraulic density envelope, pressure and transport properties fixed and changes
all inlet declarations on that branch consistently.

`fan_banks` retains declaration order, bank identity and nominal speed ratio.
Each available row reports `dtemperature_dspeed_ratio_k` for one absolute bank
speed, including its retained hydraulic flow response, capacity rate,
in-domain Reynolds coefficient derivative and downstream air response.
Series/parallel bank coordinates are independent; other bank speeds, passive
loss geometry and fan-curve data remain fixed. Unsupported coordinates retain
`status: "unsupported"` and their reason, with no partial or substitute zero.
If any consumed endpoint lacks a bank differential, its trajectory row stays
unsupported. An initial selected peak has zero boundary influence.

`control_steps` covers exactly the selected reverse prefix. At most 64 fan
banks are admitted; aggregate arrays and bounded curve scratch are included in
`control_workspace_bytes` before integration, independently of history length.
The coupled enthalpy regressions also compare ambient-only and inlet-only
perturbations, an absolute fan-speed perturbation at a nonunit nominal ratio,
latent-history inlet propagation and a joint mixed natural/inlet/fan/radiation
direction. These tests preserve the physical density envelope during inlet
edits and check explicit unsupported fan-domain endpoint reasons.

Driver 74 also returns `transient.adjoint.material_controls` from the same
complete endpoint multipliers. `conductivity_regions` scales each declared
volume region's whole effective conductivity tensor/curve at one, including
its retained phase multiplier, and reports
`dtemperature_dconductivity_multiplier_k`. The shared enthalpy owner retains
the four local conductive-load entries at the accepted temperature and mean
nodal liquid mass fraction; the existing full state tangent already owns
temperature/phase feedback. Labels aggregate element partials without volume
normalization or merging regions that share a card.

`contact_interfaces` reports `dtemperature_dresistance_multiplier_k` for
one scale on all retained face resistances of each named matching interface.
It preserves resistance maps, independent trace temperatures, consistent P1
mass integration, original card identity and face-pair count. Contact geometry,
pairing and constitutive selection remain fixed. Both coordinate families
hold their positive scales through the selected prefix, with reference value
one (logarithmic reference zero) and kelvin derivative units. Initial selected
peaks have zero material/contact influence. These are model sensitivities,
not absolute property derivatives, card mutations or material uncertainty.
One bounded active element/contact workspace and aggregate report are charged
before integration; no per-step parameter tape or extra transpose is added.

Schema v14 optionally declares a conductivity multiplier for each regional
enthalpy material. `linear-liquid-mass-fraction` interpolates between positive
solid and liquid multipliers at the tetrahedron's mean nodal liquid mass
fraction, then scales its existing sourced tensor at mean nodal temperature.
The common enthalpy owner supplies the full phase derivative in the Newton
tangent, including latent-plateau feedback. Regions without a declaration keep
their base conductivity. The law is explicitly phenomenological: its source
text and values carry Declared authority, without inferring a volume-fraction
mixture or authenticating the supplied property data. Base conductivity support
and original material-query receipts remain in force.
Each material row in the transient receipt retains its optional law, endpoint
multipliers, coordinate and source. Schema v13 migration preserves existing
storage and base conductivity, adding no phase effect. The synthetic
`cooling-enthalpy-phase-conductivity.fsim` example combines the law with
latent storage, convection, radiation and a startup pulse. Driver version 61
separates this transport model from older retained solves.

The enthalpy model also composes with declared ambient gray radiation on its
prescribed Robin patches. The shared radiative enthalpy owner holds physical
old enthalpy fixed through every trial, then reassembles the actual endpoint
residual and checks storage against nonlinear radiative heat. The original
energy tolerance and complete initial-residual target remain the acceptance
gates. Each time row and trajectory energy use source minus Neumann heat,
original convection and signed nonlinear radiation; a hotter reservoir may
supply heat. The ordinary final Robin operator and radiation receipt retain
the accepted applied boundary and its separate convection/radiation split.

Per-step radiation receipts retain physical residual/threshold/energy, applied
and nonlinear heat, source patches and all trial work. The shared owner exposes
Newton histories only for the final inner solve of each accepted time step;
these genuine rows keep that explicit scope. Separate `transient.radiation`
counters include every trial's Newton and Krylov work, and the admitted
step/radiation/Newton/Krylov/backtracking products bound the full trajectory.
Driver version 58 separates this composition from older retained solves.
`data/reference-project/cooling-enthalpy-radiation.fsim` is the synthetic
latent-storage pulse example. The model remains radiation from patch-mean
temperature to fixed reservoirs, without enclosure exchange or pointwise
integration of T(x)^4.

Declared natural convection also composes with spatial enthalpy and optional
ambient radiation. The existing sourced correlation evaluates every candidate
wall state; its original equilibrium and Rayleigh-domain refusals remain.
Every natural, radiation and Newton trial retains the same physical old
specific enthalpy. The new shared endpoint evaluator checks the complete
candidate under updated physical convection and actual nonlinear radiation,
without imposing an adjoint smoothness condition at a phase-chart corner.
Publication requires coefficient convergence, a fixed nodal residual target,
and the original absolute joule-energy gate. The target is fixed from the first
physically evaluated candidate, so a powered ambient start does not require
an undefined equilibrium coefficient. Later inner tolerances tighten.

Both grids count all natural iterations and nested radiation/Newton/Krylov
work; direct endpoint evaluations have a separate bounded allowance. Per-step
and integrated energy use actual physical convection and signed radiation.
The retained final field remains the last counted solid response; its applied
operator is kept distinct from the independently checked physical law. The
normal `natural_convection` receipt retains the actual final card evaluation.
`tests/solve/transient_enthalpy_natural.rs` supplies independent sensible and
latent energy references, both heat directions and refusal-without-field cases.
`data/reference-project/cooling-enthalpy-natural.fsim` demonstrates passive
losses during a synthetic melting run. Driver version 60 separates these
semantics from older retained solves.

The same native declaration also admits ambient gray radiation on its
prescribed Robin patches, using the existing sourced emissivity cards and
`conduction.radiation` controls. Every coupling trial holds the accepted
physical old temperature fixed. Acceptance checks the actual endpoint
area-mean radiation law, full discrete implicit residual and physical
storage-minus-input energy balance; the fixed-node reaction includes its
storage jump. Temperature-dependent conductivity remains evaluated at the
endpoint with its full tangent. Material-domain escape, exhausted coupling
work or cancellation publishes no partial trajectory. Wall-time checks are
cooperative between numerical operations; no intra-kernel wall-time guarantee
is claimed.

Each step's `radiation` object distinguishes frozen-secant applied radiation
from actual nonlinear radiation, and reports the physical residual norm,
threshold, energy residual and prescribed-temperature reaction. Existing step
and common energy fields retain their applied-operator balance. The top-level
conduction `radiation` object retains the complete final fine endpoint and
card/source evidence; `transient.radiation` aggregates both grids' work and
maximum physical energy residual. Linear work caps each solid response
(including all Newton corrections); the radiation cap bounds repeated
responses per outer airflow or natural-convection evaluation (or endpoint
without either). The checked product of outer convection, radiation and
linear caps is the disclosed full endpoint Krylov allowance. Per-step
nonlinear details describe the final inner solve, while radiation and
trajectory counters include all trials. Driver semantics version 55 includes
native natural-convection transients and keeps this behavior separate from
older retained solves.

`data/reference-project/cooling-radiative-pulse.fsim` combines a declared
startup pulse with ambient radiation using the native validate/import/solve
workflow. The model remains area-mean patch radiation to a fixed reservoir;
enclosure exchange, latent storage, coupled air storage, transient adjoints
and continuous-time peak guarantees are not inferred.

Native `airflow-convection` declarations also compose with this solid
trajectory. Existing flow-network operating points, wetted-area lowering,
convection cards and branch-major IQN-ILS exchange are reused. Every air trial
holds the accepted old solid field fixed, including all inner radiation and
conductivity corrections. The endpoint is the last counted solid response
at its retained applied references; the updated air reference does not
trigger another solve. Each branch's own heat gate remains independent.
The raw reference tolerance is additionally bounded by one percent of the
declared energy tolerance divided by `dt * sum(hA)`. An independent endpoint
joule gate replaces on-path Robin heat with actual air enthalpy gain and
retains off-path convection, nonlinear radiation, physical prescribed
reaction and storage. Radiation is excluded from the air budget. The physical
free-residual target is fixed from the initial predictor evaluated at actual
marched air references, including the new prescribed-temperature endpoint lift.
Inner solid tolerances and reference convergence reserve room under that fixed
target. Acceptance audits the counted candidate under its actual air/radiation
law without another primal solve; convergence of the applied-reference exchange
alone cannot admit the endpoint.

Per-step `conjugate` evidence retains the applied references, exchange receipt,
physical residual and its original target, actual-air coupled energy residual
and every nested numerical work count. Each air solve adds exactly two direct
endpoint audits; intermediate reference-load changes use a conservative norm
bound. Nested natural iterations accumulate the air audit count separately
from natural-owned endpoint evaluations.
`transient.conjugate` retains trajectory totals and the full product of the
airflow, radiation and solid-response work caps; `conduction.conjugate`
retains the final fine endpoint. Airflow is quasi-steady at the fixed retained
operating point, with no air storage or fan dynamics. The complete fan-cooled
pulse input is `examples/heatsink-fan/heatsink-fan-pulse.fsim`; its sourced
capacity and fan curve are explicitly synthetic.

Native `natural-convection` also composes with heat storage, endpoint k(T),
matching contact, scheduled power and ambient radiation. It reuses the
existing Churchill–Chu card and signed area-mean wall law. The initial
coefficient guess permits powered startup from ambient; each solved
candidate must satisfy the actual law and its Rayleigh domain. Exact
equilibrium and out-of-domain Rayleigh numbers retain the owner's named
refusals. Driver 64 also admits distinct natural and airflow-convection
boundaries in the same transient, with each natural update enclosing a
complete quasi-steady air exchange.

Every natural trial holds the accepted old field fixed. Acceptance requires
relative h convergence at 1e-10, the full physical free residual and the
original physical joule energy gate. The free-residual target is frozen
after the first physical candidate; later inner solves reserve room for the
outer law discrepancy. The natural loop has at most 80 evaluations. The
checked work admission includes all nested radiation/Newton/Krylov responses
and at most 160 direct endpoint audits. A successful endpoint performs
2*iterations-1 direct audits and retains its last counted solid response;
no post-convergence re-solve is hidden.

Per-step `natural` evidence separates applied and actual coefficients,
convective fluxes, physical residual/threshold, prescribed storage reaction,
radiation, energy and work. The common final energy retains the applied
operator; `conduction.natural` retains the final physical card receipt.
`transient.natural` retains both grids' work, maximum accepted physical
energy residual and fine integrated heat balance.
`data/reference-project/cooling-passive-pulse.fsim` is the complete synthetic
passive pulse input. Its final-time temperature scope remains explicit.

For mixed transient cooling, the outer endpoint evaluator reassembles the
full nodal residual under both the actual natural coefficients and the actual
air-law reference temperatures. It also gates storage against actual branch
air enthalpy, off-path reservoir heat, radiation and the prescribed-node storage
reaction. Each natural trial retains the same accepted old temperature;
coefficient or scalar heat convergence alone cannot admit a field. After the
first candidate, inner air and solid tolerances reserve room under the fixed
physical residual target. The checked work products include all nested natural,
air, radiation and nonlinear responses. The final air row describes its last
natural trial, while trajectory counters include every trial without counting
solid work twice. No air storage or fan dynamics is inferred. The explicit
transient-adjoint report uses the accepted physical target and complete mixed
boundary transpose described above.
Three G1/G3 native regressions in `tests/solve/transient_mixed.rs` compare a
nonuniform P1 endpoint with an independent implicit scalar reference, covering
both natural heat directions, serial air heating, prescribed startup storage,
k(T), radiation, deterministic replay and refusals without partial fields.

Optional `power-schedules` supply absolute delivered watts for named volume
regions, replacing their static `power` rows without applying duty factors
again. Every schedule covers the full horizon with explicit switch times.
The coarse grid lands on the union of all switches; the fine grid bisects
each coarse interval. Both count against the same total step allowance.
Short pulses are therefore integrated even when both former uniform grids
would miss them. The workload receipt retains declarations and integrated
scheduled input, while each numerical step records its actual source watts.
Unscheduled regional and surface heat inputs stay static. Final-time QoIs
and the separate temporal estimate keep their existing scope; no peak between
time nodes or probabilistic workload model is inferred.

Runtime dependencies remain Franken-only.

## Conformance tests

`tests/cli.rs` covers the grammar and all v0 verbs, the FrankenScript
executor (the tracked Journey A program reproduces the `.fsim` verbs' run id,
report content hash and package root exactly; unit-error, unbound-verb, seed,
hash-pin and capability falsifiers refuse without opening the ledger; the
JSON twin admits to the same identity; a `study.run` program on the reduced
marquee fixture matches the `study` verb's receipt modulo wall time), stable exit classes,
strict validation success, structural findings with fixes, noncanonical input
refusal, JSON escaping/line discipline, import-policy conflict/numeric
refusals, routing of both admitted import policy shapes into bounded project
I/O, the exact producer-Bead refusal for report/package, and the solve
grammar's ledger-operand requirement.

`tests/import.rs` supplies a closed reference tetrahedron and covers G0 retained
lineage, repair of deterministic duplicate/degenerate STL facets, strict
faceted-STEP decoding through topology/SDF handoff, separate nested receipts,
and exact repaired-mesh retention. G3 covers changed source identity, open-mesh
promotion refusal, mis-scaled unit refusal, dangling assignment refusal, and a
clean/dirty re-tessellation pair with identical selector statistics. It also
drills the fail-closed named-group behavior when repair removes faces. G4
covers pre-cancellation with zero publication. G5 covers content-identity
equivalence across independent ledgers and proves that changing exact STEP
sampling bits moves the frozen operation IR and retained summary. Every
recorded case runs the ledger linter.

`tests/solve.rs` covers G0 run-identity determinism and input sensitivity,
the pinned stage order and gap owners, prefix execution with a ledgered
`cli-solve-stage-gap` refusal, and the honest budget-exceeded partial with
warning fractions and refused re-resume. G3 covers missing/mismatched import
evidence, an exact-schema import-summary decoy, a caller-sealed but
non-authoritative state, a canonical higher checkpoint with one substituted
predecessor, `{}` replacements for otherwise exact import limits/policy
objects, isolated per-source and aggregate byte-cap forgeries, a genuine import
outside the project-memory solve envelope, a tiny PLY header with an unknown
zero-property element and a huge declared count, unrelated successful import
and same-run operations wider than the complete edge scan, competing valid
longest checkpoints, unknown run ids, and a corrupted retained project pin
(refused by the ledger's own read-integrity gate before the driver's identity
check — both fail closed). Identity refusals occur before progress or ledger
publication. G4 covers pre-cancellation with zero
publication, every fresh solve-owned evidence phase, a storage-row-spanning
raw artifact, an exact-4-MiB inline promotion receipt stopped after its first
shared 64-KiB controlled tile, extreme finite PLY writer spellings, incremental
UTF-8 and JSON-parser phases, named-group face-range checks, chunked retained
project comparison, resume checkpoint/project/receipt re-attestation, and the
final pre-publication boundary with exact zero publication beyond any already
durable prefix. It also proves a cancelled resume retry reproduces the same stage-receipt
identities as an uninterrupted run. G5 covers identical stage-receipt
identities as an uninterrupted run. G0 additionally lowers the committed
two-region imported fixture through its declared finite-resistance contact,
checks nonzero interface heat transfer and energy closure, and replays the same
receipts on resume. G5 covers identical stage-receipt identities across
independent fresh ledgers. The typed G4 matrices enumerate
61 fresh and 79 resume checkpoints, including entry/intermediate/completion
points around visible-ID pages, controlled candidate fields and text
conversion, operation sidecars, artifact descriptors, edge seals/pages/set
comparisons, the invocation-cached canonical project JSON and Five Explicits,
project validation/identity, entity resolution, source-indexed and stage
assignment derivation, canonical receipt construction, retained-project
parsing, state decoding, both available reconstructed stage receipts, and
order-preserving duplicate work.
Focused unit regressions cover immediate JSON phase entry, final-current
syntax priority, empty/missing evidence dispatch entry, mismatch and invalid
UTF-8 result priority, malformed PLY completion polling, named-face range
priority, exact-4-KiB/cap-plus-one labels, and the exact-1-GiB cumulative-work
boundary whose refused plus-one charge does not advance the meter. A focused
accounting regression proves accepted UTF-8 copies and direct byte comparisons
charge each processed tile once while invalid UTF-8 is not charged as accepted
copy work. Integration regressions cover multi-page unrelated history, the
8,192/8,193 visible-ID boundary, and cancellation/retry through exact-1-MiB
operation fields. A G3
import-IR regression keeps
an intentionally unsorted named-face list unchanged while refusing repeated
faces, repeated group names, and an over-ceiling label with zero stage
publication.

## No-claim boundaries

- FrankenScript executor v0 binds only the cooling project pipeline and the
  canonical study driver (`EXECUTABLE_VERBS`); the fs-ir catalog's physics
  operators (`flux.*`, `ascent.*`, …) admit but do not execute. Each stage
  driver re-reads its file exactly as the CLI verbs do; `:hash` pins are
  checked once, before execution.
- `validate` proves only canonical structural and dimensional admissibility.
  It does not prove referenced artifacts or material cards exist, a requested
  capability is installed, the project is solvable, or any physical model is
  valid.
- Geometry import binds exact raw bytes, lower-layer receipts, one promoted
  finite tessellation, assignment reports, and their lineage. The legacy FNV
  hook and caller path/label do not authenticate custody, physical/CAD
  sameness, continuum coverage, units, or topology beyond the retained
  lower-layer claims. The summary's `assignment_table` is strictly encoded and
  bounded as part of the retained bytes but is diagnostic rendering, not
  resume authority; typed assignment-report artifacts carry the retained
  assignment evidence.
- Faceted STEP support is limited to fs-io's pinned triangular root-reachable
  resource subset and estimated SDF handoff. It is not full EXPRESS/AP
  interpretation, representation/unit-context discovery, NURBS/surface
  tessellation, component nesting, self-intersection certification, or
  physical/CAD sameness. Named face groups are caller-supplied labels on the
  promoted soup, not independently certified CAD product-structure identity.
- `report`, `package`, and `compare` execute against completed runs and prove
  them by sealed evidence (see the verb section); none of them replays
  physics, and `compare` writes nothing.
- All seven stages execute; a project without an explicit conduction setup or
  without a `temperature-max` requirement refuses at the stage that needs it
  (`cli-solve-conduction-undeclared`, `cli-solve-qoi-undeclared`). The
  conduction path proves only the declared finite-mesh, static
  Dirichlet/Neumann/Robin solve, exact matching-P1 finite contact, and the
  reported algebraic residual and energy balance. It does not authenticate
  source geometry or material claims, establish mesh convergence beyond what a
  declared uniform h-ladder measured on the QoI's nodal maximum (no
  energy-norm or DWR bound), or lower nonmatching or temperature-dependent
  contact. Constant-contact admission checks source curve knots throughout
  the admitted temperature range; equal endpoints cannot hide interior
  resistance variation. Flat supported curves remain usable.
  Unsupported interface models refuse rather
  than silently assuming perfect contact. When the project declares a
  schema-v4 `(airflow-convection ...)` law, the stage performs ONE branch's
  conjugate airflow exchange (driver v10, receipt v3 `conjugate` object):
  the coefficient is derived from the named `fs-convection` card at the
  flow-network branch Reynolds number with dry-air properties frozen at
  300 K, the reference temperature is the exponential-law marched air
  temperature over the target's retained exterior area, and the solid/air
  fixed point must pass a kelvin convergence test, a scale-free watt
  balance gate, and an `fs-conduction` decomposition cross-check or the
  stage refuses by name (`cli-solve-conduction-airflow-*`). It does not
  close more than one branch, re-evaluate properties or the coefficient at
  film temperature, model recirculation, buoyancy, redistribution, or
  momentum feedback, propagate the flow bracket through the card, or claim
  experimental validation or a maturity level for the exchange.
- `material-resolve` proves that every declared region and interface resolves
  to an admitted card whose selected claim answers the required property at
  every point of the declared temperature range under the owner's continuous
  support and source-selection checks, and it retains that
  claim's replayable usage receipt. It does **not** authenticate the pack
  producer, validate a claim against any external corpus, narrow or replace a
  claim's stated uncertainty, or turn an `Unstated` uncertainty into a bound.
  A pack that decodes is admissible input, not validated data. Perfect
  contact remains refused by the binding layer, so declaring it still fails
  closed. `import-verify` proves
  canonical internal row/receipt/lineage structure, retained content hashes,
  and the documented resource-count consistency only. PLY/report parsing
  checks exact writer grammar, project source/unit/subject ordering, and
  bounded counts; it does not recompute the source FNV/parser admission,
  cross-check or replay promotion/STEP repair-root-spacing policy, authenticate
  lower-layer producers, verify mesh/named-group/request/assignment
  fingerprints, re-run named-group or selector classification, or recompute
  reported areas, volumes, and bounds. It therefore does not prove the
  geometry is watertight, meshable, or physically meaningful. `assign` binds
  declared targets to verified import evidence; it does not re-run selector
  resolution.
- The session opened by the driver uses caller-declared capability evidence
  (`open_session_declared`); no external issuer, grant authentication, or
  revocation policy is claimed. Memory consumption is charged as zero — the
  governor does not sample RSS — so governor/runtime-memory enforcement is
  currently inert while its core-seconds and wall-seconds enforcement is real.
  Separately, the solve import byte-total preflight genuinely refuses evidence
  above the declared memory value used as an admitted-input/work cap; it
  neither measures nor proves peak memory use.
- The governor is in-memory per process; exactly-once submission and durable
  session terminals (`Governor::new_durable`) are not yet wired. Budget
  continuity across resume relies on the recorded consumption in the sealed
  driver state, not on a durable governor.
- Stage-level checkpoints are driver state only. Intra-stage solver
  checkpoints (fs-exec v2 envelopes, fs-ledger solver-checkpoint receipts,
  pause acknowledgement) arrive with the physics stage producers. Migration
  of the driver state to the v2 envelope waits on fs-exec's post-restart
  authorized-resume seam, which its CONTRACT names as unfinished.
- The solve evidence cadence covers controlled visible-ID page progress,
  controlled operation/artifact tiles, guarded text conversion, prehashed
  sidecar checks, bounded descriptors and edge seals/pages, UTF-8 copy, JSON
  cursors, canonical PLY decoding and writer-token comparison, Five-Explicits
  and canonical byte comparison, project/source identity and validation,
  entity resolution, assignment filtering/derivation, named-face ranges,
  duplicate work, and exact edge-set comparison. SQL may materialize one
  guarded operation field of at most 1 MiB before its at-most-64-KiB callback
  slices, one bounded descriptor, or the bounded rows selected for a 64-ID
  page/1,024-edge read; solve claims cancellation and accounting at the API
  entry/callback/result boundaries, not inside the storage engine.
  JSON/PLY cursor advancement is checkpointed, but its already-admitted input
  bytes are not charged again per parser instruction. Standard-library
  allocation, hash-table insertion, guarded text conversion, `parse_sexpr`
  over at most 16 MiB, `findings`, canonical project/Five-Explicits rendering,
  `resolve_entities`, and each bounded receipt-row formatting call remain
  opaque between their pre/post checkpoints; derived render/receipt bytes and
  accepted UTF-8/direct-comparison passes receive the conservative charges
  described above, but no intra-call wall-time latency or peak-allocation
  claim is made. Canonical project JSON and Five Explicits render once per
  fresh invocation; resume caches them by fully attested retained-project
  source hash and fully attests any distinct source before using its own
  result. Lower-layer work outside these calls is not owned by this driver.
- The order-preserving duplicate sets use a fixed standard-library hasher.
  Their membership result is deterministic and their iteration order never
  enters retained output; the implementation does not request a randomized
  set seed.
- Resume attestation proves consistency with this driver's retained canonical
  row/receipt/lineage format. Operation-content sidecars and artifact-edge
  seals detect divergence and post-seal mutation; they are not signatures,
  external producer authentication, or authorization against an actor that can
  manufacture an entirely canonical ledger history through the public ledger
  API. That stronger authority model requires a separately authenticated
  producer capability or signed receipt scheme.
- Deterministic replay is claimed for retained stage artifacts, not for
  enforcement outcomes; kill -9 crash recovery of the ledger is fs-ledger's
  proven claim, and the driver's own kill -9 drill (resume to identical
  content root after a mid-transaction kill) is owed with the first
  physics-bearing stage.
