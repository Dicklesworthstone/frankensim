# fs-orbit — Periodic-orbit service (HB + shooting + continuation)

The first implementation of the fs-vmanifest I09 periodic-machine
orbit slot (music bead `frankensim-music-v8-root-3ez8g.11.1`). Music
is the first consumer, not the owner.

## Purpose and layer

Layer L2. General periodic-orbit machinery over the ISLAND-PLUS-
LINEAR-PORT problem shape: a nonlinear part evaluated pointwise in
time (reed law, device card, cubic spring) and a linear part
specified per harmonic in the frequency domain (`s I` for an ODE, a
TMM impedance for an acoustic load). Serves music lock questions,
flutter boundaries, and thermoacoustics through one problem-shaped
API — never instrument-shaped.

## Public types and semantics

- `OrbitProblem` — `dim`, `island(t, x, out)` (the AFT-sampled
  nonlinearity), `port(s)` (row-major `d x d` complex operator,
  default `s I`), `autonomous()`, and optional `time_domain_rhs()`.
  The latter supplies a `TimeDomainRhs<Self>` callback for the explicit
  derivative `x' = f(t, x)` in the same state and time coordinates.
  It defaults to unavailable; a plain ODE explicitly supplies
  `Some(Self::island)`, while a mass-weighted problem supplies `M^-1 g`.
- `IslandStateJacobian` / `OrbitProblem::island_state_jacobian()` optionally
  declares the complete exact state derivative of the original nonlinear island.
  Availability is separate from evaluation: a declared callback must fill every
  row-major entry and return true; refusal or nonfinite output is `TangentRefused`.
  Absent callbacks retain the original finite-difference coefficient columns.
- `HbAnchor` — `Forced { omega }`, `Autonomous { omega_guess }`
  (phase anchored by `Im X_1[0] = 0`), `Backbone { amplitude,
  omega_guess }` (conservative backbone point closed by an APPENDED
  amplitude-norm equation `|X_1[0]| = a/2`; every balance row stays
  enforced and the phase stays free, so `omega` is
  phase-independent). The authored `a` is a positive finite physical peak
  amplitude, with a representable half-amplitude Fourier coefficient; its
  relative tolerance must be finite and strictly between zero and one.
- `HbBudget` / `solve_hb` / `solve_hb_seeded` -> `HbOrbit` (omega,
  canonical coefficients, per-iteration residual trace, verified residual,
  final verification sample count, measured sampling defect).
- `HbScales` / `solve_hb_scaled_seeded` -> `ScaledHbOrbit` conditions this
  same Newton owner with separate physical state, equation and frequency
  magnitudes. The seed and returned orbit remain in original units. Acceptance
  bounds the absolute normalized Parseval RMS defect and independent sampling
  change; large static preloads never relax it. Both scaled diagnostics are
  retained alongside the original physical residual trace. This is the entry
  for mixed-unit or very small physical responses; the original unscaled
  entry points retain their existing coefficient-relative convergence rule.
- `solve_hb_backbone_coordinate_scaled_seeded` and its adaptive counterpart
  pin the fundamental peak of any chosen original physical state. A static or
  nodal first coordinate no longer excludes another modal family. A full-state
  permutation reuses the existing backbone/Newton and refinement owners; island
  rows, both axes of the exact Jacobian and frequency-dependent port, physical
  state/equation scales, and per-coordinate refinement tolerances all follow
  that same permutation. Returned coefficients retain original ordering. No
  balance row, coupling, loss or material-memory state is removed. The pinned
  fundamental amplitude is distinct from a nonlinear turning amplitude.
  Selecting coordinate zero delegates directly and retains original bits.
  Other coordinates add one prepared state buffer and allocation-free island
  mapping; returned port matrices are permuted in place. Native G1/G3 tests in
  `tests/backbone_coordinates.rs` compare nonlinear frequency to independent
  energy quadrature, preserve dense port coupling and q/v harmonics under phase
  changes, retain a static first state and a last-coordinate velocity anchor,
  refine harmonic order, and propagate declared tangent refusals. These are
  finite harmonic balances, not existence or stability certificates.
- `solve_hb_autonomous_scaled_seeded` uses a complete physical seed and a
  fixed full-waveform phase hyperplane in the normalized Parseval state chart.
  All original coefficients and balance rows remain active; the first state may
  be static or lack a fundamental, and higher harmonics can define the phase.
  Its exact phase row is separate from dynamical deflation and enters only the
  normalized residual; raw diagnostics retain original physical balance units.
  `solve_hb_autonomous_adaptive_scaled_seeded` preserves that initial reference
  across order extensions, retaining independent physical refinement checks.
  Nonautonomous problems and malformed inputs refuse before source callbacks;
  a DC-only reference refuses as `TrivialCollapse`. This local phase chart does
  not prescribe amplitude or establish primitive period, existence or stability.
  G1/G3 checks in `tests/autonomous_harmonic_phase.rs` cover a static first state,
  independent state/equation units, exact/probed tangents, deflation, a reference
  containing only a second harmonic, and phase preservation across refinement.
- `HbRefinementBudget` / `solve_hb_adaptive` /
  `solve_hb_adaptive_seeded` -> `AdaptiveHbOrbit`. These additive entry
  points retain the fixed-order Newton budget, start at its harmonic order,
  and warm-start successively larger orders up to an explicit maximum.
  Per-equation residual tolerances and per-state waveform-change tolerances
  have the corresponding physical units; frequency-change tolerance is relative.
  `HbRefinementCheck` records each solved order, projected band and sample
  count, total/retained/omitted/sampling ratios, and actual inter-order changes.
  `solve_hb_adaptive_scaled_seeded` carries explicit `HbScales` through every
  larger solve without altering the original physical refinement tolerances
  or the independent omitted-band assessment. Warm starts retain all physical
  coefficients and the unknown frequency in its newly extended slot.
  `fs-couple::phs_orbit` exposes both scaled entries on the original complete
  port-Hamiltonian ODE. Its full physical seed includes material memory;
  missing authored forcing bandwidth refuses before source evaluation. The
  original Hamiltonian, resistance, signed ports and source-refusal behavior
  remain authoritative. A native G1 integration regression checks tiny second-
  and third-harmonic displacement/momentum phasors and independent physical
  dissipation/power balance through both fixed and adaptive wrappers.
- `ShootBudget` / `solve_shooting` -> `ShootOrbit` (period, orbit
  point, Floquet multipliers from the forward-differenced monodromy
  through the fs-la complex eigensolver). Admission requires the explicit
  time-domain callback; a frequency-only problem returns
  `TimeDomainUnavailable` before evaluating its island.
- `VariationalShootControls` / `solve_shooting_variational` ->
  `VariationalShootOrbit` adds the exact derivative of the discrete shooting
  map, original physical state scales, full-state autonomous phase closure,
  bounded cancellation and the complete scaled monodromy. It requires a
  separately declared exact `time_domain_state_jacobian`; an island tangent
  does not establish the derivative of a general mass or realized port.
- `ContinuableProblem` / `ContinuationBudget` / `continue_branch` —
  TRUE pseudo-arclength (bordered corrector on coefficients and the
  parameter jointly), so parameter folds are traversed, not stalled
  at. This original entry is forced-anchor only.
- `autonomous::continue_autonomous_branch` adds complete autonomous
  coefficient/frequency/parameter continuation with a full-waveform integral
  phase condition, independent physical scales and explicit work/payload caps.
  It returns `AutonomousBranch`, with exactly two corrected seeds plus the
  declared number of accepted additional points, not a parameter endpoint.

## Forced continuation in physical units

`ContinuationScales` / `continue_branch_scaled_seeded` accept the full original
Fourier seed and separate state, equation and physical parameter magnitudes.
The starting point is solved again. Both this additive entry and the original
`continue_branch` use the same bordered corrector; the original entry retains
its amplitude seed and coefficient-relative acceptance. The scaled predictor,
secant, arclength row, Newton columns and finer-grid check use one normalized
Parseval chart. Its step budgets are dimensionless arclength, and neither DC
preload nor another coordinate's units relax absolute normalized acceptance.
The physical HB trace and sampling defect remain separate from normalized
HB/AFT/arclength diagnostics, including physical `SamplingUnresolved` fields.

Declared island tangents reuse the shared `tangent::Workspace` for coefficient
columns. `ContinuableProblem::forced_parameter_jacobian` optionally supplies
the COMPLETE original residual-parameter column in physical units. Its
`ForcedParameterPoint` retains the original parameter, actual forcing clock,
canonical coefficients and current AFT grid size. The derivative includes
parameter-dependent ports, nonlinear laws, forcing and the caller's
`omega_of(parameter)` sampling chart; a fixed-clock partial is insufficient.
The generic source author is responsible for equivalence to that same family
and clock. The corrector then applies its existing row/parameter units and
exact arclength row. Refusal, incomplete writes or nonfinite entries produce
`ParameterTangentRefused`, without a fallback probe. Without a declaration,
the original finite-difference column remains. An all-analytic corrector also
avoids replaying unmodified original source samples. No extra numerical buffer
is needed; the bounded existing probe workspace is reused. Source parameter
mutations remain owned by the supplied family; this entry does not add a
separate model or claim to restore family state after refusal. Autonomous
continuation does not consume the forced-column declaration.

`fs-couple::phs_orbit::continue_forced_frequency_response_scaled_seeded` uses
the original full pHS state and independently audits signed supply, storage and
resistance loss at every accepted frequency. Quadratic composition and linear
relaxation storage supply cached complete-state tangents through explicitly
exact `quadratic_part` metadata (borrowing exact Q when available); a partial
quadratic with a nonlinear remainder never becomes a complete tangent.
Both original pHS continuation adapters supply the exact `i*n*X_n` frequency
column in their fixed-phase forcing chart: all original state and explicit
material-memory rows remain, and DC has zero clock derivative. This declaration
is limited to their original state-independent `PeriodicPortInput` and
identity frequency map, never assumed for a general source or eliminated
frequency-dependent port. The cancellable branch polls cancellation before
its analytic column and preserves the existing source/error propagation.

Native regressions traverse both Duffing folds under independent coordinate,
equation and parameter unit changes, restart a selected higher-harmonic
waveform, verify tiny-equation AFT alias refusal/refinement, and compare exact
tangent versus residual-probe source calls. Parameter-column regressions also
retain changing material/port coefficients, DC and complex forcing with a
nonidentity clock, and refuse declared incomplete/nonfinite derivatives.
Original pHS tests retain tiny
complex forcing harmonics, Maxwell memory phasors and independently computed
cycle dissipation. These remain fixed-order numerical branches: harmonic
truncation, physical frequency validity, branch identity, precise bifurcation
location and stability need their own convergence/validation checks.

## Exact state tangents in harmonic balance

The fixed-order Newton owner and its adaptive-order callers consume the exact
island Jacobian at each current AFT sample. Compensated Fourier coefficients
of its difference from the reference-sample Jacobian are retained through order
`2N`, together with the actual discrete basis means. Real coefficient columns use
`A_(n-k) + A_(n+k)`; imaginary columns use `i(A_(n-k) - A_(n+k))`, with
conjugation for negative indices and a single real DC coefficient. Assembly
restores the constant Jacobian with those measured means and differentiates the
reference subtraction: for nonzero output order `n`, subtract
`J(0) b_k(0) mean(exp(-i n theta))`. The real DC/AC input basis values at zero
are one/two, and imaginary AC inputs are zero. The DC row restores its reference
once. The original linear port supplies its exact diagonal-in-harmonic blocks.
This retains the sampled original equation's cross-harmonic derivative, up to
floating-point arithmetic; it assumes neither a band-limited material law nor
exact numerical orthogonality of the trigonometric tables.

Fourier accumulation costs O(samples * N * d^2), followed by O(N^2 * d^2)
coefficient-block assembly. Tangent callbacks run once per sample and Newton
iterate, not once per Fourier unknown. Checked tangent buffers are reused
across Newton. Reconstruction, projection, omitted-band checks and the tangent
share one grid of M complex roots of unity. Integer modular indexing evaluates
harmonic n at sample j as root `(n*j) mod M`; the hot traversal uses modular
addition without division or floating oscillator recurrence. Conjugate roots
and cardinal phases are shared exactly, and every trigonometric argument is
bounded by pi. This avoids order-dependent argument reduction and accumulated
phase drift. Finer odd grids retain their actual sample count.

Phase storage is 16*M bytes and preparation uses at most floor(M/2)
deterministic sin/cos calls, independent of harmonic order N. Previously the
residual and analytic tangent together held 16*M*(3*N+1) bytes of phase tables.
Continuation payload admission counts the new shared allocation, including
coexisting coarse/fine grids. Direct projection and tangent summation work
retain the counts above, along with measured discrete basis means and the
constant-preserving projection. These are algorithmic/storage counts, not
whole-solver throughput or cross-ISA bit-equality claims.

High-order residual reconstruction/projection and adaptive full-equation
assessment additionally share one prepared `fs-fft::GeneralFft` per AFT grid.
Power-of-two grids select this route from order 32; odd verification grids
select it from order 256 because Bluestein needs a larger convolution. The
actual quadrature length, time samples, retained order and independent finer
grid are unchanged. Per-channel transforms reuse one data/scratch allocation;
with default FFT features there are no new per-residual allocations or source
evaluations. Small orders
retain compensated direct sums. Exact tangent projection continues to use the
original compensated owner; it is checked against independent high-order
convolution with FFT residuals active, to floating-point accuracy.

Reconstruction uses an unnormalized forward transform of conjugated Fourier
bins, avoiding an intermediate multiplication by sample count. DC is added
after AC reconstruction. Projection transforms only the already normalized
deviations from the first source sample; its DC sum remains compensated and a
constant source retains exactly zero AC. Adaptive checks cancel original static
terms before adding AC and still measure the full sampled equation and its
omitted part. A refused transform falls back to the original compensated sum
on those same samples; it cannot hide source nonfiniteness or change equations.

The additional payload is the FFT owner's admitted preparation/plan/scratch
bound plus `16*(M+d*(N+1))` bytes for its data and retained coefficients. The
root grid remains available for direct checks and tangents. Continuation
admits twice the largest FFT payload over its entire coarse/fine ladder BEFORE
source callbacks; checking only the finest odd grid would miss an active base
power-of-two plan. Forced-sensitivity preparation separately admits the sum of
its actual two FFT grids and reuses their projection for source directions.
This is a measured speed/memory tradeoff, not a claim that
FFT buffers are smaller than the root-only representation.

Bounded native measurements on the same x86_64 host with the pinned
2026-07-06 toolchain and optimized actual numerical owners supported the
separate route thresholds. At M/N=256/32, reconstruction was 19.35→6.26 us
and projection 18.87→6.79 us; at 2048/256, 1528→74.1 us and 1338→77.4 us.
At 4097/256 the odd-grid gains were smaller (2856→1916 us and 2656→1786 us).
The 1025/64 projection instead slowed from 177→442 us, so that route stays
direct. Maximum compared waveform/coefficient differences were 5.33e-15 and
3.34e-16. These are transform-kernel comparisons on a shared host, not full
Newton throughput, a quiet-machine benchmark, or cross-ISA bit equality.

The transform regression cases check every retained complex coefficient on
even/odd grids, exact-zero AC under a large constant preload, nonlinear
convolution and its tangent, refinement to an initially omitted 65th harmonic,
pre-callback continuation memory admission, and finite physical coefficients
when an FFT intermediate overflows. Fifteen orbit unit tests, 51 existing
harmonic/backbone/continuation tests and seven forced-sensitivity integration
tests passed in native builds of
the complete orbit crate with the actual math/SIMD/FFT and linear-algebra
numerical owners. Seven original exact-length FFT tests also passed. Unrelated
executor/orchestration services were excluded from dependency builds; this is
not a full Cargo/workspace gate.

Physical row/column conditioning, phase masks, the analytic backbone-amplitude
row, and the full product rule for equilibrium deflation remain in the same
Newton solve. The unknown-frequency column still probes the complete original
residual. The finer-grid residual check and adaptive omitted-band check remain
independent acceptance gates. A source tangent does not certify convergence,
continuum accuracy, derivative correctness, or a nonsmooth branch choice.

The `fs-couple::phs_orbit` consumer supplies the exact `(J-R)Q` tangent for
explicitly exact quadratic storage. Its explicit nonlinear path uses the
source owner's supplied Hamiltonian Hessian action. Original physical memory,
storage, resistance and forcing remain in the model; feedback derivatives are
not inferred. Forced continuation reuses the same tangent workspace; shooting
retains its existing differentiation path.

## Static loading and weak harmonics

The retained and independently assessed omitted-band AFT projections share one
constant-preserving owner. The first island sample is the reference; subsequent
samples are converted in place to normalized deviations. Compensated sums form
the DC mean by restoring the reference once and project AC only from deviations.
Equal finite samples therefore have exactly zero AC coefficients on every grid.
Normalization precedes accumulation, with a divided-difference fallback when
opposite finite sample extremes would overflow their raw difference. Waveform
reconstruction accumulates AC before adding DC once, including `HbOrbit::sample`.
These changes retain the original equation, budgets and sampling/order gates;
they add no sample buffer and cannot recover information already lost when a
source evaluates a weak variation beside a large offset in finite precision.

G1/G3 regressions in `tests/harmonic_balance.rs` cover large constant drives,
weak complex forcing within source-sample resolution, the omitted-band check,
and state reconstruction. `src/tangent/tests.rs` compares the centered derivative
against direct sampled products and finite differences, and checks exact
constant-Jacobian DC-to-AC decoupling. This supports
`frankensim-irp-v5-instrument-realism-pexfp.175`; it does not claim that bead's
graph/time-domain threshold acceptance.

## Conservative backbone continuation

`continue_backbone_branch` continues a conservative periodic family without
introducing a physical continuation parameter. The two complete physical
seeds each declare a positive first-harmonic peak in any original state
coordinate. Both are corrected against every original balance equation,
the shared full-waveform phase condition and that amplitude constraint.
Seed convergence metadata is not trusted. Before the second correction,
all harmonics of every state are coherently phase-rotated together.

Later points replace the amplitude row with a normalized Parseval
pseudo-arclength plane. Frequency and every Fourier coefficient are jointly
corrected, so frequency or selected amplitude may turn along the family.
The metric uses the original independent state, equation and angular-frequency
scales. Static or nodal first states, dense nonsingular equation mixing and
explicit conservative coupling remain present. No balance row, loss term,
coupling or memory coordinate is removed.

With nc real Fourier coefficients, the corrector has nc+1 unknowns and nc+2
rows: the original nc balance equations plus phase and amplitude/arclength.
The conservation-dependent equation is retained. The original fs-la
Householder QR solves this rectangular system after finite, unit-norm column
equilibration. No normal equations, pseudoinverse or diagonal jitter are
introduced. Every iteration, including an already residual-zero seed, checks
QR pivots and the authored maximum estimated 1-norm condition of its
equilibrated triangular R factor. The bounded five-iteration Hager/Higham
screen is a numerical estimate, not a rank or conditioning certificate.
Extra unresolved neutral directions refuse instead of selecting an
arbitrary family silently.

Original AFT, exact island coefficient tangents, full-waveform phase and
physical scaling remain shared owners. Declared tangent refusal cannot
become finite differences. Undeclared coefficient tangents and the frequency
column use centered normalized probes. Backtracking evaluates the full
original augmented residual. Each accepted point passes an independent
finer-grid retained-balance and sampling-defect check; physical HB diagnostics
remain separate from dimensionless phase and amplitude/arclength errors.

`BackboneContinuationBudget` bounds unknowns (also capped at 2048), retained
points, Newton iterations, step retries, complete residual evaluations and
numerical payload before source work. Payload includes simultaneous AFT/FFT
workspaces, optional analytic tangents, dense QR work and every retained
coefficient row/header and residual trace. Source callbacks are immutable and
must be deterministic. Cancellation is polled around sampled residual/tangent
passes, between columns, around QR and between steps; each native sampled or
dense factorization call remains a bounded noninterruptible unit. Success
contains exactly steps+2 points; refusals publish no partial branch.

The six new integration tests verify independent energy-period quadrature and
every velocity harmonic under original-state/row rescaling, mixed phases,
exact/probed tangents and reciprocal two-mass coupling. A 14-point
softening-hardening family crosses a genuine frequency projection turn. Its
order-21 maximum relative energy-period frequency discrepancy is 1.108e-9;
order 15 failed the unchanged physical oracle despite converged retained
balance. Additional-neutral-direction seeds, physical damping, resource
limits, cancellation and authored tangent refusal are covered. All 39 focused
native tests pass, including the original backbone, autonomous tangent,
autonomous continuation and full orbit unit suites. This compiles complete
fs-orbit with actual numerical owners; it is not a full Cargo workspace gate.

Harmonic order stays fixed. Retained-balance acceptance alone does not resolve
omitted harmonics or prove continuous-time accuracy. This entry does not
classify stability, primitive period, branch identity, bifurcations, hybrid
cycles or separatrices. The original source owns its physical bandwidth and
constitutive validity domain.

## Autonomous cycle continuation

`AutonomousScales` separately conditions original states, balance rows,
frequency and parameter. State coefficients use Parseval weights, so exact-zero
harmonic padding does not change the arclength metric. Large DC offsets do not
relax the absolute scaled residual tolerance or non-DC oscillation floor.

The corrector retains every original HB row and adds an integral phase row
from all non-DC state harmonics, plus a fixed-parameter seed row or arclength
plane. It works when the first state has no oscillatory component. Frequency
and physical parameter are independent unknowns. A phase-projected secant
predicts later points; shared exact island tangents when supplied, central
differences for the remaining columns, native LU, residual-decreasing
backtracking and bounded step halving/growth correct them. No damping,
equilibrium deflation or diagonal jitter is introduced on this path.

Original AFT buffers and trigonometric tables are reused within each corrector.
The lazily admitted tangent workspace is retained across the complete branch,
including corrector retries, and is used only while the current source declares
the derivative. A cached workspace cannot preserve a stale availability claim.
When `island_state_jacobian` is declared at the current physical parameter,
the original `tangent::Workspace` prepares every coefficient derivative from
the current waveform before any frequency or parameter perturbation. The same
physical row/column scales and exact phase/arclength rows apply. Frequency and
parameter columns retain centered original-residual probes; absent tangents
retain all original probes. Refused, nonfinite or incomplete declared tangents
propagate without fallback. Optional tangent storage, including replacement
phase tables at refinement, is admitted within the existing payload cap before
allocation. It does not increase the payload requirement of a family that
never supplies the callback. The residual limit counts actual HB evaluations;
`tangent_preparations` separately reports exact tangent preparations. Tangent
samples are bounded by the existing order, quadrature and Newton caps.
Cancellation is checked on both sides of a complete sampled tangent and at
each coefficient column, preserving the original atomic AFT/linear-algebra
boundary.

The existing independent odd-grid sampling test applies to every accepted
point. Its `HbOrbit` retains physical coefficients/frequency and unscaled HB
residual/trace/defect; separate fields retain the combined dimensionless
residual, phase/arclength rows and scaled sampling discrepancy.

Dimensions, full retained output/AFT/dense payload and total residual calls
are bounded before the corresponding work. Cancellation is polled between
residuals, Jacobian columns and steps; one native residual or LU remains atomic.
On refusal the mutable family parameter returns to the last accepted value,
or the first seed parameter before acceptance. No partial branch is published.
Callbacks must be deterministic and may not conceal evolving state.

The authored `tests/autonomous_continuation.rs` checks an analytical fold of
cycles with a static first coordinate, amplitude-dependent frequency, common
phase shifts, independent state/equation unit changes, zero harmonic padding,
pre-source resource refusals, equilibrium refusal and parameter restoration.
Its analytical family is lambda=(r^2-1)^2, omega=1+0.2*r^2. These are finite
numerical regression claims, not existence, primitive-period, branch-identity,
stability, fold-location or infinite-tail certificates. Extra-neutral
conservative families and hybrid cycles remain outside this entry's claim.

`tests/autonomous_tangent.rs` adds an independently soluble nonlinear cycle
with two retained feedback memories and a parameter-dependent port. It compares
the complete branch against residual probes and independent physical-unit
changes, checks memory phasors, dynamic tangent availability, original refusal
and restoration, and lazy payload admission. Source-call reductions are operation
counts, not timing claims.

## Invariants

- The original `solve_hb` / `solve_hb_seeded` and continuation APIs use
  fixed, disclosed harmonic truncation `N` (X-Struct).
  AFT starts with `m = (4N+4).next_power_of_two()` samples and uses
  prepared det-routed trigonometric sums. Basis tables and state/RHS
  buffers are reused through Newton and finite-difference Jacobian
  evaluations; the DC port is evaluated once per residual.
- A converged candidate is re-evaluated on a finer `2m+1` grid.
  Both its residual and the norm of the change in retained balance
  coefficients must meet the same tolerance and coefficient scaling
  (explicit absolute physical normalization for the autonomous branch above).
  Otherwise Newton resumes on the finer grid. There are at most four
  refinements; the fourth is a verification grid. Persistent disagreement
  returns `SamplingUnresolved`, and refinements also consume the existing
  Newton iteration budget.
- This sampling gate also applies to every continuation corrector.
  Accepted points retain their actual residual trace and final verified
  residual, never a fabricated zero. The original forced corrector includes
  arclength in that norm; the autonomous result separates physical HB norms
  from its normalized phase/arclength diagnostics.
  `aft_samples` records the final verification grid and
  `sampling_defect` records the measured retained-coefficient change.
- The adaptive entry points grow solved orders as `min(2N+1, max_harmonics)`.
  Each candidate's original equation is independently projected through
  harmonic `2N+1`, including harmonics above the maximum solved order when the
  cap is reached. Its complete sampled equation defect is also measured, so
  nonlinear harmonics outside the projected band cannot escape assessment
  merely because of that cutoff. The linear action is reconstructed from its
  retained coefficients; omitted state coefficients are exactly zero, so their
  equation defect is `-G_n`. No undefined or unnecessary higher-frequency port
  is evaluated, and the sampled check reuses the original island evaluations.
  The original sampled defects are saved in the no-longer-needed state buffer
  before island centering. Static linear/source terms cancel before AC is added,
  so this check preserves the constant-load and weak-harmonic arithmetic above.
  Assessment quadrature uses the same bounded finer-grid policy, starting no
  coarser than the candidate's verified sampling. Both projected coefficients
  and full/omitted sampled RMS measurements must agree between the grids.
- Adaptive acceptance requires each full measured per-equation RMS defect,
  taking the larger of its projected-band and complete sampled measurements,
  to fit its single physical tolerance. Sampled omitted energy is measured
  after subtracting the retained Fourier residual from the full equation
  defect, avoiding cancellation between nearly equal squared norms.
  Separate band diagnostics do not each receive the whole error budget.
  Quadrature disagreements must also fit the physical tolerances,
  plus waveform and frequency agreement between two actually solved orders.
  Parseval weights are one for DC and two for positive-frequency coefficients.
  Scaling is per row/state, so a large DC bias or another coordinate's units
  cannot dilute a small but relevant harmonic. The maximum must leave room
  for at least one larger solve. A retained-row failure refuses immediately
  rather than increasing order to compensate for a loose Newton tolerance.
- Warm starts retain every solved coefficient and move the solved frequency
  to its new canonical slot. Forced and autonomous phase gauges stay unchanged;
  comparison of conservative backbone shapes aligns only their free phase
  using the explicitly pinned first-state fundamental. Frequencies and
  unknown-frequency seeds must remain positive and finite.
- Canonical unknown ordering: DC block, then per-harmonic
  `Re`/`Im` blocks, state-major; the packing is the seed format.
- The explicitly scaled entry uses one row/column unit for DC and that unit
  divided by sqrt(2) for each positive-frequency real/imaginary coefficient.
  Thus Euclidean norms in the conditioned chart are real-waveform RMS norms.
  Equation units are independent of state units, as required by general
  mechanical/acoustic/electrical frequency ports. Newton differences use each
  coordinate's declared physical magnitude; the dense solve is conditioned on
  both sides, and physical steps are restored before original-model evaluation.
  The backbone's extra row is its relative fundamental-amplitude error. No
  amplitude or frequency is silently normalized in the returned waveform.
  Invalid scales, physical seeds and budgets refuse before island/port calls.
- Backbone acceptance independently bounds `|2|X_1[0]|/a - 1|` by the
  requested tolerance. Other state coordinates and DC offsets cannot relax
  this relative amplitude gate. The appended Newton row is the unsquared
  norm difference `|X_1[0]|-a/2`, with its analytic norm gradient; this avoids
  squared-unit scaling and differencing a norm with a perturbation larger
  than the physical oscillation. The dynamical residual's coefficient scale
  uses only the non-DC coefficients, with the prescribed amplitude as its
  floor instead of one SI unit, so static preloads cannot hide an incorrect
  frequency or unresolved oscillatory balance.
  Backbone state-column differences use each coordinate's retained
  oscillation scale (the target supplies a zero-coordinate fallback);
  the frequency-column and non-backbone differencing rules are unchanged.
  Deflation scales only the physical equations, consistently in the Jacobian,
  right-hand side and line-search merit; it does not scale the amplitude row.
- A converged residual below the declared tolerance is the ONLY
  accepted discrete-balance evidence. Unanchored autonomous HB solves
  refuse equilibrium by name (`TrivialCollapse`), checking the full
  non-DC spectrum rather than the fundamental alone. A backbone's independently
  resolved prescribed positive amplitude excludes equilibrium without an
  arbitrary absolute displacement floor. The original
  absolute autonomous HB collapse threshold is 1e-9 in packed state units;
  the new branch uses its explicit scaled RMS `minimum_oscillation`.
  Forced HB solves admit DC, zero, arbitrarily small, and higher-harmonic
  responses: their period is prescribed, not inferred from a cycle.
  Equilibrium deflation applies only to the original autonomous/backbone anchors.
  Shooting applies its derivative-based equilibrium check only to
  autonomous problems: a forced response can have zero velocity at an
  ordinary turning point.
- Every RK4 stage, shooting Jacobian, monodromy evaluation, and autonomous
  equilibrium check uses the admitted time-domain callback. Shooting never
  substitutes `island` for a general linear port's realized dynamics.
- In the original HB solver, conservative families carry a parity/phase
  nullspace: its Newton step regularizes a singular factorization with FIXED
  Tikhonov jitter; the converged residual gate is unchanged. Autonomous
  parameter continuation does not use that regularization.

## Consistent discrete variational shooting

The variational entry requires the complete original explicit smooth flow and
its row-major state Jacobian. Each callback must write every finite entry;
an unavailable, false or incomplete tangent refuses without a trajectory-probe
fallback. The frequency-domain island and port retain their own semantics.
For example, an island Jacobian for `M*x'=g(t,x)` does not include the mass
inverse required by `df/dx`. Realized memory coordinates cannot be omitted.

One RK4 stage owner now serves both shooting entries. The new entry propagates
the state and its full derivative at the same four source stages. An autonomous
period column differentiates the actual discrete step size at every stage,
including the stage-state chain rule. It is not replaced by the continuous
endpoint vector field. The original entry retains its forward-difference
Jacobian, original phase policy and Floquet classification, using the same
primal stages; unwritten or nonfinite time-domain rows now refuse explicitly.

Let `S=diag(state_scales)`. Newton operates on `z=S^-1*x`, with absolute
infinity-norm closure tolerance `ShootBudget::tolerance` in this chart. Every
physical state has its own positive scale: a large static preload cannot
relax another state's tolerance. The autonomous period variable is
`T/period_scale`. A fixed phase plane through the seed uses the normalized
full scaled seed vector field, and its independent residual must meet the
same tolerance. All original closure rows remain present; no first-coordinate
motion assumption or row replacement is made. An autonomous cycle must exceed
its declared normalized excursion floor, which must itself exceed closure
tolerance. An equilibrium or unresolved short monotone drift cannot be
accepted solely because its closure defect is small. Forced zero, static,
tiny and turning-point responses remain admissible at their prescribed period.

The returned `monodromy` is explicitly `S^-1*D(Phi_T)*S`, not an undisclosed
original-unit matrix. Its eigenvalues come from the existing `fs-la` complex
eigensolver. The autonomous `period_derivative` is recovered in original
state units per second. Periodic state, physical residual trace, normalized
closure and phase residuals, observed excursion and source/map counts are
reported separately. Newton uses the existing dense LU, without jitter, and
a bounded eight-trial halving line search; rejected trials count against the
period-map cap. There is no reduced state, second integration algorithm or
finite-difference monodromy in this entry.

Before callbacks or numerical allocation, controls admit the complete state
dimension (hard ceiling 256), positive scales, finite tolerances, iteration
and period-map caps, and checked payload/work envelopes. For dimension `d`,
the payload ceiling is `384*d^2 + 1024*d + 64*(max_newton+1)` bytes, excluding
source/caller storage and allocator overhead. The work score is
`maps*steps*(16*d^3+128*d^2+128*d) + (128*max_newton+8192)*d^3`; it is an
admission envelope rather than measured FLOPs or a wall-clock guarantee.
Stages reuse their storage. The shared allocation-free matrix-product kernel
is called in row tiles, with cancellation between rows and before/after
source callbacks. Dense LU/eigen calls remain separately bounded atomic stages.
Cancellation or refusal returns no partial result and never mutates the source.

This is a finite smooth discrete periodic-response and multiplier result.
Source equivalence, temporal resolution, primitive period, continuous-orbit
existence, conditioning and physical validity require separate evidence.
The new entry reports multipliers without inferring a torus, global attraction,
or hybrid saltation. Analytic tests cover mass-weighted explicit dynamics,
time-periodic nonnormal maps, cubic forced response, nonisochronous autonomous
cycles with a static first coordinate and large preload, physical-unit
covariance, exact discrete period derivatives, deterministic replay and typed
source/budget/cancellation refusals.

## Error model

Typed `OrbitError`: `BadParameter`, `NewtonStalled` (with the full
per-iteration residual trace), `SamplingUnresolved` (finest sample count,
residual and measured defect), `HarmonicsUnresolved` (last finite-band check,
including whether retained equations or omitted/order-change evidence failed),
`SingularJacobian`, `TangentRefused` (the failing AFT sample),
`ContinuationExhausted`, `TrivialCollapse`, `TimeDomainUnavailable`,
`TorusSuspected` (a
non-trivial unit-circle Floquet pair with nonzero angle:
quasi-periodic drift is a NAMED no-claim in v1, never a wrong
answer), forwarded `Eigen`.

`autonomous::AutonomousContinuationError` adds typed argument, cancellation
and residual-evaluation-limit refusals while retaining native numerical errors.

`VariationalShootError` adds invalid-control, cancellation, payload/work/map
budget, unavailable explicit flow/tangent and source-stage refusal variants,
and retains the original numerical `OrbitError` through `Numerical`.

## Determinism class

Deterministic: bounded deterministic grid refinement, fixed iteration caps and step policies,
det-routed transcendentals, no RNG, no time dependence; bitwise
repeatability is a conformance case.

## Cancellation behavior

Synchronous computation; no direct `Cx` dependency. Original HB/shooting and
forced continuation do not poll cancellation. The autonomous continuation and
variational shooting entries accept explicit cancellation callbacks; adapters
may supply their original `Cx` checkpoints. Native AFT residuals and shared
dense LU/eigen calls remain atomic; variational source stages and row-tiled
matrix products have explicit cancellation boundaries.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/orbit_conformance.rs` (ob-001..ob-006) — oracles with no
music in them: the conservative Duffing backbone vs the fs-nlmodal
pinned first-order law with the band authored at the NAMED
second-order Lindstedt coefficient `(15/256) eps^2 a^4` (measured:
the HB deviation IS that term); van der Pol at mu = 1 against TWO
independent published sources (Amore arXiv:2111.12198
T = 6.663286859323130 — matched to nine digits — and
Gasull–Giacomini–Grau arXiv:1602.00113 T ~ 6.6632866) plus the
weak-mu perturbation laws; pseudo-arclength continuation TRAVERSING
both folds of the forced-Duffing S-curve with fold frequencies
within 0.2% of the independent first-harmonic scalar law; shooting
cross-validating the same vdP orbit to 1e-12 relative period with
Floquet classification (trivial multiplier 1 + transverse 8.6e-4);
refusals by name including `TorusSuspected` on a constructed
incommensurate two-center fixture; bitwise determinism. Additional G1/G3
regressions refuse a frequency-only port before island evaluation, compare
HB and shooting for `2 x' = -x + cos(t)` against the exact phase-zero value
0.2 and Floquet multiplier `exp(-pi)`, and retain a forced cosine response
whose initial derivative vanishes.

The G1/G3 analytic regressions in `tests/harmonic_balance.rs` compare
higher-harmonic forced oscillator coefficients (including their phase and
velocity) with the closed-form transfer function, exercise small and DC
forced responses, and retain autonomous-equilibrium refusal. A seventh-order
nonlinear spring exposes the old eight-point alias (0.9968112456 amplitude
instead of the retained Galerkin answer 1.0); an independent shifted
65-point quadrature checks the corrected coefficients. Further regressions
verify truthful continuation residuals against the analytic scalar balance
and a bounded named refusal for persistently unresolved forcing.

The physically scaled path additionally recovers the exact fundamental and
generated third harmonic of a nonlinear cubic cascade under independent state,
equation and clock-unit changes, including picometre-scale coordinates. A
1e12 unrelated static coordinate cannot hide its missing third harmonic; an
insufficient Newton budget refuses. Other G1/G3 cases retain prescribed
backbone amplitude/frequency under free phase, refine the seventh-order alias
in small physical units, refuse invalid scales before source calls, and replay
the accepted coefficients and scaled residual bit-for-bit on the same toolchain.

Backbone regressions additionally correct a seed with twice the prescribed
micrometre or picometre displacement, under phase shifts and optional
deflation, without rejecting the resolved small motion as equilibrium. A
large unrelated static coordinate cannot hide an unresolved amplitude or an
incorrect frequency; adequate iterations recover both.
A nonlinear Duffing coordinate-rescaling test compares the complete retained
waveform and frequency after changing displacement scale by 1e-6; invalid
zero, negative, nonfinite and unrepresentable-half amplitudes, and relative
tolerances outside `(0,1)`, refuse before model callbacks. These verify the
discrete amplitude constraint and scale
handling, not continuous-orbit existence or a global nonlinear convergence
guarantee.

Adaptive-order regressions use the exactly soluble nonlinear cascade
`x' = -x + cos(t), y' = -y + x^3`: the third harmonic is absent from an
N=1 solution even when its retained balance is converged. The adaptive solve
checks orders 1, 3, and 7 against exact complex coefficients, and requires a
second actual solve after the missing harmonic enters. A capped N=2 solve
detects the omitted third harmonic without querying ports outside its domain;
rescaling an equation and adding DC bias do not mask that defect. Further
cases verify backbone-frequency warm starts, retained-residual refusal under
a deliberately loose Newton tolerance, admission before model callbacks, and
refusal when separately sub-budget retained/omitted bands exceed the single
combined physical RMS budget.

The ninth-harmonic cascade `x' + x = cos(t+phase) - sin(t+phase),
y' + y = T_9(x)` supplies only fundamental forcing and generates its ninth
harmonic through the Chebyshev polynomial. Orders 1 and 3 agree while their
projected residuals through 7 are zero. Independent full-equation RMS still
detects the missing ninth harmonic and refuses a cap of 3 without evaluating
an out-of-domain port. With a sufficient cap, orders 1, 3, 7, 15 and 31 recover
the exact amplitude and phase `Y_9 = exp(9i*phase)/(2*(1+9i))`, with no invented
lower output harmonics.

The scaled adaptive cubic-cascade regression changes state and equation units
independently, retains a large unrelated preload, and actually solves orders
1, 3 and 7. A cap at order 2 still checks through order 5 and refuses the
analytically known missing third harmonic. This exercises physical conditioning
and spectral resolution together rather than normalizing a solved result.

The tangent unit regressions compare complete coefficient actions against an
independent cubic convolution and an explicitly time-periodic linear matrix's
known Fourier coefficients, including DC and complex phase. Native analytic
HB regressions recover the same nonlinear fundamental/third harmonic with
fewer source evaluations, retain mixed physical units and unknown-frequency
backbones with and without deflation, and refuse false or incomplete tangent
callbacks before any substitute perturbation. These are focused G1/G3 numerical
claims; source-call counts are not wall-clock performance measurements.

## Driven sidebands about a complete periodic orbit

`sideband::PreparedPeriodicLinearization` prepares the exact original state
tangent along a supplied `HbOrbit`. For carrier `nu` and base frequency `omega`,
it solves the complete coupled equations
`sum_l [delta_kl P(i*(nu+k*omega)) - A_(k-l)] X_l = F_k`.
All original coordinates, full port matrix entries and both signs of every
retained sideband participate. Coefficients are peak complex perturbations in
`Re sum_k X_k exp(i*(nu+k*omega)*t)`, unlike the base orbit's Fourier half
coefficients. Coincident signed physical frequencies must be combined coherently
by consumers; the owner does not fold frequencies or add a damping floor.

Preparation independently checks the complete sampled base equations and
normalized tangent Fourier coefficients on two grids. It requires an exact
declared tangent and refuses unavailable, partial or nonfinite derivatives.
The cached tangent serves multiple carriers; `factor` uses the shared complex
LU owner once per carrier and each factor serves multiple original-equation
forcing vectors without further source/tangent callbacks. Independent row
contractions report retained backward error, and original-equation defects
are checked in the finite omitted bands `N < |k| <= 2*N+2`.
Exactly one-entry homogeneous equations retain their algebraically required
zero unknowns after floating-point pivoting. Up to four original-residual
corrections reuse the same factor when needed; every retained coupled row
must still meet the strict componentwise backward-error allowance. No numerical
zero threshold changes an equation or removes a response band.
Independent state/equation scales, numerical payload/work admission and
cancellation apply before publication. No physical source is mutated.

`sideband::refinement::solve_sideband_adaptive` adds bounded automatic retained
order resolution. Its input declares the complete initial signed forcing
support, which is preserved by zero-padding as orders grow `N -> 2*N+1` up to
the authored cap. Intermediate finite omitted-band defects request refinement;
singularities, retained backward-error failures and unavailable/refused exact
tangents still refuse. Acceptance requires two actually solved orders to agree
in each original state's physical complex-envelope L2 tolerance.

Low-order agreement alone is insufficient: an exact tangent containing a ninth
harmonic can produce distant sidebands that neither of two low orders retains
or checks. Every adaptive candidate therefore evaluates the FULL ORIGINAL
linearized equation envelope on independent shifted grids of `2*M+1` and
`4*M+3` samples. Original complete port matrices at every signed retained
frequency and exact tangent callbacks at the original base waveform participate;
no time-domain realization is inferred from a dispersive port. Each original
equation has its own physical envelope tolerance, independent of conditioning
scales. Full sampled defects, grid differences, finite omitted-band defects and
actual response changes are returned in a refinement trace. Complex-envelope
norms give each signed coefficient weight one, not the real-waveform AC factor
two; physical frequency coincidences remain the consumer's coherent sum.

Adaptive quadrature grows with retained order. The fixed owner's work allowance
is charged cumulatively across the actually attempted solves AND independent
assessments, before source evaluation for each step. Its payload allowance also
includes the previous solution and all live refinement scratch. Assessments
stream samples; they do not retain a samples-by-state history. Neither a partial
response nor an accepted response is published on exhaustion or cancellation.
These measurements detect generated energy outside finite Fourier guard bands,
but do not certify absence of arbitrary aliases beyond both sample grids or
replace independent base-waveform refinement.

`sideband::sensitivity` differentiates the retained coupled response with respect
to the signed carrier, at a FIXED original operating orbit, clock and phase.
`PreparedSidebandResponse::prepare_carrier_derivative` admits and caches the
declared complete `dP(i*w)/dw` at every signed shifted frequency. This derivative
includes the imaginary-axis chain factor; it is not merely `dP(s)/ds`. The exact
native/authored derivative must cover all original mass, damping, gyroscopic and
dispersive/memory terms. Unavailable, incomplete or nonfinite callback results
refuse without finite-difference substitution.

Each sensitivity query supplies both the complete physical forcing and its
carrier derivative, including any frequency-dependent source law. The existing
LU and correction owner solve `D X = F` and `D X' = F' - P' X`; no second
factorization or differentiation of solver iterations is introduced. A separate
double-double contraction checks the full differentiated original equation
`D X' + P' X - F'`, independently of the rounded formed RHS. Differentiated
omitted-band defects have their own allowance, in original-equation scale units
per angular frequency. The complete primal/derivative pair is published only
after both pass. Prepared port-derivative blocks serve later right-hand sides
without source callbacks.

Explicit additional payload/work caps include derivative preparation, the
retained primal result, differentiated solve and assessments; they add to the
original factor/solve caps and apply per query, not to an unbounded campaign.
The returned derivatives concern complex coefficients, not the time waveform or
a changing nonlinear operating branch. Finite-order backward-error evidence
does not supply a forward-sensitivity bound near an ill-conditioned resonance.
`tests/sideband_carrier_derivative.rs` checks full dispersive/gyroscopic ports and
frequency-dependent forcing against closed two-coordinate inverse derivatives,
phase-dependent periodic conversion against exact rational derivatives, and
independent derivative truncation, budget, source and cancellation refusals.

`sideband::sensitivity::refinement::solve_sideband_carrier_derivative_adaptive`
now resolves the response AND its carrier derivative. Physical equation and
state-change allowances for the derivative are independent and carry units per
angular frequency. Each actual order solves both retained systems with one LU,
then releases that factor before streaming both full original equations through
the same shifted `2*M+1` and `4*M+3` assessment grids. The derivative assessment
includes `P_nu X - F_nu` from cached ORIGINAL port-derivative entries, rather than
trusting only the rounded derivative RHS. Each tangent sample serves both
fields. Neither a low-order response plateau nor a small primal defect can
stand in for derivative convergence. The core and added derivative work caps
separately accumulate across all attempted orders; their live payload caps add.

`PreparedPeriodicLinearization::factor_at_order` factors an actual smaller
retained system from the same checked tangent, without resampling that tangent
or truncating a larger solution. The prepared owner freezes only the physical
base coefficients and clock in an explicitly budgeted snapshot, preserving the
existing constructor lifetime contract and excluding diagnostic history.
`solve_carrier_derivative_assessed_at_order` supplies candidate pairs with full
equation measurements and deferred omitted-band acceptance. Consumers must
inspect the checks and compare actual solves with `compare_with` before claiming
resolution. This supports low/high comparisons at every queried resonance
frequency with one common tangent and at most one live LU. Common quadrature
comparisons measure retained-order convergence; base and sampling adequacy
remain separately required. The finite measurements are not an infinite-tail,
forward-sensitivity or stability enclosure.

`tests/sideband_response.rs` checks a phase-shifted nonlinear cascade against
analytic converted sidebands, a parametrically driven oscillator against
independent real-time integration and order refinement, and full coupled
mass/damping/stiffness/gyroscopic/Maxwell ports against closed two-coordinate
responses. Stale base waveforms, aliases, exact singular resonance, unavailable
or refused derivatives, cancelled reuse and explicit resource refusals are
also exercised. These are focused numerical checks, not experimental evidence.
`tests/sideband_refinement.rs` additionally checks distant-band false convergence
against an independently soluble triangular periodic system, preservation of
signed complex forcing and disparate physical units, actual inter-order
agreement, and resource/source/cancellation refusal.

The source and tangent must actually repeat at the authored period. Finite
quadrature, finite checked omitted bands and small backward error establish
neither an infinite-tail bound nor a conditioning-dependent forward-error bound.
This is an infinitesimal particular response about the authored operating orbit;
it does not establish Floquet stability, startup decay or finite-probe accuracy.

## No-claim boundaries

- Quasi-periodic orbits (tori) are a NAMED refusal, not a claim.
- The original `continue_branch` remains forced-anchor. The additive autonomous
  entry solves frequency independently but is fixed-order and does not certify
  folds, stability, primitive periods, global branches or extra-neutral families.
- Collocation is described by the I09 slot but not implemented here;
  shooting is the second method.
- No a-priori truncation-error bound: the returned `N` remains disclosed
  structure. Fixed-order sampling refinement improves quadrature; the
  optional adaptive-order entry points also compare finite omitted bands,
  complete sampled equation defects and successive solved waveforms. Neither
  supplies a rigorous infinite-tail
  bound, an existence proof, or a conditioning-dependent forward-error bound.
  Agreement between finite grids is a measured convergence check, not a
  universal aliasing bound: sufficiently high frequencies can alias on both
  grids or exceed their sampling resolution. Explicit input bandwidth
  must be admitted by domain adapters. Continuation remains fixed-order;
  adaptive order does not add relative-periodic, hybrid, matrix-free, or
  resume capabilities. Shooting is an independent cross-check
  only for problems with an equivalent time-domain
  formulation. The explicit callback is a model-supplied realization, not
  a proof of that equivalence. Frequency responses are never automatically
  inverted or assigned a state-space realization by the shooting solver.
- The original shooting entry's monodromy remains forward-differenced.
  The additive variational entry differentiates its complete discrete RK4 map
  exactly, but does not supply an adaptive, hybrid or matrix-free flow tangent,
  a continuous-time error bound, or the rest of the I09 capability lattice.
- Grazing/hybrid (non-smooth) islands are out of scope in v1.

## Shared phase-grid verification

The high-order regression in `tests/analytic_harmonics.rs` checks physical
oscillator phase and amplitude through the 64th harmonic, including a drive
near resonance and an independently sampled odd verification grid. Root-grid
tests cover odd/even quadratures, conjugate symmetry and repeated periods
without phase drift. Existing nonlinear convolution, finite-difference tangent,
adaptive alias-refusal and continuation tests remain the integration checks.

A pinned-nightly native benchmark compiled the actual `basis.rs` and complete
`fs-math` source. For 512/1025 samples at order 64 and 2048/4097 samples at
order 256, its compensated reconstruction/projection comparisons against the
former tables had maximum round-trip error 6.25e-17 and maximum projection
difference 4.94e-15. Preparation and high-order projection improved locally;
some smaller warmed kernels were slower. Those shared-host measurements are
not stable release throughput guarantees. Full Cargo integration was attempted
but blocked by shared disk exhaustion and subsequent environment loss.

After execution recovered, a bounded native build compiled the complete orbit
crate against actual fs-math/SIMD and the original fs-la LU/eigenvalue/complex
LU/SVD owners plus their exact production GEMM functions. Unrelated executor,
budget and provenance surfaces of fs-la were omitted; no numerical algorithm
was substituted. The new five backbone-coordinate tests, all nine orbit unit
tests and all thirteen autonomous tangent/continuation regressions pass. This
is scoped numerical verification, not a full Cargo/workspace gate. The phase
comparison allows eight machine epsilons for independent full-angle argument
rounding. Payload-refusal tests now straddle the actual reduced allocations,
retaining their admission-before-callback checks without restoring old tables.

## Forced waveform parameter derivatives

`sensitivity::PreparedForcedSensitivity` differentiates every original DC and
complex Fourier coefficient of a fixed-order forced HB solution. It reuses the
exact real AFT coefficient Jacobian, the existing physical row/column scales and
one fs-la LU factor across physical parameter directions: `R_X X_p = -R_p`.
Complete residual declarations may include original source, frequency-port and
clock changes; the island-only entry requires a fixed clock/port and projects
complete source partials using the original centered AFT owner. Both independent
grids are reassessed, including the rounded published derivative in the original
physical linearized equations. No finite-difference tangent, coefficient clipping,
mode decoupling or regularization is introduced. Unavailable/incomplete source
tangents, stale base equations, singular or overly ill-conditioned factors,
inconsistent directions, exhausted budgets and cancellation explicitly refuse.

Returned coefficients are physical derivatives per parameter unit at fixed phase,
with the same half-amplitude convention as HbOrbit. The authored parameter scale
conditions arithmetic and absolute residual acceptance only. The cached source
must remain unchanged across directions. The condition number is an estimate,
not a certified upper bound; sampling agreement and fixed-order differentiation
supply no infinite-tail, existence, branch-uniqueness or stability guarantee.
Refine harmonic order and confirm time-domain behavior independently.

`tests/forced_sensitivity.rs` compares coupled nonlinear/relaxation force and
material derivatives with separately solved parameter perturbations, and checks
3/7/11/15-harmonic convergence against independently integrated continuous-time
variational RK4. It also covers complete phased frequency-port/clock derivatives,
physical unit changes, repeated-factor source-call counts and named refusals.

## Autonomous frequency and waveform parameter derivatives

`sensitivity::PreparedAutonomousSensitivity` extends the same exact AFT
preparation, physical scaling, LU and residual-refinement owner used by forced
sensitivities. An unknown physical angular frequency and the existing fixed
full-waveform Parseval phase hyperplane close the complete autonomous derivative.
Every original DC/real/imaginary coefficient and balance row remains present;
neither amplitude nor a selected coordinate is pinned.

The caller supplies a complete packed phase reference and an exact original
residual-frequency column at fixed coefficients and physical parameters. Both
original AFT grids independently check the base, tangent and frequency column.
Parameter directions supply the complete residual partial at fixed coefficients
and frequency, or an original island partial when the linear frequency port is
fixed. One cached bordered factor returns frequency and full-waveform derivatives
in original parameter units. Exact declarations never fall back to probes.
Invalid budgets, nonoscillatory references, incompatible base phase,
nonisolated/singular cycles, excessive conditioning, unresolved differentiated
equations and cancellation refuse without a partial result.

Four actual-source native tests compare with an independent nonlinear oscillator
with a changing mass/time port, a static first state and second-harmonic
relaxation memories affecting frequency. They verify fixed-gauge waveform
derivatives, physical unit conversions and factor reuse after refusal. All seven
existing forced-sensitivity regressions also pass through the modified shared
owner. This local fixed-order derivative establishes no existence, primitive
period, infinite-tail, uniqueness or Floquet-stability guarantee.

## Conservative backbone frequency and waveform derivatives

`sensitivity::PreparedBackboneSensitivity` differentiates a conservative family
at an authored physical first-harmonic amplitude of any original coordinate.
It shares the exact AFT Jacobians, physical scales, complete frequency-column
declaration, two-grid base/tangent checks and original residual projection of
the forced/autonomous sensitivity owner. Every original balance row remains;
a fixed full-waveform phase row and an amplitude row close the coefficients and
unknown frequency. The rectangular system reuses the original backbone
continuation's column-equilibrated Householder QR, without normal equations,
deleted balance rows, artificial damping or regularization of a neutral cycle.

`BackboneSensitivityDirection` independently states the parameter scale and
physical amplitude derivative. Zero amplitude derivative gives material,
geometry or source sensitivities at fixed first-harmonic peak; a zero original
source partial and unit amplitude derivative gives the backbone slope and all
harmonic waveform changes with amplitude. The complete residual callback holds
both coefficients and frequency fixed. The island-only entry requires the
original frequency port to remain fixed. All DC, higher-harmonic and material
states retain their original ordering and units; output uses the existing
`AutonomousOrbitSensitivity` waveform/frequency result.

One prepared factor serves every direction. Payload/work admission includes the
extra row and QR scratch before source callbacks. Existing square LU paths keep
their componentwise backward-error metric and arithmetic. The rectangular path
uses maximum row defect relative to maximum row action, avoiding division by
quadrature noise in analytically zero harmonic rows. Independently of that
relative diagnostic, BOTH original physical grids and both closure rows must
meet their absolute normalized tolerances using the rounded published derivative.
An incompatible dissipative direction therefore cannot be accepted merely
because its rectangular least-squares problem has a minimizer. Invalid anchors,
incomplete exact declarations, unresolved residuals, excessive conditioning,
resource exhaustion and cancellation publish no partial result.

The focused G1/G3 tests compare complete derivatives against independently
resolved neighboring nonlinear Duffing backbones and differentiated continuous
energy-period quadrature. They retain a static first state, dense coupled
mass/time ports, arbitrary phase, physical unit/parameter changes and higher
harmonics, and check incompatible damping, declared refusals and factor reuse.
These are local fixed-order numerical derivatives, not a proof that an arbitrary
source is conservative, a harmonic-tail bound, a primitive-period or stability
certificate. The selected first-harmonic peak is not the turning-point amplitude;
an amplitude fold or additional neutral family can make this chart refuse.

All four new regressions and all 22 original forced/autonomous sensitivity and
backbone coordinate/continuation regressions pass in a bounded native build.
It compiles the complete original `fs-orbit` crate, full `fs-math`, observability,
substrate and SIMD crates, and unchanged fs-la factor/eigen/complex-LU/SVD/GEMM
and fs-fft transform owners. Unrelated executor/provenance dependency surfaces
are excluded, with no substituted numerical methods. This scoped execution
does not replace a full Cargo/DSR workspace gate; the full dependency build was
stopped when resource pressure made the execution environment unavailable.


## Phase-fixed adaptive conservative QR points

`solve_backbone_adaptive` combines the original conservative rectangular QR
corrector with measured retained-order refinement. `BackboneRefinementControls`
bounds Newton work, complete unknown count, payload, residual-grid evaluations
and estimated conditioning. `AdaptiveBackboneOrbit` returns the final corrected
physical point, every order's physical refinement check and aggregate
evaluation/tangent counts.

The supplied complete waveform fixes one normalized Parseval phase hyperplane.
The selected original state pins `2*abs(X_1)` to its physical amplitude. That
phase and amplitude stay fixed as orders grow by `min(2H+1,max_harmonics)`; new
coefficients are zero-prolonged and then actually solved. All original balance
rows remain. Column-equilibrated Householder QR admits rank/conditioning even at
residual-zero seeds, without diagonal jitter, normal equations, removed balances
or artificial damping. Existing phase-free adaptive APIs retain their behavior.

Each candidate is assessed on independent grids using both projection through
`2H+1` and the complete sampled original time-domain equation. Acceptance requires
retained/omitted/full physical equation tolerances, sampling agreement, and
waveform/frequency agreement between two actual solves in the same fixed phase.
Per-row/state tolerances are absolute in original units; unrelated DC preload or
amplitude cannot relax them. A retained-row failure refuses immediately.
Nonlinear aliasing can prevent conservative Newton convergence before its
ordinary acceptance gate: stalled or weakly decreasing iterates therefore
request a finer-grid residual check. Only an independently measured sampling
defect above tolerance promotes that same physical iterate, within the existing
grid/Newton/global evaluation caps.

Admission reserves correction and assessment payload at every possible order
before numerical source callbacks. Caller convergence metadata is ignored
without cloning its arbitrary trace. Declared tangent refusal, rank loss,
resource exhaustion and cancellation publish no partial orbit. Finite numerical
convergence does not establish an infinite-tail bound, continuous orbit
existence, primitive period, branch identity or stability.

Four new tests cover strongly nonlinear cubic/quintic generated harmonics
against continuous energy-period quadrature and cycle energy; tiny motion under
extreme state/equation unit changes and large unrelated static preload; global
budgets/cancellation/source refusal and distinct retained versus omitted
resolution failure; and rejection of exact-residual seeds with an extra neutral
family. The 47 relevant existing continuation, coordinate, harmonic-balance,
autonomous-phase and backbone-sensitivity tests also pass in the original-source
bounded native lane; this is not a Cargo/workspace/DSR gate.

## Forced harmonic objective adjoints

`sensitivity::PreparedForcedSensitivity::prepare_adjoint` prepares one real
objective's implicit gradient with the existing forced-response LU. It retains
every original DC, real and imaginary harmonic balance coordinate, including
coupled material memory. The caller supplies the complete ordinary real-coordinate
objective gradient at fixed physical frequency. For a real state's mean square,
`q = X_0^2 + 2 sum |X_n|^2`, its entries are `2 X_0`, `4 Re X_n` and
`4 Im X_n`; no additional complex conjugation or harmonic weight is applied by
the adjoint.

With original row and coordinate units `D_r` and `D_x`, the existing conditioned
matrix is `A = D_r^-1 R_X D_x`. The new path solves
`A^T lambda_hat = D_x q_X / q_scale` using its cached transpose LU solve and
bounded iterative corrections, then returns physical residual weights
`lambda = q_scale D_r^-1 lambda_hat`. The componentwise relative solve residual
and BOTH original physical dual equations `R_X^T lambda = q_X`, independently
assembled on the coarse and fine AFT grids, must satisfy the prepared controls.
Objective scaling changes arithmetic and normalized absolute admission, not the
physical derivative.

`PreparedForcedAdjoint::gradient` accepts complete original residual parameter
columns, including material, forcing, frequency-port and clock dependence.
`gradient_island_direction` instead projects complete local-island partials
through the cached centered AFT owner and requires a fixed clock and linear
port. Both return `q_p - lambda^T R_p` in original objective/parameter units;
the explicit objective partial must include any direct clock dependence.
Both the full normalized residual-column difference and the normalized
contracted-gradient difference between the two grids are checked. Incomplete
declarations, unresolved dual equations or sampling changes, invalid units and
cancellation refuse without publishing a gradient.

The adjoint adds O(N) storage to the prepared O(N^2) factorization; it stores no
transposed dense matrix. Its bounded transpose solves and residual corrections
fit the existing sensitivity work allowance, and its peak extra vectors fit the
existing payload allowance. A complete-column parameter contraction is O(N),
apart from the caller's column construction; the island entry additionally
performs the existing AFT projection. No per-parameter waveform solve, new
factorization, baseline source/tangent evaluation or finite-difference probe is
introduced. These are operation counts, not measured performance claims.

Five focused G1/G3 regressions in `tests/forced_sensitivity.rs` check nonlinear
coupled force and material objectives against independent forward derivatives
and separately resolved neighboring orbits, frequency-port/explicit clock
dependence, and independent state/equation/objective/parameter units. They also
check incomplete declarations, cancellation, inconsistent sampling and refusal
of an original dual equation unresolved on one AFT grid.

This forced, fixed-order objective derivative supplies no existence, harmonic-tail,
branch-uniqueness or stability guarantee. The original source, objective
derivative and cached base must remain unchanged. Autonomous and conservative
phase/family closures are outside this entry.

## Cycle-mean harmonic objectives

`PreparedForcedSensitivity::mean_objective` evaluates a differentiable scalar
cycle-mean observation and its complete ordinary real coefficient gradient on
the two cached original AFT trajectories. The callback writes the value and
every physical state partial, including zeros; outputs are reset to NaN for
each sample and incomplete, refused or nonfinite callbacks publish no result.
The shared centered AFT projection retains DC and both complex quadratures;
positive-harmonic covector entries include the factor of two from real waveform
reconstruction. Its gradient feeds `prepare_adjoint` directly, without another
multiplicity, conjugation or sign change.

`ForcedMeanObjectiveControls` supplies positive objective units and separate
positive value/gradient sampling limits. Acceptance compares normalized mean
values and the full covectors scaled by the original HB column chart; matching
values cannot hide an unresolved objective derivative. The result returns both
measured defects. Existing scratch and two coefficient vectors fit the original
per-query sensitivity resource allowance, and there are no baseline source,
tangent, orbit-solve or factorization calls. Cancellation and refusal permit a
subsequent query on the unchanged cached trajectories.

These are forced, fixed-frequency/fixed-phase coefficient partials. Explicit
objective parameter, frequency and time-chart derivatives remain separate
`explicit_derivative` contributions to the adjoint; autonomous phase/amplitude
closures are outside this entry. Differentiability and harmonic refinement are
caller obligations; two finite grids do not bound the infinite harmonic tail.
G1/G3 regressions check Parseval and every real coefficient, independent state,
equation and objective units, nonlinear dissipated-power gradients against
independently resolved parameter changes, value-only and gradient-only aliasing,
complete callbacks, cancellation and refusal/retry without baseline callbacks.

## Autonomous harmonic objective adjoints

`PreparedAutonomousSensitivity::prepare_adjoint` supplies objective gradients
through an isolated self-excited cycle's complete waveform and unknown frequency.
It accepts every ordinary real coefficient partial and a separate objective
frequency partial at fixed coefficients. The original square bordered system
`[R_X R_omega; phase_X 0]` is transposed without dropping a balance equation,
freezing frequency, regularizing the neutral direction or creating another
factorization. The shared LU, bounded residual refinement and independent
physical dual checks serve both forced and autonomous objectives. Because the
raw Jacobian caches intentionally omit the dimensionless phase row, autonomous
dual verification explicitly restores its action in the original scaled chart;
the frequency dual equation includes the supplied objective frequency partial.

`PreparedAutonomousAdjoint` exposes the physical balance weights and the separate
fixed-reference phase weight. Its complete residual and local-island parameter
entries follow the forward autonomous sensitivity conventions: callbacks hold
coefficients and omega fixed, write only the original balance rows, and include
all direct material, source and frequency-port parameter terms. The fixed phase
reference has zero direct parameter partial. Returned derivatives are
`q_p - lambda_R^T R_p` in original objective/parameter units, with the solved
frequency response already included. Explicit objective parameter partials are
supplied separately and must not count the frequency chain twice. Complete
balance-direction and contracted-gradient sampling changes remain independently
bounded, including directions orthogonal to the objective weights.

The existing payload/work admission covers the added coefficient-sized storage;
no original source, state tangent, frequency-column or phase callback is repeated.
Refusal and cancellation leave the original factor and cached trajectories
reusable. Phase-sensitive objectives refer to the authored fixed phase chart.
Conservative families with an amplitude closure remain outside this square
adjoint entry; isolation, finite conditioning and harmonic order retain their
original limits, and no stability or infinite-tail certificate is inferred.
G1/G3 tests compare physical objectives with independent analytic memory cycles
and forward contractions, including changing mass/time scales, all state and
objective units, a frequency-only objective, a nonzero phase weight verified by
independently phase-aligned cycles, complete direction sampling and refusal/retry.

## Autonomous cycle means and period derivatives

`PreparedAutonomousSensitivity::mean_objective` evaluates a differentiable
cycle mean and its complete coefficient/frequency covector on the unchanged
cached autonomous trajectories. It reuses the original centered AFT observation
owner for the scalar value and every ordinary real coefficient partial. One
additional compensated scalar mean carries the declared frequency partial;
there is no second waveform reconstruction, source call or factorization.

The observation callback must write the value, all state partials, and its
complete frequency partial at FIXED phase and state on every sample. For
`f(t,x,omega)` with `t=theta/omega`, this partial is
`f_omega-(t/omega)*f_t`. A state-only mean supplies zero. Energy per cycle
from power P(x) supplies `(TAU/omega)*P`, scales its state gradient by the same
period, and supplies `-value/omega` as the frequency partial. The phase measure
of a normalized cycle mean introduces no additional period factor.
Original physical time-unit conversions remain the caller's authored chart.

The result's coefficient gradient and separate angular-frequency partial feed
the autonomous `prepare_adjoint` arguments directly. Independent sampling
admission includes the entire covector, using the original frequency scale as
well as every state scale. Agreement of values and state partials cannot hide
an unresolved frequency derivative. Incomplete/nonfinite declarations and
cancellation publish no partial result; the cached trajectories remain reusable.
The existing per-query allowance covers two coefficient vectors and the added
constant scalar storage and linear accumulation work. Explicit objective
parameter partials remain separate adjoint contributions. Differentiability,
fixed phase-reference semantics and harmonic refinement retain their original
obligations; two finite grids do not bound the infinite harmonic tail.

G1/G3 tests in `tests/autonomous_objective.rs` check energy-per-cycle derivatives
on independently known nonisochronous cycles, including changing period and
independent state, equation, time, parameter and objective units. They also
check phase shifts, exact Parseval covectors, callback counts, frequency-only
aliasing refusal, complete declarations, cancellation and deterministic reuse.


### Combined local-source and frequency-port sensitivities

`PreparedForcedSensitivity::solve_split_direction` composes a complete original
fixed-phase island partial with a complete packed linear frequency-port partial.
The cached original trajectories and centered AFT owner project the island;
the authored port column includes its coefficient, descriptor and clock changes.
A changing prescribed frequency is admitted without repeating baseline source
or tangent evaluation. A time-local source must include its fixed-phase clock
chain; no generic source-frequency law is inferred. Every port row must be
explicit and finite on both grids. The existing coupled factor and independent
physical residual checks admit the combined direction.

`PreparedForcedAdjoint::gradient_split_direction` shares that projection and
contracts the complete combined column with the original objective dual, without
a per-parameter waveform solve. Direct objective clock dependence remains in
the explicit objective partial. Existing numerical scratch/work bounds cover
the additional linear column; source callback work and storage keep their
existing exclusions. Cancellation/refusal preserves the cached base and reuse.
The split regression checks unchanged source/tangent call counts, agreement
with original island-only derivatives, incomplete and grid-inconsistent port
refusal, adjoint contraction and cancellation recovery. These remain local
finite-order derivatives, not stability or harmonic-tail convergence evidence.

### Sampled direct parameter partials of cycle-mean objectives

`PreparedForcedAdjoint::gradient_mean_split_direction` samples a complete
explicit integrand parameter partial on both original cached AFT trajectories
and composes it with the existing local-island and linear frequency-port
partials. The island-only entry holds the prescribed clock/linear port fixed;
`PreparedAutonomousAdjoint::gradient_mean_island_direction` retains the original
bordered dual and its solved-frequency chain. The scalar callback holds the
original Fourier coefficients and phase fixed, includes authored direct
parameter/time-chart dependence, and does not add the state chain again.
Autonomous direct partials also hold the unknown frequency fixed.

The shared compensated accumulation forms each physical scalar mean without a
raw sample sum. Independent explicit-partial agreement is normalized by the
declared parameter increment and the original objective units. Each grid's own
mean participates in `q_p - lambda^T R_p`; both the full residual-column change
and the complete contracted-gradient change retain their original independent
admission. Separately small objective and residual changes therefore cannot
hide an excessive combined gradient change. Every scalar is reset to NaN;
incomplete, nonfinite, refused or cancelled callbacks publish no gradient.

This adds constant scalar scratch and two cached sample passes within the
existing numerical query allowance. It introduces no baseline source/tangent
calls, waveform solve, factorization or synthetic state-gradient callback.
Source callback work/storage retain their existing exclusions, and two-grid
agreement does not certify the omitted harmonic tail. Four focused regressions
cover independently scaled physical derivatives, direct and combined aliasing,
refusal/cancellation recovery, representable extreme DC means, and an autonomous
cycle-energy direction with an independently known changing-period response.
