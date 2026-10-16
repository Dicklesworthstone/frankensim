# Bridge-informed soundboard reduction

An explicit preparation path for flat and crowned boards can project a larger
certified soundboard modal slice into the existing bounded piano engine. It addresses the
retention obstacle encountered when a useful source band contains more than
128 modes. It does not promote a new default band or establish convergence of
the underlying board mesh, material data, bridge mobility, or microphone sound.

## Selecting a model

Both piano frontends accept:

```text
--board-reduction max_modes,keep_low_modes,frequency_hz,...
--board-reduction adaptive,max_modes,keep_low_modes,tolerance,seed_hz,...
```

- `max_modes`: maximum retained board coordinates, from 1 through 128.
- `keep_low_modes`: number of original lowest modes to preserve exactly, from
  zero through `max_modes`, also no greater than the actual source count.
- Frequencies: 1 through 16 finite, positive, strictly increasing values in Hz.
  Static response is included automatically. Each frequency must lie within the
  explicitly selected source board band.
- `adaptive` adds a requested normalized complex bridge-transfer tolerance,
  finite and strictly between zero and one. The frequency list
  supplies initial response directions; the actual check also includes static
  response, every distinct source resonance and the intervening midpoints.

For example, this requests a 2600 Hz source slice and at most 96 playback
coordinates, preserving its first 24 modes:

```bash
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --equilibrate-board-mass \
  --board-band-hz 2600 --board-reduction 96,24,800,1600,2400 \
  --note 84 --velocity 1.5 --duration 3 --render piano-ritz.wav
```

This is an illustrative command, not a validated Model D preset or a
convergence result. The supplied geometry must pass the ordinary assembly and
modal certification checks. The source band is never enlarged implicitly.

For `piano_exterior`, the source ceiling comes from `board-band-hz` in the
supplied FSPE specification. The reduction flag is accepted after the ordinary
positional arguments by `response`, `admittance`, `render`, and
`render-loaded`. It composes with the supported board inertia, source string
damping, supplied polarization, hammer, pedal, and performance controls.

Without acoustic reduction, `render-loaded` accepts at most 32 board coordinates
for its passive radiation realization. Supplying `--radiation-ports 1..32`
separately selects a bounded radiating subspace while preserving all retained
board coordinates, up to 128. The projected passive load must meet the original
complete BEM matrix and resistance error limits; an insufficient acoustic rank
refuses. Every receiver still uses all structural inputs. See
`RADIATION_FEEDBACK.md` for the geometry-derived basis and full-load validation.
BEM panel, frequency-grid, fitting, conditioning, and pressure
headroom checks still apply. `response` remains pressure per supplied modal
acceleration; it is not a bridge-force response. Use `admittance` for the
force-driven string/board/radiation system.

The source solve admits at most 512 modes on this explicit path. The ordinary
unreduced path still admits a complete slice of at most 128. The string partial
budget controlled by `--modes` is independent. No path silently discards the
highest source eigenpairs to fit a budget.

Supplied crowned shells are reduced after the ordinary shell preparation. If
downbearing is supplied, the existing static solve first establishes the loaded
equilibrium and tangent stiffness; its certified tangent modes are the source
slice. Reduction preserves that equilibrium geometry and all six nodal shell
degrees of freedom. Flat-board mass and acoustic-refinement options remain
restricted to flat geometry.

In `grand_piano`, reduction requires a geometric render and excludes `--dump-board`:
the existing modal CSV cannot represent the full projected damping operator.
Geometry-only exports do not perform this reduction. A reduced modal table must
not be re-imported as though per-mode damping ratios were its complete material
law.

## Adaptive bridge-response accuracy

Fixed snapshots can describe their selected displacement directions very
accurately while missing resonances between the requested frequencies. In a
three-mode source with resonances at 100, 1000 and 1100 Hz and 1% modal damping,
a two-coordinate reduction protecting the 100 Hz mode and sampling at 100 Hz
reported a maximum snapshot residual of about 0.0011%. Nevertheless, the
actual bridge compliance errors at 1000 and 1100 Hz were approximately 95%,
with amplitudes reduced by 7.4 and 8.2 dB. This is a numerical counterexample,
not a measured piano spectrum.

The explicit adaptive form checks the actual Galerkin transfer instead:

```text
--board-reduction adaptive,96,24,0.05,800,1600,2400
```

Here 0.05 requests at most 5% normalized **complex** transfer error at
every evaluated frequency and every pair of admitted bridge force ports.
It accounts for amplitude and phase together, including transfer between
different keys and supplied string polarizations. The 96-coordinate budget, 24 exact low modes
and caller's source band remain explicit. This example is a preparation
request, not evidence that a particular 2600 Hz geometry meets it.

For force columns `b_i` and `b_j`, the source and reduced compliance are

```text
H_ij(omega)   = b_i^T [diag(lambda) - omega^2 I - i omega diag(c)]^-1 b_j
H_r,ij(omega) = (Q^T b_i)^T [K_r - omega^2 I - i omega C_r]^-1 (Q^T b_j)
error_ij      = |H_r,ij - H_ij| / (sqrt(|H_ii|) sqrt(|H_jj|))
```

The real-valued transpose here preserves the physical force/velocity port;
both evaluations use the same supplied source damping. The reduced solve uses
the **complete** projected stiffness and damping, including off-diagonal
terms. One factorization from the existing `fs_la::eigen_complex` owner is
shared by all bridge ports at each frequency.

On the diagonal, this is the original relative driving-point error. Away
from it, the source driving points provide a scale even when a cross transfer
is zero or nearly cancels. No arbitrary pressure, force or damping floor is
introduced, and independently rescaling or reversing a port leaves the
normalized error unchanged. A silent zero port must stay exactly silent.
All reciprocal pairs are checked with their signed complex response; matching
only the response at the struck key can miss motion transmitted to another
key, which matters for sympathetic excitation and the radiating board field.

A four-mode regression makes this distinction concrete: its two-coordinate
model has 5% driving-point error throughout the checked grid, satisfying a
6% request, while its normalized cross-bridge error exceeds 200% near a
158.114 Hz antiresonance. The complete transfer check refuses that reduction;
four coordinates recover both diagonal and cross responses below 1e-11 error.
Reversing and independently rescaling ports preserves the decision. These are
analytic source-model tests, not measured piano transfer data. The focused
native Ritz closure passed all 11 tests with the pinned compiler and the actual
`fs-math` and `fs-la` complex-LU implementations; a full example build was not run.

Preparation starts with the protected low modes and seed responses, evaluates
the full finite grid, and enriches the basis with real and imaginary source
responses at the worst frequency/drive-port pairs. Multiple failed receiver
ports for the same drive and frequency use one enrichment candidate. Every
successful final sweep checks the full reciprocal transfer matrix for all
admitted primary bridges and any supplied secondary directions. Cross-transfer
checks reuse the existing factorization and solves, without retaining an
all-frequency matrix of port pairs. Their source and projection work is charged
against the same explicit preparation budget.
It never drops a difficult frequency or bridge from the check to obtain
acceptance. Insufficient independent directions, retained rank or bounded
preparation work refuse the requested reduction before playback. No damping
floor or relaxed tolerance is inserted.

This is a comparison with the supplied **finite source model**, at a reported
finite set of frequencies. It does not certify response between those points,
a string-loaded resonance, receiver
pressure, modes outside the source band, or agreement with a real Model D.
The same final nodal field still supplies all bridges and acoustic outputs;
this option improves structural response selection without inventing a
separate receiver spectrum.

## Physical construction

In the admitted source-modal coordinates, mass is the identity, stiffness is
diagonal with entries `lambda_i = omega_i^2`, and viscous wood damping has
entries `c_i = 2*zeta_i*omega_i`. The complete original FE slice is checked
before any projection: count, positivity and certificate fields, finite
eigenvectors, mass normalization, and mutual mass orthogonality. An invalid
source pair cannot disappear by being omitted from the reduced result.

For every primary bridge in the admitted key set, and every supplied secondary
string direction, preparation constructs:

1. The static displacement response to its bridge-force vector.
2. Real and imaginary displacement responses at each selected harmonic
   frequency, using the supplied material damping and the
   `exp(-i*omega*t)` convention.

When `--string-polarization` is supplied, its complete card is admitted before
the source eigensolve. Its physical site, arm and two transverse directions
are projected directly from each source mode's three structural facet nodes.
The supplied primary projection must reproduce the original geometric bridge
coefficient in **every source mode**, before any direction can be discarded by
reduction. A wrong frame cannot become valid by hiding its discrepancy in a
dropped mode. Source-force snapshots and the final motion surface share the
same P1 or edge-cubic projection law, including physical rotations. There is no
512-mode runtime motion object or inferred lateral coupling.

Each nonzero real-valued snapshot is normalized in the complete source-modal
space **before** the protected low coordinates are removed. This prevents a
tiny residual tail from receiving the same weight as an entire physical
response. Protected modes enter as exact canonical vectors. The remaining
space is chosen by deterministic largest-residual pivots and two-pass modified
Gram-Schmidt, up to the requested budget and actual numerical rank. Ties follow
source snapshot order. Spare capacity is not filled with arbitrary directions.

For the resulting source-modal basis `Q`, the physical matrices are:

```text
M_r = I
K_r = Q^T diag(lambda) Q
C_r = Q^T diag(c) Q
```

Only the unprotected tail of `K_r` is diagonalized, using the existing
`fs_modal::eigh_gen_dense` owner. Protected eigenvectors, signs, eigenvalues,
frequency intervals, bridge projections, and surface shapes remain unchanged.
The final columns, including that tail rotation, are transformed into nodal FE
coordinates once. The same nodal field is then used for all bridge projections,
volume velocity, Rayleigh surface quadrature, and full-vector motion supplied
to polarization and finite-body acoustics.

The selected P1 or edge-cubic field is used consistently, including shell
translation, rotation and eccentric bridge-arm motion. Neither reduction
nor the acoustic receiver introduces a separate displacement interpolation.
Material geometry, mass, stiffness, supports, and bridge locations are supplied
by the original board.

### Full material damping

Damping is generally dense after response-space reduction, even if the source
model uses a common modal damping ratio. Retaining only its diagonal would
change the physical velocity forces and dissipated power.

`PreparedBoard::physical_damping` therefore carries the complete row-major
matrix in the returned mass-normalized board basis. Both frontends install it
through `Instrument::configure_bare_board_damping` before excitation.
Harmonic preparation installs the same matrix through
`BridgeResponse::configure_bare_board_damping`. Each bank projects it into its
actual string-loaded basis. The existing conservative dynamics and contact
compliance remain the owners of mechanical work; the existing dissipative
half-flow applies the full positive-semidefinite operator.

The per-mode `damping_ratio` fields describe diagonal entries for compatibility
with ordinary board construction. They are insufficient to reconstruct the
reduced material law. The full matrix **replaces** that diagonal loss, rather
than adding another damping channel. Lossless harmonic comparison still
validates the supplied matrix and disables material loss explicitly; it uses
the same geometrically prepared reduction space as the damped comparison.

## Reports and scope of accuracy

Preparation retains separate source and reduced metadata:

- Every original FE frequency interval remains in
  `ReductionReport::source_frequency_intervals_hz`.
- The protected prefix retains those original FE intervals.
- Mixed-tail intervals certify the **projected pencil** only. They are not
  certificates for individual eigenvalues of the unreduced FE model.
- Reports include source count, actual retained count, protected count,
  requested harmonic frequencies, and nonzero snapshot count.
- Adaptive reports additionally preserve the requested tolerance, achieved
  maximum complex driving-point error, complete evaluation frequency list,
  force-port count and enrichment-round count. Both frontends print the error
  and its scope; exterior CSV comments retain the exact frequency grid.
- The reported maximum snapshot residual is
  `||x - Q Q^T x|| / ||x||` over the normalized nonzero static, harmonic-real,
  and harmonic-imaginary displacement snapshots.

The frontends print this scope explicitly, including source and retained
frequency intervals in reduced exterior CSV comments. A small snapshot
projection residual does not bound the actual reduced transfer error away from
the sampled responses, nor the error in a loaded string/board resonance or
receiver pressure. The target set covers primary bridge forces and the actual
secondary directions of a supplied polarization card. This retains response
directions that a primary-only force may not excite. Receiver pressure is
evaluated from the resulting shared geometry; receivers are not independent
optimization targets.

An exact undamped source pole at a requested frequency refuses without adding
a damping floor. Source modes, ports, harmonic samples, retained rank, and
estimated scalar work are bounded before allocating the snapshot bank. The
leaf permits at most 176 force vectors (two directions for each of 88 courses)
and 250 million estimated source-coordinate work visits; this is an admission
budget for fixed selection, not a runtime claim. Adaptive preparation has a
separate budget of eight billion estimated complex scalar-work visits. Source
response preparation, projection, complete response sweeps and enrichment are
charged before execution. Larger source slices and ranks can exhaust it;
refusal identifies the bounded preparation work instead of returning an
unchecked basis. No real-time or wall-clock performance claim follows from
these budgets. Omitting the polarization card preserves primary-only selection.

Mesh refinement, source-band extension, retained-rank sweeps, and independent
bridge/microphone comparisons remain necessary before promoting a new default
or claiming high-band piano fidelity.

## Regression coverage and execution status

The new leaf regressions check projected static/harmonic equations, exact low
modes, mass/stiffness/damping energy forms, off-diagonal damping, deterministic
rank selection, explicit budgets, and invalid physical input. Geometry
regression compares full-source and reduced sampled bridge responses and
checks that bridge, volume, surface, and full-vector motion use one basis.
Frontend regressions cover admission and propagation into actual played and
harmonic preparation.

Native checks with the repository's pinned `nightly-2026-07-06` compiler pass:
five standalone Ritz-kernel tests, 17 focused `grand_piano` tests and 17 focused
`piano_exterior` tests. The example sets include shared regressions. They cover
loaded crown equilibrium, exact protected modes, six-DOF motion, full physical
damping and work, source-primary validation before selection, lateral-only
response retention, P1/cubic preparation, harmonic/BEM consistency, and played
stereo block invariance with felt, pedals and string stretching.

Both example test binaries were compiled from the repository's exact source
paths using a small external Cargo manifest and the pinned dependency sources.
The repository's ordinary `--locked` build refused because its Cargo lock does
not match the pinned constellation. No repository lock, dependency pin or
solver admission was changed. These focused native results are not a full
workspace/DSR acceptance, high-band convergence or performance result.

### Adaptive response extension

The adaptive extension's native exact-source owner checks pass nine tests
(five existing fixed-selection tests and four adaptive regressions). They
include the missed-resonance counterexample, multiport enrichment, preserved
energy forms and protected modes, deterministic selection, and explicit
refusal for exhausted budgets, undamped poles and nonrepresentable error
calculations.

A separate final-projection regression uses an eight-mode analytical physical
pencil with unequal nodal masses and nonzero force participation in every mode.
It reduces to seven coordinates, protects the lowest eigenpair, and compares
the source with the returned nodal modes, rotated stiffness and full material
damping at all 18 evaluation frequencies. Maximum relative complex
driving-point error is approximately 0.00181843, below the requested 0.002.
The final damping matrix has nonzero off-diagonal entries; the test checks its
physical congruence and agreement with the pre-rotation response report.

A native 40-source-mode, 88-distinct-port probe also passes, retaining 28
coordinates and checking 81 frequencies. Its last 12 source modes are
intentionally unforced; this exercises all-key admission and complete-grid
evaluation without asserting that a real soundboard has that low rank.

These extension checks use the pinned compiler and existing numerical owners.
The full piano frontends, BEM sweep, complete rendering and listening accuracy
were not validated by this extension run. The frontend admission tests have
been extended for the new syntax; the adaptive option does not promote a new
default board band or a calibrated Steinway geometry.
