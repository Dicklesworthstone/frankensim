# Bridge-informed soundboard reduction

An explicit flat-board preparation path can now project a larger certified
soundboard modal slice into the existing bounded piano engine. It addresses the
retention obstacle encountered when a useful source band contains more than
128 modes. It does not promote a new default band or establish convergence of
the underlying board mesh, material data, bridge mobility, or microphone sound.

## Selecting a model

Both piano frontends accept:

```text
--board-reduction max_modes,keep_low_modes,frequency_hz,...
```

- `max_modes`: maximum retained board coordinates, from 1 through 128.
- `keep_low_modes`: number of original lowest modes to preserve exactly, from
  zero through `max_modes`, also no greater than the actual source count.
- Frequencies: 1 through 16 finite, positive, strictly increasing values in Hz.
  Static response is included automatically. Each frequency must lie within the
  explicitly selected source board band.

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

`render-loaded` retains its independent limit of 32 board coordinates for the
passive radiation realization. Use a reduction budget no greater than 32 for
that command. BEM panel, frequency-grid, fitting, conditioning, and pressure
headroom checks still apply. `response` remains pressure per supplied modal
acceleration; it is not a bridge-force response. Use `admittance` for the
force-driven string/board/radiation system.

The source solve admits at most 512 modes on this explicit path. The ordinary
unreduced path still admits a complete slice of at most 128. The string partial
budget controlled by `--modes` is independent. No path silently discards the
highest source eigenpairs to fit a budget.

Crowned shells are currently unsupported by this reduction. In
`grand_piano`, it requires a geometric render and excludes `--dump-board`:
the existing modal CSV cannot represent the full projected damping operator.
Geometry-only exports do not perform this reduction. A reduced modal table must
not be re-imported as though per-mode damping ratios were its complete material
law.

## Physical construction

In the admitted source-modal coordinates, mass is the identity, stiffness is
diagonal with entries `lambda_i = omega_i^2`, and viscous wood damping has
entries `c_i = 2*zeta_i*omega_i`. The complete original FE slice is checked
before any projection: count, positivity and certificate fields, finite
eigenvectors, mass normalization, and mutual mass orthogonality. An invalid
source pair cannot disappear by being omitted from the reduced result.

For every primary bridge in the admitted key set, preparation constructs:

1. The static displacement response to its bridge-force vector.
2. Real and imaginary displacement responses at each selected harmonic
   frequency, using the supplied material damping and the
   `exp(-i*omega*t)` convention.

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

The selected P1 or edge-cubic field is used consistently. Neither reduction
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
- The reported maximum residual is
  `||x - Q Q^T x|| / ||x||` over the normalized nonzero static, harmonic-real,
  and harmonic-imaginary displacement snapshots.

The frontends print this scope explicitly, including source and retained
frequency intervals in reduced exterior CSV comments. A small snapshot
projection residual does not bound the actual reduced transfer error away from
the sampled responses, nor the error in a loaded string/board resonance or
receiver pressure. The target set currently covers primary bridge forces.
Secondary string directions and receivers use the same reduced geometry but
are not independent optimization targets.

An exact undamped source pole at a requested frequency refuses without adding
a damping floor. Source modes, ports, harmonic samples, retained rank, and
estimated scalar work are bounded before allocating the snapshot bank. The
leaf permits at most 88 ports and 250 million estimated source-coordinate work
visits; this is an admission budget, not a runtime claim.

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

These regressions were added without an executed result in the implementation
session: the local execution service disconnected before syntax checking,
compilation, or native tests could run. Direct GitHub publication and exact
committed-source verification succeeded; those checks are not numerical or
performance validation.
