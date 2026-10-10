# Radiation-loaded piano bridge response

`piano_exterior admittance` applies a **unit peak harmonic force at one actual
bridge station** and solves the coupled retained string/board response with the
exterior pressure reacting onto the mechanics. This is not the pressure-only
`response` command, and it does not change the one-way `render` command.

```sh
cargo run --release -p fs-couple --example piano_exterior -- \
  admittance settled.fss strings.csv acoustic-body.obj acoustic.fspe 69 bridge.csv
```

`steinway-d` may replace `strings.csv` to use all 88 raw source courses. Use the
same supplied tensions that produced any settled board. Flat or crowned/loaded
board files, full-vector skin mapping, rigid lid/cabinet components and up to 64
receivers use the existing format in `EXTERIOR_ACOUSTICS.md`. Every acoustic
part must still be explicitly mapped. No missing geometry or material is filled
in, and visual MTL properties are not interpreted as elastic constants.

## What reaches the result

All admitted speaking strings, unison members and duplex spans remain present,
including unplayed courses. The bank's moving-endpoint inertia completion,
reciprocal cross-potential, actual tension/flexural partials and complete retained
board basis are reused. The default is at most 24 partials per string;
`--modes 1..512` exposes the existing retention budget, with the same 21.6 kHz
ceiling. The structural slice is the complete
admitted slice through `board-band-hz`. No mode is silently dropped to make a
failed solve pass. High duplex partials outside the retained band are reported;
the existing endpoint stiffness/mass contributions remain.

The frequency-domain image has **no hammer or key-damper contact**. It is the
linear continuous operator underlying the moving-boundary bank, with the same
selected string loss law and physical wood damping matrix. It is
not the exact transfer of the time stepper at finite step size. Modal wood
losses and the current string spectrum are not measured Steinway calibrations.

At one angular frequency, the existing Helmholtz solver solves one batch with
a unit generalized **velocity** for each retained board coordinate. The same
solutions supply both receiver pressures and the complete modal radiation
impedance `Z = G^T A P`: `G` contains the normal-motion rows, `A` the physical
panel areas, and each column of `P` is a solved pressure field. Areas enter the
force projection once. Off-diagonal coupling and the imaginary/reactive load
are retained; no fitted damping constant or matrix symmetrization replaces them.

Under `exp(-i omega t)`, velocity is `-i omega q`. The opposing fluid load adds
`-i omega Z` to the mechanical dynamic stiffness. The many independent string
coordinates are normally eliminated into a board-sized Schur complement, then
recovered for a separate full-equation residual and power balance. Coordinates
near fixed-interface string poles instead remain in a bounded coupled border.
They are solved with the board by the existing pivoted complex LU owner, never
by dividing by an unresolved string diagonal. A fixed-interface pole is not
necessarily a singularity of the complete coupled piano: it can be a bridge
antiresonance with finite string motion.

The border retains the full cross-potential, self-stiffness, material damping
and supplied acoustic matrix. It adds neither artificial damping nor a pole
shift, pseudoinverse, modal deletion or frequency interpolation. At most 128
near-pole string coordinates may join the original board solve; exceeding this
setup budget refuses. Singular LU, nonfinite results and failed full-equation
or power checks still refuse. This does not supply an arbitrary solution at a
genuine unresolved coupled-system resonance.

## Analyze the same string directions and losses as playback

The harmonic `admittance` and `response` fronts accept the physical
`--string-polarization` and `--rt0425-string-damping` controls already used by
played output. For example:

```sh
cargo run --release -p fs-couple --example piano_exterior -- \
  admittance settled.fss steinway-d acoustic-body.obj acoustic.fspe 69 bridge.csv \
  --modes 128 --string-polarization bridge-frames.fspp --rt0425-string-damping
```

The polarization file uses the same full-vector geometric motion and complete
per-key bridge-frame projection described in `STRING_POLARIZATION.md`. Both
transverse directions of every admitted speaking string, unison member and
duplex span remain in the harmonic bank, including all unplayed keys. Their
physical endpoint inertia changes the loaded board basis before the acoustic
boundary is projected. Their signed bridge forces then contribute to the
same coupled equation. An unstruck lateral speaking string is not treated as
a duplex merely because it has no hammer contact.

The existing unit-force experiment retains its meaning: drive and `bridge`
CSV rows describe force and velocity in the primary hammer-plane direction
at the supplied key's bridge station. Lateral strings react through the shared
board; they do not add a second applied force, double the hammer, or redefine
the reported velocity as a sum over directions. A zero lateral bridge
projection preserves the original primary experiment while retaining the
uncoupled extra string coordinates.

`--rt0425-string-damping` selects the same published per-key `R_u` and `eta_u`
projection as playback, for both directions and all retained spans. It requires
the `steinway-d` source scale. The estimated common loss remains the default.
The conservative mechanics, physical string frequencies and board damping law
are unchanged by choosing the source string loss. This remains the existing
reduced transverse-string model, not the complete higher-order model from the
report or a measured-instrument calibration.

Use matching string controls, retained partial count and structural inputs
when comparing harmonic admittance with a played note. The pressure-only
`response` experiment still applies prescribed generalized board acceleration;
its columns are not bridge-force mobilities. Hammer/felt contact, key dampers
and finite-amplitude geometric string extension remain outside these linear
harmonic experiments.

The selected harmonic constructor builds the existing bank directly and reads
its actual course/member/polarization/duplex order. Regression tests compare
the recovered solution with the complete stiffness Hessian of that bank's
stored energy, compare source-loss power with its time generator, and retain
both coordinates at coincident polarization poles. These source tests do not
constitute native execution or physical validation evidence by themselves.

## Explicit conservative-structure comparison

`admittance` alone accepts `--lossless-structure`:

```sh
cargo run --release -p fs-couple --example piano_exterior -- \
  admittance settled.fss strings.csv acoustic-body.obj acoustic.fspe 69 bridge.csv \
  --modes 128 --lossless-structure
```

This deliberately sets the original wood and string material damping to zero.
The same geometry, installed tension, mass, string partials and bridge
couplings remain, and the **complete complex acoustic load stays active**.
It is a controlled lossless-structure comparison, not a claim about real wood
or string losses. The CSV header records the selection. Wood and string power
columns are zero; admitted input power balances radiation. The one-way columns
still omit the fluid reaction, so genuine undamped coupled resonances can
refuse that comparison instead of being smoothed into a finite peak.

Without the option, material damping and the original well-conditioned Schur
arithmetic remain unchanged. `response`, `render` and `render-loaded` reject
the flag; no played piano silently loses its existing physical dissipation.
Numeric retention controls keep their previous meaning. Duplicate flags and
missing option values refuse before input/output access.
Combining `--lossless-structure` with `--rt0425-string-damping` also refuses:
one requests zero material loss and the other explicitly requests a loss law.
Polarization remains compatible with the conservative-structure experiment.

Four new response tests compare exact/near/coincident partials against the
original time bank's full energy Hessian, retain normal damped behavior, and
check the complete pole-budget refusal. Two command regressions cover strict
selection and actual finite-body BEM/bridge output through a physical partial.
The existing sympathetic-string regression now probes the silent course's
actual fundamental rather than an arbitrary off-resonance frequency; its
original effect-size and reciprocity thresholds remain unchanged.

## CSV and comparisons

Each frequency has a `bridge` row for every scale key and a `receiver` row for
each microphone. Bridge values are complex velocity/force **m/s/N**; receiver
values are complex **Pa/N**. All phasors are peak, not RMS. Columns `real,imag`
are the radiation-loaded result. `one_way_real,one_way_imag` use the **same air
transfer**, but the mechanics omit its reaction. This comparison isolates
radiation loading; the receiver reference is not sound propagating in vacuum.

Repeated diagnostic columns are cycle-average watts for the one-newton peak
experiment: input power, wood loss, string loss and radiated power. The reported
residual checks `input = wood + string + radiation`. Another column reports the
backward error of the full recovered string/board equation, not just the Schur
solve. The BEM wavelength-resolution and condition-lower-bound diagnostics are
also retained. Negative power beyond the stated numerical admission, failed
balance, nonfinite data or an unresolved solve reject the whole CSV. Outputs
must be fresh paths; an OS write error can still leave a partial new file.

The requested frequency grid is used directly, without vector fitting, a
synthesized impulse response, retuned strings or an audio equalizer. A uniform
grid may miss narrow resonances. These are sampled estimates, **not** spatial,
modal-truncation or frequency-interpolation convergence certificates. Acceptance
of one force experiment is not a global passivity certificate for an arbitrary
multiport matrix.

## Remaining physical boundary

This implements radiation feedback for **harmonic bridge-force analysis**, not
for played nonlinear time-domain hammer impacts. `piano_exterior render` still
uses one-way acoustics. The soundboard is linearized about the supplied static
equilibrium; cabinet/lid components are rigid; the dynamic strings remain the
existing transverse model. Flexible rims, pin friction, axial string dynamics,
room acoustics and measured-instrument validation are not supplied by this
command. No new external Steinway mesh or factory material dataset is bundled.

Why this observable matters: Kerem Ege and Antoine Chaigne, *End conditions of
piano strings* (2011), arXiv:1101.4511, treats the bridge input admittance as the
string termination. This implementation does not reproduce their specimen or
claim agreement with their measurements.

## Complete physical board damping

The engine `Instrument::configure_bare_board_damping` and the harmonic
`BridgeResponse::configure_bare_board_damping` accept a complete symmetric
positive-semidefinite viscous matrix before excitation. It is row-major in the
original mass-normalized bare-board coordinates, with units `1/s`, and replaces
the loss derived from the individual `BoardMode::damping_ratio` values.

The bank projects this matrix through its actual string-loaded board basis.
Playback uses the existing dissipative velocity half-flow; harmonic analysis
uses that same admitted continuous matrix in its impedance and power balance.
Off-diagonal entries are retained. This supports non-proportional wood loss
and physical damping projected from a larger modal basis.

Configuration is cold and once-only. Invalid shape, nonfinite or asymmetric
data, materially indefinite damping, or a running instrument refuse without
publishing a partial change. A globally lossless bank still validates the
physical matrix, then keeps its effective damping zero.
