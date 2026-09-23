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
board files, full-vector skin mapping, rigid lid/cabinet components and one/two
receivers use the existing format in `EXTERIOR_ACOUSTICS.md`. Every acoustic
part must still be explicitly mapped. No missing geometry or material is filled
in, and visual MTL properties are not interpreted as elastic constants.

## What reaches the result

All admitted speaking strings, unison members and duplex spans remain present,
including unplayed courses. The bank's moving-endpoint inertia completion,
reciprocal cross-potential, actual tension/flexural partials and complete retained
board basis are reused. At most 24 partials per string are retained by this
wrapper, with the existing 21.6 kHz ceiling; the structural slice is the complete
admitted slice through `board-band-hz`. No mode is silently dropped to make a
failed solve pass. High duplex partials outside the retained band are reported;
the existing endpoint stiffness/mass contributions remain.

The frequency-domain image has **no hammer or key-damper contact**. It is the
linear continuous operator underlying the moving-boundary bank, with the same
authored string bending-loss spectrum and physical wood damping matrix. It is
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
coordinates are eliminated into a board-sized Schur complement, then recovered
for a separate full-equation residual and power balance. Exact unresolved
lossless string poles refuse rather than adding artificial damping or nudging
the requested frequency. The numerical solve is the existing complex LU owner.

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
