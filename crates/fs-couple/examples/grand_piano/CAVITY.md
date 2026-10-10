# Reciprocal cavity air in the played piano

`--cavity enclosure.fspc` couples a supplied sealed rectangular air volume to
the actual geometric soundboard. Its pressure changes the board motion and
therefore the strings, felt contact, bridge response and receiver pressure.
The piano keeps one mechanics clock. The cavity is prepared before excitation
and its history survives notes, pedal events and audio block boundaries.

This is an explicit enclosure model. A grand piano preset does not select it
automatically: a rigid, sealed rectangular cavity is not a reconstruction of
an open grand piano's lid, rim and openings. The existing exterior body and
receiver descriptions remain independently supplied physical inputs.

## Input

Every card requires exactly one of each row below. Coordinates and dimensions
are in metres, temperature is in kelvin, and absolute pressure is in pascals.
This example describes an **authored** one metre square enclosure:

```text
frankensim-piano-cavity-si-v1
source,estimated,authored one metre sealed enclosure
interface-origin-m,0,0,0
dimensions-m,1,1,0.3
modes,4
damping-ratio,0.02
gas,dry-air-ussa1976,293.15,101325
```

The source row requires an authority label (`estimated`, `mixed`, `published`
or `measured`) and an attribution. Importing a label does not independently
verify that attribution. The interface origin is the top face's lower x/y corner. The cavity occupies
`[x,x+Lx] × [y,y+Ly] × [z-Lz,z]`. Positive board z displacement is **away from
the gas** and increases its volume. The entire supplied soundboard must lie
in the interface plane and within that rectangle. It may cover only part of
the face; the rest is a rigid wall. A nonplanar crowned board refuses this
rectangular chart instead of being flattened or replaced.

`modes` selects the lowest 1–8 rigid-wall pressure modes, including the uniform
zero-frequency compression mode. Dimensions and the admitted gas sound speed
determine their frequencies; no resonances are entered or fitted by hand.
`damping-ratio` supplies momentum drag `d_j = 2 ζ ω_j` for each nonzero acoustic
mode. The uniform mode has no momentum coordinate and receives no drag. This
is a causal acoustic loss selection, not hysteretic loss inferred from a wall
material or a reverberation-time measurement.

Missing files, duplicate or invalid rows, unsupported geometry and exceeded
work/rate budgets fail before rendering. A rejected card cannot select an
uncoupled fallback.

## Rendering

For a supplied planar board and string scale:

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --board-geometry panel.fsb --scale strings.csv \
  --cavity enclosure.fspc --render cavity-piano.wav \
  --note 69 --velocity 1 --duration 2
```

The same option is accepted by `piano_exterior render` and
`piano_exterior render-loaded`, after their usual positional arguments. The
first computes exterior pressure from cavity-loaded motion. The second also
retains the existing passive exterior radiation reaction. Mono and stereo
receivers observe the same instrument trajectory. Exterior pressure is still
emitted by the supplied moving solid boundary; this sealed-cavity option adds
no fictitious aperture sound source.

For these exterior commands, supply a closed **outer enclosure** OBJ whose
exposed soundboard face is marked moving and whose remaining walls and bottom
are marked rigid in the acoustic card. `board-skin` and
`board-skin-continuous` alone expose both board faces and therefore refuse with
a sealed cavity. A supplied moving underside also refuses, since it would
apply exterior fluid to a face already inside the cavity. The input author
must make the cavity volume and outer enclosure describe the same object;
the two independent geometric descriptions are not automatically reconciled.

The cavity projection uses the selected P1 or edge-cubic soundboard field.
It also works after soundboard Ritz reduction: pressure overlaps are computed
in the retained geometric basis and then transformed by the string-loaded
bank's actual board basis. Full physical board damping remains intact. No
additional independent board solve or diagonal damping approximation is used.

## Force-driven harmonic response

The same card is accepted by `piano_exterior admittance`:

```sh
cargo run --release -p fs-couple --example piano_exterior -- \
  admittance panel.fsb strings.csv enclosure.obj acoustics.fspe 69 mobility.csv \
  --cavity enclosure.fspc
```

It changes the physical bridge mobility and receiver pressure per unit bridge
force. The frequency-domain model uses the same cavity springs, acoustic
inertia, momentum drag and loaded board projection as playback. Its existing
pivoted solve retains every acoustic inertia alongside string coordinates
that need explicit treatment near their own resonances. This avoids dividing
by zero at an undamped fixed-wall cavity frequency; coincident cavity/string
resonances receive the same complete coupled solve. A singular full physical
system still refuses instead of moving a pole or adding artificial damping.

Selected sweeps append a `cavity_w` column to the existing CSV. It reports
cycle-average acoustic momentum dissipation for the applied 1 N peak bridge
force. The reported work defect includes wood, string, exterior radiation and
cavity losses once each. Sweeps without a cavity retain the original column
layout. The `one_way` comparison omits exterior radiation reaction while
retaining the selected interior cavity; it does not silently change the
enclosure between the two comparisons.

`--lossless-structure` removes wood and intrinsic string losses, leaving the
card's explicit cavity drag and exterior radiation intact. Set the card's
`damping-ratio` to zero for a lossless cavity. The `response` command prescribes
modal acceleration and has no mechanical reaction solve, so it rejects
`--cavity`; use `admittance` or a render to observe cavity loading.

## Stored energy, reaction and numerical scope

For each pressure mode, the existing cavity owner supplies frequency `ω_j`,
volume norm `Λ_j` and the interface overlap
`C_rj = integral(phi_r psi_j dA)`. Define `A_j = ρ c² / Λ_j` and pressure-energy
coordinate `y_j = sqrt(A_j) C_j q + ω_j z_j`. The acoustic storage is
`sum(y_j² + p_j²)/2`, omitting the nonexistent momentum of a uniform mode.
Pressure coefficient `-sqrt(A_j) y_j` exerts the reciprocal board force
`-C_j sqrt(A_j) y_j`. Uniform compression therefore retains its physical gas
spring; it is not discarded because its modal frequency is zero.

The implementation composes the existing collective power-port exchange and
exact free-oscillator owners. It adds no piano-specific acoustic integrator.
Each isolated exchange conserves combined quadratic energy to its checked
roundoff tolerance. The symmetric composition with mechanics and exterior
radiation is second order and needs mechanical-rate convergence for accuracy;
it is not an exact full-instrument propagator or a real-time certificate.

Cavity momentum loss and exterior radiation loss are recorded separately and
included once in the combined work balance. A refused output frame restores
the cavity history together with strings, board, hammer, felt and radiation
states. Loss is not clipped to repair a work defect.

Interface integration uses bounded quadrature on the actual displacement
field and checks refinement. Its reported discrepancy is an estimate, not a
spatial convergence certificate. The existing 128-coordinate played-board
limit and independent exterior radiation budgets still apply. This feature
does not establish high-frequency piano fidelity or close the outstanding
high-band soundboard/radiation validation gap.
