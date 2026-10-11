# Elastic drum barrel

`--elastic-barrel wall.fsb` adds a vibrating cylindrical sidewall to every
`drum`, `drum-modal`, `drum-stretch`, `snare` and `snare-off` command, including
their `-wav` and `-mic` forms when the declared windows fit the audio band.
The wall exchanges energy with both heads through the existing enclosed-air
owner and contributes its actual outer-surface velocity to exterior acoustics.

The drum specification supplies the clear radius, outer radius, depth and
azimuthal divisions. Their radius difference is the wall thickness; their
mean is its midsurface radius. The barrel card supplies material, intrinsic
loss, axial divisions and a retained frequency window. `fs-plate` assembles
the shell stiffness and mass, and its eigenproblem supplies the frequencies
and mass-normalized shapes. No resonance frequencies are entered as material
parameters.

Both end rings are clamped to rigid hoops. The heads retain their fixed rims,
and the existing sticks and wire bank act on those heads. This represents
air-mediated head/wall interaction. Moving hoops, individual lugs, rimshots,
ply anisotropy and direct stick/barrel impacts require additional mechanics.

## Use

The tracked card is an **illustrative isotropic declaration**, with estimated
material and loss. Its wider structural window is intended for mechanics CSV:

```sh
cargo run --release -p fs-couple --example percussion -- \
  drum 4096 --analytic-newton --impact-substeps 8 511 \
  --elastic-barrel crates/fs-couple/examples/percussion/illustrative-barrel.fsb \
  --cavity-modes --cavity-drag-per-s 20 \
  --strike-position-m 0.06 0.01 > elastic-barrel.csv
```

Omit `--cavity-modes` to use the shared compact compression spring. Use
`snare` to retain the complete declared wire bank, or add a supplied drum
specification, moving head mute, flexible shaft or head material memory.
The barrel can use either the prepared modal owner or the joint nonlinear
owner; it does not require adding head stretching to choose its physics.

For pressure export, supply a wall card whose window lies inside 40..1640 Hz,
and whose combined head/wall modes fit the existing source budget. The
illustrative 6 kHz card intentionally does not satisfy that audio admission.
A hoop-clamped barrel may have no modes in a requested narrow window; that
refuses instead of changing its material or inventing a low resonance. A
retained oval or torsional mode may correctly have negligible uniform-volume
coupling. Distributed pressure can excite additional spatial participation.

## SI format

```text
frankensim-drum-barrel-v1
material,10000000000,0.30,700,0.005
axial_intervals,4
band_hz,40,6000
```

All three records are mandatory, with comments beginning at `#`.
`material` gives Young modulus [Pa], Poisson ratio, density [kg/m³], and a
nonnegative modal damping ratio. Each wall mode receives the declared
viscous coefficient `2 * damping_ratio * angular_frequency`. Zero means
zero intrinsic wall loss. It leaves head, contact, felt, wire, shaft and
acoustic losses independent.

The bounded preparation accepts 2..32 axial divisions, at most 512 nodes and
1,024 facets, and a complete window of 1..63 wall modes. Combined mechanical,
acoustic source and panel budgets also apply; scalar input limits do not
guarantee every combination fits. No returned mode is silently dropped.
Files are limited to 16 KiB, and every window must meet the mechanical
Nyquist guard. A sidewall neck is refused with this closed barrel because a
cutout and moving-aperture coupling are not represented.

## Pressure work and radiation

At each inner facet, three positive-area quadrature points project the
finite-thickness shell displacement, including rotations, onto the outward
normal of the enclosed gas. The compact spring receives the **negative**
integral of that row: outward wall motion expands the cavity and reduces
compression. Distributed air receives the same outward rows tested against
its pressure shapes, with reciprocal reactions in the same mechanical step.
Its uniform member replaces the compact spring, so bulk compliance is
included once.

Only the outer skin enters the closed exterior boundary. The rigid bearing
annuli stitch to its actual end-ring vertices; no inner faces are inserted
into exterior radiation. Head and wall source rows retain their actual
mechanical addresses. The normal velocity projection includes shell
rotations and contains zero participation from strikers and wire coordinates.
Existing acoustic resolution, passivity and observer-fit admission still
apply. Successful structural preparation alone does not establish a
successful pressure render.

The cylinder and its offset skins use finite planar facets. The supplied
radii describe the intended cylinder; averaged shell directors and chordal
facets introduce a geometry approximation, including small end-ring radial
offsets. The acoustic seam uses the actual offset vertices. Inner pressure
sampling must remain inside the declared cylindrical cavity. Refinement is
required to assess geometric and modal convergence for a chosen specimen.

Focused Rust regressions cover skin work, cavity excitation, shared state
addresses and exterior participation. They require native execution; source
review and independent geometry calculations do not establish solver
convergence, measured-instrument fidelity or real-time performance.
