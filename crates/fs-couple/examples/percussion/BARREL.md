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
and the wire bank acts on the resonant head. Sticks can strike the batter head
or the actual sidewall as described below. Moving hoops, individual lugs,
rimshots and layered laminate mechanics require additional operators. A
homogeneous direction-dependent material is supported as described below.

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

## Strike the sidewall

`--barrel-strike AZIMUTH_RAD HEIGHT_BELOW_BATTER_M` moves the first hand's
contact to the elastic barrel. It requires `--elastic-barrel`; the azimuth is
in `[-2*pi,2*pi]` radians and the height must lie strictly between the two
hoops. The first-hand `--strike-position-m` and a first-hand felt mallet are
refused with this selection. An optional second stick or second mallet still
plays its separately supplied batter-head position.

```sh
cargo run --release -p fs-couple --example percussion -- \
  drum 4096 --analytic-newton --impact-substeps 8 511 \
  --elastic-barrel crates/fs-couple/examples/percussion/illustrative-barrel.fsb \
  --barrel-strike 0.3 0.08 --strike-speed-m-s 1.2 \
  --cavity-modes --cavity-drag-per-s 20 > sidewall-strike.csv
```

The azimuth locates a point by intersecting a radial ray with the actual
outer-skin triangles. Strike motion then follows the **inward facet normal**,
which can differ from the analytical radial direction on a coarse mesh. A
shared-edge or vertex hit selects the lowest incident facet index and its
normal once. There is no nearest-node strike, duplicated edge impulse, or
unmodelled hoop response. The skin projection retains both shell translation
and the moment of a force applied at the finite thickness offset.

The first stick keeps its supplied speed, 0.2 mm initial exterior clearance,
inertia and force program. A flexible shaft keeps its complete basis and hand
force row, now interpreted along this local inward normal. With inward tip
motion `s` and outward skin motion `u`, contact penetration is `s + u`; the
existing contact potential supplies both reactions in the shared time step.
This contact replaces the first hand's head contact. Both heads respond through
the enclosed pressure and any independently selected second stick. Wall motion
continues to supply the existing outer acoustic velocity rows.

The hard tip retains the example's explicitly estimated Hertz contact law.
This is a normal point-contact approximation on the undeformed shell skin,
without tangential friction, finite-rotation tracking or an inferred finite
mallet face. A rigid hoop remains a boundary condition, not a radiator. These
commands illustrate the supported inputs; they are not completed render or
specimen-fidelity claims.

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

### Homogeneous orthotropic material

Version 2 supplies independent stiffness along and across the material grain,
plus its in-plane shear modulus and physical orientation:

```text
frankensim-drum-barrel-v2
orthotropic,12000000000,800000000,0.35,650000000,650,0.005,0
axial_intervals,4
band_hz,40,3000
```

The seven values after `orthotropic` are `E1` [Pa], `E2` [Pa], `nu12`, `G12`
[Pa], density [kg/m³], modal damping ratio, and grain angle [rad]. All are
required. The supplied constants describe one homogeneous plane-stress
material through the wall thickness. Reciprocity gives `nu21 = nu12*E2/E1`;
the reduced compliance must be positive definite. In particular, the
orthotropic Poisson ratio is checked through `1 - nu12*nu21 > 0`, rather
than replacing it with an isotropic Poisson bound. No missing modulus is
derived from an instrument name.

At zero angle, material axis 1 follows the barrel axis upward, toward the
batter head. Positive angles turn toward increasing cylinder azimuth; at
`pi/2`, axis 1 follows the hoop. Angles must lie in `[-pi,pi]`. The axis is
projected into each actual shell facet and then expressed in that facet's
local element frame. This matters because the two triangles of a barrel
strip have different first-edge directions: assigning one unrotated local
matrix everywhere would make the physical grain jump at mesh diagonals.

The existing `PlateSection::orthotropic_plane_stress_at_angle` owner supplies
the transformed constitutive matrix, including bending/twisting coupling.
`assemble_shell_sections` and `ShellReduction` receive the same section set.
The material therefore changes the wall eigenmodes and its response to head
pressure or a sidewall strike. Thickness, mass geometry, contact position,
cavity volume and external force programs retain their physical definitions.

The tracked `illustrative-orthotropic-barrel.fsb` is a synthetic research
declaration, not measured wood data. For example:

```sh
cargo run --release -p fs-couple --example percussion -- \
  drum 4096 --analytic-newton --impact-substeps 8 511 \
  --elastic-barrel crates/fs-couple/examples/percussion/illustrative-orthotropic-barrel.fsb \
  --barrel-strike 0.3 0.08 --strike-speed-m-s 1.2 \
  --cavity-modes > orthotropic-wall-strike.csv
```

Its 3 kHz structural window is for mechanics CSV and exceeds the current
40..1640 Hz wall-audio admission. The command has not been natively rendered.
The original version 1 isotropic files remain accepted with their original
material law. Version 1 accepts `material`; version 2 accepts `orthotropic`.
Mixing the two records or omitting a v2 constant refuses.

This image uses the homogeneous relation `A = 12*D/thickness²`. A stacked
ply laminate generally has independent membrane and bending tensors and may
also couple stretching to bending. Those tensors, inter-ply slip, glue layers,
grain variation through thickness and wood moisture dependence are not
inferred by version 2. Supplying an equivalent homogeneous material does not
establish a laminate or specimen-fidelity claim.

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
addresses and exterior participation. Sidewall-strike checks additionally
cover actual facet location and edge ties, the finite-thickness force moment,
reciprocal impact, pressure-driven head motion, flexible-shaft excitation and
failed/cancelled-step rollback. Orthotropic fixtures check principal-axis
energy in every actual facet frame, signed grain-angle conventions, the exact
isotropic limit, and changed modal frequencies and radial point compliance
with unchanged geometry. They require native execution; source
review and independent geometry calculations do not establish solver
convergence, measured-instrument fidelity or real-time performance.
