# Finite-body piano acoustics from supplied geometry

`piano_exterior response` connects the physical piano's existing modal motion
to the existing Helmholtz boundary-element solver. Unlike the ordinary Rayleigh
observer, this calculation has no infinite baffle: supply both sides and edges
of the soundboard and any rigid cabinet/lid surfaces that should scatter sound.
A lid is part of the same boundary solve, not a second sound source or an EQ.

```sh
cargo run --release -p fs-couple --example piano_exterior -- \
  response settled.fss strings.csv acoustic-body.obj acoustic.fspe response.csv
```

`steinway-d` may replace the scale CSV to select all 88 raw source Model D
courses. An explicit CSV preserves its supplied tensions; use the SAME scale
that produced a settled board. Flat and crowned/preloaded structural files are
supported. The existing source felt, shank and reciprocal string bank prepare
the loaded basis, including silent strings. No source course is discarded.

## Acoustic specification

Every singleton row below is mandatory; numbers illustrate the schema, not a
measured Steinway or a promise this bandwidth will pass on an arbitrary mesh.
Object/group labels must match the supplied OBJ exactly.

```text
frankensim-piano-exterior-si-v1
source,estimated,Explicit acoustic idealization; replace with actual attribution
obj-scale-m,1
obj-origin,0,0,0
max-skin-offset-m,0.015
medium,1.204,343
band-hz,40,800,41
board-band-hz,400
fit-order,8
min-panels-per-wavelength,6
full-scale-pa,2
receiver-m,0.6,1.0,3.0
moving,soundboard_skin
rigid,cabinet
rigid,lid
```

A second `receiver-m` row requests a second physical receiver. Positions and
medium values are coherent SI in the structural board frame. Transform OBJ
coordinates by `p_m = obj-scale-m * (p_obj - obj-origin)`; origin is in source
OBJ units. Axes must already agree with the board axes. This does not infer
orientation, centimetres/metres, acoustic impedances or mechanical material
cards from visual MTL properties. No referenced material file is opened.

Every OBJ region must match exactly one `moving`/`rigid` label; unused labels,
unmapped regions and ambiguous assignments refuse. A moving skin gets its
normal velocity from the SAME structural eigenvectors and mass-loading basis
as the string-bridge mechanics. Facet-normal projection must lie in a real
structural triangle and within the declared skin offset, at every acoustic
corner and three quadrature points. Through-thickness rotation contributes
`theta cross arm`. The crowned path retains all Cartesian motion and the
actual equilibrium geometry, not the flattened Rayleigh field.

## Geometry and numerical admission

Supply outward-oriented closed components. A finite board needs top, bottom
and edge surfaces; a floating rigid lid must also have thickness and closure.
Exact coincident OBJ seam vertices are identified, but no tolerance welding,
hole filling, face flipping or invented cabinet is performed. Duplicate faces,
inconsistent/open edge uses and nonpositive component volume refuse. These
checks are NOT a general self-intersection, nonmanifold-vertex, component-overlap
or fluid-accessibility certificate; those remain input responsibilities. Do not
supply intersecting, nested or inaccessible closed pieces as independent bodies.

The current dense BEM budget is 2048 acoustic triangles and 128 input modes,
with a bounded structural-projection product. Each requested frequency must
meet the declared minimum panels/wavelength (at least six); negative radiation
power beyond roundoff and nonfinite diagnostics refuse. The reported condition
number lower bound is not a condition-number certificate. Refine the input
mesh or choose a narrower explicitly declared band when resolution is refused;
no modes or frequencies are silently dropped.

Receivers must conservatively lie outside the body's enclosing sphere, leaving
at least two 48 kHz propagation samples and at most 0.5 seconds of guaranteed
travel time. This excludes some legitimate close microphones. Receivers below
the soundboard are legal. The BEM evaluates finite-distance pressure, including
spreading and propagation; no extra 1/r or gain is applied.

The CSV contains receiver-major/input transfer values at each frequency:
pressure [Pa] per unit generalized modal acceleration [m sqrt(kg)/s^2], with
`exp(-i omega t)` phasors. The boundary input is `v_n = i shape / omega` for
unit acceleration. All receiver responses share one boundary factorization
and modal batch at a frequency. This is a frequency-domain transfer, not yet
a real-time callback or a measured instrument response.

## Physical boundaries and assets

This is one-way linear exterior acoustics on fixed geometry about the supplied
structural equilibrium. Radiation pressure does not feed back into mechanics.
The selected rigid cabinet/lid does not flex or absorb; no room, air absorption,
string-direct radiation, nonlinear air or above-band accuracy is claimed. Mesh
resolution and the structural/modal truncation need independent convergence
and measurement comparisons before a realism claim.

A public asset lead is seavenois's CC0 Steinway D274 on BlendSwap:
https://blendswap.com/blend/7279 . The author describes it as simplified; its
render geometry is not automatically an admissible structural/acoustic mesh.
No external model bytes or unverified factory measurements are included by
this feature. A licensed surface still needs explicit acoustic part mapping,
physical structural cards, coherent scale and mesh preparation. Attribution
and measured/published/estimated labels do not themselves verify the data.
