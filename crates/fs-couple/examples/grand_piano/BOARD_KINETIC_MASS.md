# Consistent kinetic energy for the cubic soundboard

The flat-board geometry format accepts an optional row:

```text
panel-rotary-mass,consistent-cubic
```

Together with `--edge-cubic-board-mass`, this selects exact panel translation
and rotary inertia from the same declared cubic displacement used by bridge
forces and acoustic motion. Omitting the row, or selecting `lumped`, preserves
the original panel-rotary image and its floating-point path. This is an explicit
alternative for physical convergence studies, not a new default piano voicing.

For a supplied flat geometry containing the row:

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --board-geometry kinetic-board.fsb \
  --edge-cubic-board-mass --equilibrate-board-mass \
  --acoustic-refinement-levels 2 --render kinetic-piano.wav
```

The file must supply the physical board, material sections, supports and bridge
stations required by the selected scale. The ordinary `--dump-geometry` export
can be used as a starting point without treating its estimates as measurements.
The same row is consumed by exterior preparation through its shared flat-board
loader and existing edge-cubic mass selection. No new oscillator, soundboard
frequency table or output filter is introduced.

## Kinetic model

For transverse midsurface displacement `w`, the declared through-thickness
velocity field is `(-z*w_dot,x, -z*w_dot,y, w_dot)`. Thus

```text
2*T = rho*h * integral(w_dot^2 dA)
    + rho*h^3/12 * integral((w_dot,x)^2 + (w_dot,y)^2 dA).
```

The existing Bernstein field supplies `w` and its analytic physical gradients.
Their products are integrated by exact simplex polynomial moments, separately
for every supplied section thickness and density. The previous lumped panel
rotary diagonal is removed, not retained as an extra mass. Stiffener inertia
continues to follow the independent `stiffener-mass` row; none of it is removed
or inferred. Geometry, total material mass, DKT stiffness, prestress and supports
are unchanged. The new mass changes the actual generalized eigenproblem and
its mass-normalized force and motion projections, not just reported frequencies.

The selected pencil is used before the existing eigensolver, optional bridge
Ritz reduction, dense material-damping projection and loaded-string bank.
Bridge, volume, refined Rayleigh samples and full-vector exterior motion still
come from that one final physical basis. Reports identify the selected panel
inertia and the separate beam law.

## Admission and scope

`consistent-cubic` without the edge-cubic displacement selection refuses before
assembly. Unknown/duplicate rows refuse. Crowned geometry rejects this flat
row, even when its value is `lumped`, because the shell owns a different six-DOF
kinetic field. Nonrepresentable integrals or unresolved removal of the old
rotary diagonal refuse rather than adding a mass floor or clamping eigenvalues.

This is a Rayleigh/Kirchhoff **kinetic** image of the existing edge-cubic field,
not a Mindlin transverse-shear model or an exact continuum displacement field.
The DKT stiffness is unchanged; spatial shear deformation, bonding loss, rim
compliance and material calibration are not supplied by this option. A chosen
mode cap, frequency band or acoustic refinement does not become a high-band
convergence or listening-quality certificate. The high-frequency piano bead
`frankensim-5c0sj` still requires source-mesh, transfer and listening comparisons.

## Focused verification

```sh
cargo test --release -p fs-plate --lib edge_cubic
cargo test --release -p fs-couple --example grand_piano kinetic
cargo test --release -p fs-couple --example piano_exterior kinetic
```

Five new plate regressions cover exact quadratic continuum products, independent
positive Gauss integration, kinetic positivity, rigid motions, heterogeneous
sections, beam preservation and support restriction. Four new shared piano
regressions cover unchanged defaults, option refusal, the actual modal residual
and mass normalization, common bridge/surface motion, reduced material damping,
acoustic refinement and crowned-input refusal.

An executed independent NumPy/SciPy reference checked 240 randomized triangle
and section cases plus 108 analytical polynomial products. Maximum relative
rotary-matrix discrepancy against separate pointwise quadrature was 1.45e-15.
These numerical references are not native Rust execution, full-preset renders,
measured soundboard comparisons or listening validation. Native tests must pass
with the repository's actual toolchain and pinned dependencies before those
claims can be made.
