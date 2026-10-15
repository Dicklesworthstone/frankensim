# Source-bound membrane equilibrium and vibration

The same imported effective-sheet material now drives spatial deformation,
reference mass, gravity, prestress, natural frequencies and mechanical/fluid
port projections. This path does not turn fiber bulk moduli into sheet stiffness
or infer bending, heat capacity, damping or radiation from missing data.

## Run the existing textile executable

```sh
cargo run -p fs-couple --example textile_sheet -- \
  /path/to/material.pack --membrane 293.15 1.02 0.05 9.81 0.3 0.1
```

After the pack, the six numerical arguments are temperature in K, imposed
perimeter stretch ratio, downward centre-node dead force in N, downward gravity
in m/s2, material grain angle in radians, and maximum absolute engineering
Green-strain component. The last is a caller operating limit, not empirical
validation of the source. Stretch must exceed one for this tensioned example.
The original conductivity-only and through-thickness heat/plate paths remain.

The mechanical path requests areal density, thickness and the four effective
membrane properties from the existing material-state owner at the requested
temperature and fixed characterization. It uses a 0.1 m square, eight P1
triangles and prescribed perimeter displacements. This leaves the centre's three
Cartesian coordinates free; it is an integration example, not resolved cloth.

Output contains the accepted displacements and support reactions, energy,
Newton diagnostics, all three prestressed frequencies, eigenvalue residuals
and intervals, and oriented pressure/volume ports. The modal result comes from
this same accepted equilibrium and constraints, not a second unstressed solve.
Changing source mass, membrane stiffness, grain or supports changes the physical
problem; there is no target pitch or frequency retuning.

## Reuse the physics owners

`fs_solid::membrane::MembranePatch::solve_vibration` combines the original
nonlinear solve and exact material-plus-geometric tangent with
`fs_modal::eigh_gen_dense`. `VibrationSettings` limits the free-coordinate count
and relative eigenpair residual. The dense leaf has a hard cap of 128 free
coordinates; cancellation is checked at tangent columns and around that leaf.
Nonpositive or unresolved eigenvalues refuse rather than being clipped, dropped
from a frequency window or converted to artificial pitched oscillators.

`MembraneVibration` retains a borrow of the original source-bound patch, the
accepted equilibrium, free-DOF map and mass-normalized modes. Generalized
coordinates have units m sqrt(kg). `project_nodal_force` and `reconstruct` use the
same shapes: generalized force and velocity preserve nodal virtual work.
`volume_coupling` gives B in m2/sqrt(kg), with swept volume rate Q = B dot qdot
and pressure force f = B p. The oriented area vectors are integrated on the
actual equilibrium geometry. These linearized ports supply neither a sealed
cavity nor a radiation law; those must come from their physical owners.

The executable regression feeds the resulting modes and projected force into
the existing `ModalAcousticTimeModel` and reconstructs nodal displacement. It
selects zero damping and zero pressure transfer explicitly, so this is a
lossless small-vibration check, **not synthesized microphone audio**.

## Focused native checks

```sh
cargo test -p fs-solid --test membrane -- --test-threads=1
cargo test -p fs-couple --example textile_sheet -- --test-threads=1
```

Native compilation and these tests were not executable in the authoring host
because Cargo/rustc were unavailable. Independent symbolic and NumPy/SciPy
references passed, including the loaded-pyramid frequency of 129.76741009 Hz,
but do not substitute for these native checks. Scope is fixed-state quadratic
sheet elasticity, dead-load prestress and linearized vibration. No wrinkle
relaxation, bending/contact, time-evolving material state, fluid loading,
finite-amplitude vibration, continuum stability or empirical fabric claim.
