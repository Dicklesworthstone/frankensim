# Material-driven magnetic current, reaction and thermal flow

`ElectricSetup::solve_with_cell_motion` extends the existing voxel current solver
with `J = sigma(-grad(phi) + u x B)`. It accepts actual SI velocities and an imposed
magnetic flux-density field through `ElectricMotionInput`; the normal cell material
callback supplies resistivity, or explicit insulation. The same half-cell geometry
used for electrical resistance integrates motional voltage and scatters its
transpose magnetic force. Graded grids, heterogeneous fluid/solid assignments,
voltage/current electrodes and the existing electrical contact partitions remain
supported. Occupancy never changes just because a fluid conducts electricity.

The returned `MovingConductorSolution` separates conductor-frame Joule heat from
mechanical work. Both identities are checked at the supplied frozen state:

```
electrode power = sum(I^2 R_bulk + I^2 R_contact) + sum(force dot velocity)
sum(force dot velocity) = -sum(face current * motional voltage)
```

Heat is never replaced by `J dot E`: doing so would mistake mechanical exchange
for irreversible heating. A short-circuited moving conductor can therefore heat
while delivering zero electrode power; its Lorentz work is negative (braking).
Contacts have zero thickness and no motional voltage of their own. Ideal electrode
cells have no bulk half resistance or magnetic-force volume; these electrical
constraints do not model the mechanical or electrochemical electrode bodies.
Every conducting component must have a voltage reference. Floating conductors
are explicitly refused rather than silently dropping their possible loop currents.

## Live fluid integration

`MaterialSlab::step_magnetothermal` uses the accepted slab's local temperature and
Guo-corrected velocity (including its retained previous-step force convention).
It solves the current, adds the resulting Lorentz force density to the caller's
external forces, and supplies the spatial bulk/contact Joule powers once to the
existing `step_with_power_density`. That path retains buoyancy and independently
measured thermal wall exchange. No alternate material database, current solver,
flow solver or heat discretization is introduced.

All three vector components participate in 3D. In 2D only an out-of-plane Bz is
admitted: the method refuses an in-plane magnetic field rather than discarding
unresolved out-of-plane current or force. The same explicit extrusion depth as
the electrical/thermal adapter controls current and represented energy, not the
2D velocity/temperature evolution. Electrical outer faces are insulated even
where the fluid is laterally periodic; this is not a periodic electrical solve.
Thermal wall cells remain insulating and cannot contain electrodes or Joule heat.

The electrical sources are evaluated once at the accepted initial state and held
over the step. They are recomputed on the next call. `MagneticHeatStep` retains
actual current/heat/force fields, electrode energy, Joule energy and independent
thermal exchange. `lagged_mechanical_work_j` is initial Lorentz power times dt,
**not** the measured kinetic-energy change. Its `coupling_defect_j` checks the
source/heat transfer with that lagged work; it does not establish conservation of
total energy across the split fluid time integrator. Strong magnetic damping
needs time refinement and can make explicit coupling unsuitable.

Material, geometry, field, temperature, power, transfer and cancellation failures
leave the original slab intact. Input and candidate velocities must pass the
existing necessary lattice-Mach limit, but this is not a general stability proof.
Scratch includes the original, trial and thermal-owner tentative states. The
caller must supply a physically admissible imposed field and establish the
low-induction regime. There is no magnetic induction solve, source-field dynamics,
magnetization/Hall force, moving mesh, latent heat or empirical calibration claim.
Thermal coefficients remain the slab's retained material snapshot.

## Data-only runnable example

The two files differ only in resistivity. B is a separate environmental input;
changing its sign reverses the driven motion, not the initial resistive heat.
The same seven-value material format is used by `material_joule_flow`.

```sh
cargo run -p fs-lbm --example material_magnetic_flow -- \
  crates/fs-lbm/examples/material_magnetic_flow/conducting.txt 2 100
cargo run -p fs-lbm --example material_magnetic_flow -- \
  crates/fs-lbm/examples/material_magnetic_flow/resistive.txt 2 100
cargo test -p fs-lbm --test electric_motion --test material_magnetic_flow
```

CSV contains time, electrode current, Joule power, instantaneous Lorentz power,
mean velocity, mean Kelvin temperature and the coupling-accounting residual.
Eleven Rust regressions are authored across the two increments. Native compilation,
tests and example execution remain unverified in the implementation environment:
Cargo/rustc are absent. Independent incidence-circuit and thermal-LBM equation
checks are numerical reference evidence, not execution or validation of Rust.
