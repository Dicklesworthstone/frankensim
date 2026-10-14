# Electrode-free conducting bodies and fluid braking

`ElectricSetup::solve_isolated_cell_motion` admits one or more unwired conducting
components with an empty electrode list. It uses the existing material callback,
voxel current matrix, half-cell motional voltage, contact-resistance heat partition
and reciprocal Lorentz force. The original grounded methods are unchanged.
Neither a fluid/solid label nor a material name selects electrical conductivity.
Missing data remain an error; `None` explicitly declares an insulating cell.

## A gauge is not an electrode

An isolated potential is determined only up to one additive constant per
connected conductor. The lowest-index conductor in each component is assigned
zero potential by symmetric row/column elimination. This is only a numerical
coordinate choice. It does not turn the cell into a perfect conductor or remove
its resistance, motional voltage, Joule heat, contacts, mass or magnetic force.
Disconnected conducting cells also have a gauge and no artificial leakage.

Every face current is independently accumulated at both endpoints, including
the gauge cells. Each current residual must satisfy the caller's explicit
`current_tolerance_a`. The existing independent power and transpose-work checks
remain active. There are no electrode currents or electrical supply power:

    conductor-frame resistive heat + magnetic mechanical power = 0

This frozen-state identity is not a certificate for time-integrated total energy.
Uniform translation in uniform B produces a cancelling potential gradient, not
invented drag. Nonuniform motional voltage can drive circulating currents. The
regressions include both cases, graded cells, explicit interfaces and separate
bodies. Potentials are gauge-relative, not predicted absolute charging voltages.

## Existing thermal flow, including within-step feedback

`MaterialSlab::step_magnetothermal_isolated` and
`MaterialSlab::step_magnetothermal_midpoint_isolated` expose this same source in
2D and 3D. Both use `ElectricHeatInput` with an empty electrode list and an extra
ampere tolerance. The source-producing currents, Joule heat and force enter the
existing material-bound flow/thermal step, not a new discretization. Gauges are
recomputed from actual material conductivity on every source evaluation.

The explicit method holds initial-state sources for one step. The midpoint
method reuses the same bounded temperature/velocity iteration, final local heat,
impulse/work checks, and rollback as the wired method. Thermal walls insulate
electrically; the 2D path requires out-of-plane Bz and explicit extrusion depth.
Electrical exterior faces insulate even where fluid boundaries are periodic.
Failed material, current, field, thermal, convergence or cancellation admission
leaves the accepted slab unchanged; captured query receipts remain caller-owned.

The imposed magnetic field is stationary within a solve. The caller must establish
its physical validity and the low-induction regime. There is no induction equation,
Hall effect, electrochemistry, magnetic-source dynamics or changing mesh. Thermal
coefficients remain the retained material snapshot. Mechanical source work is
not the measured kinetic-energy change. Viscous/numerical losses are not silently
converted into Joule heat. The path does not solve a mixture of floating bodies
and externally driven electrode components; use the existing grounded path for
wired circuits. No empirical calibration or unconditional stability is claimed.

## Data-only executable

The existing seven-value synthetic material files are reused unchanged. Field
values are separate inputs. The example starts the fluid at 0.02 m/s with fixed
thermal walls and a spatially varying, divergence-free imposed Bz. Viscous wall
losses also slow the fluid; the magnetic contribution is separately reported.

```sh
cargo run -p fs-lbm --example material_isolated_braking -- \
  crates/fs-lbm/examples/material_magnetic_flow/conducting.txt 20 60 20
cargo run -p fs-lbm --example material_isolated_braking -- \
  crates/fs-lbm/examples/material_magnetic_flow/resistive.txt 20 60 20
cargo test -p fs-lbm --test electric_isolated_motion --test material_isolated_magnetic_flow
```

Seven core and four coupled Rust regressions are authored. Native tests and
example execution remain unverified in the implementation environment because
Cargo/rustc are absent. Independent NumPy incidence and full singular-system
references check 200 heterogeneous cases, gauge invariance, continuity, braking
and heat. Independent coupled lattice equations also check the fluid transfer;
these are reference equations, not Rust execution or experimental validation.
