# Coupled magnetic-current, force and heat evolution

`MaterialSlab::step_magnetothermal_midpoint` closes temperature and velocity
feedback within a single step of the existing 2D/3D thermal-flow solvers:

```
(T_next, u_next) = existing_flow_heat_step(
    accepted_populations,
    Joule(T_mid, u_mid), Lorentz(T_mid, u_mid) + external_force)
T_mid = (T_old + T_next) / 2
u_mid = (u_old + u_next) / 2
```

Both endpoints use the existing Guo-corrected velocity convention, including
its retained force field. This is a definition of the coupled discrete map,
not a claim of second-order temporal accuracy or exact kinetic-energy balance.
The force-discretization convention must not be independently changed in this
adapter. Background: Guo, Zheng and Shi, Physical Review E 65, 046308 (2002),
DOI 10.1103/PhysRevE.65.046308.

Every trial restarts from the same accepted populations. The material law,
current matrix, contact heat partition, reciprocal Lorentz force, heat source
conversion, buoyancy, and thermal wall exchange retain their existing owners.
The source-evaluation and power-accounting helpers are shared with the original
explicit path. Iterating does not repeatedly advance time or add heat again.

## Acceptance and return value

`material::electrothermal::magnetic::midpoint::MagneticMidpointOptions` provides
separate coordinate tolerances in K and m/s, Joule-source tolerance in J,
force-impulse tolerance in N s, and mechanical-source-work tolerance in J.
The coordinate residual alone cannot accept a step. A fresh current solve at
the actual candidate midpoint checks L1 local heat, impulse-component and work
differences. Opposite regions cannot cancel a local source error. Under-relaxation
changes the iteration only; a small relaxed update is not a convergence test.

The returned `MagneticMidpointStep` contains the source-producing electrical,
force and heat fields, their evaluation temperature and velocity, all residuals,
iteration/solve counts, actual thermal wall exchange, and source-energy accounts.
The final consistency solve does not replace the sources actually applied.

`midpoint_mechanical_work_j` is evaluated magnetic power times dt, NOT measured
kinetic-energy change. The coupling defect checks source transfer, not total
energy of the split hydrodynamic evolution. Density, heat capacity and thermal
transport remain the retained material snapshot during a step. B, external
forces, contact data and electrode drives are fixed during each nonlinear call.

All unsuccessful trials leave the accepted slab unchanged, including force
history. Captured material-query receipts remain caller-owned. There are at
most 128 trials and two current solves per trial. Candidate temperature and
lattice-Mach bounds remain active; a midpoint iteration can still fail or be
unsuitable. A nonconvergence error does not prove physical runaway. No magnetic
induction, Hall response, phase change, empirical calibration, total-energy or
unconditional stability claim is added. Cancellation is polled between trials,
in field reductions, and by the current solver; the existing thermal step is
still the cancellation granularity for thermal evolution. Peak memory includes
original, trial and thermal-owner tentative grids plus two electrical results.

## Existing executable, unchanged material input

Append `--midpoint` to opt in. The seven-value material files, geometry, drivers,
and original explicit invocation are unchanged. Midpoint CSV output appends
trial count and the five residuals; the explicit CSV columns remain unchanged.

```
cargo run -p fs-lbm --example material_magnetic_flow -- \
  crates/fs-lbm/examples/material_magnetic_flow/conducting.txt 20 100 --midpoint
cargo test -p fs-lbm --test material_magnetic_flow --test material_magnetic_midpoint
```

Eight new Rust regressions cover actual source closure, zero-field compatibility,
bounded magnetic braking, independent force/work gates, rollback/cancellation,
2D depth and reversal, material substitution and input refusal. They are authored
but unexecuted in the implementation environment (Cargo/rustc unavailable).
Independent Python circuit plus lattice equations and SciPy nonlinear roots are
reference checks only; they do not replace native Rust or experimental testing.
