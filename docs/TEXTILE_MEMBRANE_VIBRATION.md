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

## Material damping and exact coupled modal motion

`fs_solid::membrane::viscosity::MembraneViscosity` adds an explicitly
selected, objective Kelvin--Voigt loss law to the original sheet geometry.
Resolve the six reduced plane-stress viscosity coefficients together with
the sheet properties in one complete material-state bundle, then derive
`ResolvedSheet` and `PlaneStressViscosity` from that same bundle. The selected
viscosity plane is `Axes12`, matching the sheet's ordered warp/weft axes.
Reference thickness comes from the sheet; neither density nor elastic
stiffness is used to invent a loss coefficient.

`evaluate` returns the full resisting nodal force and total/per-triangle
irreversible power. It uses engineering Green-strain rate and integrates
against reference area times reference thickness. This is an explicit
reference-configuration constitutive choice; it does not establish that
small-strain source data are calibrated at finite strain. Admitted natural
references use the hyperelastic owner's same transformation:
for G=diag(g1,g2,g3), D=diag(g1^-2,g2^-2,(g1*g2)^-1), the installed
viscosity is det(G) D V D. Original mesh area, thickness and mass remain
unchanged; the normal reference length contributes through det(G). An identity
reference preserves the original arithmetic.

For a thermal natural reference, use
`MembraneViscosity::new_with_thermal_endpoints` with the original per-triangle
viscosity laws and original material cards. It first verifies each complete
reference sheet bundle, then re-evaluates those same six pinned viscosity claims
at the hyperelastic model's exact endpoint material context. Endpoint passivity,
axes and the complete source frequency band are checked again. A different hot
claim cannot replace an expired original source. Reference thickness and areal
mass retain their original source receipts and are never queried at the hot
endpoint. `reference_laws()` retains the original laws; `laws()` exposes the
actual endpoint coefficients used for force and heat.

The ordinary `new` constructor still refuses a thermal natural reference without
this explicit binding. Endpoint binding prepares a fixed-temperature mechanical
model; it does not evolve temperature or silently rebind an accepted trajectory.

`project_modes` uses the exact geometry retained by `MembraneModes` and
preserves the complete matrix `C = Phi^T C_spatial Phi`, including
off-diagonal terms. Its source square-root rows preserve directional
zero-loss channels. Every retained eigenfrequency enclosure must fit the
declared source band. That admission does not prove that arbitrary forcing
or a finite-amplitude trajectory stays within the band's physical validity.

A flat, prestretched membrane has zero first-order in-plane strain rate
under purely transverse motion. This mechanism therefore leaves those
linear modes undamped even when their prestress stiffness is positive.
Curvature and in-plane motion produce the corresponding material loss.
Bending, air loading, support friction and material memory require their
own models; none is inserted as a damping floor.

`viscosity::modal::PreparedMembraneDamping` advances the full retained
equation `q'' + C q' + Lambda q = f` through the existing
`fs_phs::linear_exact::PreparedExactZoh`. `from_vibration` binds the original
accepted equilibrium; `new` also accepts an existing `MembraneModes`
basis, including a gas-loaded basis supplied by `fs-couple`. The original
gas stiffness must remain in that basis without being applied a second time.

State order is `[q_0, ..., q_(n-1), qdot_0, ..., qdot_(n-1)]`, with unit
modal mass. Inputs are held generalized forces from `project_nodal_force`;
outputs are interval-mean modal velocities. The returned record separates
mechanical energy change, applied work and independently integrated material
heat. `step_into` reuses prepared work arrays and introduces no new
exponential or owner allocation on a successful warm call. A late
cancellation or refused step preserves both caller output arrays.

Preparation retains all modes and all damping couplings. With two states and
one force port per membrane mode, the existing exact owner's augmented
dimension limit permits at most 42 modes. A larger basis refuses instead
of being silently truncated. `LinearSplitBudget.max_modes` counts scalar
state coordinates (two per structural mode); projection also takes its own
explicit work/storage limit. This path provides linear motion around the
retained configuration, not nonlinear damped membrane time integration.

## Finite-amplitude material damping

`membrane::dynamics::viscosity::DampedMembraneMotion` consumes an existing
`MembraneMotion` and its same `MembraneViscosity`. The original motion
selects lumped, dense consistent or sparse consistent reference inertia;
its constraints, nonlinear constitutive law and prepared factors are retained.

At each held endpoint geometry the material Gram rows give
`C(u) = sum s_r s_r^T`. Each rank-one viscous channel advances by its exact
exponential, with independently integrated positive heat. Forward and reverse
channel sweeps form a symmetric damping stage. Half a damping stage, the
original full nonlinear Verlet step, and another half damping stage form a
second-order composition. Noncommuting channels and evolving geometry still
have splitting error: use time refinement for motion and heat. The existing
complete modal exponential remains the exact choice for its admitted linear
model.

The new path has no modal truncation or 42-mode ceiling. Lumped damping touches
at most nine local coordinates per source row and costs O(triangles) per sweep.
Consistent damping uses the original prepared mass inverse for its nonlocal
velocity response. All modes of inertia retain stationary-support impulses,
including the prescribed/free consistent-mass contribution. Warm steps reuse
the wrapper's state and work buffers; original constitutive owners keep their
own costs.

`MembraneDampingBudget` bounds source-row storage, consistent inverse actions,
positive step duration and absolute signed energy defect. The record reports
physical material heat separately from mechanical energy, held-load work and
numerical error. Last-step triangle heat and accumulated heat are available
for physical accounting; neither introduces a temperature evolution law.
Energy closure alone does not establish phase accuracy or stability.

A flat prestretched transverse launch initially has zero material loss, but
finite transverse deformation creates Green-strain rate and activates the
corresponding damping. Source null directions remain intact. Admitted natural
references, including explicitly bound thermal endpoints, use the same
transformed Gram rows described above. Trials that violate the material domain, step/work budget, energy
budget or cancellation leave accepted motion and heat unchanged.

## Focused native checks

```sh
cargo test -p fs-solid --test membrane -- --test-threads=1
cargo test -p fs-solid --test membrane_viscosity -- --test-threads=1
cargo test -p fs-couple --example textile_sheet -- --test-threads=1
```

The damping sources received API/physics review and rustfmt parsing; new
runtime and test modules were formatted. Native Cargo compilation and tests
remain unexecuted: the restored host has a full shared filesystem, a stale
checkout without these commits, and a refused sibling-dependency preflight.
Independent numerical references checked natural-frame work/heat, exact
channel flows and second-order refinement. These do not substitute for the
native checks above.

The example remains a fixed-state elastic membrane and lossless vibration
demonstration. The opt-in viscosity APIs add complete linear modal evolution
and nonlinear damped motion when the source supplies the required coefficients.
The nonlinear wrapper retains the original dead-load/stationary-support scope;
it supplies no follower pressure, gas feedback, bending/contact, evolving
material state, continuum stability certificate or empirical fabric validation.
