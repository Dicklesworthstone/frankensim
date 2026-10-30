
### Native finite-gas heat programs in plate-valve performances

Native plate-valve files accept `finite_gas`, `gas_energy`, `gas_feedback`,
`gas_heat` and signed `heat SAMPLE POWER_W` records after the complete DG input.
They retain original source packs and derive inventory volume from the native
geometry. WAV metadata records the initial inventory, original physical clock,
power assignments and solve controls before type erasure. This is closed,
fixed-volume gas heating with the existing internal pressure observer. It does
not supply mean mass/enthalpy transport, phase change, wall heating, an
independent upstream/branch inventory, or a computed exterior mean-pressure
jump. Such inputs refuse instead of changing their physical boundary.



### Heterogeneous specimens in one finite moisture atmosphere

`FiniteMoistureAtmosphere::advance_shared` solves one endpoint vapor balance
against all original `MoistureSlab` diffusion owners. Each specimen retains its
own material, geometry, source sorption, optional finite surface-transfer law
and transport controls. Every trial restarts all accepted fields and evaluates
all exposed faces at the same RH. A single gas endpoint receives the aggregate
water transfer; simultaneous uptake and drying do not duplicate the reservoir
or acquire sequential specimen-order bias. Stable IDs define deterministic
assembly order, and the accepted clocks, temperature and total-cell budget must
agree before any native solve. The original single-specimen `advance` uses this
same implementation.

Acceptance checks the represented vapor, every noncondensing species, total
mass and the sum of original solid-water changes. A failed later specimen,
source/domain exit, root failure or cancellation publishes no partial gas or
solid state. Returned slabs remain the original source-bound inputs to swelling,
wet inertia and structural stiffness consumers. Gas composition, pressure,
transport and sound speed come from the original thermochemical owner.

The physical scope remains isothermal, well mixed, fixed gas volume and dilute
bound-water diffusion on fixed dry skeletons. Gas internal-energy change is
reported separately and is not sorption heat or a complete thermostat duty.
There is no condensation, gas momentum or automatic downstream mechanical
publication. G1/G3/G4 regressions cover an independent heterogeneous one-cell
balance over continuing intervals, opposite-sign uptake/drying, canonical input
permutation, a companion material changing an unchanged specimen's native plate
response, and whole-candidate refusal/cancellation.


### Shared material for immersed flow, additional inertia and sphere films

`material_immersed::RetainedFluidSpheres` binds one original `MaterialScaling`
and its `FluidTransport` source to the existing D3Q19 moving-sphere assembly.
The same frozen transport snapshot determines the bulk BGK relaxation,
reference-density unit conversion and missing-film dynamic viscosity. Film
cutoff/regularization are explicit geometry/calibration inputs; there is no
independent viscosity knob. Each `AdditionalMaterialSphere` retains its
original `ResolvedBody`, requiring a matching uniform spherical material
region. Source density supplies mass and inertia without material-name cases.

Held SI forces/couples become native impulses using the same dx/dt and density
scales. Actual midpoint gaps drive the original cutoff-subtracted normal film
and optional reciprocal tangent/spin modes. One native solve publishes all
body poses, momenta and fluid populations together. The accepted-step counter
advances only after native acceptance; cancellation/refusal leaves physical
state unchanged. Native receipts retain separate resolved-fluid, film, wall,
external-load and rotational-discretization accounts. Checked conversions
expose their mass, inertia, impulse, angular impulse and energy in SI.

The model deliberately retains the native INTERIOR fluid. Sourced carrier
inertia is ADDITIONAL to it; this is not a displaced-fluid rigid-solid closure,
and subtracting a guessed interior mass is not supported. This lane is frozen,
incompressible Newtonian flow with declared spherical geometry and stationary
or periodic domain faces. It adds no finite heating, phase transition, contact,
moving-wall model, compressible/rarefied film, measured cutoff calibration or
trajectory accuracy certificate. Low-Mach and native geometric admission are
necessary numerical conditions, not proof of those continuum assumptions.

`tests/material_immersed.rs` supplies three focused G1/G3/G4 groups: original
material substitutions through real bulk/film/body evolution; an independent
SI missing-film impulse/work formula and actual population-plus-carrier
momentum; physical load conversion, continuation and cancellation/refusal
rollback; and source/geometry/allocation-domain admission. Fixtures use
synthetic source cards and do not imply experimental validation.
