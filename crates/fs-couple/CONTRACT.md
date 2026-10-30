
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
