

## Compact borrowed free-surface impulse staging

`FreeSurface3::prepare_impulse_stencil` retains a mutable borrow of the original native owner while sampling a bounded, strictly increasing wet-cell stencil. Native populations, located retained liquid mass, physical Guo-frame velocity and original cell centers feed the shared population impulse kernel. Gas, inactive/static solid, piston, unresolved interface inventory, excessive speed and unrepresentable states refuse. Staging preserves original populations until the final commit; dropping a stage/candidate or final cancellation leaves the original unchanged. Candidate readback exposes each original mass/position and actual endpoint velocity so a finite receiving owner can independently gate reciprocal angular momentum and work before joint publication. No fluid clock, held force, fill/classification, geometry or unlocated carry is changed.

`MaterialFreeSurface3::prepare_native_impulse_stencil` reuses original physical impulse admission and scenario speed, while explicitly returning native-unit staging for a composing owner’s SI transaction. It does not silently reinterpret lattice readback as physical units.

The native stencil's three focused tests execute actual retained populations, force-frame stress, readback and transactional refusal. A further native refinement regression checks retained solid keys outside the physical domain and beyond representable child Morton coordinates: refinement refuses before allocation, preserving source mass, geometry and populations. All six native refinement cases pass after this guard. The prior three SI receipt checks remain unchanged; full streaming/material integration was not executed.

## Temperature-activated molecular wetting

`MolecularKineticLaw::try_new_arrhenius` optionally evolves the equilibrium jump
frequency as `K(T)=K_ref exp[(E_a/k_B)(1/T_ref-1/T)]` inside the original calibrated
temperature domain. Activation energy is explicit joules per molecular event;
the frequency in `MolecularKineticControls` is its value at the declared reference
temperature. This constant-prefactor Arrhenius closure adds no inferred adhesion
measurement, Eyring prefactor, extrapolation or activation-heat source. The
original constructor still selects constant frequency. Reference-temperature and
zero-activation responses preserve that original model exactly.

The same law feeds the existing dense, sparse and enthalpy wall-angle owners,
retaining original interface claims and actual local temperature. Accepted samples
report the local jump frequency separately from static material measurements.
The rate and velocity scales must be finite and positive across both temperature
bounds; thermal domain refusal retains native flow, enthalpy and wetting history.
Per-length molecular loss remains diagnostic, without claiming a resolved line
measure, a global surface-energy balance or heat deposition. G1/G3 unit cases
check independent forward kinetics, low-speed friction and constant-model limits;
the existing thermal integration case also exercises activated rates and rollback.

## Material-flow enthalpy reconstruction

`MaterialFlowAdvection::with_reconstruction` attaches the original thermal owner's shared `LimitedEnthalpyReconstruction` to solved material face fluxes. Selected contiguous axial cells reconstruct total specific enthalpy, including a declared passive equal-density latent inventory; paired face transfers, source-law inversion, accepted history and the combined conduction/advection timestep restriction remain with `fs-thermal`. Unselected cells and external streams retain donor enthalpy. Empty plans retain the donor-cell arithmetic. This is an explicit first-order-in-time stationary-control-volume model; it does not infer multidimensional junction slopes, phase mobility, density-change flow, pressure work or moving interfaces. The `advection()` accessor exposes the underlying graph and actual transport audit; implicit clients must attach the same `reconstruction()` plan to their own native implicit transport owner.

## Hydrostatic rest under material density-dependent gravity

`MaterialFreeSurfaceScaling::free_surface3_hydrostatic` initializes a planar
column using the original material-bound Guo force, collision, streaming and
gas-pressure reconstruction. Gravity points toward one declared Cartesian
supporting wall; zero gravity is also admitted. Complete liquid layers retain
the original stationary solid tiles and at least one layer of gas. The caller
supplies a positive absolute pressure reference and a density-excursion bound.
Capillarity and wetting must be explicitly absent.

For lattice acceleration magnitude `a`, set `h = a/(2*cs^2)`. The nearest wet
layer has `rho = 1/(1-h)`; successive deeper layers multiply by
`(1+h)/(1-h)`. This surface factor follows the native missing-gas-link rule.
A continuum pressure profile evaluated at cell centers alone does not supply
the same discrete rest state. Layer densities are computed once and reused
while the original constructor owns liquid mass, sparse activation and the
physical-zero-velocity Guo frame. Ordinary `step()` retains the same gas datum.

This balances the weakly compressible numerical operator. It is not a liquid
equation of state or an incompressible mass projection: native mass contains
the retained density perturbation. Pressure and submerged traction approach
their incompressible hydrostatic limits as `dt/dx` decreases with the physical
configuration fixed; refining `dx` at fixed `dt/dx` does not remove that
compressibility error. The constructor adds no moving-solid or partially
wetted-body force closure.

G1/G3/G4 regressions in `tests/material_free_surface3/hydrostatic.rs` exercise
all six gravity directions, zero-gravity equality, pressure-datum invariance,
actual reflected-link traction and torque with decreasing compressibility,
mass retention and transactional refusal. The declared density, source-domain,
geometry and cancellation bounds remain with the original owners.
