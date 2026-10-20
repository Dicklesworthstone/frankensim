# Material viscosity in the actual momentum operator

`fs_lbm::conjugate::simple_flow_with_viscosity` accepts one dynamic molecular
viscosity per voxel and uses the existing laminar staggered MAC/SIMPLEC owner.
Density remains uniform. Geometry, boundaries, fans, convection, porous drag,
pressure correction, numerical budgets and cancellation keep their original
owners. The legacy constant-property entry keeps its original arithmetic.

The variable-viscosity entry evaluates the complete Newtonian stress

```
tau = mu * (grad(u) + transpose(grad(u)))
rho * div(u tensor u) = -grad(p) + div(tau) + f
div(u) = 0
```

Changing the coefficient of a component Laplacian alone would miss reciprocal
shear forces when viscosity varies. The new entry uses `2*mu*du_a/dx_a` on
normal faces and both velocity gradients on each shear edge. Reciprocal
components and opposite control-volume sides use the same positive harmonic
edge coefficient, accumulated in fixed cell order. Half-cell wall patches use
their actual parallel areas. Existing Darcy terms use the same local molecular
viscosity; Forchheimer terms keep the unchanged reference density. Cell Reynolds
numbers use the applied local viscosity.

This is a scalar, Newtonian, constant-density, steady laminar model. Viscosity
must be positive and representable; shape mismatches and LVEL are refused. It
does not imply variable-density continuity, an evolving free surface,
phase-dependent mobility, non-Newtonian rheology, turbulent transport, kinetic
storage or viscous heating. Source selection and a viscosity field's physical
temperature/composition epoch belong to its material caller. Staircase wall and
open-boundary approximations remain those of the original discretization;
resolution and convergence must be established for each scenario.

Focused regressions in `conjugate/simple/material_viscosity_tests.rs` exercise
piecewise-viscous Couette traction, a manufactured reciprocal shear force under
a spatial viscosity gradient, variable-viscosity Darcy flow, reciprocal edge
identity and unsupported input/closure refusal. These use the actual momentum
owner, not a replacement reference solver. Merely authoring them is not a claim
that the native test suite was executed in an environment lacking dependencies.

## Material density drives natural convection on the same heat timeline

`LiquidFlowEnthalpyInput::feedback_boussinesq` connects the original pinned
liquid and caloric sources to the existing `MaterialFlowFeedback` owner.
Heating can now initiate flow without a fan: each actual thermal substep uses
the same material evaluation's density, viscosity and conductivity for buoyancy,
full Newtonian stress and optionally assigned fluid-face conduction.

The reused `BoussinesqMaterial` projection retains the source response and
receipts, and keeps reference density in inertia, cell masses and native face
mass rates. The body force is `(rho(T)-rho0)*gravity` in N/m³. Pressure openings
are relative to the reference hydrostatic state. No separately authored
expansion coefficient or second buoyancy force is added. The native public
`simple_flow_with_viscosity_and_force` entry applies that force to the original
staggered momentum volumes, with unchanged pressure correction, resistance,
fan and boundary owners. Its force vectors must be finite and zero in solids.

The generic `MaterialFlowFeedback::try_new_boussinesq` entry accepts any existing
concrete fluid-transport source. The liquid convenience entry additionally
requires the original Cp/density selections, non-temperature context and fully
liquid support. It allows a temperature-dependent source density only under
explicit Boussinesq departure bounds. The registered borrowed caloric curve
still owns h/T and heat capacity; the projected reference Cp does not replace
that caloric law. Caloric density metadata cannot change the fixed cell mass.

Bounds apply at initial/source epochs and proposed thermal endpoints. Invalid
gravity, source support, excessive density/Cp departure, a changed caloric
binding or a refused inner history discards the whole candidate. Applied
epochs retain the actual voxel buoyancy field alongside source samples and flow.

This admits small density departures in quasisteady natural/mixed convection;
the user must establish that approximation and thermal timestep refinement.
The heat audit covers original heat sources, conduction and enthalpy transport.
It does not include gravitational/kinetic storage, pressure/viscous work, moving
boundaries, phase mobility or free surfaces. Bounds are admission checks, not
experimental validation or a closed total-energy certificate.

G1/G3 regressions in `tests/liquid_flow_enthalpy/feedback.rs` change only the
source density curve in a heated porous column with two pressure openings and
no fan. They compare native flow and heat with independent Darcy resistance,
conduction and donor-cell balances, including same-state source pins and fixed
mass. G4 endpoint-departure and binding refusals preserve accepted heat/history.
