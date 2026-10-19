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
