

### Finite-inertia guided walls on original free-surface liquid momentum

`free_surface_wall::GuidedFreeSurfaceWall` retains the complete original
material-bound D3Q19 free-surface liquid and a native finite-mass rigid state.
One fully flooded planar wall translates along an explicit unit tangent guide,
with no normal displacement or spin. Each accepted step solves the shared
midpoint velocity against the ACTUAL native wall impulse and an optional held
guide force. Safeguarded scalar secant/bisection reuses the native liquid owner;
no fitted drag, added-mass coefficient or replacement fluid algorithm is used.

All trials start from the same accepted fluid. A new immutable native wall
candidate API shares exactly the existing moving-wall admission and receipts,
avoiding an additional full-population clone per trial. Complete domain cells,
aggregate trial-cell work and trial count are bounded before the respective
scan or allocation. Only the accepted candidate advances the liquid, rigid
momentum, midpoint displacement, elapsed time and accepted-step count together.

The native rigid impulse receipt, actual midpoint mismatch and liquid boundary
work must meet separate physical tolerances. Small representable impulse/work
defects are retained without subtracting a large kinetic-energy baseline.
The raw guide reaction closes native numerical momentum; pressure-perturbation
and explicitly pressure-referenced physical support loads remain separate.
An absolute load is not invented when the source has no pressure reference.

This supplies finite translational wall feedback, not moving occupancy, wall
rotation or a complete angular-momentum audit: native wall link moments are
unavailable. Original isothermal VOF viscosity, forcing, capillarity and phase
conversion keep their numerical/physical limits; boundary work balance does
not turn that discretization into an energy-conserving thermal fluid. Five
integration regressions cover independent D3Q19/finite-inertia momentum and
work oracles, bitwise native replay and continuation, pressure-datum shifts,
bounded geometry/work/trials, cancellation and late native liquid refusal.
