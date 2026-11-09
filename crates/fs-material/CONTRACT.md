

## Conservative gas-energy coordinates below a temperature quantum

`GasCaloricBranch::sample_density_internal_energy_coordinate` is an explicit
numeric-stage view for an ODE owner retaining its own conservative energy
coordinate. An otherwise refused nonzero request can sample the unchanged
reference temperature only if the exact original Cv integral to the directional
adjacent float proves it is the closest represented temperature, the neighbor
is inside the admitted branch and the caller's absolute energy tolerance is
met. The unresolved energy is returned in the sample residual and MUST remain
in that separate coordinate. Ordinary converged samples preserve their bits.

This does not relax standalone gas heating or accepted parcel reconstruction:
those strict APIs still refuse a nonzero increment that would disappear into
an unchanged state. Caloric bubble and hydraulic pocket stages retain the
energy coordinate; their original-source endpoint and closed-window energy
checks remain strict. Source failures, bounds and cancellation propagate.

All nine focused original prepared/floating-point-floor regressions passed,
including the retained-coordinate test and the existing strict heat-increment
refusals. The native cavitation source subset passed 82 tests. These are
focused numerical checks, not full workspace/DSR or experimental validation.
