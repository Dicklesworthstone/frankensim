
## Source-integrated Fourier profiles at gray conducting gaps

`thermal::conduction::GrayLayerGap::equal_area_full_view` explicitly declares
two equal-area, opaque diffuse-gray faces with unit mutual view factor and a
fixed nonradiative gap HTC. `TwoLayerConduction::evaluate_gray_gap` resolves
their two material-face temperatures using the original k(T) Kirchhoff
integrals and inverses. The ordinary contact and gray-gap paths share the
existing supported scalar flux bracket and safeguarded iteration. No extra
temperature profile, thermal timestepper, wall storage or source query is added.

The interface relation is `j = h_g*(Ua-Ub) + epsilon_eff*sigma*(Ua^4-Ub^4)`;
each half cell satisfies its own integrated Fourier equation at that same j.
The full analytic bulk-temperature derivatives include the ratios of bulk
to surface conductivity. Individual source spans need cover only their own
material profiles, including when the gap bridges disjoint temperature spans.
Disabling radiation delegates to the original contact-resistance owner and
retains its exact result and near-perfect-contact representation gate.

`ImplicitLayerLink::with_gray_gap` binds this complete response to the original
chain/graph enthalpy solve, with the original caloric dT/dh chain rule. One
equal/opposite heat transfer enters its energy and entropy accounts; reported
conductive and radiative components are descriptions of that same transfer.
The attachment requires zero authored contact resistance, no separate contact
law, equal original covered areas and planar geometry. Source layers retain
their explicit pure-phase restriction. Radial radiation geometry, nonunit view
factors, phase fronts inside a half cell, evolving emissivity and gap storage
remain outside this closure.

Regression definitions compare independent nonlinear material integrals and
surface balances, the original constant-k gray link, exact radiation-disabled
contact recovery, reversed ordered materials, full flux tangents, and a
manufactured original-source chain/graph endpoint with conserved heat. Source,
ownership, nonlinear work and cancellation refusals retain the original state.
