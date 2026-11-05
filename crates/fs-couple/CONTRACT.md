
### Resolved electrostatic geometry with a finite material circuit

The default-off `resolved-electrostatic-circuit-motion` feature exposes
`electrostatic::motion::circuit::PlanarCapacitorCircuitMotion`. It joins the
existing source-bound planar dielectric path, native material RLC circuit,
native articulated body impulse and bounded displacement feedback owners.
The actual field is rebuilt and solved at every changed trial geometry, so
heterogeneous permittivity and fringing inside the authored planar domain
change capacitance, charge/current and mechanical motion together.

For `H(q,s)=q²/(2C(s))`, the circuit uses the harmonic-mean endpoint capacitance.
The symmetric charge/geometry discrete gradient supplies mechanical work
`-(q0²+q1²) delta(1/C)/4`. Its zero-displacement limit uses the original field's
analytic shape derivative. This closes the variable-charge identity as well as
the original fixed-charge limit; replacing the charge-square mean with q0*q1
would violate that identity. Native source work, bulk resistor and contact heat,
unrelated capacitor/mutual-inductor storage and actual body kinetic work remain
separate observations. The physical electrical endpoint account includes the
resolved field once. The native circuit subproblem's `delta_h` uses its auxiliary
constant capacitance and is not the physical changing-geometry storage change.

Construction supplies an initial joint assembly, replacing one unbound scalar
capacitor placeholder while retaining its original charge coordinate. All other
material snapshots, graph topology, contacts, storage laws and inertia remain
fixed during evolution. A competing source-owned capacitor, nonlinear law or
coupled-capacitance block on the active branch refuses. Native constraints must
admit unchanged charge/flux at every coefficient binding; incompatible ideal
clamps refuse without projection. The existing held-voltage motion path serves
direct ideal voltage control. Cancellation, domain exits and failed independent
electrical/mechanical/total work checks return no joint candidate.

The original fixed-thickness, explicitly isotropic, strain-independent dielectric
mesh path and single guided translating body restrictions remain. Heat is reported
for explicit receiving owners, not silently deposited into a thermal state.
This does not add electrostriction, dielectric charge transport, relaxation,
breakdown, general remeshing, 3-D motion, contact events or thermal feedback.
Work admission is not a time/spatial error bound or experimental device validation.

Focused G1/G3/G4 regressions are in `tests/electrostatic_motion/circuit.rs`:
isolated-control parity, independently reduced driven RC/RLC plate equations,
separate passive storage/heat, original branch/material identity, continuing
tapered heterogeneous-field feedback, competing source and clamp refusal, and
late cancellation/retry. Run with
`cargo test -p fs-couple --features resolved-electrostatic-circuit-motion --test electrostatic_motion`.
