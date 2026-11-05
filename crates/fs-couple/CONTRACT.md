
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

## Complete-field finite-opening interior/exterior acoustics

`cavity_modal::interior::exterior::RadiatingInterior` couples the complete
original pressure FE field to the native finite-opening exterior BEM matrix.
It constructs `InteriorMultiport::new_full` (or its original boundary-layer
variant) from the original FE mesh, signed P1 boundary modes and coherent source
gas. Every original pressure DOF remains, subject to the caller's dense cap and
the FE hard ceiling. Original-coordinate field recovery uses the existing
Robin-loaded factorization without an eigenmode cutoff or unloaded-pole
impedance inversion.

Opening triangles must match the retained outward closed exterior exactly and
have planar support. All modes on radiating faces must be included; driven
interior feeds cannot share those faces. Different signed radiation modes may
overlap, including zero-net-flow modes. An exterior P0 mesh that erases an entire
supplied mode refuses rather than silently dropping it. Native BEM velocity
shapes remain dimensionless P1 panel means, so the unchanged source artifact
retains its generalized force/velocity units [N s/m]. The FE reaction matrix is
`conj(Z_native[r,s])/(A_r A_s)` [Pa s/m^3], using TWO geometric support areas,
never signed net-flow gains. All self and mutual entries are retained, and a
resolved native passivity deficit refuses without clipping or matrix repair.

With FE flow positive inward and exterior flow outward, the interface is
`P + Z Q_inward = 0`. Positive real power references only set numerical
coordinates. Source changes reuse the full pressure factor, radiation matrix
and small port closure. Returned amplitudes retain original FE coordinate
order. Source work pairs solved pressure with PRESCRIBED feed flow; radiated
work and native viscous/thermal wall loss are independently evaluated. Original
feed-flow and opening-pressure constraints are checked in power coordinates,
with default relative ceiling 1e-8 and an explicit tolerance entry point.
Port and source-minus-radiation-minus-wall-power residuals remain observable;
no residual is assigned to a physical heat channel.

This is a bounded harmonic model with one coherent homogeneous gas, planar
openings, a rigid remaining exterior and the selected BEM route's existing
mesh/band admission. Pressure, trace and exterior convergence are independent
obligations. Geometry/support, source identity, finite data, explicit work caps,
physical residuals and cooperative cancellation produce refusals. There is no
new dependency direction or external numerical runtime. A causal time-domain
radiation realization and the broader pexfp.99/.103/.104 in-loop acceptance
requirements remain open.

Four native regressions compare against an independently assembled original
CSR FE/BEM mixed block and separately computed BEM matrix; they exercise full
P2 coordinates, mutual loading, live zero-net modes, reference invariance,
finite-radiation versus ideal-release behavior, native wall-loss balance and
binding/refusal paths. Deliberately corrupting the cached closure must fail
original physical constraint admission. Pinned Rust parser/format checks and
independent three-port closure/unit/work algebra passed. Native execution is
still pending: the current environment kills the unchanged upstream dependency
build before these crates. Source checks are not native numerical validation.
