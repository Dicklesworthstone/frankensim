# fs-ports

The L3 domain-neutral owner of the existing PortSchema v2, StreamPort v1,
closed-window accounting, relation admission, and scalar coupling primitives.
The implementation is extracted unchanged from fs-couple; that crate reexports
these exact types and remains the owner of assembled physical systems.

Only fs-iface and fs-qty are direct runtime dependencies. Geometry and physics
operators may use the common types without importing acoustic, fluid or thermal
solvers. No new quantities, wire versions, source authority, material laws,
cancellation behavior or whole-system passivity claims are introduced.

The existing fs-couple/tests/couple.rs battery continues to exercise this owner
through its public reexports. fs-couple/tests/port_identity.rs pins type and
error identity across both entry points. FEEC's terminal-relative and integral
batteries consume the direct fs-ports dependency.

## Port-schema adoption contract

`PortedModule` binds a stable port order to the pre-integrated transfers of the
last accepted step, positive into the module. Rejection preserves the preceding
record; before the first accepted step it is `None`. Logical ticks never replace
physical-time integration. Schema admission alone grants no energy closure.

`JunctionSchema` describes 2..=4097 conjugate endpoints using common effort,
common flow, or an exactly skew finite Dirac matrix. Sign pullbacks apply to
both power factors. Pairwise metadata admission is the existing
`ConservativeJunction::new` contract, including orientation, clocks and roles.
The ordinary two-port common-effort instance retains its exact acceptance and
errors. Mixed domains require explicit power-preserving coordinate charts;
dimension relabelling is not a chart. No new solver or wire version is added.

G0 unit tests in `adoption::tests` cover N=2/3/64, cap and cap+1, exact two-port
admission parity, and non-skew/sign/identity refusal. They establish structural
power cancellation, not whole-system passivity or component audit closure.

Coordinate charts accompany canonical schemas. A mass-normalized chart retains
the raw N/sqrt(kg), m sqrt(kg)/s units and an explicit reference-mass conversion;
it does not claim a physical point force or measured modal mass. The signed
volume projection, when present, is retained separately. `ModuleWindowAudit`
distinguishes a real report, no accepted step, and a named unavailable physical
account. Observation errors are visible through `last_window_transfer_error`.
