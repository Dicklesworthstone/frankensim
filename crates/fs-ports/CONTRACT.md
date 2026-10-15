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
