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

Closed-window first-law, mass, element, charge and entropy residuals accumulate
the two stored inventories and every signed boundary term in one compensated,
fixed-order reduction. Entropy mechanisms remain separate through that final
reduction and through exchange grouping. A weak heat or diffusive contribution
must not disappear merely because two separately rounded large totals agree.
The existing coupling battery covers one-unit deficits beside 1e16-unit
inventories/streams and cancellation within and across entropy rows. This is a
floating-point conservation check, not an interval bound or a recovery of
precision already lost by the producer of an individual inventory.

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

## Moving control-surface streams

`StreamPort::try_from_uniform_surface` evaluates the original source state on
an oriented uniform face or quadrature point. With material velocity u and cut
velocity w, mass uses `rho*A*(u-w).n`, momentum uses
`mdot*u+A*(p*n-tau*n)`, and the existing enthalpy chart includes
`mdot*(h+|u|²/2+gz)+p*A*w.n-A*(tau*n).u`. The swept pressure term is
separate from deviatoric work and is present at zero-mass moving walls.
Constituents and entropy use the same signed relative mass flow. Reusing the
retained chart refuses replacement of its geometry-derived mass or momentum
rates.

`SurfaceTransport` exposes relative/swept volume flow, advective and stress
momentum, material-side stress power, work on the moving cut, and angular flux
about an explicit origin. Open-stream convection is not silently interpreted
as a body force; slipping-wall work is not inferred heat. Source EOS/stress
validity, heat and species diffusion, spatial quadrature, mesh/GCL advancement
and finite-body evolution remain with their original owners.

G0/G3 native tests cover impermeable-wall impulse/work, Galilean boosts,
uniform-flow geometric conservation, off-origin angular flux and rigid-motion
power, orientation reversal, stationary-chart parity and invalid/rebound
inputs. These are local transport identities, not a whole-system stability
or physical-validation claim.
