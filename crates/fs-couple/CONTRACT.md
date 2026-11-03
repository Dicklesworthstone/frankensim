

### Prescribed axial rotation in material film surfaces

`ResolvedImmersedBodies::with_spinning_prescribed_planes` accepts one signed
angular speed in rad/s per unchanged `PlaneFilmSi` source descriptor. Rotation
is about each plane normal through its translating source point. Rates convert
with `Omega_lattice = Omega_SI * timestep_s`, with no length factor; local rigid
wall velocity uses the original cell-size/time conversion. Source SI rates are
retained separately from accepted native plane state. Read-only observations
expose current SI angular speed and local velocity without mutating the owner.

All film resistance, equal/opposite linear and angular impulse, and translational
plus rotational wall work stay inside the original simultaneous native
fluid/body solve. Translation, material rebind and cancellation/retry retain the
same plane spin and accepted clock. Zero rates follow the original prescribed
plane path. Existing `PlaneFilmSi` construction remains source-compatible, and
the spatial fluid-force forwarding API is preserved. No separate force update,
post-step impulse, wall trajectory or finite carrier is introduced.

The resolved fluid boundary must independently match the prescribed motion.
The native positive-gap, local-speed, missing-resistance and asymptotic rules
remain authoritative. This correction does not supply finite axial torsional
resistance, tilting planes or solid contact. The native six-test production-source
suite exercises the coupled resistance/work model. Three material integration
tests cover SI conversion, actuator work and momentum, shared-clock rebind,
zero-rate continuation and refusal/retry; they parsed but were not executed in
the delivery environment. No full material/grid integration pass is claimed.
