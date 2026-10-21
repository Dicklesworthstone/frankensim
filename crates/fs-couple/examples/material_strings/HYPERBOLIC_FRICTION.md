# Hyperbolic friction in material-string performances

The existing `material_strings` executable can select the shared `fs-tribo`
hyperbolic sliding law. It uses the same implicit sticking/sliding solve,
physical surface velocity (including selected torsion), coupled body, source
checks, work ledger and pressure renderer as other dry friction laws.

For an explicitly declared research input, the record is:

```text
friction hyperbolic MU_STATIC MU_HIGH_SPEED WEAKENING_SPEED_M_S
```

All three values are required. The coefficients satisfy
`0 <= MU_HIGH_SPEED <= MU_STATIC`; the weakening speed is finite and positive.
The kinetic law is
`mu(v) = MU_HIGH_SPEED + (MU_STATIC-MU_HIGH_SPEED)/(1+abs(v)/WEAKENING_SPEED_M_S)`.
At zero slip the solver determines a static reaction up to `MU_STATIC*N`;
it does not substitute a smoothed force or an invented velocity floor.
The physical weakening speed is independent of any ODE regularization.

For an ordered source-backed interface, use one declaration per bow track:

```text
friction interfaces 1
bow-interface bow contact.fsintpk T_K AREA_M2 QUERY_P_PA QUERY_V_M_S hyperbolic STATIC_PROPERTY HIGH_SPEED_PROPERTY SPEED_PROPERTY
```

The last three fields select distinct properties from that interface pack:
two dimensionless coefficients and a positive speed in m/s. Each must be a
source-declared constant with explicit common T/p/v validity. The original
ordered surfaces, history, pack and query receipts remain retained. Surface A
must match the actual exposed string material, including a bonded sleeve.
Missing data, wrong units, unsupported state, or inconsistent surface identity
refuse; the runtime checks actual nominal pressure and solved slip speed.
No body or bulk-string coefficient is changed by this selection.

This is an opt-in, memoryless steady-sliding benchmark, not a thermal rosin,
wear, hair-layer or contact-aging model. It has no stored frictional energy;
accepted traction times actual slip supplies dissipation. No existing violin,
cello, bass or other instrument preset is retuned by selecting this new API.
Measured instrument accuracy requires independent physical validation and
modal/time-resolution studies; this law alone does not establish it.

Focused native checks (not executed in the editing environment):

```sh
cargo test -p fs-tribo --test hyperbolic
cargo test -p fs-couple --lib g1_trait_stick_slip_matches_modal_friction_owner_bitwise
cargo test -p fs-couple --example material_strings hyperbolic
```
