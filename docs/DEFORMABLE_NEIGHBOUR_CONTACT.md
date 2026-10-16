# Deformable neighbours and contact-dependent heat exchange

`fs_solid::bridge::neighbours` advances independently deformable deposited
roads under temperature histories, gravity, supplied line loads and explicit
frictionless contact pairs. This is a bounded part of issue #1 R5, not a
whole-print predictor or a new material model.

Supply one `Body` for each road and one `Contact` per reference station pair.
Both roads must already exist and share the same accepted time. The gap is
`reference_gap + w_upper - w_lower`. One solved reaction acts upward on the
upper road and downward on the lower. Neither is replaced by a rigid stop;
a middle road in a stack participates in both adjacent pairs.

The calculation reuses the implicit `LayeredStrip::bending_response` and
`fs_solver::orthant` owner. Signed curvature factors across all roads form
one relative flexibility Gram matrix. Temperature-dependent material state
is probed without aging it, then actual material points advance once. Removal
or lift-off does not erase earlier branch strain or residual stress. Both
road meshes and material/temperature fields may differ.

Bodies, contact rows and total curvature columns share one work allowance.
Actual final gaps and force-times-gap residuals are independently checked.
External dead-load work and internal contact work remain separate and are
reconciled with material storage, thermal/modulus work and physical/numerical
dissipation. Internal endpoint-quadrature work need not vanish on first
engagement or an imposed gap change; it is not heat or a passivity certificate.
`advance_controlled` and executor-scoped `advance_to` return a complete staged
assembly. Late errors and cancellation return no partially updated body set.

## Coupled example

```sh
cargo test -p fs-solid bridge::neighbours
cargo test -p fs-couple --example neighbour_contact_cooling
cargo run -p fs-couple --example neighbour_contact_cooling
```

The example initializes a hot, compliant road over a colder, stiffer road.
It advances actual cell enthalpies, evaluates those temperatures in both
Maxwell beam histories, and solves equal-and-opposite neighbour contact.
A GIVEN finite material patch exchanges heat while the reaction exceeds the
supplied force threshold. An imposed upward line load subsequently causes
lift-off and removes that thermal link, restoring ambient-facing area.

Patch area, total conductance and force threshold are explicit synthetic
inputs. They are not derived from a point reaction. Both conducting and
insulating cases retain the same supplied patch model; zero conductance still
covers its area. No adhesion or healing is inferred. This reduced patch model
requires independent justification and data before use for a real print.

Each call stages heat, the entire mechanical assembly, and thermal-contact
replacement before publication. Only accepted contact affects the next thermal
interval. Thermal and mechanical energy accounts stay separate: no viscous or
impact heat is silently inserted. Eight seconds uses 80 bounded mechanical
steps with thermal subcycling; this is endpoint-staggered integration, not
within-step contact-time localization or an error-bounded coupling algorithm.

## Limits and verification

Contact acts at supplied material-reference stations in parallel x/z frames.
Independent axial motion means this is not general common-world contact
pairing. Point constraints do not exclude penetration between stations, create
finite-area pressures, or discover current contact footprints. There is no
friction, adhesion, large deformation, impact or moving-front deposition.
Built-in pin/clamp restraints remain bilateral. No-pair updates recover
independent beam evolution; clone is an in-memory checkpoint, not a codec.

Six mechanics regression functions and two example regressions are registered.
Native compilation, tests and the example remain unverified here: DSR, its
absolute-path fallback, RCH and Cargo could not start. Independent numerical
references passed 150 constrained systems and the synthetic force/thermal
coupling cases. Those references do not execute Rust or qualify material data.
