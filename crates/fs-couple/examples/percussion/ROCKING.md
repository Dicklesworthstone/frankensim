# Small-angle cymbal rocking

`--cymbal-rocking MAXIMUM_ANGLE_RAD` adds two rigid tilts to the existing vertical
translation and complete elastic mode window. It works with `splash`, `hihat`,
and their `-wav` and `-mic` forms. The supplied cap must be finite, positive, and
at most **0.1 radian**. It is an explicit model restriction, not a certified
finite-rotation error bound. The original vertical-only image remains the default.

```sh
cargo run --release -p fs-couple --example percussion -- \
  splash 4096 --analytic-newton --impact-substeps 8 511 \
  --cymbal-rocking 0.03 --strike-position-m 0.068 0.01
```

The same flag enables independent rocking of both hi-hat shells. Their existing
pedal carriage, stand felts, inter-shell contact stations, optional gas storage,
material memory, mallets, flexible sticks, and exterior mute retain their own
physical ports. Each hand can select its target using `--strike-target` or
`--second-stick-target`; see [HIHAT.md](HIHAT.md). A lower-shell strike is measured
in its inward direction, while reported rotation vectors use world axes.

## Physical model

The reduction constructs geometric rigid velocity fields from the supplied
midsurface coordinates. Translational and rotary nodal masses come from the
same shell thickness, facet density, and mass matrix used by the elastic model.
Two fixed-order mass-orthogonalization passes preserve the requested heave/tilt
subspace, including translation/rotation cross inertia on asymmetric specimens.
No ideal-disk inertia or hand-selected rocking frequency is inserted.

The reference origin is the origin of the supplied shell coordinates. It is
an explicit geometric datum; it does not infer the center of mass or a stand
pivot. The image retains constant inertia and first-order rigid kinematics,
`u = translation + rotation × reference_position`. It omits lateral sliding,
finite rotation, centrifugal terms, and gyroscopic coupling. Contact locations,
normal directions, and acoustic boundary charts remain at their declared
reference geometry. This is a bounded approximation for small motion.

A linearized tilt passed directly into Green–Lagrange strain would create an
unphysical quartic restoring spring. The optional reduction first verifies the
rigid geometry, mass orthogonality, free-shell stiffness nullity, and complete
small-signal projection. It then excludes the rigid coordinates from internal
membrane strain and bending energy. Elastic nonlinear strain is retained in
full. Hereditary bending memory uses the same strain-free rigid directions.

Rigid motion still participates in physical point and finite-skin displacement,
contact forces, stand-pad compression, and reciprocal radiation ports. Actual
contact and felt loading therefore resist the tilt. The six stand contacts per
shell and their independent felt histories remain in the joint implicit solve.
The existing acoustic boundary and fitting restrictions continue to apply;
this option does not establish a converged moving-boundary radiation model.

## Accepted motion and output

The time owner checks the norm of each shell's **combined physical rotation
vector** at initial admission and before every accepted endpoint. A state that
exceeds its cap is refused before state, felt history, or time is committed.
Reference, prepared, analytic, and substepped execution share that gate. A
refused output tick retains its scheduled player force for retry. Substeps do
not relax the angle restriction or convert a model-domain refusal into motion.

Mechanical CSV output adds `cymbal_rigid_rotation_x_rad`, `_y_rad`, and `_z_rad`
for a single shell. Hi-hat output uses `upper_rigid_rotation_*_rad` and
`lower_rigid_rotation_*_rad`. These are world-axis rotation-vector components,
not Euler angles. The lower shell uses the same proper half-turn about world x
as its physical placement. Acoustic outputs consume the retained rigid source
coordinates through the existing boundary projection.

The full elastic frequency window is preserved. Three rigid coordinates plus
that window must fit the existing **32 modes per shell** and facet-work limits.
Oversized images refuse; no elastic modes are silently dropped. Native execution,
mesh/time/band convergence, and measured stand and instrument responses are
still required for sound-fidelity or real-time claims.

## Focused regression coverage

`fs-plate/tests/shell_rocking.rs` checks physical inertia on an asymmetric,
nonuniform shell; absence of rigid restoring energy and bending loss; preserved
elastic energy; finite-skin lever arms; contact/radiation work reciprocity; and
refusal of excessive combined angle, forged rigid columns, or a grounded pencil.
The percussion tests cover initial and endpoint refusal across all nonlinear
execution images, retryable player inputs, off-center stand loading, and the
paired hi-hat composition.

These Rust tests were authored but not executed in the editing environment,
which has no working Cargo/rustc or remote Rust runner. Independent numerical
inertia and geometry checks supplement source review; they do not establish a
successful Rust build or native playback.
