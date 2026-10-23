# Felt mounting with heave and two backing tilts

`--mount-spec` can now use one backing with vertical translation and two
principal-axis rotations. Every supplied felt station acts on that same body.
Its signed lever arms return both force and moment through the existing
nonlinear ImpactSystem; material conditioning and Kelvin recovery are not
separately stepped. This is a mechanical connection, not an audio filter.

The shared `felt::mount::modal::ModalMountedFelt` accepts an existing
mass-normalized ImpactBody and complete support-point compression rows. It also
supports caller-supplied elastic support modes without inventing a stand beam,
inertia tensor, material, or separate constitutive law. The example adapter
below is the smaller rigid-body specialization, not a complete stand model.

## Input and units

Keep the existing `frankensim-felt-mount-v1` file and add one optional record:

```text
rocking,com_x_m,com_y_m,Ixx_kg_m2,Iyy_kg_m2,Kx_Nm_rad,Ky_Nm_rad,Rx_Nms_rad,Ry_Nms_rad
```

The eight values are required when the record is present. Inertias are positive
principal moments about the backing assembly's centre of mass; their axes must
be parallel to world X and Y. These are not contact-effective masses or inertia
about the instrument's centre. The existing `carrier` row still supplies the
assembly mass, heave stiffness/resistance, and initial vertical motion. Ground
stiffness and resistance are diagonal and nonnegative. Rotations start at the
reference attitude with zero angular velocity; no hidden angular launch is
introduced. COM horizontal translation and yaw remain constrained.

For a patch at `(x,y)`, the linearized backing height change is

```text
u_z + (y-com_y)*theta_x - (x-com_x)*theta_y
```

The normalized compression row includes those signed lever arms divided by
`sqrt(Ixx)` and `sqrt(Iyy)`. Upper and lower pads have opposite compression
signs. The transposed row supplies the reciprocal reaction and moment; do not
normalize the point row. Kinetic energy is
`(mass*v_z^2 + Ixx*omega_x^2 + Iyy*omega_y^2)/2`.

Reference gaps still come from the actual upper/lower thickness-offset shell
skins at each supplied station. Adding the record does not replace those skins
with a flat midsurface, change loaded areas, duplicate washer creep coefficients,
or reset material history. Omitting the record retains the original translating
mount. Every support coordinate is appended after the existing sticks, shell,
mute and flexible shafts; source, force and observer addresses are preserved.

## Example

`estimated-rocking-mount.fsmount` is an **illustrative input**, not measured
Zildjian hardware or an identified stand. From the workspace root:

```sh
cargo run --release -p fs-couple --example percussion -- --shell-profile crates/fs-couple/examples/percussion/estimated-splash.profile \
  splash 4096 --analytic-newton --impact-substeps 8 511 \
  --cymbal-rocking 0.03 \
  --mount-spec crates/fs-couple/examples/percussion/estimated-rocking-mount.fsmount
```

The same selection is accepted by `splash-wav` and `splash-mic`, including
supplied shell profiles/meshes and existing prepared, substepped and acoustic
feedback paths. The backing affects pressure through shell mechanics; it is not
silently added as a BEM radiating surface. Audio export remains offline.

Mechanical CSV appends physical backing displacement/velocity and, when selected,
`mount_rotation_x_rad`, `mount_angular_velocity_x_rad_s`,
`mount_rotation_y_rad`, and `mount_angular_velocity_y_rad_s`. The observations
read only the appended mechanical prefix, not following felt or acoustic memory.

## Domain and verification

This adapter is a **small-rotation linearization** with fixed vertical patch
projections and reference contact stations. It has no finite-rotation kinematics,
sliding, sleeve/hole collision, full flexible stand, or solved gravity/preload
equilibrium. There is currently **no enforced backing-angle cutoff**. The
`--cymbal-rocking` limit guards the cymbal reduction, not this backing. Inspect
backing angles and restrict use to small excursions; a physical large-angle
support needs a different kinematic owner, not a clamped angular readout.

Five shared-owner regressions cover translational equivalence, force/moment
virtual work, modal motion and energy, rollback, and invalid inputs. Five example
regressions cover the optional specification, principal inertia/lever units,
source-address preservation, torque-driven shell response, and CSV/state layout.
At implementation time these Rust tests were authored but **not compiled or run
locally**, because the environment had no Rust toolchain. Independent Python
physical-equation checks are not a substitute for those native tests, a full
shell/BEM render, specimen identification, listening evaluation or deadline
benchmarks. No calibration, full-band acoustic convergence or real-time claim
is made.
