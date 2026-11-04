# Moving points in rotating carriers

Both `fs_mbd::PreparedRigidPointMotion` and
`fs_mbd::articulated::PreparedLinkPointMotion` expose `relative_point_motion`.
Supply the point's current local coordinates and their first and second time
derivatives in those rotating coordinates, all at the same accepted instant.
The standalone body uses its centre-of-mass origin; the articulated query uses
its original link origin. An existing child link already owns its joint motion:
do not add that motion again as an extra relative input.

For `x = X + R r`, the native fixed-point query supplies translation, angular
acceleration and centripetal acceleration. The relative query adds
`R r_dot` to velocity and `2 omega cross (R r_dot) + R r_ddot` to acceleration.
It neither solves the relative motion nor changes mass, reaction forces,
gravity or the accepted mechanical state. Zero relative derivatives preserve
the original fixed-point result. Queries allocate no memory.

## Migrating acoustic sources and microphones

`fs_couple::pcm_wav::observation::moving::rigid::RelativeAcousticMotion` holds
local displacement from the history's original reference attachment plus its
local velocity and acceleration. Use `sample_relative` / `sample_link_relative`
for individual events and `append_relative` / `append_link_relative` for the
existing paired, preallocated histories. Microphones can request only
`kinematics_relative` / `kinematics_link_relative`, or record zero source strength.

For a coherent multi-link event, the existing `LinkAcousticUpdate` bindings and
scratch buffer are accepted by `append_link_frame_with_relative_motion`. Supply
one relative-motion record per binding, with zero records for fixed points.
The original stage-then-publish implementation owns both paths: malformed
counts, invalid motion, source failures, exhausted histories and cancellation
publish no partial frame. Retrying overwrites every scratch slot before use.

The returned motion and strength slices enter the same trajectory preparation,
retarded-time solver, Doppler checks and pressure/PCM renderer as fixed
attachments. No separate integration or resampling clock is introduced. Fixed
and moving sources share the original force-rotation calculation; point motion
is queried once per accepted sample, not once per source component.

Source forces remain explicit forces **on the fluid**, expressed in the carrier
axes. The supplied local force derivative must be the **total derivative along
the source path**, including any spatial variation encountered by migration.
The adapter adds carrier-axis rotation, not an inferred force gradient or a new
source orientation. A migrating geometric contact location is not the velocity
of either material surface: do not substitute it for material slip in friction
or contact work. Finite-source compactness, motion/strength sampling convergence,
source bandwidth, smoothness and radiation reaction retain their existing
owners. Never interpolate through an unresolved impact or topology change.
