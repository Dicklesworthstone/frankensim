# Articulated mechanisms on the physical acoustic clock

`fs_couple::pcm_wav::observation::moving::rigid` accepts both standalone rigid
motion and the original `fs_mbd::articulated::PreparedLinkPointMotion`.
The explicit link methods are `RigidAcousticAttachment::kinematics_link`,
`sample_link`, and `RigidAcousticHistory::append_link`.

An attachment used with these methods is measured from the **link origin**, not
its centre of mass. Prepare the link from its matching accepted dynamics and
kinematics receipts. Inverse/prescribed-base forward receipts take their original
gravity vector in `prepare_point_motion`; free-floating receipts already retain
physical gravity and take no additional gravity argument. Reuse one prepared
snapshot for every source and microphone on that link.

The native owner supplies position, velocity and physical Cartesian acceleration,
including joint, Coriolis, centrifugal and angular-acceleration contributions.
The adapter transports explicitly supplied force **on the fluid** with
`F_world = R F_link` and
`Fdot_world = R Fdot_link + omega_world cross F_world`.
The local derivative is the time derivative of the link-coordinate components.
A constant local force can radiate when the link rotates. Volume injection and
its derivative remain scalars; neither gravity nor a motor/support reaction is
automatically an acoustic source.

## Record a coherent accepted frame

The `rigid::articulated` module exports `LinkAcousticUpdate` and
`append_link_frame`. Reserve each `RigidAcousticHistory` before the mechanical
loop, and reserve the scratch once. For two histories, the recording step is:

```rust,ignore
let mut scratch = [None; 2]; // Allocate once outside the accepted-event loop.
// Inside the loop, root and slider are native prepared link snapshots.
append_link_frame(
    accepted_time_s,
    &mut [
        LinkAcousticUpdate {
            history: &mut source_history,
            motion: &slider,
            strength: force_on_fluid_in_slider_axes,
        },
        LinkAcousticUpdate {
            history: &mut microphone_history,
            motion: &root,
            strength: BodyFrameCompactSourceStrength::default(),
        },
    ],
    &mut scratch,
    &mut cancelled,
)?;
```

All paired motion/strength endpoints publish together or none do. This includes
late source failures, exhausted capacity and cancellation after staging. Scratch
may change on failure; every retry overwrites it before publication. Work is
linear in the explicit attachment count, with no append-time allocation,
history cloning or repeated dynamics solve. The final non-interruptible commit
is one paired push per attachment, so callers bound its latency by frame size.

Compile each history's `kinematics()` with the existing
`PreparedKinematicTrajectory`, prepare sources against the same microphone using
`PreparedMovingReceiverSource`, and pass `sources()` through
`CompactSourceHistory` to `MovingCompactRenderer::with_moving_receiver`.
The original coherent sum, Doppler admission, causal coverage, callback clock,
optional decimator and PCM stream remain the only observation pipeline.
Laboratory histories in uniform wind use `PreparedKinematicTrajectory::in_uniform_flow`
with one shared frame; the link adapter does not apply an implicit wind correction.

The source owner still supplies a physical radiation model and a bandwidth that
includes force rotation. Recording endpoints does not establish interpolation
convergence, finite-source compactness or measured SPL. Do not interpolate
across unresolved impacts, topology changes or mismatched mechanical events.
`tests/articulated_acoustic_history.rs` and its `render` submodule exercise native
link motion, transactional recording and the existing pressure/PCM path.
