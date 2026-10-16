# CONTRACT: fs-otto-wasm

## Purpose and layer

Layer: **L6 WASM boundary**. This crate exposes the source-bounded US 194,047
mechanism composition owned by `fs-mbd::otto` to browser consumers. It performs
admission and stable JSON-envelope serialization only. Generic revolute and
prismatic joints plus slider-crank closure remain in `fs-mbd`.

## Public boundary

`otto_topology_step(crank_angle_rad, crank_radius, connecting_rod_length,
engine_rpm) -> String` returns the complete connected pose: joint axes,
four-stroke phase, crank/wrist pin coordinates, rod angle and closure span,
half-speed side-shaft angle, admission-slide coordinate, exhaust lift, and a
normalized display governor pose.

Success returns `{"ok":{...}}`. Refusal returns
`{"refusal":{"code", "message", "repairs"}}` and never a partial success.

## Invariants and refusal

- All inputs must be finite; geometry must be positive and the rod must be
  longer than the crank radius.
- Engine speed is admitted only in `[0, 600]` rpm, the museum control domain.
- The crank-to-wrist distance is recomputed on every call and equals the
  caller-supplied rod length within binary64 roundoff.
- One independent crank drive determines all eight scalar joint coordinates.

## Cancellation and resources

The call is synchronous, constant-time, performs no I/O, starts no threads or
tasks, and has no partial commit. Cancellation is the call boundary.

## Unsafe boundary and dependencies

No unsafe code. The crate depends only on `fs-mbd` plus target-specific
`wasm-bindgen`; it owns no duplicate mechanism law.

## No-claim boundaries

The patent prints no construction dimensions, mass, inertia, speed, load,
pressure, valve lift, governor setting, or efficiency. Caller dimensions are
display-scale geometry, governor spread is normalized presentation state, and
this boundary does not claim combustion CFD, bearing contact, friction,
vibration, stress, wear, or a historically measured operating point.
