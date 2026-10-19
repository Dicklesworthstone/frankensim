# Normal support motion during plucking and bowing

The Rust `MultiBowPerformance::with_stop_motion` extension moves the existing
compliant stops along their admitted fixed normals during a performance. It
uses the same material string/body/cavity state, bow interfaces, load controls,
contact potential, mechanical solver and pressure output. It does not add a
standalone instrument engine or change a string's speaking length or frequency.

```rust,ignore
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::motion::{
    ObstacleMotionLimits, StopMotionKnot, StopMotionProgram,
};

// `performance` already has one admitted material-string stop. At 48 kHz,
// this presses it by 0.1 mm in 10 ms, holds for 10 ms, then withdraws in 10 ms.
let performance = performance.with_stop_motion(StopMotionProgram {
    knots: vec![
        StopMotionKnot { sample: 480, offsets_m: vec![-0.0001] },
        StopMotionKnot { sample: 960, offsets_m: vec![-0.0001] },
        StopMotionKnot { sample: 1440, offsets_m: vec![0.0] },
    ],
    limits: ObstacleMotionLimits {
        maximum_abs_offset_m: 0.001,
        maximum_speed_m_s: 0.1,
    },
    max_knots: 16,
    max_coordinate_values: 256,
}, gate)?;
```

The performance horizon must include every knot. Each knot has one offset for
EVERY stop in source order. Positive offsets open clearance; negative offsets
press inward. `sample` is elapsed physical samples: sample 1 is the endpoint of
the first mechanical step. At sample zero all translations are zero relative
to the original source geometry. Initial indentation belongs in that geometry
and the existing contact-aware preload, not an unaccounted instantaneous move.

Paths are continuous and piecewise linear between their integer-clock knots.
After the last knot the final position is held, not returned automatically to
zero. There is no expansion to a samples-by-stops buffer. Bow and external-force
controls retain their original semantics; the normal endpoint joins their one
joint mechanical step. A callback boundary does not reset any cursor or motion.

The normal contact's relative closure is `N^T q - offset`. Its original potential
secant and internal-loss law see both endpoints of that RELATIVE coordinate.
The applied reaction contributes `-R * delta_offset` of actuator work; positive
work is supplied to the instrument and negative work is extracted. This work is
reported separately from bow work and contact loss, and participates in the
complete energy check before accepted geometry or pressure is published.

`MultiBowStringBank::with_stop_motion` and `step_with_stop_offsets` provide the
same mechanism for direct stepping. Existing `step` calls hold the last accepted
positions. The generic owner lives in `friction::multiple::obstacles::motion`,
so the mechanism is not restricted to musical strings. Failed pressure, contact,
speed, energy or cancellation checks leave accepted positions and motion intact.

The source station, normal and receiver projection remain fixed. A receiver's
own solved normal displacement and the prescribed offset are additive; the
receiver still receives the opposite mechanical reaction. The prescribed offset
is an actuator input, not a solved finger, bow-hair ribbon, contact-area evolution,
station-slide model or changing-material state. Bow normal loads remain separately
prescribed. Contact time/spatial refinement is necessary; conservation alone is
not experimental accuracy or a nonlinear anti-aliasing guarantee.

`StopMotionProgram::validate` exposes the same preflight for scene/file loaders.
It does not evolve the instrument or replace per-step endpoint admission.

This increment is a Rust API. The existing `--stop-profile` file grammar still
specifies static reference geometry; it does not yet accept motion knots.
