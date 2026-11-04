//! Scalar condensation and bounded source-domain trajectory checks.
use super::{
    CancelGate, JkrBranch, JkrNormalMotion, JkrNormalMotionError, JkrNormalMotionState, elastic,
    finite, poll,
};

pub(super) struct Segment {
    pub duration: f64,
    pub indentation: f64,
    pub velocity: f64,
    pub force: f64,
    pub support_force: f64,
    pub support_peak: f64,
    pub travel: f64,
    pub event: bool,
}

pub(super) fn next(
    law: &JkrNormalMotion<'_>,
    state: JkrNormalMotionState,
    load: f64,
    h: f64,
    barrier: Option<f64>,
    gate: &CancelGate,
) -> Result<Segment, JkrNormalMotionError> {
    poll(gate)?;
    let from = state.contact.indentation_m();
    let velocity = state.velocity_m_s;
    if state.contact.branch() == JkrBranch::Detached {
        let touch_acceleration = finite(
            (load - law.support.mean(from, 0.0)?) / law.mass_kg,
            "touch acceleration",
        )?;
        if let Some(duration) = arrival(-from, velocity, touch_acceleration, h, 1.0)? {
            return segment(law, state, load, duration, 0.0, 0.0, true, barrier);
        }
        let to = if law.support.stiffness == 0.0 {
            finite(
                from + h * velocity + 0.5 * h * (h * (load / law.mass_kg)),
                "free approach",
            )?
        } else {
            let acceleration = finite(
                (load - law.support.force(from)?) / law.mass_kg,
                "supported flight acceleration",
            )?;
            let denominator = finite(
                1.0 + 0.25 * h * (h * (law.support.stiffness / law.mass_kg)),
                "supported flight denominator",
            )?;
            finite(
                from + (h * velocity + 0.5 * h * (h * acceleration)) / denominator,
                "supported flight approach",
            )?
        };
        return segment(law, state, load, h, to, 0.0, false, barrier);
    }
    let fold = law.normal.minimum_attached_indentation_m();
    // Only a monotone outward segment may jump directly to this boundary.
    // A segment that first reverses is resolved before testing the fold again.
    if velocity <= 0.0 {
        let force = law.normal.attached_path_average(from, fold)?;
        let total_force = finite(force + law.support.mean(from, fold)?, "fold total reaction")?;
        let acceleration = finite((load - total_force) / law.mass_kg, "fold acceleration")?;
        if let Some(duration) = arrival(fold - from, velocity, acceleration, h, -1.0)? {
            let derivative =
                law.normal.attached_path_end_derivative(from, fold)? + 0.5 * law.support.stiffness;
            let slope = finite(
                1.0 + 0.5 * duration * (duration * (derivative / law.mass_kg)),
                "fold root slope",
            )?;
            if slope <= 0.0 {
                return Err(JkrNormalMotionError::Unresolved(
                    "fold step needs refinement",
                ));
            }
            return segment(law, state, load, duration, fold, force, true, barrier);
        }
    }
    let force0 = law.normal.evaluate_attached(from)?.normal_force_n;
    let total_force0 = finite(force0 + law.support.force(from)?, "initial normal reaction")?;
    let residual = |to: f64| -> Result<(f64, f64), JkrNormalMotionError> {
        let force = law.normal.attached_path_average(from, to)?;
        let total_force = finite(force + law.support.mean(from, to)?, "mean normal reaction")?;
        Ok((
            finite(
                (to - from) - h * velocity + 0.5 * h * (h * ((total_force - load) / law.mass_kg)),
                "normal root residual",
            )?,
            force,
        ))
    };
    let initial = finite(
        -h * velocity + 0.5 * h * (h * ((total_force0 - load) / law.mass_kg)),
        "initial normal residual",
    )?;
    let initial_derivative =
        law.normal.attached_path_end_derivative(from, from)? + 0.5 * law.support.stiffness;
    let initial_slope = finite(
        1.0 + 0.5 * h * (h * (initial_derivative / law.mass_kg)),
        "initial normal root slope",
    )?;
    if initial_slope <= 0.0 {
        return Err(JkrNormalMotionError::Unresolved(
            "normal step needs refinement",
        ));
    }
    if initial.abs() <= root_tolerance(law, from.abs()) {
        return segment(law, state, load, h, from, force0, false, barrier);
    }
    // F'(delta) increases along the original JKR branch. The residual is
    // convex; these opposite-sign brackets select its upper, positive-slope
    // root connected to the small-step solution. An unresolved bracket halves.
    let (mut lo, mut hi) = if initial < 0.0 {
        let hi = law.normal.envelope().maximum_indentation_m;
        if residual(hi)?.0 < 0.0 {
            return Err(JkrNormalMotionError::Unresolved("normal upper bracket"));
        }
        (from, hi)
    } else {
        if residual(fold)?.0 > 0.0 {
            return Err(JkrNormalMotionError::Unresolved("normal lower bracket"));
        }
        (fold, from)
    };
    let predictor = from - initial;
    let mut to = if predictor > lo && predictor < hi {
        predictor
    } else {
        lo.midpoint(hi)
    };
    for _ in 0..law.budget.maximum_root_iterations {
        poll(gate)?;
        let (value, force) = residual(to)?;
        let tolerance = root_tolerance(law, from.abs().max(to.abs()).max((h * velocity).abs()));
        let derivative =
            law.normal.attached_path_end_derivative(from, to)? + 0.5 * law.support.stiffness;
        let slope = finite(
            1.0 + 0.5 * h * (h * (derivative / law.mass_kg)),
            "normal root slope",
        )?;
        if value.abs() <= tolerance && slope > 0.0 {
            return segment(law, state, load, h, to, force, to == fold, barrier);
        }
        if value < 0.0 {
            lo = to;
        } else {
            hi = to;
        }
        let newton = to - value / slope;
        let next = if slope > 0.0 && newton > lo && newton < hi {
            newton
        } else {
            lo.midpoint(hi)
        };
        if next == to {
            return Err(JkrNormalMotionError::Unresolved("normal root precision"));
        }
        to = next;
    }
    Err(JkrNormalMotionError::Unresolved("normal root iterations"))
}

// Resolve a local root more tightly than the whole-call coordinate allowance:
// many accepted substeps contribute to the independent native moment check.
// The user tolerance remains a hard upper bound, including below roundoff.
fn root_tolerance(law: &JkrNormalMotion<'_>, scale: f64) -> f64 {
    law.coordinate_tolerance(scale)
        .min(8.0 * f64::EPSILON * scale.max(law.normal.indentation_scale_m()))
}

// Solve the constant-acceleration boundary equation in a unit time interval.
// Scaling all three displacement coefficients avoids squaring a dimensional
// velocity. The q-form retains the small root without cancellation.
fn arrival(
    distance: f64,
    velocity: f64,
    acceleration: f64,
    horizon: f64,
    direction: f64,
) -> Result<Option<f64>, JkrNormalMotionError> {
    let linear = finite(horizon * velocity, "event flight")?;
    let quadratic = finite(0.5 * horizon * (horizon * acceleration), "event curvature")?;
    let scale = linear.abs().max(quadratic.abs()).max(distance.abs());
    if scale == 0.0 {
        return Ok(None);
    }
    let (a, b, c) = (quadratic / scale, linear / scale, -distance / scale);
    let roots = if a == 0.0 {
        if b == 0.0 {
            return Ok(None);
        }
        [-c / b, f64::INFINITY]
    } else {
        let discriminant = b * b - 4.0 * a * c;
        if discriminant < 0.0 {
            return Ok(None);
        }
        let q = -0.5 * (b + discriminant.sqrt().copysign(b));
        if q == 0.0 {
            return Ok(None);
        }
        [q / a, c / q]
    };
    let mut earliest = None;
    for fraction in roots {
        if fraction > 0.0 && fraction <= 1.0 {
            let time = horizon * fraction;
            if time == 0.0 {
                return Err(JkrNormalMotionError::Unresolved("event duration underflow"));
            }
            let end_velocity = finite(velocity + time * acceleration, "event velocity")?;
            if direction * end_velocity >= 0.0 && earliest.is_none_or(|old| time < old) {
                earliest = Some(time);
            }
        }
    }
    Ok(earliest)
}

#[allow(clippy::too_many_arguments)]
fn segment(
    law: &JkrNormalMotion<'_>,
    state: JkrNormalMotionState,
    load: f64,
    h: f64,
    to: f64,
    force: f64,
    event: bool,
    barrier: Option<f64>,
) -> Result<Segment, JkrNormalMotionError> {
    let from = state.contact.indentation_m();
    let support_force = law.support.mean(from, to)?;
    let reaction = finite(force + support_force, "total normal reaction")?;
    let acceleration = finite((load - reaction) / law.mass_kg, "normal acceleration")?;
    let velocity = finite(state.velocity_m_s + h * acceleration, "normal velocity")?;
    if velocity.abs() > law.budget.maximum_normal_speed_m_s {
        return Err(JkrNormalMotionError::Range("normal speed"));
    }
    let residual = finite(
        (to - from) - 0.5 * h * (state.velocity_m_s + velocity),
        "normal kinematic residual",
    )?;
    if residual.abs() > law.coordinate_tolerance(from.abs().max(to.abs())) {
        return Err(JkrNormalMotionError::Unresolved(
            "normal kinematic precision",
        ));
    }
    let (mut low, mut high) = (from.min(to), from.max(to));
    let mut travel = (to - from).abs();
    if acceleration != 0.0 {
        let turn = -state.velocity_m_s / acceleration;
        if turn > 0.0 && turn < h {
            let extremum = finite(
                from + 0.5 * turn * state.velocity_m_s,
                "normal turning approach",
            )?;
            low = low.min(extremum);
            high = high.max(extremum);
            travel = finite(
                (extremum - from).abs() + (to - extremum).abs(),
                "normal path length",
            )?;
        }
    }
    if low < law.budget.minimum_indentation_m {
        return Err(JkrNormalMotionError::Range("minimum normal approach"));
    }
    let support_peak = law
        .support
        .force(low)?
        .abs()
        .max(law.support.force(high)?.abs());
    if support_peak > law.budget.maximum_force_n {
        return Err(JkrNormalMotionError::Range("support force along path"));
    }
    if state.contact.branch() == JkrBranch::Detached {
        if high > 0.0 || (!event && high == 0.0) {
            return Err(JkrNormalMotionError::Unresolved("unresolved contact entry"));
        }
    } else {
        let fold = law.normal.minimum_attached_indentation_m();
        if low < fold || (!event && low == fold) {
            return Err(JkrNormalMotionError::Unresolved("unresolved separation"));
        }
        let low_force = law.normal.evaluate_attached(low)?.normal_force_n;
        let high_force = law.normal.evaluate_attached(high)?.normal_force_n;
        let x = 3.0_f64.cbrt();
        let force_minimum = law.normal.displacement_scale_m * (x * x * x * x - 4.0 * x);
        let peak = low_force.abs().max(high_force.abs()).max(
            if low <= force_minimum && high >= force_minimum {
                law.normal.force_controlled_pull_off_n().abs()
            } else {
                0.0
            },
        );
        if peak > law.budget.maximum_force_n {
            return Err(JkrNormalMotionError::Range("normal force along path"));
        }
        elastic::check_barrier(law, state, load, low, high, barrier)?;
    }
    Ok(Segment {
        duration: h,
        indentation: to,
        velocity,
        force,
        support_force,
        support_peak,
        travel,
        event,
    })
}
