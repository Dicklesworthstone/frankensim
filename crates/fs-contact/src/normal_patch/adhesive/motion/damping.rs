//! Passive midpoint damping and optional sealed quadrature history.
use super::{
    CancelGate, JkrNormalMotion, JkrNormalMotionError, JkrNormalMotionState, JkrNormalMotionStep,
    advance, finite,
};

/// One original positive-loss midpoint interval, with its integration weights.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrNormalDissipationSegment {
    offset_s: f64,
    duration_s: f64,
    midpoint_velocity_m_s: f64,
    work_j: f64,
}
impl JkrNormalDissipationSegment {
    pub const fn start_offset_s(&self) -> f64 {
        self.offset_s
    }
    pub fn end_offset_s(&self) -> f64 {
        self.offset_s + self.duration_s
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub const fn midpoint_velocity_m_s(&self) -> f64 {
        self.midpoint_velocity_m_s
    }
    pub const fn dissipated_work_j(&self) -> f64 {
        self.work_j
    }
    /// Original discrete power; downstream transfer still audits its own dt.
    pub fn power_w(&self) -> f64 {
        self.work_j / self.duration_s
    }
}

/// Optional sealed original-step and positive-loss history. Storage is bounded
/// by maximum_substeps (at most 4096 records), allocated fallibly only as used.
/// Existing `advance` retains its allocation-free, Copy step result.
#[derive(Debug, PartialEq)]
pub struct JkrNormalMotionTrace {
    step: JkrNormalMotionStep,
    segments: Vec<JkrNormalDissipationSegment>,
}
impl JkrNormalMotionTrace {
    pub const fn step(&self) -> &JkrNormalMotionStep {
        &self.step
    }
    pub fn segments(&self) -> &[JkrNormalDissipationSegment] {
        &self.segments
    }
}
impl JkrNormalMotion<'_> {
    /// Same mechanics and acceptance gates as `advance`, with sealed positive
    /// viscous-loss intervals for a separately authorized downstream consumer.
    /// Refusal returns neither a step nor any partial dissipation authority.
    pub fn advance_with_dissipation_trace(
        &self,
        previous: JkrNormalMotionState,
        external_force_n: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<JkrNormalMotionTrace, JkrNormalMotionError> {
        let mut segments = Vec::new();
        let step = advance::run(
            self,
            previous,
            external_force_n,
            dt_s,
            gate,
            Some(&mut segments),
        )?;
        Ok(JkrNormalMotionTrace { step, segments })
    }
}

pub(super) fn record(
    segments: &mut Vec<JkrNormalDissipationSegment>,
    maximum: usize,
    offset: f64,
    duration: f64,
    velocity: f64,
    work: f64,
) -> Result<(), JkrNormalMotionError> {
    if work == 0.0 {
        return Ok(());
    }
    if work < 0.0 || duration <= 0.0 || offset + duration <= offset {
        return Err(JkrNormalMotionError::Range("dissipation interval"));
    }
    finite(work / duration, "dissipation source power")?;
    if segments.len() >= maximum {
        return Err(JkrNormalMotionError::TrialBudget);
    }
    if segments.len() == segments.capacity() {
        let additional = segments.capacity().max(8).min(maximum - segments.len());
        segments
            .try_reserve_exact(additional)
            .map_err(|_| JkrNormalMotionError::Range("dissipation trace allocation"))?;
    }
    segments.push(JkrNormalDissipationSegment {
        offset_s: offset,
        duration_s: duration,
        midpoint_velocity_m_s: velocity,
        work_j: work,
    });
    Ok(())
}

pub(super) fn factor(law: &JkrNormalMotion<'_>, h: f64) -> Result<f64, JkrNormalMotionError> {
    if law.damping_n_s_per_m == 0.0 {
        return Ok(0.0);
    }
    finite(
        0.5 * h * (law.damping_n_s_per_m / law.mass_kg),
        "damping time factor",
    )
}
pub(super) fn force(law: &JkrNormalMotion<'_>, velocity: f64) -> Result<f64, JkrNormalMotionError> {
    if law.damping_n_s_per_m == 0.0 {
        return Ok(0.0);
    }
    finite(law.damping_n_s_per_m * velocity, "viscous support force")
}

/// Constitutive force bounds along the admitted quadratic geometry/linear
/// velocity envelope, separate from the held midpoint force used for impulse.
pub(super) fn bounds(
    law: &JkrNormalMotion<'_>,
    from: f64,
    to: f64,
    v0: f64,
    v1: f64,
    h: f64,
    spring_peak: f64,
) -> Result<(f64, f64), JkrNormalMotionError> {
    if law.damping_n_s_per_m == 0.0 {
        return Ok((0.0, spring_peak));
    }
    let viscous0 = force(law, v0)?;
    let viscous1 = force(law, v1)?;
    let viscous_peak = viscous0.abs().max(viscous1.abs());
    let mut combined = finite(
        law.support.force(from)? + viscous0,
        "initial support reaction",
    )?
    .abs()
    .max(finite(law.support.force(to)? + viscous1, "final support reaction")?.abs());
    if law.support.stiffness > 0.0 && v0 != v1 {
        // In unit time u, r(u)=r0+B*u+C*u². The derivative has at most one
        // interior zero, and is checked without numerical differentiation.
        let twice_c = finite(
            law.support.stiffness * (h * (v1 - v0)),
            "support reaction curvature",
        )?;
        let b = finite(
            law.support.stiffness * (h * v0) + law.damping_n_s_per_m * (v1 - v0),
            "support reaction slope",
        )?;
        if twice_c != 0.0 {
            let u = -b / twice_c;
            if u > 0.0 && u < 1.0 {
                let d = finite(
                    from + (u * h) * (v0 + 0.5 * u * (v1 - v0)),
                    "support reaction extremum",
                )?;
                let v = finite(v0 + u * (v1 - v0), "support reaction extremum velocity")?;
                combined = combined.max(
                    finite(
                        law.support.force(d)? + force(law, v)?,
                        "interior support reaction",
                    )?
                    .abs(),
                );
            }
        }
    }
    if viscous_peak.max(combined) > law.budget.maximum_force_n {
        return Err(JkrNormalMotionError::Range(
            "viscous support force along path",
        ));
    }
    Ok((viscous_peak, combined))
}
