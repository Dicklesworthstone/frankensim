//! Continuous coupled dry friction under held physical forces.
//!
//! The declared law is `dv/dt = W*(drive + friction)`, with
//! `friction_i = -capacity_i * sign(v_i)` while moving and the closed interval
//! `[-capacity_i, capacity_i]` at rest. Static and kinetic bounds are equal.
//! For constant positive-definite `W`, drives and capacities, its unique
//! velocity path is piecewise affine. Each segment ends at the first actual
//! zero crossing; all resting ports resolve their coupled reactions together.
//! No time discretization or velocity capture band replaces those events.
//!
//! Travel, source work and dissipation are integrated over every segment.
//! An integrated port impulse paired with the two outer endpoint velocities
//! does not give that port's continuous friction work when other modes couple
//! to it. Only the assembled endpoint pairing reconstructs the kinetic change.
//! Geometry, normal loading, unequal static/kinetic history and heat partition
//! remain with their physical owners.

/// Distinct static/kinetic bounds with explicit coupled mode admissibility.
#[path = "relaxation/stiction.rs"]
pub mod stiction;

use super::{
    CoupledFrictionImpulseControls, CoupledFrictionImpulseError, CoupledFrictionImpulseInput,
    Factor, MAX_FRICTION_IMPULSE_PORTS, Vector, normalized, positive,
    solve_coupled_friction_impulses,
};
use fs_math::dd::Dd;

type Accurate = [Dd; MAX_FRICTION_IMPULSE_PORTS];

/// Work-conjugate scalar ports of one unchanged physical response.
#[derive(Clone, Copy, Debug)]
pub struct CoupledFrictionFlowInput<'a> {
    /// Initial velocities [m/s]. A nonzero value is moving, however small.
    pub velocity_before_m_s: &'a [f64],
    /// Complete row-major `n*n` positive-definite mobility [1/kg].
    pub mobility_per_kg: &'a [f64],
    /// Held, nonnegative dry-friction force capacities [N].
    pub capacity_n: &'a [f64],
    /// Held signed applied forces, in the same port coordinates [N].
    pub drive_n: &'a [f64],
    /// Nonnegative physical duration [s]; zero is an admitted identity.
    pub duration_s: f64,
}

/// Bounded event resolution and original-system physical acceptance gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledFrictionFlowControls {
    /// Existing coupled projection controls. Matrix admission and iteration
    /// bounds retain their meanings. Velocity, impulse and energy tolerances
    /// additionally gate the complete original continuous endpoint; the impulse
    /// budget bounds actual integrated drive, friction and resultant impulses.
    pub projection: CoupledFrictionImpulseControls,
    /// Maximum affine segments, in `1..=4096`. A refusal never enlarges a step.
    pub maximum_segments: usize,
    /// Original force-balance allowance on a declared resting coordinate [m/s²].
    /// This is not a threshold for capturing a moving coordinate.
    pub acceleration_tolerance_m_s2: f64,
    /// Maximum actual absolute drive, friction and resultant port force [N].
    /// The original physical friction intervals are never reduced to this bound.
    pub maximum_force_n: f64,
}

impl CoupledFrictionFlowControls {
    pub fn validate(self) -> Result<(), CoupledFrictionFlowError> {
        self.projection.validate()?;
        if !(1..=4096).contains(&self.maximum_segments)
            || !positive(self.acceleration_tolerance_m_s2)
            || !self.maximum_force_n.is_finite()
            || self.maximum_force_n < 0.0
        {
            return Err(CoupledFrictionFlowError::Invalid("flow controls"));
        }
        Ok(())
    }
}

/// Continuous path candidate. Array tails beyond `count` remain zero.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct CoupledFrictionFlowReport {
    pub count: usize,
    /// Integrated actual friction force [N s].
    pub friction_impulse_ns: Vector,
    /// Integrated held drive [N s].
    pub drive_impulse_ns: Vector,
    /// Actual `v_before + W*(friction_impulse + drive_impulse)` [m/s].
    /// Native mechanics should apply these same impulses and check this endpoint.
    pub velocity_after_m_s: Vector,
    /// Piecewise-affine path endpoint [m/s]. Only originally resting ports,
    /// or ports reaching an isolated stopping event, acquire an exact zero.
    /// Floating native impulse response remains separately visible above.
    pub flow_velocity_after_m_s: Vector,
    /// Original impulse-response endpoint minus the analytical path [m/s].
    pub velocity_residual_m_s: Vector,
    /// Independent `W^-1*(v_after-v_before)` minus actual total impulse [N s].
    pub impulse_residual_ns: Vector,
    /// Positive rounding excess above `capacity*duration`, otherwise zero [N s].
    /// Actual impulses are reported unchanged and this residual is gated.
    pub impulse_capacity_residual_ns: Vector,
    /// Signed integral of the actual piecewise-affine port velocity [m].
    pub displacement_m: Vector,
    /// Integral of the absolute port velocity, including every reversal [m].
    pub distance_m: Vector,
    /// Smallest absolute speed on the whole affine path, including every stop [m/s].
    pub minimum_speed_m_s: Vector,
    /// Largest absolute speed on the whole affine path, including internal events [m/s].
    pub maximum_speed_m_s: Vector,
    /// Nonpositive continuous friction work, separately integrated per port [J].
    pub friction_work_j: Vector,
    /// Signed continuous source work, separately integrated per port [J].
    pub drive_work_j: Vector,
    /// Nonnegative total loss `-sum(friction_work_j)` [J]; no thermal allocation.
    pub dissipated_j: f64,
    pub total_drive_work_j: f64,
    /// Independently reconstructed in the original inverse-mobility metric [J].
    pub kinetic_energy_change_j: f64,
    /// Kinetic change minus continuous drive and friction work [J].
    pub energy_residual_j: f64,
    pub segments: usize,
    /// Number of individual port stopping events, including simultaneous ones.
    pub stops: usize,
    pub projection_iterations: usize,
    pub factorizations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum CoupledFrictionFlowError {
    Invalid(&'static str),
    Projection(CoupledFrictionImpulseError),
    ForceBudget {
        index: usize,
        force_n: f64,
        maximum_n: f64,
    },
    ImpulseBudget {
        index: usize,
        impulse_ns: f64,
        maximum_ns: f64,
    },
    SegmentLimit {
        segments: usize,
        elapsed_s: f64,
    },
    Unresolved(&'static str),
    EnergyBalance {
        dissipated_j: f64,
        residual_j: f64,
    },
    Cancelled,
}

impl From<CoupledFrictionImpulseError> for CoupledFrictionFlowError {
    fn from(error: CoupledFrictionImpulseError) -> Self {
        match error {
            CoupledFrictionImpulseError::Cancelled => Self::Cancelled,
            _ => Self::Projection(error),
        }
    }
}
impl core::fmt::Display for CoupledFrictionFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "continuous coupled friction: {self:?}")
    }
}
impl std::error::Error for CoupledFrictionFlowError {}

/// Advance the original coupled dry-friction inclusion for a physical duration.
///
/// Moving ports retain their exact opposing force. At rest, the same bounded
/// coupled owner resolves the simultaneous force balance in the principal
/// mobility. Its work-conjugate impulse problem uses the actual remaining
/// duration, not a fictitious mass or a diagonal replacement. Resting ports
/// whose capacity cannot hold leave in the acceleration selected by that solve.
/// The first exact zero crossing ends each affine segment.
///
/// All arrays, factors and path sums have bounded stack storage. Cancellation,
/// singular/dependent mobility, unrepresentable events and exhausted work or
/// physical budgets refuse without changing borrowed input or returning a
/// partial path. Success retains actual original-matrix endpoint and kinetic
/// audits, including numerical resting/event residuals.
pub fn advance_coupled_friction(
    input: &CoupledFrictionFlowInput<'_>,
    controls: CoupledFrictionFlowControls,
    cancelled: impl FnMut() -> bool,
) -> Result<CoupledFrictionFlowReport, CoupledFrictionFlowError> {
    advance_with_law(input, controls, &mut EqualBounds, cancelled)
}

// One event/path owner for both physical force laws. The equal-bound branch
// retains its original acceleration selection and all arithmetic unchanged.
struct PhaseForces {
    friction: Vector,
    // None selects the original equal-bound projection rule. Some declares the
    // already resolved resting coordinates; their actual residual is still gated.
    holding_mask: Option<u32>,
}
trait PhaseLaw {
    type Error: From<CoupledFrictionFlowError> + From<CoupledFrictionImpulseError>;
    #[allow(clippy::too_many_arguments)]
    fn resolve<C: FnMut() -> bool>(
        &mut self,
        input: &CoupledFrictionFlowInput<'_>,
        controls: CoupledFrictionFlowControls,
        velocity: &Accurate,
        scales: &Vector,
        remaining: Dd,
        sums: &mut Path,
        cancelled: &mut C,
    ) -> Result<PhaseForces, Self::Error>;
}
struct EqualBounds;
impl PhaseLaw for EqualBounds {
    type Error = CoupledFrictionFlowError;
    fn resolve<C: FnMut() -> bool>(
        &mut self,
        input: &CoupledFrictionFlowInput<'_>,
        controls: CoupledFrictionFlowControls,
        velocity: &Accurate,
        _scales: &Vector,
        remaining: Dd,
        sums: &mut Path,
        cancelled: &mut C,
    ) -> Result<PhaseForces, Self::Error> {
        Ok(PhaseForces {
            friction: forces(input, controls, velocity, remaining, sums, cancelled)?,
            holding_mask: None,
        })
    }
}

#[allow(clippy::too_many_lines)]
fn advance_with_law<L: PhaseLaw>(
    input: &CoupledFrictionFlowInput<'_>,
    controls: CoupledFrictionFlowControls,
    law: &mut L,
    mut cancelled: impl FnMut() -> bool,
) -> Result<CoupledFrictionFlowReport, L::Error> {
    poll(&mut cancelled)?;
    controls.validate()?;
    let n = input.velocity_before_m_s.len();
    if n > MAX_FRICTION_IMPULSE_PORTS
        || input.mobility_per_kg.len() != n * n
        || input.capacity_n.len() != n
        || input.drive_n.len() != n
        || !input.duration_s.is_finite()
        || input.duration_s < 0.0
    {
        return Err((CoupledFrictionFlowError::Invalid("flow dimensions or duration")).into());
    }
    let mut scales = [0.0; MAX_FRICTION_IMPULSE_PORTS];
    let mut velocity = [Dd::ZERO; MAX_FRICTION_IMPULSE_PORTS];
    for i in 0..n {
        poll(&mut cancelled)?;
        if !input.velocity_before_m_s[i].is_finite()
            || !input.drive_n[i].is_finite()
            || !input.capacity_n[i].is_finite()
            || input.capacity_n[i] < 0.0
            || !positive(input.mobility_per_kg[i * n + i])
        {
            return Err((CoupledFrictionFlowError::Invalid("flow port data")).into());
        }
        scales[i] = input.mobility_per_kg[i * n + i].sqrt();
        velocity[i] = Dd::from_f64(input.velocity_before_m_s[i]);
        force_budget(i, input.drive_n[i], controls)?;
    }
    for i in 0..n {
        poll(&mut cancelled)?;
        for j in 0..i {
            let a = normalized(input.mobility_per_kg[i * n + j], scales[i], scales[j])?;
            let b = normalized(input.mobility_per_kg[j * n + i], scales[j], scales[i])?;
            if (a - b).abs() > controls.projection.relative_matrix_tolerance {
                return Err((CoupledFrictionFlowError::Invalid("mobility symmetry")).into());
            }
        }
    }
    let zeros = [0.0; MAX_FRICTION_IMPULSE_PORTS];
    let original = CoupledFrictionImpulseInput {
        velocity_before_m_s: input.velocity_before_m_s,
        mobility_per_kg: input.mobility_per_kg,
        capacity_ns: &zeros[..n],
    };
    let all = if n == 32 { u32::MAX } else { (1_u32 << n) - 1 };
    let mut factor = Factor::empty();
    factor.prepare(&original, &scales, all, controls.projection, &mut cancelled)?;
    let mut sums = Path::new(n);
    for i in 0..n {
        sums.minimum_speed[i] = input.velocity_before_m_s[i].abs();
        sums.maximum_speed[i] = input.velocity_before_m_s[i].abs();
    }
    let duration = Dd::from_f64(input.duration_s);
    let mut elapsed = Dd::ZERO;
    while n != 0 && elapsed.lt(duration) {
        poll(&mut cancelled)?;
        if sums.segments == controls.maximum_segments {
            return Err((CoupledFrictionFlowError::SegmentLimit {
                segments: sums.segments,
                elapsed_s: elapsed.to_f64(),
            })
            .into());
        }
        let remaining = checked(duration - elapsed, "remaining time")?;
        let phase = law.resolve(
            input,
            controls,
            &velocity,
            &scales,
            remaining,
            &mut sums,
            &mut cancelled,
        )?;
        let force = phase.friction;
        let mut acceleration = [Dd::ZERO; MAX_FRICTION_IMPULSE_PORTS];
        for i in 0..n {
            poll(&mut cancelled)?;
            let raw = response_row(input, i, &force)?;
            acceleration[i] = if phase
                .holding_mask
                .is_some_and(|mask| mask & (1_u32 << i) != 0)
            {
                if velocity[i] != Dd::ZERO
                    || raw.abs().to_f64() > controls.acceleration_tolerance_m_s2
                {
                    return Err(CoupledFrictionFlowError::Unresolved(
                        "declared resting force balance",
                    )
                    .into());
                }
                Dd::ZERO
            } else if phase.holding_mask.is_none()
                && velocity[i] == Dd::ZERO
                && input.capacity_n[i] != 0.0
            {
                // Rest is an original state or an isolated mathematical event.
                // The interior force face imposes zero acceleration, with its
                // actual original force-balance error retained by final audits.
                let cap = input.capacity_n[i];
                let outward = (force[i] == -cap && raw.to_f64() >= 0.0)
                    || (force[i] == cap && raw.to_f64() <= 0.0);
                if outward {
                    raw
                } else if raw.abs().to_f64() <= controls.acceleration_tolerance_m_s2 {
                    Dd::ZERO
                } else {
                    return Err(
                        (CoupledFrictionFlowError::Unresolved("resting force balance")).into(),
                    );
                }
            } else {
                raw
            };
        }
        let mut step = remaining;
        let mut stopping = 0_u32;
        for i in 0..n {
            if velocity[i] == Dd::ZERO
                || acceleration[i] == Dd::ZERO
                || velocity[i].to_f64().is_sign_positive()
                    == acceleration[i].to_f64().is_sign_positive()
            {
                continue;
            }
            // An overflowing positive crossing time lies beyond every admitted
            // finite duration; do not manufacture an earlier numerical event.
            if (velocity[i].to_f64() / acceleration[i].to_f64()).is_infinite() {
                continue;
            }
            let crossing = checked(-velocity[i] / acceleration[i], "stopping time")?;
            if crossing.to_f64() <= 0.0 {
                return Err(
                    (CoupledFrictionFlowError::Unresolved("representable stopping time")).into(),
                );
            }
            if crossing.lt(step) {
                step = crossing;
                stopping = 1_u32 << i;
            } else if crossing == step {
                stopping |= 1_u32 << i;
            }
        }
        let next_time = if step == remaining {
            duration
        } else {
            checked(elapsed + step, "elapsed time")?
        };
        if step.to_f64() <= 0.0 || !elapsed.lt(next_time) || duration.lt(next_time) {
            return Err((CoupledFrictionFlowError::Unresolved("physical time progress")).into());
        }
        for i in 0..n {
            poll(&mut cancelled)?;
            let raw_after = checked(velocity[i] + acceleration[i] * step, "segment velocity")?;
            let after = if stopping & (1_u32 << i) != 0 {
                if raw_after.abs().to_f64() > controls.projection.velocity_tolerance_m_s {
                    return Err(
                        (CoupledFrictionFlowError::Unresolved("stopping event residual")).into(),
                    );
                }
                Dd::ZERO
            } else {
                if velocity[i] != Dd::ZERO
                    && raw_after != Dd::ZERO
                    && velocity[i].to_f64().is_sign_positive()
                        != raw_after.to_f64().is_sign_positive()
                {
                    return Err(
                        (CoupledFrictionFlowError::Unresolved("stopping event order")).into(),
                    );
                }
                raw_after
            };
            let midpoint = velocity[i] * Dd::from_f64(0.5) + after * Dd::from_f64(0.5);
            let travel = checked(midpoint * step, "segment travel")?;
            let friction_work = checked(Dd::from_f64(force[i]) * travel, "friction work")?;
            if friction_work.to_f64() > 0.0 {
                return Err((CoupledFrictionFlowError::Unresolved(
                    "continuous friction passivity",
                ))
                .into());
            }
            sums.displacement[i] = checked(sums.displacement[i] + travel, "signed path")?;
            sums.distance[i] = checked(sums.distance[i] + travel.abs(), "absolute path")?;
            sums.friction_impulse[i] = checked(
                sums.friction_impulse[i] + Dd::from_f64(force[i]) * step,
                "friction impulse",
            )?;
            sums.friction_work[i] =
                checked(sums.friction_work[i] + friction_work, "friction loss")?;
            let endpoint_speed = after.abs().to_f64();
            sums.minimum_speed[i] = sums.minimum_speed[i].min(endpoint_speed);
            sums.maximum_speed[i] = sums.maximum_speed[i].max(endpoint_speed);
            velocity[i] = after;
        }
        sums.segments += 1;
        sums.stops += stopping.count_ones() as usize;
        elapsed = next_time;
    }
    finish(
        input,
        controls,
        &factor,
        &scales,
        &velocity,
        &sums,
        &mut cancelled,
    )
    .map_err(Into::into)
}

struct Path {
    displacement: Accurate,
    distance: Accurate,
    friction_impulse: Accurate,
    friction_work: Accurate,
    minimum_speed: Vector,
    maximum_speed: Vector,
    segments: usize,
    stops: usize,
    projection_iterations: usize,
    factorizations: usize,
}
impl Path {
    fn new(count: usize) -> Self {
        Self {
            displacement: [Dd::ZERO; MAX_FRICTION_IMPULSE_PORTS],
            distance: [Dd::ZERO; MAX_FRICTION_IMPULSE_PORTS],
            friction_impulse: [Dd::ZERO; MAX_FRICTION_IMPULSE_PORTS],
            friction_work: [Dd::ZERO; MAX_FRICTION_IMPULSE_PORTS],
            minimum_speed: [0.0; MAX_FRICTION_IMPULSE_PORTS],
            maximum_speed: [0.0; MAX_FRICTION_IMPULSE_PORTS],
            segments: 0,
            stops: 0,
            projection_iterations: 0,
            factorizations: usize::from(count != 0),
        }
    }
}

fn forces(
    input: &CoupledFrictionFlowInput<'_>,
    controls: CoupledFrictionFlowControls,
    velocity: &Accurate,
    remaining: Dd,
    sums: &mut Path,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vector, CoupledFrictionFlowError> {
    let n = input.velocity_before_m_s.len();
    let mut force = [0.0; MAX_FRICTION_IMPULSE_PORTS];
    let mut resting = [0; MAX_FRICTION_IMPULSE_PORTS];
    let mut count = 0;
    for i in 0..n {
        if velocity[i] == Dd::ZERO {
            if input.capacity_n[i] != 0.0 {
                resting[count] = i;
                count += 1;
            }
        } else {
            force[i] = if velocity[i].to_f64() > 0.0 {
                -input.capacity_n[i]
            } else {
                input.capacity_n[i]
            };
        }
    }
    let mut free_acceleration = [Dd::ZERO; MAX_FRICTION_IMPULSE_PORTS];
    let mut needs_projection = false;
    for a in 0..count {
        poll(cancelled)?;
        free_acceleration[a] = response_row(input, resting[a], &force)?;
        needs_projection |= free_acceleration[a] != Dd::ZERO;
    }
    if needs_projection {
        let mut matrix = [0.0; MAX_FRICTION_IMPULSE_PORTS * MAX_FRICTION_IMPULSE_PORTS];
        let mut free = [0.0; MAX_FRICTION_IMPULSE_PORTS];
        let mut capacity = [0.0; MAX_FRICTION_IMPULSE_PORTS];
        for a in 0..count {
            poll(cancelled)?;
            let i = resting[a];
            let raw = free_acceleration[a];
            free[a] = checked(raw * remaining, "resting free velocity")?.to_f64();
            capacity[a] = checked(
                Dd::from_f64(input.capacity_n[i]) * remaining,
                "resting capacity",
            )?
            .to_f64();
            if capacity[a] == 0.0 || (raw != Dd::ZERO && free[a] == 0.0) {
                return Err(CoupledFrictionFlowError::Unresolved(
                    "representable resting problem",
                ));
            }
            for b in 0..count {
                matrix[a * count + b] = input.mobility_per_kg[i * n + resting[b]];
            }
        }
        let mut projection = controls.projection;
        projection.velocity_tolerance_m_s = projection.velocity_tolerance_m_s.min(
            (controls.acceleration_tolerance_m_s2 * remaining.to_f64()).max(f64::from_bits(1)),
        );
        // This auxiliary problem covers the remaining time. Its actual force
        // budget is checked again after recovering the original physical force.
        // The public impulse budget governs the integrated path, not a longer
        // trial segment which a subsequent stop may shorten.
        projection.maximum_impulse_ns = (Dd::from_f64(controls.maximum_force_n) * remaining)
            .to_f64()
            .min(f64::MAX);
        let resolved = solve_coupled_friction_impulses(
            &CoupledFrictionImpulseInput {
                velocity_before_m_s: &free[..count],
                mobility_per_kg: &matrix[..count * count],
                capacity_ns: &capacity[..count],
            },
            projection,
            &mut *cancelled,
        )?;
        sums.projection_iterations += resolved.iterations;
        sums.factorizations += resolved.factorizations;
        for a in 0..count {
            let i = resting[a];
            let impulse = resolved.impulse_ns[a];
            force[i] = if impulse == -capacity[a] {
                -input.capacity_n[i]
            } else if impulse == capacity[a] {
                input.capacity_n[i]
            } else {
                checked(Dd::from_f64(impulse) / remaining, "resting force")?.to_f64()
            };
            if force[i].abs() > input.capacity_n[i] || (impulse != 0.0 && force[i] == 0.0) {
                return Err(CoupledFrictionFlowError::Unresolved(
                    "original force capacity",
                ));
            }
        }
    }
    for i in 0..n {
        poll(cancelled)?;
        force_budget(i, force[i], controls)?;
        force_budget(
            i,
            checked(
                Dd::from_f64(force[i]) + Dd::from_f64(input.drive_n[i]),
                "resultant force",
            )?
            .to_f64(),
            controls,
        )?;
    }
    Ok(force)
}

fn response_row(
    input: &CoupledFrictionFlowInput<'_>,
    i: usize,
    friction: &Vector,
) -> Result<Dd, CoupledFrictionFlowError> {
    let n = input.velocity_before_m_s.len();
    let mut value = Dd::ZERO;
    for j in 0..n {
        let total = Dd::from_f64(input.drive_n[j]) + Dd::from_f64(friction[j]);
        value = value + Dd::from_f64(input.mobility_per_kg[i * n + j]) * total;
    }
    checked(value, "original acceleration")
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn finish(
    input: &CoupledFrictionFlowInput<'_>,
    controls: CoupledFrictionFlowControls,
    factor: &Factor,
    scales: &Vector,
    velocity: &Accurate,
    sums: &Path,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<CoupledFrictionFlowReport, CoupledFrictionFlowError> {
    let n = input.velocity_before_m_s.len();
    let mut report = CoupledFrictionFlowReport {
        count: n,
        friction_impulse_ns: [0.0; MAX_FRICTION_IMPULSE_PORTS],
        drive_impulse_ns: [0.0; MAX_FRICTION_IMPULSE_PORTS],
        velocity_after_m_s: [0.0; MAX_FRICTION_IMPULSE_PORTS],
        flow_velocity_after_m_s: [0.0; MAX_FRICTION_IMPULSE_PORTS],
        velocity_residual_m_s: [0.0; MAX_FRICTION_IMPULSE_PORTS],
        impulse_residual_ns: [0.0; MAX_FRICTION_IMPULSE_PORTS],
        impulse_capacity_residual_ns: [0.0; MAX_FRICTION_IMPULSE_PORTS],
        displacement_m: [0.0; MAX_FRICTION_IMPULSE_PORTS],
        distance_m: [0.0; MAX_FRICTION_IMPULSE_PORTS],
        minimum_speed_m_s: sums.minimum_speed,
        maximum_speed_m_s: sums.maximum_speed,
        friction_work_j: [0.0; MAX_FRICTION_IMPULSE_PORTS],
        drive_work_j: [0.0; MAX_FRICTION_IMPULSE_PORTS],
        dissipated_j: 0.0,
        total_drive_work_j: 0.0,
        kinetic_energy_change_j: 0.0,
        energy_residual_j: 0.0,
        segments: sums.segments,
        stops: sums.stops,
        projection_iterations: sums.projection_iterations,
        factorizations: sums.factorizations,
    };
    let mut total_impulse = [0.0; MAX_FRICTION_IMPULSE_PORTS];
    let mut friction_work = Dd::ZERO;
    let mut drive_work = Dd::ZERO;
    for i in 0..n {
        poll(cancelled)?;
        report.friction_impulse_ns[i] = sums.friction_impulse[i].to_f64();
        report.drive_impulse_ns[i] = checked(
            Dd::from_f64(input.drive_n[i]) * Dd::from_f64(input.duration_s),
            "drive impulse",
        )?
        .to_f64();
        total_impulse[i] = checked(
            Dd::from_f64(report.friction_impulse_ns[i]) + Dd::from_f64(report.drive_impulse_ns[i]),
            "total impulse",
        )?
        .to_f64();
        for impulse in [
            report.friction_impulse_ns[i],
            report.drive_impulse_ns[i],
            total_impulse[i],
        ] {
            if impulse.abs() > controls.projection.maximum_impulse_ns {
                return Err(CoupledFrictionFlowError::ImpulseBudget {
                    index: i,
                    impulse_ns: impulse,
                    maximum_ns: controls.projection.maximum_impulse_ns,
                });
            }
        }
        let capacity = checked(
            Dd::from_f64(input.capacity_n[i]) * Dd::from_f64(input.duration_s),
            "integrated capacity",
        )?;
        report.impulse_capacity_residual_ns[i] =
            (Dd::from_f64(report.friction_impulse_ns[i].abs()) - capacity)
                .to_f64()
                .max(0.0);
        if report.impulse_capacity_residual_ns[i] > controls.projection.impulse_tolerance_ns {
            return Err(CoupledFrictionFlowError::Unresolved(
                "integrated force capacity",
            ));
        }
        report.flow_velocity_after_m_s[i] = velocity[i].to_f64();
        report.displacement_m[i] = sums.displacement[i].to_f64();
        report.distance_m[i] = sums.distance[i].to_f64();
        report.friction_work_j[i] = sums.friction_work[i].to_f64();
        let work = checked(
            Dd::from_f64(input.drive_n[i]) * sums.displacement[i],
            "drive work",
        )?;
        report.drive_work_j[i] = work.to_f64();
        friction_work = checked(friction_work + sums.friction_work[i], "total friction work")?;
        drive_work = checked(drive_work + work, "total drive work")?;
    }
    let mut change = [0.0; MAX_FRICTION_IMPULSE_PORTS];
    for i in 0..n {
        poll(cancelled)?;
        let mut actual = Dd::from_f64(input.velocity_before_m_s[i]);
        for (j, &impulse) in total_impulse.iter().enumerate().take(n) {
            actual =
                actual + Dd::from_f64(input.mobility_per_kg[i * n + j]) * Dd::from_f64(impulse);
        }
        report.velocity_after_m_s[i] = checked(actual, "original endpoint")?.to_f64();
        report.velocity_residual_m_s[i] = checked(
            Dd::from_f64(report.velocity_after_m_s[i]) - velocity[i],
            "path endpoint residual",
        )?
        .to_f64();
        change[i] = checked(
            Dd::from_f64(report.velocity_after_m_s[i]) - Dd::from_f64(input.velocity_before_m_s[i]),
            "velocity change",
        )?
        .to_f64();
        if report.velocity_residual_m_s[i].abs() > controls.projection.velocity_tolerance_m_s {
            return Err(CoupledFrictionFlowError::Unresolved(
                "original path endpoint",
            ));
        }
    }
    let reconstructed = factor.solve(&change, scales, cancelled)?;
    let mut kinetic = Dd::ZERO;
    for i in 0..n {
        poll(cancelled)?;
        report.impulse_residual_ns[i] = checked(
            Dd::from_f64(reconstructed[i]) - Dd::from_f64(total_impulse[i]),
            "impulse reconstruction",
        )?
        .to_f64();
        if report.impulse_residual_ns[i].abs() > controls.projection.impulse_tolerance_ns {
            return Err(CoupledFrictionFlowError::Unresolved(
                "original impulse reconstruction",
            ));
        }
        let midpoint = Dd::from_f64(input.velocity_before_m_s[i]) * Dd::from_f64(0.5)
            + Dd::from_f64(report.velocity_after_m_s[i]) * Dd::from_f64(0.5);
        kinetic = checked(
            kinetic + Dd::from_f64(reconstructed[i]) * midpoint,
            "kinetic change",
        )?;
    }
    report.dissipated_j = (-friction_work).to_f64();
    report.total_drive_work_j = drive_work.to_f64();
    report.kinetic_energy_change_j = kinetic.to_f64();
    report.energy_residual_j = checked(
        kinetic - drive_work - friction_work,
        "continuous energy closure",
    )?
    .to_f64();
    if report.dissipated_j < 0.0
        || report.energy_residual_j.abs() > controls.projection.energy_tolerance_j
    {
        return Err(CoupledFrictionFlowError::EnergyBalance {
            dissipated_j: report.dissipated_j,
            residual_j: report.energy_residual_j,
        });
    }
    poll(cancelled)?;
    Ok(report)
}

fn force_budget(
    index: usize,
    force: f64,
    controls: CoupledFrictionFlowControls,
) -> Result<(), CoupledFrictionFlowError> {
    if !force.is_finite() || force.abs() > controls.maximum_force_n {
        Err(CoupledFrictionFlowError::ForceBudget {
            index,
            force_n: force,
            maximum_n: controls.maximum_force_n,
        })
    } else {
        Ok(())
    }
}
fn checked(value: Dd, label: &'static str) -> Result<Dd, CoupledFrictionFlowError> {
    if value.hi.is_finite() && value.lo.is_finite() {
        Ok(value)
    } else {
        Err(CoupledFrictionImpulseError::NonFinite(label).into())
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CoupledFrictionFlowError> {
    if cancelled() {
        Err(CoupledFrictionFlowError::Cancelled)
    } else {
        Ok(())
    }
}
