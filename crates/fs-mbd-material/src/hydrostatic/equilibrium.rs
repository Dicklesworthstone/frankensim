//! Fixed-attitude heave initialization on the original hydrostatic load owner.
use super::{
    DisplacementMesh, HydrostaticError, HydrostaticLoad, HydrostaticReservoir, finite,
    gravity_direction, poll,
};
use crate::principal::PrincipalFrame;
use fs_mbd::{Pose, RigidBodyState, Vec3};

/// Hard bound on reservoir load queries in one heave solve. Each query retains
/// the reservoir's own geometry, waterline and cancellation bounds.
pub const MAX_HEAVE_EVALUATIONS: usize = 128;

/// A caller-declared heave bracket, absolute force allowance and work cap.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeaveEquilibriumControls {
    /// Lower translation from the reference pose, opposite gravity [m].
    /// Negative values move the body down; lateral position and attitude stay fixed.
    pub lower_offset_m: f64,
    /// Upper translation from the same reference pose [m], greater than the lower.
    pub upper_offset_m: f64,
    /// Positive absolute allowance for (buoyancy + weight) along heave [N].
    /// This numerical residual is not a height or source-uncertainty bound.
    pub maximum_force_defect_n: f64,
    /// Includes bracket endpoints; must lie in 2..=128.
    pub maximum_evaluations: usize,
}

/// A held-attitude, zero-momentum initialization and its actual reservoir load.
/// This is not a roll/pitch equilibrium or a floating-stability certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeaveEquilibrium {
    /// Candidate on the original principal frame, with both momenta zero.
    pub state: RigidBodyState,
    /// Actual load at the returned pose, including any solved tank interfaces.
    pub load: HydrostaticLoad,
    /// Translation from the reference pose, positive opposite gravity [m].
    pub offset_m: f64,
    /// Actual buoyancy plus weight in the world frame [N].
    pub net_force_world_n: Vec3,
    /// Signed net force projected opposite gravity [N]; positive means upward.
    pub heave_force_residual_n: f64,
    /// External torque needed to hold this attitude, in principal axes [N m].
    /// Gravity acts at the centre of mass, so this is minus the buoyancy torque.
    /// It is reported, not applied by the free-body integrator.
    pub holding_torque_body_nm: Vec3,
    /// Actual load queries, including any evaluated bracket endpoints.
    pub evaluations: usize,
}

/// A failed initialization publishes no candidate and never changes the input.
#[derive(Clone, Debug, PartialEq)]
pub enum HeaveEquilibriumError {
    /// Invalid controls, arithmetic, geometry, reservoir or cancellation refusal.
    Physics(HydrostaticError),
    /// Neither endpoint meets the force allowance and their residuals have the
    /// same sign. This refuses this bracket, not every possible equilibrium.
    NotBracketed {
        /// Lower requested offset [m].
        lower_offset_m: f64,
        /// Upper requested offset [m].
        upper_offset_m: f64,
        /// Signed upward force at the lower endpoint [N].
        lower_residual_n: f64,
        /// Signed upward force at the upper endpoint [N].
        upper_residual_n: f64,
    },
    /// The evaluation budget or representable bracket was exhausted without
    /// satisfying the absolute force allowance. No approximate state is returned.
    NotConverged {
        /// Smallest observed signed force residual [N].
        best_residual_n: f64,
        /// Requested absolute force allowance [N].
        maximum_n: f64,
        /// Actual number of load evaluations.
        evaluations: usize,
    },
}
impl core::fmt::Display for HeaveEquilibriumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "hydrostatic heave equilibrium: {self:?}")
    }
}
impl std::error::Error for HeaveEquilibriumError {}
impl From<HydrostaticError> for HeaveEquilibriumError {
    fn from(error: HydrostaticError) -> Self {
        Self::Physics(error)
    }
}

/// Initialize heave by balancing source-derived buoyancy and body weight at
/// fixed attitude and lateral position. Works through the same `load_frame`
/// boundary for prescribed uniform/layered/continuous reservoirs and finite
/// tanks; their frozen source snapshots are not resolved again by this solve.
/// Finite tanks solve their conserved inventories at every candidate pose.
///
/// The bracket is a translation along `-gravity / |gravity|` relative to
/// `reference.pose()`. ITP steps (a secant point truncated toward the midpoint
/// and projected so the bracket shrinks at least as fast as bisection, with
/// one step of slack) retain an opposite-sign bracket; an endpoint already
/// meeting the allowance may be returned directly. No uniqueness, roll/pitch
/// balance or stability is inferred.
/// Reported holding torque must be supplied by the caller's actual constraint
/// if the body is to remain at the chosen attitude after initialization.
///
/// This initializes a new state with zero momentum. It is not a physical
/// damping operation or an energy/momentum-conserving live-state replacement.
/// The caller commits only a successful result. Unsupported poses, overflow,
/// failed inner waterline solves and cancellation propagate without changing
/// the body, reference state or reservoir. The scalar solve allocates no heap
/// memory beyond work already performed by the original reservoir load query.
pub fn solve_heave_equilibrium<R: HydrostaticReservoir + ?Sized, C: FnMut() -> bool>(
    reservoir: &R,
    mesh: &DisplacementMesh,
    body: &PrincipalFrame,
    reference: RigidBodyState,
    controls: HeaveEquilibriumControls,
    mut cancelled: C,
) -> Result<HeaveEquilibrium, HeaveEquilibriumError> {
    poll(&mut cancelled)?;
    let mut low = controls.lower_offset_m;
    let mut high = controls.upper_offset_m;
    if !low.is_finite()
        || !high.is_finite()
        || !((high - low).is_finite() && high > low)
        || !controls.maximum_force_defect_n.is_finite()
        || controls.maximum_force_defect_n <= 0.0
        || !(2..=MAX_HEAVE_EVALUATIONS).contains(&controls.maximum_evaluations)
    {
        return Err(HydrostaticError::Invalid(
            "heave requires an ordered finite bracket, positive force tolerance and bounded evaluations",
        )
        .into());
    }
    let gravity = reservoir.gravity();
    let (down, _) = gravity_direction(gravity)?;
    let up = Vec3::new(-down.x, -down.y, -down.z);
    let weight = gravity.acceleration_world().scale(body.properties().mass());
    if !weight.is_finite() {
        return Err(HydrostaticError::Invalid("body weight is unrepresentable").into());
    }
    let reference_pose = reference.pose();
    let evaluate = |offset: f64, evaluations: usize, cancel: &mut C| {
        poll(cancel)?;
        let pose = Pose::new(
            reference_pose.position_world().add(up.scale(offset)),
            reference_pose.orientation(),
        )?;
        let state = RigidBodyState::new(pose, Vec3::ZERO, Vec3::ZERO)?;
        let load = reservoir.load_frame(mesh, body, state, &mut *cancel)?;
        let net_force_world_n = load.wrench().force_world.add(weight);
        let holding_torque_body_nm = load.wrench().torque_body.scale(-1.0);
        if !net_force_world_n.is_finite() || !holding_torque_body_nm.is_finite() {
            return Err(HydrostaticError::Invalid("heave load is unrepresentable"));
        }
        let heave_force_residual_n = finite(net_force_world_n.dot(up))?;
        poll(cancel)?;
        Ok::<_, HydrostaticError>(HeaveEquilibrium {
            state,
            load,
            offset_m: offset,
            net_force_world_n,
            heave_force_residual_n,
            holding_torque_body_nm,
            evaluations,
        })
    };
    let lower = evaluate(low, 1, &mut cancelled)?;
    let mut f_low = lower.heave_force_residual_n;
    if f_low.abs() <= controls.maximum_force_defect_n {
        return Ok(lower);
    }
    let upper = evaluate(high, 2, &mut cancelled)?;
    let mut f_high = upper.heave_force_residual_n;
    if f_high.abs() <= controls.maximum_force_defect_n {
        return Ok(upper);
    }
    if f_low.is_sign_positive() == f_high.is_sign_positive() {
        return Err(HeaveEquilibriumError::NotBracketed {
            lower_offset_m: low,
            upper_offset_m: high,
            lower_residual_n: f_low,
            upper_residual_n: f_high,
        });
    }
    let mut best_residual_n = if f_low.abs() <= f_high.abs() {
        f_low
    } else {
        f_high
    };
    let mut evaluations = 2;
    // ITP steps (Oliveira and Takahashi, ACM TOMS 47(1), 2020). After the k-th
    // interior evaluation the bracket is at most `initial_span * 2^(1 - k)`,
    // bisection's width one evaluation earlier. At a slope kink on the root
    // (a face crossing a waterline or density interface) plain secant steps
    // keep replacing one end; the projection bounds that. The first step is
    // the untruncated secant point, exact for an affine residual. Truncation
    // uses the paper's kappa1 = 0.2 / initial_span and kappa2 = 2.
    let initial_span = high - low;
    let mut bracket_bound = initial_span;
    for evaluation in 3..=controls.maximum_evaluations {
        poll(&mut cancelled)?;
        let span = high - low;
        // Scale before adding opposite-sign force magnitudes; direct secant
        // subtraction could overflow despite finite endpoint loads.
        let scale = f_low.abs().max(f_high.abs());
        let left = f_low.abs() / scale;
        let right = f_high.abs() / scale;
        let secant = low + left / (left + right) * span;
        let midpoint = low + 0.5 * span;
        let toward_midpoint = if midpoint >= secant { 1.0 } else { -1.0 };
        let truncation = if evaluation == 3 {
            0.0
        } else {
            0.2 * span * (span / initial_span)
        };
        let truncated = if truncation <= (midpoint - secant).abs() {
            secant + toward_midpoint * truncation
        } else {
            midpoint
        };
        let reach = (bracket_bound - 0.5 * span).max(0.0);
        bracket_bound *= 0.5;
        let mut offset = if (truncated - midpoint).abs() <= reach {
            truncated
        } else {
            midpoint - toward_midpoint * reach
        };
        if offset <= low || offset >= high {
            offset = midpoint;
        }
        if offset <= low || offset >= high {
            break;
        }
        let candidate = evaluate(offset, evaluation, &mut cancelled)?;
        evaluations = evaluation;
        let residual = candidate.heave_force_residual_n;
        if residual.abs() <= controls.maximum_force_defect_n {
            return Ok(candidate);
        }
        if residual.abs() < best_residual_n.abs() {
            best_residual_n = residual;
        }
        if residual.is_sign_positive() == f_low.is_sign_positive() {
            low = offset;
            f_low = residual;
        } else {
            high = offset;
            f_high = residual;
        }
    }
    Err(HeaveEquilibriumError::NotConverged {
        best_residual_n,
        maximum_n: controls.maximum_force_defect_n,
        evaluations,
    })
}
