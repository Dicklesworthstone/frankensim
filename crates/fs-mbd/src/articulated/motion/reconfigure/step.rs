//! Continuous prescribed tensor changes through the original free-floating owner.
//!
//! Two momentum-preserving parameter transfers surround one native smooth step
//! at the supplied midpoint inertia. The parameter work remains separate from
//! force/actuator work and numerical integration error.

use super::super::{
    ArticulatedModel, FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
    FreeFloatingMotionStep, step_free_floating, step_free_floating_with_effort,
};
use super::{
    InertiaReconfiguration, ReconfigurationControls, ReconfigurationError, transfer_inertia,
    validate_compatible_models, validate_controls,
};
use core::fmt;
use fs_ga::Vec3;

/// One accepted changing-inertia step, retaining the original physical receipts.
#[derive(Debug)]
pub struct InertiaMotionStep {
    /// Fixed-configuration transfer from the initial to midpoint inertia.
    pub first_transfer: InertiaReconfiguration,
    /// Original smooth evolution at the midpoint model. Its endpoint precedes
    /// the final parameter transfer; its forces and work keep native semantics.
    pub motion: FreeFloatingMotionStep,
    /// Fixed-configuration transfer at the actual smooth endpoint. This owns
    /// the accepted state under the final model.
    pub final_transfer: InertiaReconfiguration,
    /// Signed sum of the two parameter-work contributions, J. Negative work
    /// removes resolved mechanical kinetic energy; no heat allocation is inferred.
    pub net_parameter_work_j: f64,
    /// Sum of ABSOLUTE parameter work, J. Opposite-signed transfers cannot hide
    /// work from the whole-window `maximum_parameter_work_j` budget.
    pub absolute_parameter_work_j: f64,
    /// Actual final-minus-initial kinetic energy, less parameter work and all
    /// original smooth work. This is numerical error, not physical dissipation.
    /// Lab-frame endpoint subtraction can lose small signals under a large
    /// common boost; the separate parameter-work receipts remain conditioned.
    pub kinetic_work_residual: f64,
    /// Actual world linear-momentum change less original external impulse, N s.
    pub linear_momentum_residual: Vec3,
    /// Actual world-origin angular-momentum change less external impulse, N m s.
    pub angular_momentum_residual: Vec3,
}

impl InertiaMotionStep {
    /// Borrow the actual accepted state under the final inertia model.
    #[must_use]
    pub fn state_after(&self) -> &FreeFloatingMotionState {
        &self.final_transfer.state_after
    }

    /// Consume the receipt and take the accepted state without cloning it.
    #[must_use]
    pub fn into_state_after(self) -> FreeFloatingMotionState {
        self.final_transfer.state_after
    }
}

/// Atomic refusal of either parameter transfer or the original smooth step.
#[derive(Debug)]
pub enum InertiaMotionError {
    /// Parameter, model compatibility, state, momentum or work admission failed.
    Reconfiguration(ReconfigurationError),
    /// The original continuous dynamics or effort callback refused.
    Motion(FreeFloatingMotionError),
    /// The actual gross parameter work exceeds the complete-window cap.
    ParameterWorkBudget {
        /// Requested sum of absolute parameter work, J.
        requested: f64,
        /// Caller-declared complete-window cap, J.
        maximum: f64,
    },
    /// A complete-window work or momentum quantity cannot be represented.
    InvalidAccounting(&'static str),
}

impl fmt::Display for InertiaMotionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Reconfiguration(source) => write!(f, "inertia transfer refused: {source}"),
            Self::Motion(source) => write!(f, "changing-inertia motion refused: {source}"),
            Self::ParameterWorkBudget { requested, maximum } => write!(
                f,
                "gross inertia parameter work {requested} J exceeds window budget {maximum} J",
            ),
            Self::InvalidAccounting(what) => {
                write!(f, "changing-inertia accounting refused: {what}")
            }
        }
    }
}
impl std::error::Error for InertiaMotionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Reconfiguration(source) => Some(source),
            Self::Motion(source) => Some(source),
            _ => None,
        }
    }
}
impl From<ReconfigurationError> for InertiaMotionError {
    fn from(value: ReconfigurationError) -> Self {
        Self::Reconfiguration(value)
    }
}
impl From<FreeFloatingMotionError> for InertiaMotionError {
    fn from(value: FreeFloatingMotionError) -> Self {
        Self::Motion(value)
    }
}

/// Advance with prescribed initial, midpoint and final full link inertia tensors.
///
/// Models must retain the same per-link mass, COM, frames, joint definitions,
/// limits and topology. They supply the inertia samples at elapsed times zero,
/// `duration/2` and `duration`. This reduced model assumes no additional intrinsic
/// angular/shape momentum, mass flux or moving-COM degree of freedom. Parameter
/// work measures exchange with the resolved mechanical motion; hidden actuator
/// storage and dissipation require their own physical model.
///
/// At fixed configuration, the first transfer preserves every canonical momentum
/// while changing to the midpoint tensor. The original native smooth owner then
/// advances the entire duration at that tensor. The final transfer preserves all
/// momenta at the actual accepted endpoint while changing to the final tensor.
/// No transfer advances configuration or applies a continuous-load impulse.
///
/// Second-order refinement requires a smooth prescribed inertia history, a
/// midpoint sample accurate to O(duration²), smooth resolved loads and timestep
/// resolution of all mechanical modes. The transfer placement is symmetric;
/// the original explicit midpoint integrator is not time-reversible, variational
/// or energy-exact. Discontinuous tensor changes should use `transfer_inertia`
/// as explicit events instead of claiming continuous second-order accuracy.
///
/// Held link-frame wrenches rotate with the links; uniform gravity is world-fixed.
/// The work cap applies to the SUM OF ABSOLUTE work over both transfers. The
/// second transfer receives only the remaining budget. Every input is borrowed
/// unchanged; all admissions and final cancellation must succeed for publication.
///
/// # Errors
/// Incompatible models, invalid controls/duration, either transfer's physical
/// refusal, original smooth-step refusal, aggregate work exhaustion, nonfinite
/// accounting, or cancellation. No intermediate result is published.
#[allow(clippy::too_many_arguments)]
pub fn advance(
    previous: &ArticulatedModel,
    midpoint: &ArticulatedModel,
    next: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: FreeFloatingMotionLoads<'_>,
    duration: f64,
    controls: ReconfigurationControls<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<InertiaMotionStep, InertiaMotionError> {
    advance_using(
        previous,
        midpoint,
        next,
        state,
        duration,
        controls,
        |transferred, cancel| step_free_floating(midpoint, transferred, loads, duration, cancel),
        cancelled,
    )
}

/// Changing-inertia evolution with the original native stage-dependent effort law.
///
/// The callback runs exactly at the native elapsed stages `0` and `duration/2`.
/// Its stage-zero state is the VIRTUAL, momentum-transferred state under the
/// midpoint tensor; its velocities generally differ from the supplied old-model
/// velocities. The second call receives the native reconstructed midpoint state.
/// Both callbacks use the same midpoint model and a scratch effort slice reset
/// to `loads.generalized_force`. Capture an absolute start time if needed.
///
/// There are no effort callbacks during parameter transfers. Smooth state/time
/// dependent laws retain the refinement scope of [`advance`]: the O(duration)
/// change of the virtual start velocity changes the explicit midpoint predictor
/// only at O(duration²). Callback-owned storage/loss accounting and side effects
/// remain the caller's responsibility; side effects cannot be rolled back.
///
/// # Errors
/// The refusals of [`advance`] and any error returned by the effort callback.
#[allow(clippy::too_many_arguments)]
pub fn advance_with_effort(
    previous: &ArticulatedModel,
    midpoint: &ArticulatedModel,
    next: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: FreeFloatingMotionLoads<'_>,
    duration: f64,
    controls: ReconfigurationControls<'_>,
    effort: impl FnMut(f64, &FreeFloatingMotionState, &mut [f64]) -> Result<(), FreeFloatingMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<InertiaMotionStep, InertiaMotionError> {
    advance_using(
        previous,
        midpoint,
        next,
        state,
        duration,
        controls,
        |transferred, cancel| {
            step_free_floating_with_effort(midpoint, transferred, loads, duration, effort, cancel)
        },
        cancelled,
    )
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), InertiaMotionError> {
    if cancelled() {
        Err(FreeFloatingMotionError::Cancelled.into())
    } else {
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn advance_using<C: FnMut() -> bool>(
    previous: &ArticulatedModel,
    midpoint: &ArticulatedModel,
    next: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    duration: f64,
    controls: ReconfigurationControls<'_>,
    evolve: impl FnOnce(
        &FreeFloatingMotionState,
        &mut C,
    ) -> Result<FreeFloatingMotionStep, FreeFloatingMotionError>,
    mut cancelled: C,
) -> Result<InertiaMotionStep, InertiaMotionError> {
    poll(&mut cancelled)?;
    if !duration.is_finite() || duration <= 0.0 {
        return Err(FreeFloatingMotionError::InvalidDuration.into());
    }
    validate_controls(previous, controls)?;
    validate_compatible_models(previous, midpoint)?;
    poll(&mut cancelled)?;
    validate_compatible_models(midpoint, next)?;
    poll(&mut cancelled)?;
    let first_transfer = transfer_inertia(previous, midpoint, state, controls, &mut cancelled)?;
    let first_work = first_transfer.parameter_work_j;
    if !first_work.is_finite() {
        return Err(InertiaMotionError::InvalidAccounting(
            "first parameter work",
        ));
    }
    let first_absolute_work = first_work.abs();
    if first_absolute_work > controls.maximum_parameter_work_j {
        return Err(InertiaMotionError::ParameterWorkBudget {
            requested: first_absolute_work,
            maximum: controls.maximum_parameter_work_j,
        });
    }
    let remaining = controls.maximum_parameter_work_j - first_absolute_work;
    let motion = evolve(&first_transfer.state_after, &mut cancelled)?;
    poll(&mut cancelled)?;
    let final_transfer = transfer_inertia(
        midpoint,
        next,
        &motion.state_after,
        ReconfigurationControls {
            maximum_parameter_work_j: remaining,
            ..controls
        },
        &mut cancelled,
    )?;
    let net_parameter_work_j = first_work + final_transfer.parameter_work_j;
    let absolute_parameter_work_j = first_absolute_work + final_transfer.parameter_work_j.abs();
    if !net_parameter_work_j.is_finite() || !absolute_parameter_work_j.is_finite() {
        return Err(InertiaMotionError::InvalidAccounting(
            "total parameter work",
        ));
    }
    // Even an ulp of representational overshoot refuses the whole window.
    if absolute_parameter_work_j > controls.maximum_parameter_work_j {
        return Err(InertiaMotionError::ParameterWorkBudget {
            requested: absolute_parameter_work_j,
            maximum: controls.maximum_parameter_work_j,
        });
    }
    let before = &first_transfer.momentum_before;
    let after = &final_transfer.momentum_after;
    let kinetic_work_residual = after.kinetic_energy
        - before.kinetic_energy
        - net_parameter_work_j
        - motion.actuator_work
        - motion.external_work
        - motion.internal_body_work;
    let linear_momentum_residual = after.linear_momentum_world
        - before.linear_momentum_world
        - motion.external_linear_impulse_world;
    let angular_momentum_residual = after.angular_momentum_world
        - before.angular_momentum_world
        - motion.external_angular_impulse_world;
    if ![
        kinetic_work_residual,
        linear_momentum_residual.x,
        linear_momentum_residual.y,
        linear_momentum_residual.z,
        angular_momentum_residual.x,
        angular_momentum_residual.y,
        angular_momentum_residual.z,
    ]
    .iter()
    .all(|x| x.is_finite())
    {
        return Err(InertiaMotionError::InvalidAccounting(
            "endpoint work/momentum residual",
        ));
    }
    poll(&mut cancelled)?;
    Ok(InertiaMotionStep {
        first_transfer,
        motion,
        final_transfer,
        net_parameter_work_j,
        absolute_parameter_work_j,
        kinetic_work_residual,
        linear_momentum_residual,
        angular_momentum_residual,
    })
}
