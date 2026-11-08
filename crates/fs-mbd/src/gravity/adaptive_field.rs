//! Adaptive intervals in a prescribed field, using the original rigid step.

use super::adaptive::{NativeStep, advance_native};
use super::{
    GravityAdvanceError, GravityAdvanceOptions, GravityAdvanceReceipt, RigidGravityField, poll,
    step_in_field_with_loads,
};
use crate::{MassProperties, RigidBodyState, Wrench};

/// Advance a complete interval under a stationary conservative gravity field.
///
/// This is the bounded step-doubling controller used by mutual gravity, applied
/// to the original [`step_in_field_with_loads`] map. It works with uniform,
/// central, quadrupole and superposed fields, including custom trait objects.
/// The accepted result is two native half-steps, without extrapolating a pose or
/// projecting momenta. The single final state is `receipt.states_after[0]`;
/// fixed-size state arrays avoid heap allocation in this adapter/controller.
///
/// `other_wrench` has the same held-world-force/held-body-torque semantics on
/// every substep as on the original step. Its work is not included in the
/// conservative potential change. Only accepted substeps contribute to that
/// change; rejected trials do not alter the returned energy account.
///
/// Fields and callbacks must be repeatable and side-effect-free. Cancellation,
/// an exhausted trial budget, or an unresolved minimum step returns no state.
/// The estimator controls local state error, not a rigorous trajectory bound.
/// Changing steps according to state does not preserve fixed-step symplecticity
/// or time reversal. Endpoint field admission still does not certify source or
/// model-boundary clearance between samples; required event admission remains
/// caller-owned. Malformed field/native arithmetic errors are never retried.
pub fn advance_in_field(
    field: &(impl RigidGravityField + ?Sized),
    frame: u64,
    state: RigidBodyState,
    properties: MassProperties,
    other_wrench: Wrench,
    duration_seconds: f64,
    options: GravityAdvanceOptions,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GravityAdvanceReceipt<[RigidBodyState; 1]>, GravityAdvanceError> {
    options.validate(duration_seconds)?;
    state.validate().map_err(super::GravityError::from)?;
    properties.validate().map_err(super::GravityError::from)?;
    other_wrench.validate().map_err(super::GravityError::from)?;
    poll(&mut cancelled)?;
    field.load(frame, state.pose(), properties)?.validate()?;
    poll(&mut cancelled)?;
    advance_native(
        [state],
        &[properties],
        duration_seconds,
        options,
        |states, duration, cancelled| {
            let result = step_in_field_with_loads(
                field,
                frame,
                states[0],
                properties,
                other_wrench,
                duration,
                cancelled,
            )?;
            Ok(NativeStep {
                states_after: [result.step.state_after],
                diagnostics_before: result.step.diagnostics_before,
                diagnostics_after: result.step.diagnostics_after,
                potential_energy_change_j: result.potential_energy_change_j,
            })
        },
        cancelled,
    )
}
