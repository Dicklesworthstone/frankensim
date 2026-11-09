//! Finite-inductance actuation through the original sustained-contact midpoint.
//!
//! Contact mobility includes the same electrical joint impedance as the free
//! response. The accepted constrained rate therefore sets back EMF, current,
//! motor effort and work together. Ground reaction is not added after a free
//! motor step, and no load callback is repeated to reconstruct the winding.

use super::{
    ArticulatedLoads, ArticulatedModel, ArticulatedState, GeneralizedKinematics,
    JointImpedance, JointVoltageDrive, JointVoltageReceipt, PreparedDrives,
    Refusal, VoltageDrivenState, cancelled, invalid,
};
use crate::partitioned::articulated::contact::{
    self as native, ArticulatedFrictionalContact, ContactPath, ContactStep,
    ContactStepConfig,
};

/// One accepted native contact step and its continuing winding/source accounts.
#[derive(Clone, Debug)]
pub struct VoltageDrivenContactStep<P> {
    /// Original full dry/fluid motion, actual contact impulses, friction/support
    /// work and caller's complete-path admission. Auxiliary contact projection
    /// energy is not physical heat and is not included in electrical accounts.
    pub mechanical: ContactStep<P>,
    /// The original winding receipts at that SAME accepted constrained midpoint.
    pub drives: Vec<JointVoltageReceipt>,
}

impl<P> VoltageDrivenContactStep<P> {
    /// Carry accepted body/joint/aero-memory state and winding currents forward.
    #[must_use]
    pub fn into_state(self) -> VoltageDrivenState {
        VoltageDrivenState {
            mechanical: self.mechanical.step.state,
            currents_a: self.drives.into_iter().map(|d| d.current_after_a).collect(),
        }
    }
}

/// Advance original sustained Coulomb contact and inductive motors together.
///
/// Electrical feedback enters the original effective mass BEFORE contact
/// mobility and reactions are solved. Aerodynamic hinge loads, dry/fluid bias,
/// passive joint impedance and electrical torque share the original native
/// midpoint. ContactPath supplies its actual accepted rates and endpoint;
/// this owner neither differences positions nor evaluates aerodynamics again.
///
/// Original contact mode, incoming-impact, friction-turning, path-coverage and
/// cancellation requirements remain in force. Electrical capability/balance
/// admission precedes the caller's required geometry admission; failure returns
/// neither current, motion, geometry result nor source debit. Native cancellation
/// polls bracket this bounded, at-most-eight-drive admission and the path owner.
///
/// Contact friction/support work remains in the mechanical receipt, separate
/// from winding/bridge heat and source work. Geometry and general nonlinear
/// energy accuracy retain the original contact owner's refinement limits.
/// This fixed-feature sustained lane does not locate joint stops or impacts,
/// compose contact across multiple trees, or provide a shared electrical bus.
#[allow(clippy::too_many_arguments)]
pub fn step_with_contacts<F, P>(
    model: &ArticulatedModel,
    impedance: &[JointImpedance],
    state: &VoltageDrivenState,
    memory_inputs: &[f64],
    drives: &[JointVoltageDrive],
    t_s: f64,
    dt_s: f64,
    contacts: &[ArticulatedFrictionalContact],
    config: ContactStepConfig<'_>,
    evaluate: F,
    admit_path: impl FnOnce(&ContactPath<'_>) -> Result<P, Refusal>,
    mut is_cancelled: impl FnMut() -> bool,
) -> Result<VoltageDrivenContactStep<P>, Refusal>
where
    F: Fn(f64, &ArticulatedState, &GeneralizedKinematics) -> Result<ArticulatedLoads, Refusal>,
{
    let prepared = PreparedDrives::new(model, impedance, state, drives, dt_s, &mut is_cancelled)?;
    let mut receipts = None;
    let mechanical = native::step_with_contacts(
        model,
        &prepared.impedance,
        &state.mechanical,
        memory_inputs,
        t_s,
        dt_s,
        contacts,
        config,
        |time, stage, map| prepared.apply(evaluate(time, stage, map)?),
        |path| {
            let accepted = prepared.drives.iter().copied().map(|drive| {
                drive.accepted_receipt(
                    path.endpoint,
                    &path.midpoint.motion.joint_velocity,
                    dt_s,
                )
            }).collect::<Result<Vec<_>, _>>()?;
            let admission = admit_path(path)?;
            receipts = Some(accepted);
            Ok(admission)
        },
        &mut is_cancelled,
    )?;
    let drives = receipts.ok_or_else(|| invalid("missing accepted constrained winding state"))?;
    if is_cancelled() {
        return Err(cancelled());
    }
    Ok(VoltageDrivenContactStep { mechanical, drives })
}
