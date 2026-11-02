//! Native adhesive recoil after a prescribed thermal change in the normal law.
//!
//! Bulk temperatures change the actual elastic response; interface temperature
//! changes work of adhesion and junction strengths. The new radius determines
//! both the constitutive torque and all scaled native mobility entries. Normal
//! potential change is signed parameter work, separate from actual frictional
//! heat. The normal approach and body poses are held during this split: this is
//! not a normal-motion solve, thermal expansion or autonomous heat transport.

use fs_contact::interface_binding::adhesive_shear::{AdhesiveShearError, SourceAdhesiveShear};
use fs_contact::normal_patch::adhesive::{
    JkrResponse, JkrState, SourceJkrSphereContact,
    temperature::{AdhesiveThermalControl, SourceJkrThermalRefresh},
};
use fs_exec::CancelGate;
use fs_tribo::partial_slip::TangentFrame;

use super::{
    RigidAdhesiveContact, RigidAdhesiveControls, RigidAdhesiveError, RigidAdhesiveReport,
    prepare_rigid_adhesive,
};
use crate::tangential_friction::{RigidContactBinding, poll};
use crate::{MassProperties, RigidBodyState};

/// Actual contact storage and native work, with no inferred thermal partition.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveThermalNormalReport {
    /// Native impulses, velocities, momenta and physical junction loss.
    pub mechanical: RigidAdhesiveReport,
    /// Original normal response reconstructed from its accepted state token.
    pub normal_before: JkrResponse,
    /// Actual normal equilibrium after the prescribed temperature change.
    pub normal_after: JkrResponse,
    /// Signed change in the represented elastic-plus-surface potential, J.
    pub normal_potential_change_j: f64,
    /// Signed supply associated with changing the normal-law parameters, J.
    /// This is neither friction heat nor a thermodynamic entropy calculation.
    pub parameter_work_j: f64,
    /// Independently differenced native kinetic plus normal potential energy, J.
    pub represented_storage_change_j: f64,
    /// Storage change + actual junction loss - parameter work, J.
    pub energy_residual_j: f64,
}

/// Complete body and source proposal. Retrying never changes accepted owners.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidAdhesiveThermalNormalCandidate {
    /// Actual endpoints in the original binding's A/B order.
    pub bodies_after: [RigidBodyState; 2],
    /// New normal owner, retaining original cards, claims and geometry.
    pub normal_source: SourceJkrSphereContact,
    /// Successor branch token at the unchanged physical approach.
    pub normal_state: JkrState,
    /// Same original junction claims at the new interface temperature.
    pub shear_source: SourceAdhesiveShear,
    /// Independent storage and native mechanics accounts.
    pub report: RigidAdhesiveThermalNormalReport,
}

/// Source, geometry, native impulse or complete storage refusal.
#[derive(Debug)]
pub enum RigidAdhesiveThermalNormalError {
    Invalid(&'static str),
    Shear(AdhesiveShearError),
    Mechanical(RigidAdhesiveError),
    Energy { residual_j: f64 },
}
impl From<AdhesiveShearError> for RigidAdhesiveThermalNormalError {
    fn from(error: AdhesiveShearError) -> Self {
        Self::Shear(error)
    }
}
impl From<RigidAdhesiveError> for RigidAdhesiveThermalNormalError {
    fn from(error: RigidAdhesiveError) -> Self {
        Self::Mechanical(error)
    }
}
impl core::fmt::Display for RigidAdhesiveThermalNormalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid adhesive thermal normal update: {self:?}")
    }
}
impl std::error::Error for RigidAdhesiveThermalNormalError {}

/// Recompute the real disk and native recoil from one held-approach update.
///
/// The accepted branch token, old normal model and old source-owned disk must
/// all agree. A held-load update with nonzero normal motion refuses, since it
/// needs its own normal kinematics and loading-device work. Lost equilibrium
/// and branch events are already refused by the normal thermal owner.
///
/// The existing native solver independently checks actual impulses, complete
/// body kinetic work, source-law friction heat and momentum closure. The final
/// combined check also includes signed normal potential/parameter work; it
/// never turns thermal parameter work or numerical residuals into friction heat.
#[allow(clippy::too_many_arguments)]
pub fn prepare_rigid_adhesive_with_normal_refresh(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidAdhesiveContact<'_>,
    normal: &SourceJkrSphereContact,
    state: JkrState,
    refreshed: SourceJkrThermalRefresh,
    dt_s: f64,
    controls: RigidAdhesiveControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveThermalNormalCandidate, RigidAdhesiveThermalNormalError> {
    poll(gate).map_err(RigidAdhesiveError::from)?;
    let before = *refreshed.before();
    let after = *refreshed.response();
    if refreshed.original_model_identity() != normal.identity()
        || refreshed.control() != AdhesiveThermalControl::HeldIndentation
        || state != before.state()
        || contact.source.normal_state().identity() != normal.normal_state().identity()
        || contact.source.contact_area_m2().to_bits() != before.contact_area_m2.to_bits()
        || contact
            .source
            .patch()
            .parameters()
            .contact_radius_m
            .to_bits()
            != before.contact_radius_m.to_bits()
        || refreshed.mechanical_work_j() != 0.0
        || before.indentation_m.to_bits() != after.indentation_m.to_bits()
    {
        return Err(RigidAdhesiveThermalNormalError::Invalid(
            "same accepted normal branch and actual disk at held physical approach required",
        ));
    }
    let shear = contact
        .source
        .shear_source()
        .rebind_jkr_temperature(normal, &refreshed)?;
    let source = shear.sliding_patch_jkr(refreshed.source().normal_state(), &after)?;
    let mechanical = prepare_rigid_adhesive(
        bodies,
        properties,
        RigidAdhesiveContact { source, ..contact },
        dt_s,
        controls,
        gate,
    )?;
    let potential_change = after.reversible_energy_j - before.reversible_energy_j;
    let parameter_work = refreshed.parameter_work_j();
    let total_before = mechanical.report.kinetic_before_j + before.reversible_energy_j;
    let total_after = mechanical.report.kinetic_after_j + after.reversible_energy_j;
    let storage_change = total_after - total_before;
    let residual = storage_change + mechanical.report.dissipated_energy_j - parameter_work;
    for balance in [
        potential_change - refreshed.potential_change_j(),
        potential_change - parameter_work,
        residual,
        mechanical.report.energy_residual_j + potential_change - parameter_work,
    ] {
        if !balance.is_finite() || balance.abs() > controls.energy_tolerance_j {
            return Err(RigidAdhesiveThermalNormalError::Energy {
                residual_j: balance,
            });
        }
    }
    poll(gate).map_err(RigidAdhesiveError::from)?;
    Ok(RigidAdhesiveThermalNormalCandidate {
        bodies_after: mechanical.bodies_after,
        normal_source: refreshed.into_source(),
        normal_state: after.state(),
        shear_source: shear,
        report: RigidAdhesiveThermalNormalReport {
            mechanical: mechanical.report,
            normal_before: before,
            normal_after: after,
            normal_potential_change_j: potential_change,
            parameter_work_j: parameter_work,
            represented_storage_change_j: storage_change,
            energy_residual_j: residual,
        },
    })
}

/// Publish both native momenta, the normal branch and both source owners together.
///
/// Source, geometry, work, cancellation and stopping refusals leave every input
/// unchanged. Temperature has been explicitly proposed by the normal owner;
/// this entry supplies no implicit caloric solve or temperature-driven expansion.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_adhesive_with_normal_refresh(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    binding: RigidContactBinding,
    frame: TangentFrame,
    normal: &mut SourceJkrSphereContact,
    state: &mut JkrState,
    shear: &mut SourceAdhesiveShear,
    refreshed: SourceJkrThermalRefresh,
    dt_s: f64,
    controls: RigidAdhesiveControls,
    gate: &CancelGate,
) -> Result<RigidAdhesiveThermalNormalReport, RigidAdhesiveThermalNormalError> {
    let source = shear.sliding_patch_jkr(normal.normal_state(), refreshed.before())?;
    let candidate = prepare_rigid_adhesive_with_normal_refresh(
        bodies,
        properties,
        RigidAdhesiveContact {
            binding,
            frame,
            source,
        },
        normal,
        *state,
        refreshed,
        dt_s,
        controls,
        gate,
    )?;
    poll(gate).map_err(RigidAdhesiveError::from)?;
    bodies[binding.body_a] = candidate.bodies_after[0];
    bodies[binding.body_b] = candidate.bodies_after[1];
    *normal = candidate.normal_source;
    *state = candidate.normal_state;
    *shear = candidate.shear_source;
    Ok(candidate.report)
}
