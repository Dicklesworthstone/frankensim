//! Original-source thermal parameter changes with retained bonded twist.
//!
//! The original normal refresh supplies independently queried bulk moduli,
//! interface adhesion and attached area at held approach. The signed normal
//! potential and torsional spring changes remain separate parameter work, not
//! plastic loss or heat. The independently calibrated critical twist must stay
//! valid, and the original mode-III/adhesion ceiling is checked again. This is
//! an open endpoint update, not a closed thermoelastic or mixed-mode solution.

use super::{
    SourceJkrTorsionContact, SourceJkrTorsionError, SourceJkrTorsionState, admit_release, identity,
};
use crate::normal_patch::adhesive::temperature::{AdhesiveThermalControl, SourceJkrThermalRefresh};
use crate::normal_patch::adhesive::{JkrBranch, JkrResponse};
use fs_blake3::ContentHash;
use fs_math::dd::Dd;
use fs_tribo::adhesive::torsion::AdhesiveTorsionError;
use std::sync::Arc;

/// Signed storage changes at held approach, elastic twist and total angle.
/// The original normal potential is never merged into torsional spring storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceJkrTorsionTemperatureReport {
    /// Complete original torsion owner, including source and calibration.
    pub original_model_identity: ContentHash,
    /// Complete refreshed torsion owner with the same calibration and limits.
    pub refreshed_model_identity: ContentHash,
    pub before: JkrResponse,
    pub after: JkrResponse,
    pub normal_storage_before_j: f64,
    pub normal_storage_after_j: f64,
    pub torsion_storage_before_j: f64,
    pub torsion_storage_after_j: f64,
    /// Signed original normal potential change; normal actuator work is zero.
    pub normal_parameter_work_j: f64,
    /// Signed spring-coefficient work at the exact retained elastic twist.
    pub torsion_parameter_work_j: f64,
    /// Sum of normal and torsion parameter work, never dissipated heat.
    pub parameter_work_j: f64,
}

/// Sealed original history and refreshed source/history pair. Acceptance must
/// publish source and state together; editing diagnostics cannot rebind memory.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceJkrTorsionTemperatureStep {
    previous: SourceJkrTorsionState,
    source: SourceJkrTorsionContact,
    state: SourceJkrTorsionState,
    report: SourceJkrTorsionTemperatureReport,
}

impl SourceJkrTorsionTemperatureStep {
    pub const fn previous_state(&self) -> SourceJkrTorsionState {
        self.previous
    }
    pub const fn source(&self) -> &SourceJkrTorsionContact {
        &self.source
    }
    pub const fn state(&self) -> SourceJkrTorsionState {
        self.state
    }
    pub const fn report(&self) -> &SourceJkrTorsionTemperatureReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        SourceJkrTorsionContact,
        SourceJkrTorsionState,
        SourceJkrTorsionTemperatureReport,
    ) {
        (self.source, self.state, self.report)
    }
}

impl SourceJkrTorsionContact {
    /// Continue loaded bonded twist through an authenticated original normal
    /// temperature refresh at exactly held indentation and angular coordinate.
    ///
    /// Original ordered cards, geometry, envelope and critical twist stay fixed.
    /// The caller establishes the unchanged yield calibration's validity at the
    /// new temperatures. Both bulk shear moduli are refreshed independently;
    /// the original maximum mode-III release/adhesion screen then runs using
    /// the new moduli, area and work of adhesion. An inadmissible contact refuses.
    ///
    /// Double-double elastic/plastic memory survives without a twist increment,
    /// plastic flow, snap release or manufactured heat. This operation does not
    /// certify an intervening temperature path or stress-free recruitment of
    /// new annuli. The refreshed source bundle is moved without a second clone;
    /// every refusal leaves the accepted source and state unchanged.
    pub fn continue_at_temperature(
        &self,
        previous: &SourceJkrTorsionState,
        refresh: SourceJkrThermalRefresh,
    ) -> Result<SourceJkrTorsionTemperatureStep, SourceJkrTorsionError> {
        self.check(previous)?;
        let normal = refresh.source();
        let old_elastic = self.normal.normal_state().elastic();
        let elastic = normal.normal_state().elastic();
        // The retained endpoint is already the original canonical JKR response;
        // exact comparison authenticates it without another normal root solve.
        if refresh.original_model_identity() != self.normal.identity()
            || refresh.before() != &self.response
            || refresh.control() != AdhesiveThermalControl::HeldIndentation
            || refresh.mechanical_work_j() != 0.0
            || refresh.before().branch != JkrBranch::Attached
            || refresh.response().branch != JkrBranch::Attached
            || refresh.before().indentation_m.to_bits() != self.response.indentation_m.to_bits()
            || refresh.response().indentation_m.to_bits() != self.response.indentation_m.to_bits()
            || refresh.response().normal_state_identity() != normal.normal_state().identity()
            || normal.effective_radius_m().to_bits() != self.normal.effective_radius_m().to_bits()
            || normal.interaction_range_m().to_bits() != self.normal.interaction_range_m().to_bits()
            || normal.envelope() != self.normal.envelope()
            || elastic.interface().card_identity() != old_elastic.interface().card_identity()
            || [elastic.surface_a_state(), elastic.surface_b_state()]
                .map(|state| state.resolved().card_identity())
                != [old_elastic.surface_a_state(), old_elastic.surface_b_state()]
                    .map(|state| state.resolved().card_identity())
        {
            return Err(SourceJkrTorsionError::StateMismatch);
        }
        let after = *refresh.response();
        let torsion = self.law.refresh_material_response(
            &previous.constitutive,
            after.contact_radius_m,
            [elastic.surface_a_state(), elastic.surface_b_state()]
                .map(|state| state.elastic_state().shear_modulus_pa()),
        )?;
        let maximum_mode_iii_release = admit_release(normal, torsion.law(), self.controls)?;
        let normal_parameter_work_j = refresh.parameter_work_j();
        let torsion_parameter_work_j = torsion.storage_change_j();
        let parameter_work_j = (Dd::from_f64(normal_parameter_work_j)
            + Dd::from_f64(torsion_parameter_work_j))
        .to_f64();
        if !parameter_work_j.is_finite() {
            return Err(AdhesiveTorsionError::Unrepresentable {
                field: "normal and torsion thermal parameter work",
            }
            .into());
        }
        let identity = identity(normal, after.indentation_m, self.controls);
        let report = SourceJkrTorsionTemperatureReport {
            original_model_identity: self.identity,
            refreshed_model_identity: identity,
            before: self.response,
            after,
            normal_storage_before_j: self.response.reversible_energy_j,
            normal_storage_after_j: after.reversible_energy_j,
            torsion_storage_before_j: torsion.storage_before_j(),
            torsion_storage_after_j: torsion.storage_after_j(),
            normal_parameter_work_j,
            torsion_parameter_work_j,
            parameter_work_j,
        };
        let (law, constitutive) = torsion.into_parts();
        Ok(SourceJkrTorsionTemperatureStep {
            previous: *previous,
            source: Self {
                normal: Arc::new(refresh.into_source()),
                response: after,
                law,
                identity,
                maximum_mode_iii_release,
                controls: self.controls,
            },
            state: SourceJkrTorsionState {
                identity,
                constitutive,
            },
            report,
        })
    }
}
