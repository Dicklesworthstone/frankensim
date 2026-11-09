//! Original-source thermal parameter changes with retained rolling memory.
//!
//! The held approach and scalar angle do no mechanical work. Independently
//! refreshed bulk and interface temperatures change the JKR potential and
//! rolling spring coefficient; both signed changes remain explicit parameter
//! work. The independently calibrated critical shift must remain valid at the
//! new temperatures. This endpoint operation supplies neither a temperature
//! path stability certificate nor closed thermoelastic feedback or heat.

use super::{SourceJkrRollingContact, SourceJkrRollingError, SourceJkrRollingState, identity};
use crate::normal_patch::adhesive::temperature::{AdhesiveThermalControl, SourceJkrThermalRefresh};
use crate::normal_patch::adhesive::{JkrBranch, JkrResponse};
use fs_blake3::ContentHash;
use fs_math::dd::Dd;
use fs_tribo::rolling_loss::adhesive::AdhesiveRollingError;
use std::sync::Arc;

/// Signed storage changes at held approach, elastic shift and total angle.
/// Normal adhesion potential remains separate from rolling spring storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceJkrRollingTemperatureReport {
    /// Complete original rolling owner, including source and calibration.
    pub original_model_identity: ContentHash,
    /// Complete refreshed rolling owner with the same calibration and limits.
    pub refreshed_model_identity: ContentHash,
    pub before: JkrResponse,
    pub after: JkrResponse,
    pub normal_storage_before_j: f64,
    pub normal_storage_after_j: f64,
    pub rolling_storage_before_j: f64,
    pub rolling_storage_after_j: f64,
    /// Signed original normal potential change; normal actuator work is zero.
    pub normal_parameter_work_j: f64,
    /// Signed spring-coefficient work at the exact retained elastic shift.
    pub rolling_parameter_work_j: f64,
    /// Sum of normal and rolling parameter work, never dissipated heat.
    pub parameter_work_j: f64,
}

/// Sealed original history and refreshed source/history pair. Acceptance must
/// publish the source and state together; editing a report cannot alter either.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceJkrRollingTemperatureStep {
    previous: SourceJkrRollingState,
    source: SourceJkrRollingContact,
    state: SourceJkrRollingState,
    report: SourceJkrRollingTemperatureReport,
}

impl SourceJkrRollingTemperatureStep {
    pub const fn previous_state(&self) -> SourceJkrRollingState {
        self.previous
    }
    pub const fn source(&self) -> &SourceJkrRollingContact {
        &self.source
    }
    pub const fn state(&self) -> SourceJkrRollingState {
        self.state
    }
    pub const fn report(&self) -> &SourceJkrRollingTemperatureReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        SourceJkrRollingContact,
        SourceJkrRollingState,
        SourceJkrRollingTemperatureReport,
    ) {
        (self.source, self.state, self.report)
    }
}

impl SourceJkrRollingContact {
    /// Continue a loaded contact through an authenticated original normal
    /// temperature refresh at exactly held indentation and rolling angle.
    ///
    /// The original ordered cards, geometry, approximation envelope, critical
    /// displacement and small-shift ceiling remain fixed. The caller must
    /// establish validity of that unchanged yield calibration across the
    /// temperature change. The original normal and rolling applicability gates
    /// run again; branch loss refuses rather than discarding rolling memory.
    ///
    /// Complete double-double elastic/plastic history survives. This creates
    /// no rolling increment, plastic loss, snap release or caloric heat. The
    /// supplied refresh owns the material query and is moved into the successor
    /// without copying its source bundle. Failure leaves accepted state intact.
    pub fn continue_at_temperature(
        &self,
        previous: &SourceJkrRollingState,
        refresh: SourceJkrThermalRefresh,
    ) -> Result<SourceJkrRollingTemperatureStep, SourceJkrRollingError> {
        self.check(previous)?;
        let normal = refresh.source();
        let old_elastic = self.normal.normal_state().elastic();
        let elastic = normal.normal_state().elastic();
        // Both responses are sealed in their owners. Comparing the retained
        // canonical original endpoint avoids an unnecessary JKR root solve.
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
            return Err(SourceJkrRollingError::StateMismatch);
        }
        let after = *refresh.response();
        let rolling = self.law.refresh_material_response(
            &previous.constitutive,
            after.contact_radius_m,
            normal.normal_state().adhesion_energy_j_per_m2(),
            elastic.reduced_modulus_pa(),
        )?;
        let normal_parameter_work_j = refresh.parameter_work_j();
        let rolling_parameter_work_j = rolling.storage_change_j();
        let parameter_work_j = (Dd::from_f64(normal_parameter_work_j)
            + Dd::from_f64(rolling_parameter_work_j))
        .to_f64();
        if !parameter_work_j.is_finite() {
            return Err(AdhesiveRollingError::Unrepresentable {
                field: "normal and rolling thermal parameter work",
            }
            .into());
        }
        let identity = identity(normal, after.indentation_m, self.controls);
        let report = SourceJkrRollingTemperatureReport {
            original_model_identity: self.identity,
            refreshed_model_identity: identity,
            before: self.response,
            after,
            normal_storage_before_j: self.response.reversible_energy_j,
            normal_storage_after_j: after.reversible_energy_j,
            rolling_storage_before_j: rolling.storage_before_j(),
            rolling_storage_after_j: rolling.storage_after_j(),
            normal_parameter_work_j,
            rolling_parameter_work_j,
            parameter_work_j,
        };
        let (law, constitutive) = rolling.into_parts();
        Ok(SourceJkrRollingTemperatureStep {
            previous: *previous,
            source: Self {
                normal: Arc::new(refresh.into_source()),
                response: after,
                law,
                identity,
                controls: self.controls,
            },
            state: SourceJkrRollingState {
                identity,
                constitutive,
            },
            report,
        })
    }
}
