//! Exact thermal continuation of an unchanged original Mindlin constitutive law.
//!
//! Original bulk and interface sources may be re-resolved at independent actual
//! temperatures without erasing retained shear memory. This owner proves that
//! the selected constitutive coefficients are constant throughout each changed
//! temperature interval. It does not infer this from equal endpoint samples:
//! every changed source bundle must use the same scalar ConstantWithinValidity
//! claims, whose admitted validity boxes contain the whole intervening segment.
//!
//! Normal load/geometry, elasticity, friction, material roles, non-temperature
//! conditions and all source selections remain fixed. Under those conditions
//! the original analytic spring history and stored energy are unchanged. A
//! temperature-dependent coefficient needs a physical thermoelastic/plastic
//! history-transfer owner and refuses here, even when its endpoints are equal.

use super::{SourceMindlinContact, SourceMindlinError, SourceMindlinState};
use fs_matdb::EvaluationDecision;
use fs_material::state_point::{
    ResolvedInterfaceStatePoint, ResolvedMaterialStatePoint, ResolvedScalarProperty,
};
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{QuantitySpec, Temperature};
use fs_tribo::TangentialLaw;
use std::collections::BTreeMap;

/// Source continuation either preserves the complete checkpoint or publishes none.
#[derive(Clone, Debug, PartialEq)]
pub enum SourceMindlinTemperatureError {
    /// Original source/history, frozen geometry or thermal-path admission failed.
    Source(SourceMindlinError),
    /// Cancellation observed before returning a complete transferred history.
    Cancelled,
}
impl core::fmt::Display for SourceMindlinTemperatureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source Mindlin thermal continuation: {self:?}")
    }
}
impl std::error::Error for SourceMindlinTemperatureError {}

impl SourceMindlinContact {
    /// Preserve this accepted history under an already resolved thermal endpoint.
    ///
    /// Both bulk temperatures and the interface temperature may differ and move
    /// independently. The caller resolves `refreshed` through the original
    /// material owners; every original full card, selected property claim,
    /// quantity convention, non-T coordinate, normal model and operating envelope
    /// must survive unchanged. A changed temperature admits only scalar claims
    /// explicitly constant throughout their validity boxes. Equal endpoints of
    /// a curve or surface cannot prove that the intervening physics was constant.
    ///
    /// The native law must be exactly unchanged, including N, patch radius, both
    /// G/nu pairs, mu and history capacity. The candidate retains the complete
    /// original analytic memory with the refreshed source identities. This is a
    /// zero-motion, zero-work continuation: no contact solve, normal update,
    /// thermal expansion, coefficient-work correction or heat deposit occurs.
    /// Cancellation/refusal leaves both sources and the accepted history intact.
    pub fn continue_at_temperature(
        &self,
        state: &SourceMindlinState,
        refreshed: &Self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceMindlinState, SourceMindlinTemperatureError> {
        poll(cancelled)?;
        if !self.stored_energy_j(state).is_finite() {
            return Err(invalid(
                "thermal continuation requires the original Mindlin history",
            ));
        }
        let old_normal = self.normal_source();
        let new_normal = refreshed.normal_source();
        let old_model = old_normal.model();
        let new_model = new_normal.model();
        let old = old_model.normal_state();
        let new = new_model.normal_state();
        if self.law() != refreshed.law()
            || self.normal_response().approach_m.to_bits()
                != refreshed.normal_response().approach_m.to_bits()
            || self.normal_response().reversible_energy_j.to_bits()
                != refreshed.normal_response().reversible_energy_j.to_bits()
            || old_normal.geometry() != new_normal.geometry()
            || old_normal.envelope() != new_normal.envelope()
            || old_model.model_card() != new_model.model_card()
            || old_model.law() != new_model.law()
            || old.elastic().interface() != new.elastic().interface()
        {
            return Err(invalid(
                "thermal continuation must preserve original geometry, load, law and interface",
            ));
        }
        for (a, b) in [
            (
                old.elastic().surface_a_state().resolved(),
                new.elastic().surface_a_state().resolved(),
            ),
            (
                old.elastic().surface_b_state().resolved(),
                new.elastic().surface_b_state().resolved(),
            ),
        ] {
            material(a, b, cancelled)?;
        }
        interface(old.interface_state(), new.interface_state(), cancelled)?;
        interface(
            self.friction_source(),
            refreshed.friction_source(),
            cancelled,
        )?;
        let candidate = SourceMindlinState {
            normal_model_identity: new_model.identity(),
            friction_state_identity: refreshed.friction_source().identity(),
            constitutive: state.constitutive.clone(),
        };
        poll(cancelled)?;
        Ok(candidate)
    }
}

fn material(
    before: &ResolvedMaterialStatePoint,
    after: &ResolvedMaterialStatePoint,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), SourceMindlinTemperatureError> {
    if before.card_identity() != after.card_identity() || before.material() != after.material() {
        return Err(invalid(
            "thermal continuation requires each original ordered full bulk card",
        ));
    }
    let changes = coordinates(
        before.query_point(),
        after.query_point(),
        before.axis_quantities(),
        after.axis_quantities(),
    )?;
    properties(before.properties(), after.properties(), changes, cancelled)
}

fn interface(
    before: &ResolvedInterfaceStatePoint,
    after: &ResolvedInterfaceStatePoint,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), SourceMindlinTemperatureError> {
    if before.card_identity() != after.card_identity()
        || before.surface_materials() != after.surface_materials()
    {
        return Err(invalid(
            "thermal continuation requires the original ordered interface card",
        ));
    }
    let changes = coordinates(
        before.query_point(),
        after.query_point(),
        before.axis_quantities(),
        after.axis_quantities(),
    )?;
    properties(before.properties(), after.properties(), changes, cancelled)
}

fn coordinates(
    before: &[(String, f64)],
    after: &[(String, f64)],
    old_axes: &BTreeMap<String, QuantitySpec>,
    new_axes: &BTreeMap<String, QuantitySpec>,
) -> Result<bool, SourceMindlinTemperatureError> {
    if before.len() != after.len()
        || old_axes != new_axes
        || old_axes.get("T").is_some_and(|q| {
            q.dims() != Temperature::DIMS
                || q.semantic_type().is_some_and(|t| {
                    t.kind() != QuantityKind::AbsoluteTemperature || t.form() != ValueForm::Static
                })
        })
    {
        return Err(invalid(
            "thermal continuation requires the original kelvin coordinate schemas",
        ));
    }
    let mut temperature = None;
    for ((axis, a), (other, b)) in before.iter().zip(after) {
        if axis != other
            || !a.is_finite()
            || !b.is_finite()
            || (axis != "T" && a.to_bits() != b.to_bits())
        {
            return Err(invalid(
                "thermal continuation cannot change an unowned source coordinate",
            ));
        }
        if axis == "T" {
            if *a <= 0.0 || *b <= 0.0 {
                return Err(invalid(
                    "thermal continuation requires positive absolute temperatures",
                ));
            }
            temperature = Some(a.to_bits() != b.to_bits());
        }
    }
    temperature.ok_or_else(|| invalid("thermal continuation requires every original T coordinate"))
}

fn properties(
    before: &[ResolvedScalarProperty],
    after: &[ResolvedScalarProperty],
    temperature_changes: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), SourceMindlinTemperatureError> {
    if before.len() != after.len() {
        return Err(invalid(
            "thermal continuation cannot drop original property claims",
        ));
    }
    for (a, b) in before.iter().zip(after) {
        poll(cancelled)?;
        if a.requirement() != b.requirement()
            || a.answer().receipt.selected != b.answer().receipt.selected
            || a.value_si().to_bits() != b.value_si().to_bits()
        {
            return Err(invalid(
                "thermal continuation requires unchanged original selected properties",
            ));
        }
        // ConstantWithinValidity is an explicit scalar claim on a convex box,
        // not the scalar value obtained by evaluating a nonlinear source. Both
        // endpoints were admitted by their original selected source owners.
        if temperature_changes
            && [a, b].into_iter().any(|p| {
                !p.source_axes().is_empty()
                    || p.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity
            })
        {
            return Err(invalid(
                "thermal continuation needs source-constant coefficients over the full path",
            ));
        }
    }
    Ok(())
}

fn invalid(reason: &'static str) -> SourceMindlinTemperatureError {
    SourceMindlinTemperatureError::Source(SourceMindlinError::Invalid(reason))
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SourceMindlinTemperatureError> {
    if cancelled() {
        Err(SourceMindlinTemperatureError::Cancelled)
    } else {
        Ok(())
    }
}
