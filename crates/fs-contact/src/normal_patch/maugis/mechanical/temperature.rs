//! Original-source temperature changes at unchanged admitted mechanical geometry.

use super::{MaugisAttachedPotential, MaugisContinuationError};
use crate::normal_patch::adhesive::temperature::{
    AdhesiveTemperatureError, AdhesiveThermalControl, SourceMaugisThermalRefresh,
};

impl MaugisAttachedPotential {
    /// Stage a sealed original-source refresh at the actual held indentation.
    ///
    /// The original work/energy budgets and explicit contact/noncontact policy
    /// remain unchanged. A selected noncontact continuation is rebuilt from the
    /// new finite-zone equations; losing that branch refuses even if a contact
    /// equilibrium still exists. This never infers attachment or detachment.
    ///
    /// The receipt's predecessor must be the canonical original response at the
    /// supplied geometry. Original source equilibria retain their resolved
    /// radius, whose reported indentation can differ by roundoff from that
    /// geometry. The successor is evaluated at the unchanged actual geometry,
    /// and its potential change must agree with the original signed parameter
    /// work under the existing contact-sized energy budget.
    ///
    /// This returns a candidate without mutating either source. Its signed
    /// parameter work belongs to the prescribed temperature driver; no caloric
    /// heat, mechanical impulse or loss is inferred.
    pub fn with_held_indentation_refresh(
        &self,
        indentation_m: f64,
        refreshed: &SourceMaugisThermalRefresh,
    ) -> Result<Self, AdhesiveTemperatureError> {
        let source = refreshed.source();
        if refreshed.original_model_identity() != self.identity()
            || refreshed.control() != AdhesiveThermalControl::HeldIndentation
            || refreshed.mechanical_work_j() != 0.0
            || source.effective_radius_m() != self.source().effective_radius_m()
            || source.interaction_range_m() != self.source().interaction_range_m()
            || source.envelope() != self.source().envelope()
        {
            return Err(AdhesiveTemperatureError::Invalid(
                "Maugis potential refresh must retain original source, held control and geometry",
            ));
        }
        let before = self.evaluate(indentation_m)?;
        if *refreshed.before() != before
            || refreshed.response().phase() != before.phase()
        {
            return Err(AdhesiveTemperatureError::Invalid(
                "Maugis potential refresh requires the actual predecessor and selected phase",
            ));
        }
        // Do not reinvert a rounded equilibrium coordinate on a no-op. The
        // original receipt intentionally preserves its exact accepted radius.
        if source.identity() == self.identity() {
            if *refreshed.response() != before
                || refreshed.parameter_work_j() != 0.0
                || refreshed.potential_change_j() != 0.0
            {
                return Err(AdhesiveTemperatureError::Invalid(
                    "unchanged Maugis source cannot alter equilibrium or parameter work",
                ));
            }
            return Ok(self.clone());
        }

        let mut candidate = Self::new(source.clone(), self.budget())?;
        if self.noncontact_continuation().is_some() {
            candidate = candidate.with_noncontact_continuation()?;
        }
        let canonical = candidate.evaluate(before.indentation_m())?;
        let after = candidate.evaluate(indentation_m)?;
        if canonical != *refreshed.response() || after.phase() != before.phase() {
            return Err(AdhesiveTemperatureError::Invalid(
                "Maugis refresh must retain its original successor and selected physical phase",
            ));
        }
        let work = refreshed.parameter_work_j();
        let residual = (after.potential_j() - before.potential_j()) - work;
        let scale = before
            .potential_j()
            .abs()
            .max(after.potential_j().abs())
            .max(work.abs());
        let tolerance = self.budget().energy_absolute_tolerance_j
            + self.budget().energy_relative_tolerance * scale;
        if !work.is_finite()
            || !residual.is_finite()
            || !tolerance.is_finite()
            || residual.abs() > tolerance
            || (refreshed.potential_change_j() - work).abs() > tolerance
        {
            return Err(MaugisContinuationError::EnergyClosure {
                residual_j: residual,
                tolerance_j: tolerance,
            }
            .into());
        }
        Ok(candidate)
    }
}
