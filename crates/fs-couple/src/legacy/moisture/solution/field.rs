//! Actual native specimen endpoints for the original finite-solution owner.
//!
//! A heterogeneous body has several original material bindings and storage
//! measures. This private interface carries its physical inventories directly;
//! it does not reinterpret that body as one homogeneous MoistureField.

use super::*;
use crate::moisture::reservoir::water_sum;

pub(super) trait SolutionParticipant {
    type Field;
    type Report;

    fn id(&self) -> u64;
    fn time_s(&self) -> f64;
    fn temperature_kelvin(&self) -> f64;
    fn cells(&self) -> usize;
    fn sealed(&self) -> bool;

    /// Source-specific admission that must precede all scalar volume trials.
    fn admit(&self, _cx: &Cx<'_>) -> Result<(), SolutionMoistureError> {
        Ok(())
    }

    /// Every trial starts from this original accepted body. Its selected
    /// native owner retains all physical boundary/interface/source laws.
    fn advance(
        &self,
        dt_s: f64,
        temperature: f64,
        activity: f64,
        cx: &Cx<'_>,
    ) -> Result<(Self::Field, Self::Report), SolutionMoistureError>;

    /// Compensated local physical water changes; opposing exchanges survive.
    fn water_change(&self, next: &Self::Field) -> f64;

    /// The original or final actual stored water, independently of the change.
    fn water_mass_kg(&self, endpoint: Option<&Self::Field>) -> f64;
}

impl<F: ReservoirField<N>, const N: usize> SolutionParticipant for ReservoirParticipant<'_, F, N> {
    type Field = F;
    type Report = F::Report;

    fn id(&self) -> u64 {
        self.id
    }
    fn time_s(&self) -> f64 {
        self.field.inventory().time_s
    }
    fn temperature_kelvin(&self) -> f64 {
        self.field.inventory().temperature_kelvin()
    }
    fn cells(&self) -> usize {
        self.field.inventory().moisture.len()
    }
    fn sealed(&self) -> bool {
        self.faces.iter().all(Option::is_none)
    }
    fn advance(
        &self,
        dt_s: f64,
        temperature: f64,
        activity: f64,
        cx: &Cx<'_>,
    ) -> Result<(F, F::Report), SolutionMoistureError> {
        Ok(self.field.advance_activity(
            dt_s,
            self.faces,
            temperature,
            activity,
            self.transport,
            cx,
        )?)
    }
    fn water_change(&self, next: &F) -> f64 {
        self.field.water_change_to(next)
    }
    fn water_mass_kg(&self, endpoint: Option<&F>) -> f64 {
        let inventory = endpoint.unwrap_or(self.field).inventory();
        water_sum(
            inventory
                .moisture
                .iter()
                .enumerate()
                .map(|(i, m)| m * inventory.dry_mass_at(i)),
        )
    }
}
