//! On-stream Bond estimates for an admitted fixed-card grinding circuit.
//!
//! Bond uses fresh dry-solids throughput and the fresh F80 / combined external
//! overflow-plus-purge P80. Circulating mill tonnage is not the energy basis.
//! See GMG01-MP-2021, *Determining the Bond Efficiency of Industrial Grinding
//! Circuits*, sections 4–5, as cited by [`crate::comminution`]. No measured plant
//! efficiency, power-draw law, or transferability of the caller's Wi is inferred.
//!
//! All returned rates are steady on-stream rates. A reporting-period utilization
//! must not reduce the running throughput used to admit a source's flow domain.

use super::{CircuitError, CircuitSettings, CircuitSolution, ClosedGrindingCircuit};
use crate::comminution::{BondGrindingCard, ComminutionError, SpecificGrindingEnergy};
use crate::separation::{ParticleStream, Partition, SeparationError};
use fs_matdb::QueryPoint;
use fs_qty::{Dimensionless, Length, MassFlowRate, Power};

/// Original-owner refusal or an unrepresentable power/throughput request.
#[derive(Clone, Debug, PartialEq)]
pub enum CircuitPowerError {
    /// The original circuit refused its prototype or scaled operating state.
    Circuit(CircuitError),
    /// The original work-index card refused the circuit sizes or power request.
    Bond(ComminutionError),
    /// Invalid or unrepresentable scaling, or a changed discrete size quantile.
    Invalid(&'static str),
}
impl core::fmt::Display for CircuitPowerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Bond circuit power: {self:?}")
    }
}
impl std::error::Error for CircuitPowerError {}
impl From<CircuitError> for CircuitPowerError {
    fn from(value: CircuitError) -> Self {
        Self::Circuit(value)
    }
}
impl From<ComminutionError> for CircuitPowerError {
    fn from(value: ComminutionError) -> Self {
        Self::Bond(value)
    }
}
impl From<SeparationError> for CircuitPowerError {
    fn from(value: SeparationError) -> Self {
        Self::Circuit(CircuitError::Separation(value))
    }
}
type Result<T> = core::result::Result<T, CircuitPowerError>;

/// Immutable admitted stream shape, source cards and declared circuit conditions.
///
/// Construction solves the original circuit; a publicly editable solution is
/// never accepted as proof of admission. Every scaled operating point solves
/// again, including the actual mill dry-flow and slurry-concentration domains.
#[derive(Clone, Debug, PartialEq)]
pub struct BondCircuitPrototype {
    circuit: ClosedGrindingCircuit,
    bond: BondGrindingCard,
    reference: CircuitSolution,
    mill_conditions: QueryPoint,
    classifier_conditions: QueryPoint,
    settings: CircuitSettings,
    fresh_f80: Length,
    external_p80: Length,
    specific_energy: SpecificGrindingEnergy,
}

/// Original evidence and the independently replayed operating point.
#[derive(Clone, Debug, PartialEq)]
pub struct BondCircuitPoint {
    /// Unchanged cards, reference admission, conditions and numerical settings.
    pub prototype: BondCircuitPrototype,
    /// Actual on-stream fresh, circulating, overflow and purge material accounts.
    pub solution: CircuitSolution,
    /// Estimated net power at the grinding drive pinion, excluding auxiliaries.
    pub net_pinion_power: Power,
    /// Original Bond card evaluated at this circuit's external size reduction.
    pub specific_energy: SpecificGrindingEnergy,
    /// Discrete 80%-passing size of the fresh feed; not circulating mill feed.
    pub fresh_f80: Length,
    /// Discrete 80%-passing size of combined external overflow and purge.
    pub external_p80: Length,
}

impl BondCircuitPrototype {
    /// Admit the reference state and caller-selected work-index application.
    ///
    /// The return partition is retained verbatim, including an independently
    /// supplied tiny purge. Size quantiles use representative class diameters,
    /// without interpolation or an inferred continuous particle distribution.
    pub fn new(
        circuit: &ClosedGrindingCircuit,
        bond: &BondGrindingCard,
        fresh: &ParticleStream,
        mill_conditions: &QueryPoint,
        classifier_conditions: &QueryPoint,
        settings: CircuitSettings,
        underflow_recycle: Partition,
    ) -> Result<Self> {
        let reference = circuit.solve_with_underflow_recycle(
            fresh,
            mill_conditions,
            classifier_conditions,
            settings,
            underflow_recycle,
        )?;
        let (fresh_f80, external_p80) = circuit_sizes(&reference)?;
        let specific_energy = bond.specific_energy(fresh_f80, external_p80)?;
        Ok(Self {
            circuit: circuit.clone(),
            bond: bond.clone(),
            reference,
            mill_conditions: mill_conditions.clone(),
            classifier_conditions: classifier_conditions.clone(),
            settings,
            fresh_f80,
            external_p80,
            specific_energy,
        })
    }

    /// Original mill/classifier cards with their evidence and applicability.
    pub const fn circuit(&self) -> &ClosedGrindingCircuit {
        &self.circuit
    }
    /// Original empirical work index, test context and size domain.
    pub const fn bond(&self) -> &BondGrindingCard {
        &self.bond
    }
    /// Independently admitted reference, including the exact return/purge split.
    pub const fn reference(&self) -> &CircuitSolution {
        &self.reference
    }
    /// Explicit caller conditions, in mill then classifier order.
    pub const fn conditions(&self) -> (&QueryPoint, &QueryPoint) {
        (&self.mill_conditions, &self.classifier_conditions)
    }
    /// Original closure allowance and circulating-load limit.
    pub const fn settings(&self) -> CircuitSettings {
        self.settings
    }

    /// Scale to the fresh on-stream capacity estimated from net pinion power.
    ///
    /// The original Bond owner computes throughput. All class masses, assayed
    /// elements and carrier water scale together; both process owners then
    /// evaluate their actual operating state before a result can be returned.
    pub fn at_net_pinion_power(&self, power: Power) -> Result<BondCircuitPoint> {
        let fresh = self.bond.dry_throughput(
            self.fresh_f80,
            self.external_p80,
            power,
            Dimensionless::new(1.0),
        )?;
        self.at_operating_point(fresh, power)
    }

    /// Estimate required net pinion power for positive fresh on-stream throughput.
    ///
    /// Power equals fresh dry flow times the original Bond specific energy.
    /// Neither circulating tonnage nor slurry water is charged as fresh ore.
    pub fn at_fresh_throughput(&self, fresh: MassFlowRate) -> Result<BondCircuitPoint> {
        positive(
            fresh.value(),
            "fresh throughput must be finite and positive",
        )?;
        let power = positive(
            self.specific_energy.value() * fresh.value(),
            "required pinion power is unrepresentable",
        )?;
        self.at_operating_point(fresh, Power::new(power))
    }

    fn at_operating_point(&self, fresh: MassFlowRate, power: Power) -> Result<BondCircuitPoint> {
        let factor = positive(
            fresh.value() / self.reference.fresh_feed.dry_solids()?.value(),
            "fresh-feed scale is unrepresentable",
        )?;
        let scaled = self
            .reference
            .fresh_feed
            .scaled(Dimensionless::new(factor))?;
        let solution = self.circuit.solve_with_underflow_recycle(
            &scaled,
            &self.mill_conditions,
            &self.classifier_conditions,
            self.settings,
            self.reference.underflow_recycle,
        )?;
        let (fresh_f80, external_p80) = circuit_sizes(&solution)?;
        if (fresh_f80, external_p80) != (self.fresh_f80, self.external_p80) {
            return Err(CircuitPowerError::Invalid(
                "scaling changed a discrete size quantile; Bond energy basis is not preserved",
            ));
        }
        let specific_energy = self.bond.specific_energy(fresh_f80, external_p80)?;
        Ok(BondCircuitPoint {
            prototype: self.clone(),
            solution,
            net_pinion_power: power,
            specific_energy,
            fresh_f80,
            external_p80,
        })
    }
}

fn circuit_sizes(solution: &CircuitSolution) -> Result<(Length, Length)> {
    let external = ParticleStream::mix(&[&solution.product, &solution.purge])?;
    let fraction = Dimensionless::new(0.8);
    Ok((
        solution.fresh_feed.discrete_passing_size(fraction)?,
        external.discrete_passing_size(fraction)?,
    ))
}
fn positive(value: f64, reason: &'static str) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(CircuitPowerError::Invalid(reason))
    }
}

#[cfg(test)]
#[path = "power/tests.rs"]
mod tests;
