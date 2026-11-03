//! Bond power estimates for independently admitted continuous grinding states.
//!
//! Each call solves the original continuous circuit before evaluating its fresh
//! feed F80 and combined external overflow-plus-purge P80. Specific energy comes
//! from the original [`BondGrindingCard`]; power uses fresh dry tonnage, without
//! charging circulating solids or carrier water as new ore. See GMG01-MP-2021,
//! sections 4–5, cited by [`crate::comminution`]. The returned rates are on-stream
//! rates; reporting-period utilization does not change source admission.
//!
//! Fixed dry hold-up can change residence and the discrete product quantile as
//! throughput changes. No invariant-PSD scaling or inverse power/capacity solve
//! is used. A discrete P80 jump need not give a smooth monotone power curve.
//! Selection frequencies, daughter distributions and work-index test support
//! remain separate original inputs. This adapter supplies no rate–power
//! calibration, drive efficiency, equipment power-draw law or auxiliary load.

use super::holdup::{HoldupError, HoldupSolution, HoldupSolveSettings};
use super::{ContinuousCircuitError, ContinuousCircuitSolution, ContinuousGrindingCircuit};
use crate::comminution::circuit::CircuitSettings;
use crate::comminution::{BondGrindingCard, ComminutionError, SpecificGrindingEnergy};
use crate::separation::{ParticleStream, Partition, SeparationError};
use fs_matdb::QueryPoint;
use fs_qty::{Dimensionless, Length, Mass, Power, Time};

/// Original physical/source refusal or unrepresentable estimated power.
#[derive(Clone, Debug, PartialEq)]
pub enum ContinuousBondPowerError {
    /// The original prescribed-residence circuit refused the requested state.
    Circuit(ContinuousCircuitError),
    /// The original fixed-hold-up solve refused its state or numerical settings.
    Holdup(HoldupError),
    /// The original Bond card refused the actual size reduction or source domain.
    Bond(ComminutionError),
    /// The original particle-stream owner refused a mixture or size quantile.
    Stream(SeparationError),
    /// A positive finite power cannot be represented from the admitted inputs.
    Arithmetic(&'static str),
}
impl core::fmt::Display for ContinuousBondPowerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "continuous grinding Bond power: {self:?}")
    }
}
impl std::error::Error for ContinuousBondPowerError {}
impl From<ContinuousCircuitError> for ContinuousBondPowerError {
    fn from(value: ContinuousCircuitError) -> Self {
        Self::Circuit(value)
    }
}
impl From<HoldupError> for ContinuousBondPowerError {
    fn from(value: HoldupError) -> Self {
        Self::Holdup(value)
    }
}
impl From<ComminutionError> for ContinuousBondPowerError {
    fn from(value: ComminutionError) -> Self {
        Self::Bond(value)
    }
}
impl From<SeparationError> for ContinuousBondPowerError {
    fn from(value: SeparationError) -> Self {
        Self::Stream(value)
    }
}
type Result<T> = core::result::Result<T, ContinuousBondPowerError>;

/// Original continuous circuit and independently declared work-index application.
#[derive(Clone, Debug, PartialEq)]
pub struct ContinuousBondCircuit {
    circuit: ContinuousGrindingCircuit,
    bond: BondGrindingCard,
}

#[derive(Clone, Debug, PartialEq)]
enum Response {
    Residence(ContinuousCircuitSolution),
    DryHoldup(HoldupSolution),
}
impl Response {
    fn solution(&self) -> &ContinuousCircuitSolution {
        match self {
            Self::Residence(solution) => solution,
            Self::DryHoldup(result) => &result.solution,
        }
    }
}

/// Energy estimate paired with its original, independently solved operating state.
#[derive(Clone, Debug, PartialEq)]
pub struct ContinuousBondCircuitPoint {
    model: ContinuousBondCircuit,
    response: Response,
    /// Estimated on-stream net pinion power, excluding auxiliaries.
    pub net_pinion_power: Power,
    /// Original Bond energy in joules per kilogram of fresh dry solids.
    pub specific_energy: SpecificGrindingEnergy,
    /// Actual discrete 80%-passing diameter of fresh feed, not circulating feed.
    pub fresh_f80: Length,
    /// Actual discrete 80%-passing diameter of overflow and purge combined.
    pub external_p80: Length,
}
impl ContinuousBondCircuitPoint {
    /// Original physical models and work-index evidence.
    pub fn model(&self) -> &ContinuousBondCircuit {
        &self.model
    }
    /// Original admitted streams, class hold-ups, source query and closure result.
    pub fn solution(&self) -> &ContinuousCircuitSolution {
        self.response.solution()
    }
    /// Complete original fixed-hold-up diagnostics when that control was requested.
    pub fn holdup_solution(&self) -> Option<&HoldupSolution> {
        match &self.response {
            Response::DryHoldup(result) => Some(result),
            Response::Residence(_) => None,
        }
    }
}

impl ContinuousBondCircuit {
    /// Retain original models. Their actual operating point is admitted per call.
    /// No public, caller-edited circuit result is accepted as admission evidence.
    pub const fn new(circuit: ContinuousGrindingCircuit, bond: BondGrindingCard) -> Self {
        Self { circuit, bond }
    }
    /// Original continuous selection/breakage and classifier sources.
    pub const fn circuit(&self) -> &ContinuousGrindingCircuit {
        &self.circuit
    }
    /// Original work index, test preparation and size-support declaration.
    pub const fn bond(&self) -> &BondGrindingCard {
        &self.bond
    }

    /// Solve an original prescribed-residence state and estimate its Bond power.
    pub fn at_residence_time(
        &self,
        fresh: &ParticleStream,
        residence: Time,
        mill_conditions: &QueryPoint,
        classifier_conditions: &QueryPoint,
        settings: CircuitSettings,
        underflow_recycle: Partition,
    ) -> Result<ContinuousBondCircuitPoint> {
        let solution = self.circuit.at_residence_time(
            fresh,
            residence,
            mill_conditions,
            classifier_conditions,
            settings,
            underflow_recycle,
        )?;
        self.estimate(Response::Residence(solution))
    }

    /// Solve residence at the requested dry inventory and recompute actual F80/P80.
    /// Original bracket, mass residual, work and source diagnostics are retained.
    pub fn at_dry_holdup(
        &self,
        fresh: &ParticleStream,
        target: Mass,
        mill_conditions: &QueryPoint,
        classifier_conditions: &QueryPoint,
        circuit_settings: CircuitSettings,
        underflow_recycle: Partition,
        settings: HoldupSolveSettings,
    ) -> Result<ContinuousBondCircuitPoint> {
        let solution = self.circuit.at_dry_holdup(
            fresh,
            target,
            mill_conditions,
            classifier_conditions,
            circuit_settings,
            underflow_recycle,
            settings,
        )?;
        self.estimate(Response::DryHoldup(solution))
    }

    fn estimate(&self, response: Response) -> Result<ContinuousBondCircuitPoint> {
        let flows = &response.solution().circuit;
        let external = ParticleStream::mix(&[&flows.product, &flows.purge])?;
        let fraction = Dimensionless::new(0.8);
        let fresh_f80 = flows.fresh_feed.discrete_passing_size(fraction)?;
        let external_p80 = external.discrete_passing_size(fraction)?;
        let specific_energy = self.bond.specific_energy(fresh_f80, external_p80)?;
        let power = flows.fresh_feed.dry_solids()?.value() * specific_energy.value();
        if !power.is_finite() || power <= 0.0 {
            return Err(ContinuousBondPowerError::Arithmetic(
                "fresh dry flow times Bond energy is unrepresentable",
            ));
        }
        Ok(ContinuousBondCircuitPoint {
            model: self.clone(),
            response,
            net_pinion_power: Power::new(power),
            specific_energy,
            fresh_f80,
            external_p80,
        })
    }
}

#[cfg(test)]
#[path = "power/tests.rs"]
mod tests;
