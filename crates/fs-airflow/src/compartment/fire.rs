//! Prescribed external heat release coupled to original-source gas inventories.
//!
//! This helper transfers actual lower-layer species and enthalpy, supplies
//! convective power to the upper layer, and returns radiant power separately.
//! It infers no fuel mass, reaction, product yields or artificial upper seed.
//! It is not the inventory-limited chemistry owner: adding that owner's released
//! heat again here would double-count chemical energy. Heskestad's unconfined
//! correlation remains nominal; no weak-plume penetration, wall/corner or
//! oxygen-limitation correction is implemented. Rates are frozen at evaluation.

use super::{Layer, LayerRates, TwoZone, ZoneError};
use fs_convection::fire::{
    FireAmbientV1, FireCorrelationErrorV1, FireEstimateV1, HeskestadEntrainmentV1, HeskestadFireV1,
    HeskestadRegionV1, SpecificHeatCapacityV1,
};
use fs_qty::{Acceleration, Density, Length, MassFlowRate, Power};

/// Explicit external HRR and geometry; no default convective fraction is used.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrescribedFire {
    /// Total externally supplied heat release [W], nonnegative.
    pub total_power: Power,
    /// Convective part [W], between zero and total_power.
    pub convective_power: Power,
    /// Positive equivalent source diameter [m].
    pub diameter: Length,
    /// Source base above the floor [m], in [0, room height).
    pub base_height: Length,
}

/// Executed interlayer source branch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PrescribedFireRegime {
    /// Total and convective HRR are exactly zero.
    Dormant,
    /// The source is at/above the interface or no lower layer exists.
    InUpperLayer,
    /// Actual lower-to-upper transfer with the original Heskestad branch.
    Entrainment(HeskestadRegionV1),
}

/// Conservative gas rates and separately accounted radiant power.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrescribedFireReport {
    /// Canonical lower/upper rates; their species sums are exactly zero.
    pub rates: [LayerRates; 2],
    /// Total minus convective power [W], not deposited into gas by this helper.
    pub radiant_power: Power,
    /// Entrained lower-layer gas [kg/s], excluding any external fuel source.
    pub plume_mass_rate: MassFlowRate,
    /// Explicit executed branch, including absent-layer handling.
    pub regime: PrescribedFireRegime,
    /// Original nominal estimate and source card, only if evaluated.
    pub plume_estimate: Option<FireEstimateV1<HeskestadEntrainmentV1>>,
    /// Actual original-source lower-layer properties used by that estimate.
    pub plume_ambient: Option<FireAmbientV1>,
}

/// A refusal returns no source increment and never changes the room.
#[derive(Clone, Debug, PartialEq)]
pub enum PrescribedFireError {
    /// Invalid power partition, geometry or gravity.
    Invalid(&'static str),
    /// The original fire correlation refused an input/regime.
    Correlation(FireCorrelationErrorV1),
    /// The original thermochemical or extensive gas owner refused a state.
    Gas(ZoneError),
    /// Cancellation before the complete report was returned.
    Cancelled,
}
impl core::fmt::Display for PrescribedFireError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "prescribed compartment fire: {self:?}")
    }
}
impl std::error::Error for PrescribedFireError {}
impl From<ZoneError> for PrescribedFireError {
    fn from(value: ZoneError) -> Self {
        Self::Gas(value)
    }
}
impl From<FireCorrelationErrorV1> for PrescribedFireError {
    fn from(value: FireCorrelationErrorV1) -> Self {
        Self::Correlation(value)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PrescribedFireError> {
    if cancelled() {
        Err(PrescribedFireError::Cancelled)
    } else {
        Ok(())
    }
}

impl PrescribedFire {
    /// Evaluate against current NASA calorics and actual lower-layer density.
    ///
    /// Entrainment is removed and added using one identical carried enthalpy
    /// rate; TwoZone subsequently supplies the equal/opposite interface work.
    /// An absent upper layer may form from this real transfer. No absent lower
    /// layer is seeded or queried for density. Receiver density is not used:
    /// the CFAST weak-plume penetration constraint is outside this closure.
    pub fn evaluate(
        &self,
        room: &TwoZone,
        gravity: Acceleration,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PrescribedFireReport, PrescribedFireError> {
        poll(&mut cancelled)?;
        let (total, convective, diameter, base, g) = (
            self.total_power.value(),
            self.convective_power.value(),
            self.diameter.value(),
            self.base_height.value(),
            gravity.value(),
        );
        if !total.is_finite()
            || total < 0.0
            || !convective.is_finite()
            || convective < 0.0
            || convective > total
        {
            return Err(PrescribedFireError::Invalid(
                "require finite 0 <= convective <= total power",
            ));
        }
        if !diameter.is_finite()
            || diameter <= 0.0
            || !g.is_finite()
            || g <= 0.0
            || !base.is_finite()
            || base < 0.0
            || base >= room.height().value()
        {
            return Err(PrescribedFireError::Invalid(
                "source diameter, base height or gravity",
            ));
        }
        let mut report = PrescribedFireReport {
            rates: [LayerRates::default(); 2],
            radiant_power: Power::new(total - convective),
            plume_mass_rate: MassFlowRate::new(0.0),
            regime: PrescribedFireRegime::Dormant,
            plume_estimate: None,
            plume_ambient: None,
        };
        if total > 0.0 {
            report.rates[Layer::Upper.index()].add_heat(self.convective_power)?;
            let lower = room.layer(Layer::Lower);
            let mass = lower.mass_kg();
            let relative_height = room.interface_height().value() - base;
            if mass == 0.0 || relative_height <= 0.0 {
                report.regime = PrescribedFireRegime::InUpperLayer;
            } else {
                let density = lower
                    .density_kg_m3()
                    .ok_or(PrescribedFireError::Invalid("absent plume donor"))?;
                let (mut cp, mut correction) = (0.0, 0.0);
                for (mass_i, model) in lower
                    .species_masses()
                    .iter()
                    .zip(room.gas().source().models())
                {
                    poll(&mut cancelled)?;
                    if *mass_i == 0.0 {
                        continue;
                    }
                    let properties = model
                        .evaluate_properties(lower.temperature())
                        .map_err(|e| PrescribedFireError::Gas(ZoneError::Gas(e.into())))?;
                    let term =
                        (*mass_i / mass) * properties.cp().value() / model.molar_mass().value();
                    let adjusted = term - correction;
                    let next = cp + adjusted;
                    correction = (next - cp) - adjusted;
                    cp = next;
                }
                let ambient = FireAmbientV1::new(
                    lower.temperature(),
                    Density::new(density),
                    SpecificHeatCapacityV1::new(cp),
                    gravity,
                )?;
                let fire = HeskestadFireV1::new(
                    self.total_power,
                    self.convective_power,
                    self.diameter,
                    ambient,
                )?;
                let estimate = fire.entrainment(Length::new(relative_height))?;
                let entrained = estimate.value();
                let carried = room.carried_rate(Layer::Lower, entrained.mass_rate)?;
                report.rates[Layer::Lower.index()] = carried.opposite();
                report.rates[Layer::Upper.index()].add(&carried)?;
                report.plume_mass_rate = entrained.mass_rate;
                report.regime = PrescribedFireRegime::Entrainment(entrained.region);
                report.plume_estimate = Some(estimate);
                report.plume_ambient = Some(ambient);
            }
        }
        poll(&mut cancelled)?;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
