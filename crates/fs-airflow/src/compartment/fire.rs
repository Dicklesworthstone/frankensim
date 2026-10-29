//! Prescribed external heat release coupled to original-source gas inventories.
//!
//! This helper transfers actual lower-layer species and enthalpy, supplies
//! convective power to the upper layer, and returns radiant power separately.
//! It infers no fuel mass, reaction, product yields or artificial upper seed.
//! It is not the inventory-limited chemistry owner: adding that owner's released
//! heat again here would double-count chemical energy. The optional weak-plume
//! limit uses the actual compositions and original NASA enthalpy differences
//! to bound entrainment at neutral buoyancy. Wall/corner and oxygen-limitation
//! corrections are not inferred. Rates are frozen at evaluation.

use super::{Layer, LayerRates, TwoZone, ZoneError};
use fs_convection::fire::{
    FireAmbientV1, FireCorrelationErrorV1, FireEstimateV1, HeskestadEntrainmentV1, HeskestadFireV1,
    HeskestadRegionV1, SpecificHeatCapacityV1,
};
use fs_qty::{Acceleration, Density, Length, MassFlowRate, Power, Temperature};
use fs_thermochem::MassSpecificEnergyQuantityV1;

/// Physical origin of the optional weak-plume constraint. The variable-caloric,
/// unequal-composition extension is derived here, not claimed as a CFAST fit.
pub const WEAK_PLUME_SOURCE: &str = "NIST TN 1889v1 (October 2015), section 3.3.3, equations (3.18)-(3.19); https://doi.org/10.6028/NIST.TN.1889v1";

/// Explicit choice of interlayer penetration model.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlumePenetration {
    /// The original unobstructed Heskestad estimate, without a receiver limit.
    Unconfined,
    /// Cap entrainment so the heated donor mixture is no denser than the upper
    /// layer at the shared pressure. This is a reduced steady plume closure;
    /// momentum overshoot, detrainment and a resolved stalled plume are absent.
    NeutralBuoyancy,
}

/// Finite thermal barrier to penetration of a present upper layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlumePenetrationLimit {
    /// Donor-composition temperature at the receiver density [K]. For equal
    /// compositions this is exactly the upper-layer temperature.
    pub neutral_temperature: Temperature,
    /// Original frozen-donor enthalpy rise to that temperature [J/kg].
    pub specific_enthalpy_rise: MassSpecificEnergyQuantityV1,
    /// Convective power divided by the required enthalpy rise [kg/s].
    pub maximum_mass_rate: MassFlowRate,
    /// True exactly when this bound reduced the nominal Heskestad transfer.
    pub active: bool,
}

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
    /// Positive radiant power, with exactly zero convective power and transfer.
    RadiativeOnly,
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
    /// The caller's explicit choice; the original `evaluate` uses Unconfined.
    pub penetration: PlumePenetration,
    /// Present only when a present upper layer requires positive donor heating.
    /// Absence also covers an absent upper layer or an already buoyant donor.
    pub penetration_limit: Option<PlumePenetrationLimit>,
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
    /// layer is seeded or queried for density. This entry point retains the
    /// unconfined correlation; use `evaluate_with_penetration` to opt into the
    /// actual-receiver neutral-buoyancy constraint.
    pub fn evaluate(
        &self,
        room: &TwoZone,
        gravity: Acceleration,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PrescribedFireReport, PrescribedFireError> {
        self.evaluate_with_penetration(room, gravity, PlumePenetration::Unconfined, cancelled)
    }

    /// Evaluate a nominal plume and an explicitly selected penetration limit.
    ///
    /// NIST TN 1889v1 section 3.3.3 bounds entrainment by Qc/[Cp(Tu-Tl)].
    /// At the original common ideal-gas pressure, unequal compositions instead
    /// require `T_neutral = Tu * R_upper/R_lower`. The heat needed to reach
    /// that density is the original donor-mixture enthalpy difference, with no
    /// constant-Cp approximation or formation-energy subtraction. Thus the
    /// admitted rate is `min(m_Heskestad, Qc/Delta_h)` when Delta_h is positive.
    ///
    /// The neutral state must lie inside the retained caloric source domain;
    /// extrapolation refuses. An absent upper layer imposes no barrier. A
    /// radiation-only source returns no gas heat or species transfer and does
    /// not evaluate an irrelevant flame-height fit. No fuel/products or plume
    /// momentum are inferred. As in CFAST's reduced closure, all supplied Qc
    /// still enters the upper layer; this is not a stalled-plume simulation.
    pub fn evaluate_with_penetration(
        &self,
        room: &TwoZone,
        gravity: Acceleration,
        penetration: PlumePenetration,
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
            penetration,
            penetration_limit: None,
        };
        if total > 0.0 && convective == 0.0 {
            report.regime = PrescribedFireRegime::RadiativeOnly;
        } else if convective > 0.0 {
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
                let mut actual_mass_rate = entrained.mass_rate;
                if penetration == PlumePenetration::NeutralBuoyancy {
                    if let Some(mut limit) =
                        neutral_buoyancy_limit(room, convective, &mut cancelled)?
                    {
                        limit.active = limit.maximum_mass_rate.value() < actual_mass_rate.value();
                        if limit.active {
                            actual_mass_rate = limit.maximum_mass_rate;
                        }
                        report.penetration_limit = Some(limit);
                    }
                }
                let carried = room.carried_rate(Layer::Lower, actual_mass_rate)?;
                report.rates[Layer::Lower.index()] = carried.opposite();
                report.rates[Layer::Upper.index()].add(&carried)?;
                report.plume_mass_rate = actual_mass_rate;
                report.regime = PrescribedFireRegime::Entrainment(entrained.region);
                report.plume_estimate = Some(estimate);
                report.plume_ambient = Some(ambient);
            }
        }
        poll(&mut cancelled)?;
        Ok(report)
    }
}

fn accumulate(sum: &mut f64, correction: &mut f64, term: f64) {
    let adjusted = term - *correction;
    let next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

fn neutral_buoyancy_limit(
    room: &TwoZone,
    convective_power: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<PlumePenetrationLimit>, PrescribedFireError> {
    let lower = room.layer(Layer::Lower);
    let upper = room.layer(Layer::Upper);
    let lower_mass = lower.mass_kg();
    let upper_mass = upper.mass_kg();
    if upper_mass == 0.0 {
        return Ok(None);
    }
    let mut inverse_masses = [0.0; 2];
    let mut corrections = [0.0; 2];
    for ((&lower_i, &upper_i), model) in lower
        .species_masses()
        .iter()
        .zip(upper.species_masses())
        .zip(room.gas().source().models())
    {
        poll(cancelled)?;
        for (i, fraction) in [lower_i / lower_mass, upper_i / upper_mass]
            .into_iter()
            .enumerate()
        {
            accumulate(
                &mut inverse_masses[i],
                &mut corrections[i],
                fraction / model.molar_mass().value(),
            );
        }
    }
    let neutral_temperature = upper.temperature().value() * (inverse_masses[1] / inverse_masses[0]);
    if !neutral_temperature.is_finite() || neutral_temperature <= 0.0 {
        return Err(PrescribedFireError::Invalid(
            "unrepresentable neutral plume temperature",
        ));
    }
    if neutral_temperature <= lower.temperature().value() {
        return Ok(None);
    }
    let domain = room.gas().domain();
    if !(domain.minimum_temperature..=domain.maximum_temperature).contains(&neutral_temperature) {
        return Err(PrescribedFireError::Gas(ZoneError::CaloricDomain));
    }
    let neutral_temperature = Temperature::new(neutral_temperature);
    let (mut rise, mut correction) = (0.0, 0.0);
    for (&species_mass, model) in lower
        .species_masses()
        .iter()
        .zip(room.gas().source().models())
    {
        poll(cancelled)?;
        if species_mass == 0.0 {
            continue;
        }
        let change = model
            .molar_enthalpy_change(lower.temperature(), neutral_temperature)
            .map_err(|e| PrescribedFireError::Gas(ZoneError::Gas(e.into())))?
            .value();
        accumulate(
            &mut rise,
            &mut correction,
            (species_mass / lower_mass) * change / model.molar_mass().value(),
        );
    }
    if !rise.is_finite() || rise <= 0.0 {
        return Err(PrescribedFireError::Invalid(
            "nonpositive or unrepresentable neutral plume enthalpy rise",
        ));
    }
    let maximum_mass_rate = convective_power / rise;
    if !maximum_mass_rate.is_finite() || maximum_mass_rate <= 0.0 {
        return Err(PrescribedFireError::Invalid(
            "unrepresentable weak-plume mass limit",
        ));
    }
    Ok(Some(PlumePenetrationLimit {
        neutral_temperature,
        specific_enthalpy_rise: MassSpecificEnergyQuantityV1::new(rise),
        maximum_mass_rate: MassFlowRate::new(maximum_mass_rate),
        active: false,
    }))
}

#[cfg(test)]
mod tests;
