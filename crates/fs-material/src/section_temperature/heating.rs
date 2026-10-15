//! Live section heating through the shared conservative enthalpy timeline.
//!
//! One fixed-geometry, uniform-temperature control volume represents each
//! annulus. The original section supplies volume and invariant reference mass;
//! source-backed caloric charts supply energy storage. The same temperature
//! evaluation drives Ohmic heat, thermal links and the returned elastic section.
//! No replacement heat integrator, caloric formula or electrical law is used.
//!
//! Layer lumping and resistance geometry are caller modelling assumptions.
//! Density must be constant on this fixed-volume rung. Thermal strain/work,
//! radial temperature profiles within a layer, skin effects and melting are
//! outside it. Source stiffness is not an explicit-stability bound: refine the
//! supplied maximum timestep for electrical/thermal feedback independently.

use std::collections::BTreeMap;
use fs_thermal::enthalpy::{Boundary, CellGeometry, CellId, DepositedCell,
    DepositionBudget, DepositionError, DepositionObserver, DepositionSimulation,
    EnergyAudit, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent};
use crate::conductor::{ConductorError, OhmicDrive};
use crate::phase::{PhaseStateError, ResolvedHeatCapacityCurve};
use crate::section::{ConcentricElasticSection, transport::ElectrothermalSectionResponse};
use crate::state_point::{DENSITY_PROPERTY, MaterialStatePointError,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY};
use super::ConcentricSectionTemperatureLaw;

/// Explicit resistance geometry for adjacent annuli `i` and `i + 1`.
/// The total resistance is `shape[0]/k_i + contact + shape[1]/k_(i+1)`.
/// Shapes have units 1/m; for a planar half path they are distance/area.
/// Curved paths must supply their actual resistance geometry, not that planar
/// approximation. Covered area controls exposure, independently of resistance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SectionThermalInterface {
    /// Area covered on each incident layer [m2].
    pub area_m2: f64,
    /// Positive half-path resistance geometry [1/m], in ordered layer roles.
    pub resistance_geometry_per_m: [f64; 2],
    /// Additional interface resistance [K/W], including explicit zero.
    pub contact_resistance_k_w: f64,
}

/// Material, topology, numerical or conservation refusal; no state is changed.
#[derive(Debug)]
pub enum SectionHeatingError {
    /// Existing state-point/material admission.
    Material(MaterialStatePointError),
    /// Existing conductor admission.
    Electrical(ConductorError),
    /// Existing source-backed caloric chart admission.
    Caloric(PhaseStateError),
    /// Existing thermal solver, energy or cancellation refusal.
    Thermal(DepositionError),
    /// Incompatible physical inputs or derived quantities.
    Invalid(&'static str),
}
impl core::fmt::Display for SectionHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "section heating refusal: {self:?}")
    }
}
impl std::error::Error for SectionHeatingError {}
impl From<MaterialStatePointError> for SectionHeatingError {
    fn from(e: MaterialStatePointError) -> Self { Self::Material(e) }
}
impl From<ConductorError> for SectionHeatingError {
    fn from(e: ConductorError) -> Self { Self::Electrical(e) }
}
impl From<PhaseStateError> for SectionHeatingError {
    fn from(e: PhaseStateError) -> Self { Self::Caloric(e) }
}
impl From<DepositionError> for SectionHeatingError {
    fn from(e: DepositionError) -> Self { Self::Thermal(e) }
}

#[derive(Clone)]
struct HeatingObserver<'a> {
    law: ConcentricSectionTemperatureLaw<'a>,
    section: ConcentricElasticSection,
    interfaces: Vec<SectionThermalInterface>,
    length_m: f64,
    drive: OhmicDrive,
    applied_power_w: f64,
    electrical_energy_j: f64,
}

/// Coupled thermal history and source-resolved mechanical section.
/// All state is private; only a successful `advance_to` returns a new candidate.
/// Source references include the exact Cp/density receipts, not just a chart ID.
#[derive(Clone)]
pub struct HeatedConcentricSection<'a> {
    thermal: DepositionSimulation<'a>,
    caloric: Vec<&'a ResolvedHeatCapacityCurve>,
    observer: HeatingObserver<'a>,
}

/// Per-advance energy account and the actual final electrical constitutive state.
#[derive(Debug)]
pub struct SectionHeatingReport {
    /// Shared thermal audit; generated electrical heat is already counted here.
    pub thermal: EnergyAudit,
    /// Integrated layer Joule heat, using each actual thermal substep [J].
    pub electrical_energy_j: f64,
    /// Endpoint response, not the time-integrated or start-of-step power.
    pub endpoint: ElectrothermalSectionResponse,
}

impl<'a> HeatedConcentricSection<'a> {
    /// Build one existing thermal cell per section layer. Supply source-backed
    /// caloric charts and boundary data in the original spatial order, plus
    /// exactly one interface between each neighboring pair. Cell/material IDs
    /// are their zero-based spatial positions; no sorting by material identity.
    ///
    /// Cp and density claims, fixed context, coordinate conventions and card
    /// identity must match the mechanical bundle. No missing property becomes
    /// insulation. A nonconstant density or nonsolid chart refuses: this rung
    /// cannot conserve mass by silently changing fixed geometry.
    pub fn try_new(
        law: ConcentricSectionTemperatureLaw<'a>,
        caloric: &[&'a ResolvedHeatCapacityCurve],
        length_m: f64,
        boundaries: &[Boundary],
        interfaces: Vec<SectionThermalInterface>,
        budget: DepositionBudget,
    ) -> Result<Self, SectionHeatingError> {
        positive(length_m)?;
        let reference = law.reference();
        let count = reference.layers().len();
        if caloric.len() != count || boundaries.len() != count || interfaces.len() != count - 1 {
            return Err(SectionHeatingError::Invalid("caloric/boundary/interface counts do not match section layers"));
        }
        for interface in &interfaces {
            positive(interface.area_m2)?;
            for &shape in &interface.resistance_geometry_per_m { positive(shape)?; }
            nonnegative(interface.contact_resistance_k_w)?;
        }
        let mut schedule = Vec::with_capacity(count + interfaces.len());
        let mut materials: Vec<(MaterialId, &'a dyn ThermalLaw)> = Vec::with_capacity(count);
        let mut inner = reference.inner_radius_m();
        for (index, ((layer, integrals), &caloric)) in reference.layers().iter()
            .zip(reference.layer_integrals()).zip(caloric).enumerate()
        {
            let binding = &law.bindings[index];
            if caloric.curve().material_card_identity() != layer.material.card_identity() {
                return Err(SectionHeatingError::Invalid("caloric and section layer cards differ"));
            }
            let rho = layer.material.property(DENSITY_PROPERTY)
                .ok_or(SectionHeatingError::Invalid("layer lacks density"))?.value_si();
            if caloric.curve().knots().iter()
                .any(|k| k.liquid_mass_fraction != 0.0 || k.bulk_density_kg_m3 != rho)
            {
                return Err(SectionHeatingError::Invalid("fixed section heating requires a solid constant-density caloric span"));
            }
            for pair in caloric.receipts().chunks_exact(2) {
                for (receipt, name) in pair.iter().zip([SPECIFIC_HEAT_CAPACITY_PROPERTY, DENSITY_PROPERTY]) {
                    let property = layer.material.property(name)
                        .ok_or(SectionHeatingError::Invalid("layer lacks a caloric property"))?;
                    if receipt.selected != property.answer().receipt.selected
                        || receipt.axis_quantities != *layer.material.axis_quantities()
                        || !receipt.query_point.iter().filter(|(key, _)| key != &binding.axis)
                            .eq(layer.material.query_point().iter().filter(|(key, _)| key != &binding.axis))
                    {
                        return Err(SectionHeatingError::Invalid("caloric and mechanical source claims or conditions differ"));
                    }
                }
            }
            let temperature = binding.point.axes()[&binding.axis];
            let state = caloric.state_at_temperature(temperature)?;
            let volume = positive(integrals.area_m2() * length_m)?;
            let mass = positive(integrals.line_density_kg_m() * length_m)?;
            let surface = positive(2.0 * core::f64::consts::PI * (inner + layer.outer_radius_m)
                * length_m + 2.0 * integrals.area_m2())?;
            materials.push((MaterialId(index as u64), caloric.curve()));
            schedule.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
                id: CellId(index as u64), segment_id: 0, body_id: 0,
                material_id: MaterialId(index as u64), mass_kg: mass,
                specific_h_j_kg: state.specific_enthalpy_j_kg(),
                geometry: CellGeometry { start_m: [0.0; 3], end_m: [0.0, 0.0, length_m],
                    volume_m3: volume, surface_area_m2: surface },
                boundary: boundaries[index],
            }) });
            inner = layer.outer_radius_m;
        }
        let observer = HeatingObserver { section: reference.clone(), law, interfaces, length_m,
            drive: OhmicDrive::Current(0.0), applied_power_w: 0.0, electrical_energy_j: 0.0 };
        // Admit electrical/thermal data before allocating the thermal timeline.
        observer.electrical()?;
        for (index, &conductance) in observer.link_conductances()?.iter().enumerate() {
            schedule.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
                id: index as u64, a: CellId(index as u64), b: CellId((index + 1) as u64),
                area_m2: observer.interfaces[index].area_m2, conductance_w_k: conductance,
            }) });
        }
        let (thermal, _) = DepositionSimulation::new(schedule, &materials, budget)?
            .advance_to(0.0, StepConfig { max_dt_s: 1.0, energy_tolerance_j: 0.0 }, || Ok(()))?;
        Ok(Self { thermal, caloric: caloric.to_vec(), observer })
    }

    /// Current mechanical coefficients and all material receipts at the SAME
    /// temperatures as the accepted thermal history. No stale room-temperature E.
    pub const fn section(&self) -> &ConcentricElasticSection { &self.observer.section }
    /// Thermal history, reference masses, topology and declared boundaries.
    pub const fn thermal(&self) -> &DepositionSimulation<'a> { &self.thermal }
    /// Original caloric sources and their complete source receipts.
    pub fn caloric(&self) -> &[&'a ResolvedHeatCapacityCurve] { &self.caloric }

    /// Advance with a constant ideal terminal drive over this interval.
    /// Currents, per-layer Joule heat AND thermal-link conductances are refreshed
    /// at each shared solver substep. The final material section is re-resolved
    /// before publication. Any late source/domain/energy/budget/cancellation
    /// refusal discards both histories; `self` is never mutated.
    ///
    /// The existing explicit conduction CFL and caller max_dt are retained.
    /// They do not prove stability for a nonlinear temperature-dependent drive.
    pub fn advance_to(
        &self, target_s: f64, drive: OhmicDrive, config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, SectionHeatingReport), SectionHeatingError> {
        checkpoint()?;
        let mut observer = self.observer.clone();
        observer.drive = drive;
        observer.electrical()?;
        let (thermal, observer, audit) = self.thermal.advance_observed(
            &observer, target_s, config, &mut checkpoint,
        )?;
        let report = SectionHeatingReport { thermal: audit,
            electrical_energy_j: observer.electrical_energy_j, endpoint: observer.electrical()? };
        checkpoint()?;
        Ok((Self { thermal, observer, caloric: self.caloric.clone() }, report))
    }
}

impl HeatingObserver<'_> {
    fn refresh(&mut self, thermal: &DepositionSimulation<'_>) -> Result<(), SectionHeatingError> {
        let temperatures = (0..self.law.bindings.len()).map(|i| thermal.temperature_k(CellId(i as u64)))
            .collect::<Result<Vec<_>, _>>()?;
        let section = self.law.at_temperatures_k(&temperatures)?;
        for (current, reference) in section.layer_integrals().iter().zip(self.law.reference().layer_integrals()) {
            if current.line_density_kg_m() != reference.line_density_kg_m() {
                return Err(SectionHeatingError::Invalid("temperature changed density on a fixed-volume constant-mass section"));
            }
        }
        self.section = section;
        Ok(())
    }
    fn electrical(&self) -> Result<ElectrothermalSectionResponse, SectionHeatingError> {
        Ok(self.section.thermal_transport()?.with_ohmic_transport(self.length_m)?.evaluate(self.drive, 0.0)?)
    }
    fn link_conductances(&self) -> Result<Vec<f64>, SectionHeatingError> {
        let conductivity = self.section.layers().iter().map(|layer| {
            let k = layer.material.property(THERMAL_CONDUCTIVITY_PROPERTY)
                .ok_or(SectionHeatingError::Invalid("layer lacks thermal conductivity"))?.value_si();
            nonnegative(k)
        }).collect::<Result<Vec<_>, _>>()?;
        self.interfaces.iter().enumerate().map(|(i, interface)| {
            let (a, b) = (conductivity[i], conductivity[i + 1]);
            if a == 0.0 || b == 0.0 { return Ok(0.0); }
            let r = positive(interface.resistance_geometry_per_m[0] / a
                + interface.resistance_geometry_per_m[1] / b + interface.contact_resistance_k_w)?;
            positive(1.0 / r)
        }).collect()
    }
}

impl DepositionObserver for HeatingObserver<'_> {
    type Error = SectionHeatingError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.electrical_energy_j = 0.0;
        self.applied_power_w = 0.0;
        Ok(())
    }
    fn sources(&mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error>
    {
        checkpoint()?;
        self.refresh(thermal)?;
        let response = self.electrical()?;
        self.applied_power_w = response.joule_power_w;
        checkpoint()?;
        Ok(response.layers.iter().enumerate().map(|(i, layer)| (CellId(i as u64), layer.joule_power_w)).collect())
    }
    fn conductances(&mut self, _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<u64, f64>, Self::Error>
    {
        checkpoint()?;
        Ok(self.link_conductances()?.into_iter().enumerate().map(|(i, g)| (i as u64, g)).collect())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Err(SectionHeatingError::Invalid("section geometry is fixed during thermal evolution"))
    }
    fn step(&mut self, _: f64, duration_s: f64, _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>, checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<(), Self::Error>
    {
        checkpoint()?;
        let joules = nonnegative(self.applied_power_w * duration_s)?;
        if self.applied_power_w != 0.0 && joules == 0.0 {
            return Err(SectionHeatingError::Invalid("nonzero section heating energy underflowed"));
        }
        self.electrical_energy_j = nonnegative(self.electrical_energy_j + joules)?;
        self.refresh(thermal)?;
        self.electrical()?;
        self.link_conductances()?;
        checkpoint()?;
        Ok(())
    }
}

fn positive(value: f64) -> Result<f64, SectionHeatingError> {
    if value.is_finite() && value > 0.0 { Ok(value) }
    else { Err(SectionHeatingError::Invalid("positive finite section thermal quantity required")) }
}
fn nonnegative(value: f64) -> Result<f64, SectionHeatingError> {
    if value.is_finite() && value >= 0.0 { Ok(value) }
    else { Err(SectionHeatingError::Invalid("nonnegative finite section thermal quantity required")) }
}

#[cfg(test)]
mod tests;
