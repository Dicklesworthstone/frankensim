//! Pressure-driven series flow from the existing material resistance law.
//!
//! At one frozen material state per section, the SAME conserved mass rate m
//! obeys dp = A*m + B*m*|m|, where A=sum(a*mu/rho), B=sum(b/rho).
//! Solve this scalar reduction without a flow search, then re-evaluate EVERY
//! original section for Reynolds support, signed pressure, tangent and heat.
//! No temperature equation, pump curve, interpolation or fluid default is added.

use super::{FluidTransport, MaterialFlowResistance, ResistanceError,
    ResistanceEvaluation};
use fs_exec::CancelGate;
use fs_math::det;
use fs_qty::MassFlowRate;

/// The same solved flow and material snapshot in the existing heat-transfer owner.
pub mod convection;

/// One physical series section and its ORIGINAL material owner. Geometry does
/// not absorb viscosity or density; distinct local states retain their receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct PressureSection<F> {
    pub resistance: MaterialFlowResistance,
    pub medium: F,
}

/// Cold-work bounds and an independently checked pressure closure. Neither
/// tolerance measures material uncertainty or the constant-state approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SeriesPressureBudget {
    pub maximum_sections: usize,
    pub absolute_pressure_tolerance_pa: f64,
    pub relative_pressure_tolerance: f64,
}

#[derive(Debug)]
pub enum SeriesPressureError {
    Invalid(&'static str),
    Section { index: usize, source: ResistanceError },
    Closure { residual_pa: f64, tolerance_pa: f64 },
    Cancelled,
}
impl core::fmt::Display for SeriesPressureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material series pressure: {self:?}")
    }
}
impl std::error::Error for SeriesPressureError {}
type Result<T> = core::result::Result<T, SeriesPressureError>;
fn invalid(s: &'static str) -> SeriesPressureError { SeriesPressureError::Invalid(s) }
fn poll(gate: &CancelGate) -> Result<()> {
    if gate.is_requested() { Err(SeriesPressureError::Cancelled) } else { Ok(()) }
}
fn finite(v: f64) -> Result<f64> {
    if v.is_finite() { Ok(v) } else { Err(invalid("unrepresentable series pressure arithmetic")) }
}

/// The original material objects, actual section responses and one conserved
/// signed mass flow. Private fields keep a stale pressure/heat tuple from being
/// mistaken for an evaluation. Stored objects are not re-queried by accessors.
#[derive(Debug)]
pub struct SeriesPressureFlow<F> {
    sections: Vec<PressureSection<F>>,
    evaluations: Vec<ResistanceEvaluation>,
    driving_pressure_pa: f64,
    mass_flow_kg_s: f64,
    pressure_residual_pa: f64,
    dissipated_power_w: f64,
    dpressure_dmass_flow: f64,
}

impl<F: FluidTransport> SeriesPressureFlow<F> {
    /// Solve fixed-state series hydraulics, including reverse and zero drive.
    ///
    /// Samples each original material EXACTLY ONCE. Coefficients are assembled
    /// in mass-flow coordinates, so unequal local densities do not incorrectly
    /// enforce identical local volume flows. Every original Reynolds band is
    /// checked at the resulting mass rate; even zero flow may be unsupported.
    /// Source failure never falls back to another medium or loss mechanism.
    ///
    /// The rationalized inverse scales by single-term flow estimates. It does
    /// not square dimensional pressure/loss coefficients or subtract similar
    /// roots. No iteration bracket, Newton tolerance or solver workspace is
    /// required. Work/storage are O(number of sections), under an explicit cap.
    ///
    /// This is a prescribed differential pressure, NOT a source pressure update
    /// or pump model. Local states/geometry stay frozen. Sum(dp_i*Q_i) is actual
    /// hydraulic dissipation and must be deposited or explicitly accounted for
    /// by a thermal consumer. It is not generally dp_total times a reference Q
    /// when densities differ. No pressure work, compressible acceleration,
    /// cavitation, buoyancy or transient momentum model is inferred.
    pub fn solve(
        sections: Vec<PressureSection<F>>,
        driving_pressure_pa: f64,
        budget: SeriesPressureBudget,
        gate: &CancelGate,
    ) -> Result<Self> {
        poll(gate)?;
        if sections.is_empty() || sections.len() > budget.maximum_sections
            || !driving_pressure_pa.is_finite()
            || !budget.absolute_pressure_tolerance_pa.is_finite()
            || budget.absolute_pressure_tolerance_pa < 0.0
            || !budget.relative_pressure_tolerance.is_finite()
            || !(0.0..1.0).contains(&budget.relative_pressure_tolerance)
        { return Err(invalid("nonempty bounded sections, finite drive and pressure tolerances required")); }
        let mut snapshots = Vec::with_capacity(sections.len());
        let (mut linear, mut quadratic) = (Sum::default(), Sum::default());
        for (index, section) in sections.iter().enumerate() {
            poll(gate)?;
            let state = section.medium.transport_state().map_err(|source|
                SeriesPressureError::Section { index, source: ResistanceError::Material(source) })?;
            let a = finite((section.resistance.viscous_geometry * state.dynamic_viscosity_pa_s())
                / state.density_kg_m3())?;
            let b = finite(section.resistance.inertial_geometry / state.density_kg_m3())?;
            if (section.resistance.viscous_geometry > 0.0 && a == 0.0)
                || (section.resistance.inertial_geometry > 0.0 && b == 0.0) {
                return Err(invalid("active series pressure mechanism underflows"));
            }
            linear.add(a)?; quadratic.add(b)?; snapshots.push(state);
        }
        let mass_flow_kg_s = invert(driving_pressure_pa, linear.total()?, quadratic.total()?)?;
        let mut evaluations = Vec::with_capacity(sections.len());
        let (mut pressure, mut heat, mut tangent) = (Sum::default(), Sum::default(), Sum::default());
        for (index, (section, snapshot)) in sections.iter().zip(snapshots).enumerate() {
            poll(gate)?;
            // Evaluate the ORIGINAL law, not another copy of its pressure,
            // Reynolds or dissipation formula. No second material query.
            let evaluation = section.resistance.evaluate(&snapshot, MassFlowRate::new(mass_flow_kg_s))
                .map_err(|source| SeriesPressureError::Section { index, source })?;
            pressure.add(evaluation.pressure_drop_pa)?;
            heat.add(evaluation.dissipated_power_w)?;
            tangent.add(evaluation.dpressure_dmass_flow)?;
            evaluations.push(evaluation);
        }
        let actual = pressure.total()?;
        let pressure_residual_pa = finite(actual - driving_pressure_pa)?;
        let tolerance = budget.absolute_pressure_tolerance_pa.max(finite(
            budget.relative_pressure_tolerance * actual.abs().max(driving_pressure_pa.abs()))?);
        if pressure_residual_pa.abs() > tolerance {
            return Err(SeriesPressureError::Closure { residual_pa: pressure_residual_pa, tolerance_pa: tolerance });
        }
        let result = Self { sections, evaluations, driving_pressure_pa, mass_flow_kg_s,
            pressure_residual_pa, dissipated_power_w: heat.total()?, dpressure_dmass_flow: tangent.total()? };
        poll(gate)?;
        Ok(result)
    }
}
impl<F> SeriesPressureFlow<F> {
    pub fn sections(&self) -> &[PressureSection<F>] { &self.sections }
    pub fn evaluations(&self) -> &[ResistanceEvaluation] { &self.evaluations }
    pub const fn driving_pressure_pa(&self) -> f64 { self.driving_pressure_pa }
    pub const fn mass_flow_kg_s(&self) -> f64 { self.mass_flow_kg_s }
    pub const fn pressure_residual_pa(&self) -> f64 { self.pressure_residual_pa }
    pub const fn dissipated_power_w(&self) -> f64 { self.dissipated_power_w }
    /// Partial derivative at the FROZEN materials, not a thermal-feedback tangent.
    pub const fn dpressure_dmass_flow(&self) -> f64 { self.dpressure_dmass_flow }
}

fn invert(pressure: f64, a: f64, b: f64) -> Result<f64> {
    if pressure == 0.0 { return Ok(0.0); }
    if a < 0.0 || b < 0.0 || (a == 0.0 && b == 0.0) {
        return Err(invalid("positive series resistance required"));
    }
    let p = pressure.abs();
    // Same scaled algebra as the native passive graph's mixed branch inverse.
    let magnitude = if a == 0.0 { det::sqrt(p) / det::sqrt(b) }
        else if b == 0.0 { p / a }
        else {
            let linear = p / a;
            let quadratic = det::sqrt(p) / det::sqrt(b);
            if linear <= quadratic {
                let ratio = linear / quadratic;
                linear / (0.5 + 0.5 * det::sqrt(1.0 + 4.0 * ratio * ratio))
            } else {
                let ratio = quadratic / linear;
                quadratic / (0.5 * ratio + det::sqrt(1.0 + 0.25 * ratio * ratio))
            }
        };
    if !magnitude.is_finite() || magnitude <= 0.0 {
        return Err(invalid("series mass flow is unrepresentable"));
    }
    Ok(magnitude.copysign(pressure))
}

#[derive(Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, x: f64) -> Result<()> {
        let next = finite(self.value + x)?;
        let lost = if self.value.abs() >= x.abs() { (self.value-next)+x }
            else { (x-next)+self.value };
        self.correction = finite(self.correction + lost)?; self.value = next; Ok(())
    }
    fn total(self) -> Result<f64> { finite(self.value + self.correction) }
}

#[cfg(test)]
mod tests;
