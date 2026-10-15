//! Uniform source-bound DC heating at the body's ACTUAL thermal/phase midpoint.
//! R, ideal current/voltage control and I^2 R remain with fs-material::conductor.
//! This adapter supplies no circuit, current-distribution or phase mixture law.

use super::PowerSample;
use super::super::layers::MaterialConductionLayers;
use fs_matdb::PropertyKey;
use super::super::source::{ThermalMaterialError, ThermalMaterialPoint};
use fs_material::conductor::{ConductorError, OhmicDrive, ResolvedConductor,
    ELECTRICAL_RESISTIVITY_PROPERTY, resolve_uniform_conductor};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::MaterialPropertySelection;

/// One whole-body uniform conductor. Geometry stays fixed over the thermal
/// step; phase-dependent properties may change but expansion is not inferred.
#[derive(Clone, Debug)]
pub struct MaterialJouleHeating<'a> {
    /// Same source/coordinates as the finite body's enthalpy and heat transport.
    pub material: ThermalMaterialPoint<'a>,
    /// Exact existing resistivity claim-selection/pinning policy.
    pub selection: MaterialPropertySelection,
    /// Physical conducting length [m], not a reference for implicit expansion.
    pub length_m: f64,
    /// Physical current-carrying area [m2], distinct from heat-transfer area.
    pub area_m2: f64,
    /// Ideal signed current or voltage, explicitly held for this physical step.
    pub drive: OhmicDrive,
}

/// Original source-bound electrical state and actually delivered electrical work rate.
#[derive(Clone, Debug, PartialEq)]
pub struct JouleHeatingSample {
    /// Original resistivity query, geometry and usage receipt.
    pub conductor: ResolvedConductor,
    /// Signed current from the existing ideal-drive owner [A].
    pub current_a: f64,
    /// Signed voltage across that same conductor [V].
    pub voltage_v: f64,
    /// Actual body midpoint, including phase and source-chart identity.
    pub body: EquilibriumPhaseState,
    /// V*I minus the conductor owner's heat rate [W]; roundoff check only.
    pub electrical_power_residual_w: f64,
}

/// Neither a missing electrical datum nor a failed source law becomes insulation.
#[derive(Clone, Debug, PartialEq)]
pub enum JouleHeatingError {
    /// Same thermal-material coordinate and static-law refusal as conduction.
    Thermal(ThermalMaterialError),
    /// Original conductor/drive/property failure.
    Conductor(ConductorError),
    /// Geometry or independently checked power is unrepresentable/inconsistent.
    Invalid(&'static str),
}
impl From<ThermalMaterialError> for JouleHeatingError {
    fn from(e: ThermalMaterialError) -> Self { Self::Thermal(e) }
}
impl From<ConductorError> for JouleHeatingError {
    fn from(e: ConductorError) -> Self { Self::Conductor(e) }
}
impl core::fmt::Display for JouleHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Thermal(e) => write!(f, "material Joule heating: {e}"),
            Self::Conductor(e) => write!(f, "material Joule heating: {e}"),
            Self::Invalid(s) => write!(f, "material Joule heating: {s}"),
        }
    }
}
impl std::error::Error for JouleHeatingError {}

/// A whole-body material conductor with NO competing prescribed current or
/// voltage. Ideal drives and dynamic circuits lower this SAME source binding.
#[derive(Clone, Debug)]
pub struct MaterialOhmicBody<'a> {
    /// Original enthalpy-card identity and thermal/phase coordinates.
    pub material: ThermalMaterialPoint<'a>,
    /// Existing explicit resistivity claim selection.
    pub selection: MaterialPropertySelection,
    /// Fixed physical current path length [m].
    pub length_m: f64,
    /// Fixed physical current-carrying area [m2].
    pub area_m2: f64,
}
impl MaterialOhmicBody<'_> {
    /// Resolve the original conductor at the ACTUAL enthalpy midpoint. This is
    /// the same geometry, source support and static-law admission as ideal Joule
    /// heating, without inventing a zero-current drive for a dynamic circuit.
    pub fn resolve_conductor(&self, body: EquilibriumPhaseState, body_volume_m3: f64)
        -> Result<ResolvedConductor, JouleHeatingError> {
        let volume = self.length_m*self.area_m2;
        if !volume.is_finite() || volume <= 0.0 || !body_volume_m3.is_finite() || body_volume_m3 <= 0.0
            || (volume/body_volume_m3-1.0).abs() > 64.0*f64::EPSILON {
            return Err(JouleHeatingError::Invalid("conductor geometry must represent the whole thermal body"));
        }
        let point = self.material.at(body)?;
        let conductor = resolve_uniform_conductor(self.material.material, &point, self.selection.clone(),
            self.length_m, self.area_m2)?;
        let property = conductor.material().property(ELECTRICAL_RESISTIVITY_PROPERTY)
            .ok_or(JouleHeatingError::Invalid("resolved electrical resistivity is absent"))?;
        self.material.admit_static(property)?;
        Ok(conductor)
    }
}

impl<'source> MaterialJouleHeating<'source> {
    /// The exact same conductor binding without an ideal-drive assumption.
    pub fn ohmic_body(&self) -> MaterialOhmicBody<'source> {
        MaterialOhmicBody { material:self.material, selection:self.selection.clone(),
            length_m:self.length_m, area_m2:self.area_m2 }
    }

    /// Bind heat transport to this SAME card, condition point and phase axes.
    /// Thermal and electrical path geometries remain distinct physical inputs.
    /// Each property retains its own explicit source-selection policy.
    pub fn conduction_layers<'a>(
        &'a self, conductivity: &'a PropertyKey, selection: MaterialPropertySelection,
        area_m2: f64, gas_distance_m: f64, body_distance_m: f64,
    ) -> MaterialConductionLayers<'a> {
        MaterialConductionLayers { material: self.material.material, context: self.material.context,
            temperature_axis: self.material.temperature_axis, liquid_fraction_axis: self.material.liquid_fraction_axis,
            conductivity, selection, area_m2, gas_distance_m, body_distance_m }
    }

    /// Evaluate the existing Ohmic law at T/phase from the same body enthalpy
    /// chart. The conductor must cover this whole body's volume (`length*area`).
    /// Geometric agreement is a floating-point consistency check, not evidence
    /// of uniform current or negligible internal temperature gradients.
    ///
    /// Scalar and explicitly declared temperature/fraction source curves are
    /// admitted. Zero drive still checks geometry, source support and identity;
    /// missing/zero resistivity does not silently become zero power. Pressure
    /// and all other material conditions stay explicit on the retained context.
    /// Use as the power callback of `GasPiston::powered_thermal_candidate` so
    /// the electrical heat and both thermal/mechanical owners accept together.
    pub fn resolve(&self, body: EquilibriumPhaseState, body_volume_m3: f64)
        -> Result<PowerSample<JouleHeatingSample>, JouleHeatingError> {
        let conductor = self.ohmic_body().resolve_conductor(body, body_volume_m3)?;
        let current = conductor.current_a_for_drive(self.drive)?;
        let power = conductor.joule_power_w(current)?;
        let voltage = match self.drive {
            OhmicDrive::Voltage(v) => v,
            OhmicDrive::Current(_) => current*conductor.resistance_ohm(),
        };
        let delivered = voltage*current;
        let defect = delivered-power;
        if !voltage.is_finite() || !delivered.is_finite() || !defect.is_finite()
            || (current != 0.0 && (voltage == 0.0 || delivered == 0.0))
            || defect.abs() > 64.0*f64::EPSILON*power {
            return Err(JouleHeatingError::Invalid("electrical work rate and thermal power do not agree"));
        }
        Ok(PowerSample { power_w: power, source: JouleHeatingSample { conductor, current_a: current,
            voltage_v: voltage, body, electrical_power_residual_w: defect } })
    }
}

#[cfg(test)]
mod tests;
