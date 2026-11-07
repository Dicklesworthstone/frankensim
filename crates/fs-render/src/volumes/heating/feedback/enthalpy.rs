//! Direct-beam absorption on the original native specific-enthalpy coordinate.
//!
//! Temperature is not a complete coordinate on a sharp latent plateau. The
//! existing source-bound enthalpy tables retain distinct optical states there;
//! their original caloric identity and storage-basis checks remain authoritative.
//! This module selects that coordinate in the same cached beam, pulse and
//! refinement owner. It introduces no phase fraction or mixing law.

use super::*;
use crate::volumes::heating::temperature_field::enthalpy::CompiledEnthalpyAbsorption;
pub use crate::volumes::heating::temperature_field::enthalpy::EnthalpyAbsorptionTable;

/// A path cell bound to its exact native material, caloric datum and receiver.
/// The selected table must name that material's registered caloric identity.
#[derive(Clone, Copy, Debug)]
pub struct EnthalpyOpticalCell {
    pub grid_cell: usize,
    pub thermal_cell: CellId,
    pub material_id: MaterialId,
    pub absorption_table: usize,
}

impl<const N: usize> TemperatureOpticalHeating<N> {
    /// Prepare direct-beam absorption versus the original native h [J/kg].
    ///
    /// Every non-vacuum path cell requires an explicit original material and
    /// caloric identity. Both the h interval and the independent kelvin support
    /// are checked at preparation, each source state, events and final endpoints.
    /// Zero/negative h is meaningful only relative to the retained source datum.
    ///
    /// The original spectral compiler, DDA cache, positive attenuation kernel,
    /// schedules and complete-trajectory refinement are shared with temperature
    /// feedback. The existing budget's knot ceiling counts h knots here.
    /// Geometry/density remain fixed; no phase mixing or hysteresis is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_enthalpy(
        cx: &Cx<'_>,
        grid: &AbsorbingGrid<'_>,
        beam: IncidentBeam<N>,
        tables: &[EnthalpyAbsorptionTable],
        bindings: &[EnthalpyOpticalCell],
        thermal: &DepositionSimulation<'_>,
        budget: OpticalFeedbackBudget,
        power_tolerance_w: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        Self::compile_beam(
            cx,
            grid,
            beam,
            tables.iter().map(SourceTable::Enthalpy),
            bindings.iter().map(|binding| PreparedBinding {
                cell: ThermalOpticalCell {
                    grid_cell: binding.grid_cell,
                    thermal_cell: binding.thermal_cell,
                    absorption_table: binding.absorption_table,
                },
                material_id: Some(binding.material_id),
            }),
            thermal,
            budget,
            power_tolerance_w,
            "org.frankensim.render.enthalpy-beam.v1",
            &mut checkpoint,
        )
    }
}

impl TemperatureOpticalHeating<1> {
    /// Compose authored linear kappa(h) with the existing monochromatic beta(T).
    ///
    /// Beta remains the original bulk m/W law on its declared temperature and
    /// irradiance support. No enthalpy or phase dependence of beta is inferred.
    /// Combined h/T knot and table counts obey the existing resource limits.
    /// The pulse owner integrates each affine power waveform with both
    /// coefficients frozen at the same native start state.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_enthalpy_two_photon(
        cx: &Cx<'_>,
        grid: &AbsorbingGrid<'_>,
        beam: IncidentBeam<1>,
        linear_tables: &[EnthalpyAbsorptionTable],
        bindings: &[EnthalpyOpticalCell],
        nonlinear: nonlinear::TwoPhotonFeedback<'_>,
        thermal: &DepositionSimulation<'_>,
        budget: OpticalFeedbackBudget,
        power_tolerance_w: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        checkpoint()?;
        nonlinear::admit(
            linear_tables.len(),
            linear_tables
                .iter()
                .map(EnthalpyAbsorptionTable::knot_count),
            nonlinear,
            budget,
        )?;
        let result = Self::try_new_enthalpy(
            cx,
            grid,
            beam,
            linear_tables,
            bindings,
            thermal,
            budget,
            power_tolerance_w,
            &mut checkpoint,
        )?;
        Self::attach_two_photon(result, nonlinear, thermal, &mut checkpoint)
    }
}

pub(super) enum SourceTable<'a> {
    Temperature(&'a TemperatureAbsorptionTable),
    Enthalpy(&'a EnthalpyAbsorptionTable),
}

impl SourceTable<'_> {
    pub(super) fn knot_count(&self) -> usize {
        match self {
            Self::Temperature(table) => table.knot_count(),
            Self::Enthalpy(table) => table.knot_count(),
        }
    }

    pub(super) fn wavelength_bounds_nm(&self) -> (f64, f64) {
        match self {
            Self::Temperature(table) => table.wavelength_bounds_nm(),
            Self::Enthalpy(table) => table.wavelength_bounds_nm(),
        }
    }

    pub(super) fn compile<const N: usize>(
        &self,
        wavelengths_nm: &[f64; N],
        hash: &mut DomainHasher,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<CompiledAbsorption<N>, OpticalFeedbackError> {
        match self {
            Self::Temperature(table) => table
                .compile(wavelengths_nm, hash, checkpoint)
                .map(CompiledAbsorption::Temperature),
            Self::Enthalpy(table) => table
                .compile_beam(wavelengths_nm, hash, checkpoint)
                .map(CompiledAbsorption::Enthalpy),
        }
    }
}

#[derive(Clone, Copy)]
pub(super) struct PreparedBinding {
    pub(super) cell: ThermalOpticalCell,
    pub(super) material_id: Option<MaterialId>,
}

#[derive(Debug)]
pub(super) enum CompiledAbsorption<const N: usize> {
    Temperature(CompiledTable<N>),
    Enthalpy(CompiledEnthalpyAbsorption<N>),
}

impl<const N: usize> CompiledAbsorption<N> {
    /// Return the actual native h as well as its coefficient so the optical
    /// snapshot identity distinguishes states on a constant-temperature plateau.
    pub(super) fn at(
        &self,
        thermal: &DepositionSimulation<'_>,
        segment: &BoundSegment,
        temperature_k: f64,
    ) -> Result<([f64; N], Option<f64>), OpticalFeedbackError> {
        match self {
            Self::Temperature(table) => Ok((table.at(segment.thermal_cell, temperature_k)?, None)),
            Self::Enthalpy(table) => {
                let (h, values) = table.at(
                    thermal,
                    segment.thermal_cell,
                    segment
                        .material_id
                        .ok_or_else(|| invalid("native material binding required"))?,
                    temperature_k,
                )?;
                Ok((values, Some(h)))
            }
        }
    }

    pub(super) fn check(
        &self,
        thermal: &DepositionSimulation<'_>,
        segment: &BoundSegment,
        temperature_k: f64,
    ) -> Result<(), OpticalFeedbackError> {
        match self {
            Self::Temperature(table) => {
                table.check_temperature(segment.thermal_cell, temperature_k)
            }
            Self::Enthalpy(table) => table
                .native_value(
                    thermal,
                    segment.thermal_cell,
                    segment
                        .material_id
                        .ok_or_else(|| invalid("native material binding required"))?,
                    temperature_k,
                )
                .map(|_| ()),
        }
    }
}

#[cfg(test)]
mod tests;
