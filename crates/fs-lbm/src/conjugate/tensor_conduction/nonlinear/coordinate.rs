//! Cell coordinates for the shared original-law tensor Newton solve.

use super::super::{ChtError, PseudoStep, StorageRate, finite, finite_positive};

/// Cell storage is expressed in its original coordinate; all face unknowns
/// remain Kelvin. Static dispatch preserves the original sensible fast path.
pub(super) trait CellCoordinate {
    const IS_TEMPERATURE: bool;

    fn reference_temperatures(&self) -> Option<&[f64]>;

    fn normalize(&self, _cell: usize, value: f64) -> Result<f64, ChtError> {
        finite("thermal.tensor_conduction.cell_coordinate", value)?;
        Ok(value)
    }

    fn temperature_and_slope(&self, cell: usize, value: f64) -> Result<[f64; 2], ChtError>;

    fn prescribed_value(&self, cell: usize, temperature_k: f64) -> Result<f64, ChtError>;

    fn storage_rate(&self, cell: usize, value: f64) -> Result<StorageRate, ChtError>;

    fn storage_tangent(&self, cell: usize) -> f64;

    /// Constitutive energy is absolute in the original chart, never the
    /// centered Newton increment. Kelvin-only coordinates provide no energy.
    #[cfg(feature = "enthalpy-coupling")]
    fn conductivity_point(
        &self,
        _cell: usize,
        _value: f64,
        _temperature_k: f64,
    ) -> Result<Option<super::material::TensorCaloricConductivityPoint>, ChtError> {
        Ok(None)
    }
}

pub(super) struct KelvinCoordinate<'a, 'b> {
    pub step: Option<&'a PseudoStep<'b>>,
}

impl CellCoordinate for KelvinCoordinate<'_, '_> {
    const IS_TEMPERATURE: bool = true;

    fn reference_temperatures(&self) -> Option<&[f64]> {
        self.step.map(|step| step.previous)
    }

    fn normalize(&self, _cell: usize, value: f64) -> Result<f64, ChtError> {
        finite_positive("thermal.tensor_conduction.initial_temperature_k", value)?;
        Ok(value)
    }

    fn temperature_and_slope(&self, _cell: usize, value: f64) -> Result<[f64; 2], ChtError> {
        finite_positive("thermal.tensor_conduction.temperature_k", value)?;
        Ok([value, 1.0])
    }

    fn prescribed_value(&self, _cell: usize, temperature_k: f64) -> Result<f64, ChtError> {
        finite_positive("thermal.fixed_temperature", temperature_k)?;
        Ok(temperature_k)
    }

    fn storage_rate(&self, cell: usize, value: f64) -> Result<StorageRate, ChtError> {
        Ok(self
            .step
            .map_or(StorageRate::default(), |step| StorageRate {
                power_w: step.coefficient_w_k[cell] * (value - step.previous[cell]),
                arithmetic_w: step.coefficient_w_k[cell]
                    * (value.abs() + step.previous[cell].abs()),
            }))
    }

    fn storage_tangent(&self, cell: usize) -> f64 {
        self.step.map_or(0.0, |step| step.coefficient_w_k[cell])
    }
}
