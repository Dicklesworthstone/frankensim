//! Explicitly isotropic use of one original conductivity/enthalpy source.
//!
//! The prepared material owner retains every source and caloric knot. This
//! view borrows that owner, preserves its energy datum and restricted support,
//! and performs no resampling, source selection, phase mixing or table copying.

use super::super::caloric::TensorCaloricLaw;
use super::{
    ConductivityTensor, TensorCaloricConductivity, TensorCaloricConductivityPoint,
    TensorConductivityLaw,
};
use crate::conjugate::{ChtError, poll};
use fs_conduction::enthalpy::{DepositedCell, DepositionError, ThermalLaw, ThermalStorageBasis};
use fs_exec::CancelGate;
use fs_material::thermal::conduction::phase::{
    PhaseConductivityError, PhaseConductivitySample, PreparedPhaseConductivity,
};

/// One original phase-conductivity source explicitly declared isotropic.
///
/// The scalar material source alone does not establish physical isotropy;
/// constructing this view declares that modeling choice. The native tensor
/// owner uses its cellwise K(h) and analytic dK/dh, including a latent plateau.
/// It does not reinterpret that path as a single-valued K(T).
///
/// Supply this same view to both conductivity and caloric storage fields to
/// retain the source's restricted energy domain and selected one-sided
/// temperature derivative at its boundaries. The underlying caloric identity,
/// reference enthalpy and original entropy relation remain unchanged.
#[derive(Clone, Copy, Debug)]
pub struct IsotropicPhaseConductivity<'a, 'm> {
    source: &'a PreparedPhaseConductivity<'m>,
}

impl<'a, 'm> IsotropicPhaseConductivity<'a, 'm> {
    /// Borrow a prepared original source with positive conductivity throughout.
    ///
    /// The source owner has already admitted continuous piecewise support.
    /// Checking its retained knots admits positivity on every original segment.
    /// Explicit zero-conductivity paths cannot enter the native SPD face owner.
    ///
    /// # Errors
    /// Cancellation or a nonpositive original source conductivity is refused.
    pub fn new(
        source: &'a PreparedPhaseConductivity<'m>,
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        poll(gate)?;
        for (index, &[_, conductivity]) in source.knots().iter().enumerate() {
            if index.is_multiple_of(256) {
                poll(gate)?;
            }
            if !conductivity.is_finite() || conductivity <= 0.0 {
                return Err(refusal(
                    "the native tensor owner requires strictly positive original phase conductivity",
                ));
            }
        }
        poll(gate)?;
        Ok(Self { source })
    }

    /// The original immutable owner, including its selected claim and paths.
    pub const fn source(&self) -> &'a PreparedPhaseConductivity<'m> {
        self.source
    }

    fn sample(&self, energy: f64) -> Result<PhaseConductivitySample, DepositionError> {
        self.source
            .state_at_specific_enthalpy(energy)
            .map_err(thermal_error)
    }
}

impl TensorConductivityLaw for IsotropicPhaseConductivity<'_, '_> {
    fn conductivity_at_temperature(
        &self,
        _temperature_k: f64,
    ) -> Result<ConductivityTensor, ChtError> {
        Err(refusal(
            "phase conductivity requires the original specific enthalpy",
        ))
    }

    fn mean_conductivity(&self, _from_k: f64, _to_k: f64) -> Result<ConductivityTensor, ChtError> {
        Err(refusal(
            "a Kelvin interval cannot select original latent phase inventory",
        ))
    }

    fn caloric_conductivity(
        &self,
        point: TensorCaloricConductivityPoint,
    ) -> Result<Option<TensorCaloricConductivity>, ChtError> {
        if point.storage_basis != ThermalStorageBasis::Enthalpy
            || point.caloric_identity != self.source.caloric_curve().caloric_identity()
        {
            return Err(refusal(
                "phase conductivity and cell storage must retain the same original enthalpy chart",
            ));
        }
        let sample = self
            .source
            .state_at_specific_enthalpy(point.specific_energy_j_kg)
            .map_err(source_error)?;
        let original_temperature = sample.temperature_k();
        if !point.temperature_k.is_finite()
            || point.temperature_k <= 0.0
            || (point.temperature_k - original_temperature).abs()
                > 64.0 * f64::EPSILON * point.temperature_k.abs().max(original_temperature.abs())
        {
            return Err(refusal(
                "cell temperature disagrees with the original phase-conductivity chart",
            ));
        }
        let conductivity = sample.conductivity_w_m_k();
        let derivative = sample.conductivity_enthalpy_derivative_kg_w_m_k_j();
        Ok(Some(TensorCaloricConductivity::new(
            ConductivityTensor::new([
                [conductivity, 0.0, 0.0],
                [0.0, conductivity, 0.0],
                [0.0, 0.0, conductivity],
            ])?,
            [
                [derivative, 0.0, 0.0],
                [0.0, derivative, 0.0],
                [0.0, 0.0, derivative],
            ],
        )?))
    }
}

impl ThermalLaw for IsotropicPhaseConductivity<'_, '_> {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.source.caloric_curve().caloric_identity()
    }

    fn validate_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        self.sample(cell.specific_h_j_kg)?;
        self.source.caloric_curve().validate_cell(cell)
    }

    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        Ok(self.sample(energy)?.temperature_k())
    }

    fn bulk_density_kg_m3(&self, energy: f64) -> Result<Option<f64>, DepositionError> {
        self.sample(energy)?;
        self.source.caloric_curve().bulk_density_kg_m3(energy)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        // The unchanged full chart's bound is conservative on this restriction.
        self.source.caloric_curve().maximum_temperature_slope()
    }

    fn specific_entropy_change_j_kg_k(
        &self,
        opening: f64,
        endpoint: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.sample(opening)?;
        self.sample(endpoint)?;
        ThermalLaw::specific_entropy_change_j_kg_k(self.source.caloric_curve(), opening, endpoint)
    }
}

impl TensorCaloricLaw for IsotropicPhaseConductivity<'_, '_> {
    fn specific_energy_bounds_j_kg(&self) -> [f64; 2] {
        self.source.enthalpy_bounds_j_kg()
    }

    fn temperature_slope_kg_k_j(&self, energy: f64) -> Result<f64, DepositionError> {
        Ok(self
            .sample(energy)?
            .temperature_enthalpy_derivative_kg_k_j())
    }

    fn temperature_and_slope_at_specific_energy(
        &self,
        energy: f64,
    ) -> Result<[f64; 2], DepositionError> {
        let sample = self.sample(energy)?;
        Ok([
            sample.temperature_k(),
            sample.temperature_enthalpy_derivative_kg_k_j(),
        ])
    }
}

fn thermal_error(error: PhaseConductivityError) -> DepositionError {
    match error {
        PhaseConductivityError::OutsideEnthalpyDomain { .. } => DepositionError::MaterialDomain,
        PhaseConductivityError::Thermal(error) => error,
        _ => DepositionError::InvalidInput("original phase-conductivity source evaluation failed"),
    }
}

fn source_error(error: PhaseConductivityError) -> ChtError {
    match error {
        PhaseConductivityError::Thermal(DepositionError::Cancelled) => ChtError::Cancelled,
        PhaseConductivityError::OutsideEnthalpyDomain { .. } => ChtError::InvalidInput {
            field: "thermal.tensor_conduction.phase_source_domain",
            reason: error.to_string(),
        },
        _ => refusal(&error.to_string()),
    }
}

fn refusal(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "thermal.tensor_caloric.original_source_refusal",
        reason: reason.into(),
    }
}
