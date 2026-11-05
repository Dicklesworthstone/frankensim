//! Thermal refraction from one original single-component gas and temperature field.
//!
//! The original gas evaluator supplies rho(T,p), with its original validity
//! limits, at each optical stage. At fixed pressure/composition its exact
//! derivatives are rho_T=-rho/T and rho_TT=2rho/T². These chain into the actual
//! thermal gradient/Hessian, then into the existing molecular Lorentz--Lorenz
//! and native ray/bundle owners. No renderer gas constant, density fit, finite
//! difference or per-stage allocation is introduced.
//!
//! This is a frozen isobaric single-component ideal gas, including an explicitly
//! declared pseudocomponent. The original gas model does not certify gas-phase
//! existence, measured accuracy, or its empirical transport extrapolations.
//! Composition/pressure gradients, flow advection, condensation and contact-face
//! traversal are outside this adapter. Original thermal/provider refusals remain
//! authoritative. Molecular scattering is separate from true absorption.

use super::{
    GasDensityGradientField, GasDensityGradientSample, GasDensityJet, GasDensitySecondOrderField,
    GasRefractiveFieldError, MaterialGasRefractiveField, MaterialGasScatteringBinding,
};
use crate::gradient_index::GradientIndexError;
use crate::gradient_index::thermal::TemperatureField;
use crate::material_state::thermo_optic::{
    TemperatureGradientField, TemperatureGradientSample, TemperatureSecondOrderField,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_geom::Point3;
use fs_material::gas::GasState;
use fs_material::gas_scattering::{GasScatteringError, MolecularRayleighMixture};
use fs_qty::{Pressure, Temperature};

/// Borrowed thermal/source data with optical coefficients prepared once.
/// The explicit reference temperature prepares molecular coefficients only;
/// every spatial sample uses the original field's actual temperature.
pub struct IsobaricGasField<'a, F: TemperatureGradientField + ?Sized> {
    temperature: &'a F,
    model: &'a MolecularRayleighMixture,
    reference: MaterialGasScatteringBinding,
    pressure: Pressure,
    snapshot_identity: ContentHash,
    identity: ContentHash,
}

impl<'a, F: TemperatureGradientField + ?Sized> IsobaricGasField<'a, F> {
    pub fn try_new(
        temperature: &'a F,
        model: &'a MolecularRayleighMixture,
        pressure: Pressure,
        reference_temperature: Temperature,
        snapshot_identity: ContentHash,
    ) -> Result<Self, IsobaricGasError> {
        if model.component_specs().len() != 1 {
            return Err(IsobaricGasError::UnsupportedComposition);
        }
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(IsobaricGasError::MissingSnapshotIdentity);
        }
        let original = model
            .evaluate(reference_temperature.value(), pressure.value())
            .map_err(IsobaricGasError::Source)?;
        // The spec is the original source input, not a numerically inferred
        // replacement. Retain the mixture's water label as well as its EOS.
        let mut scalar = GasState::try_new(
            &model.component_specs()[0],
            reference_temperature.value(),
            pressure.value(),
        )
        .map_err(IsobaricGasError::Thermodynamics)?;
        scalar.water_mole_fraction = original.gas().water_mole_fraction;
        if scalar != *original.gas() {
            return Err(IsobaricGasError::Source(GasScatteringError::StateMismatch));
        }
        let mut hash = DomainHasher::new("org.frankensim.render.isobaric-gas-field.v1");
        hash.update(original.content_identity().as_bytes());
        hash.update(snapshot_identity.as_bytes());
        let reference =
            MaterialGasScatteringBinding::try_new(original).map_err(IsobaricGasError::Optical)?;
        Ok(Self {
            temperature,
            model,
            reference,
            pressure,
            snapshot_identity,
            identity: hash.finalize(),
        })
    }

    pub const fn temperature_source(&self) -> &'a F {
        self.temperature
    }
    pub const fn model(&self) -> &'a MolecularRayleighMixture {
        self.model
    }
    pub const fn reference_binding(&self) -> &MaterialGasScatteringBinding {
        &self.reference
    }
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    pub const fn snapshot_identity(&self) -> ContentHash {
        self.snapshot_identity
    }
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Reuse the original molecular phase/group law and nonlinear density
    /// Hessian, with zero true absorption and the native dilute-index limit.
    pub fn bind_refraction(
        &self,
    ) -> Result<MaterialGasRefractiveField<'_, Self>, GasRefractiveFieldError> {
        self.reference.bind_refractive_field(self, self.identity)
    }

    pub fn sample_state(&self, position_m: Point3) -> Result<GasState, IsobaricGasError> {
        admit_position(position_m)?;
        let temperature = self
            .temperature
            .temperature_k(position_m)
            .map_err(IsobaricGasError::Field)?;
        self.state_at_temperature(temperature)
    }

    fn state_at_temperature(&self, temperature: f64) -> Result<GasState, IsobaricGasError> {
        let mut state = GasState::try_new(
            &self.model.component_specs()[0],
            temperature,
            self.pressure.value(),
        )
        .map_err(IsobaricGasError::Thermodynamics)?;
        state.water_mole_fraction = self.reference.state().gas().water_mole_fraction;
        state
            .validate_coherence()
            .map_err(IsobaricGasError::Thermodynamics)?;
        Ok(state)
    }

    pub fn sample_gradient(
        &self,
        position_m: Point3,
    ) -> Result<GasDensityGradientSample, IsobaricGasError> {
        admit_position(position_m)?;
        let temperature = self
            .temperature
            .temperature_gradient(position_m)
            .map_err(IsobaricGasError::Field)?;
        self.gradient(temperature)
    }

    fn gradient(
        &self,
        temperature: TemperatureGradientSample,
    ) -> Result<GasDensityGradientSample, IsobaricGasError> {
        let g = temperature.gradient_k_per_m;
        if ![g.x, g.y, g.z].into_iter().all(f64::is_finite) {
            return Err(IsobaricGasError::Field(GradientIndexError::InvalidSample));
        }
        let state = self.state_at_temperature(temperature.temperature_k)?;
        let gradient = g.scale(-state.density / state.temperature);
        if ![gradient.x, gradient.y, gradient.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(IsobaricGasError::Unrepresentable);
        }
        Ok(GasDensityGradientSample {
            density_kg_m3: state.density,
            gradient_kg_m4: gradient,
        })
    }
}

impl<F: TemperatureSecondOrderField + ?Sized> IsobaricGasField<'_, F> {
    pub fn sample_jet(&self, position_m: Point3) -> Result<GasDensityJet, IsobaricGasError> {
        admit_position(position_m)?;
        let temperature = self
            .temperature
            .temperature_jet(position_m)
            .map_err(IsobaricGasError::Field)?;
        if !temperature
            .diagonal_k_per_m2
            .into_iter()
            .chain(temperature.off_diagonal_k_per_m2)
            .all(f64::is_finite)
        {
            return Err(IsobaricGasError::Field(GradientIndexError::InvalidSample));
        }
        let sample = self.gradient(temperature.sample)?;
        let t = temperature.sample.temperature_k;
        let first = -sample.density_kg_m3 / t;
        let second = -2.0 * first / t;
        let v = temperature.sample.gradient_k_per_m;
        let g = [v.x, v.y, v.z];
        let diagonal = core::array::from_fn(|i| {
            first.mul_add(temperature.diagonal_k_per_m2[i], (second * g[i]) * g[i])
        });
        let pairs = [(0, 1), (0, 2), (1, 2)];
        let off_diagonal = core::array::from_fn(|i| {
            let (a, b) = pairs[i];
            first.mul_add(temperature.off_diagonal_k_per_m2[i], (second * g[a]) * g[b])
        });
        if !diagonal.into_iter().chain(off_diagonal).all(f64::is_finite) {
            return Err(IsobaricGasError::Unrepresentable);
        }
        Ok(GasDensityJet {
            sample,
            hessian_diagonal_kg_m5: diagonal,
            hessian_off_diagonal_kg_m5: off_diagonal,
        })
    }
}

impl<F: TemperatureGradientField + ?Sized> GasDensityGradientField for IsobaricGasField<'_, F> {
    fn sample_density(&self, p: Point3) -> Result<GasDensityGradientSample, GradientIndexError> {
        self.sample_gradient(p)
            .map_err(IsobaricGasError::into_gradient_error)
    }
}
impl<F: TemperatureSecondOrderField + ?Sized> GasDensitySecondOrderField
    for IsobaricGasField<'_, F>
{
    fn sample_density_second_order(&self, p: Point3) -> Result<GasDensityJet, GradientIndexError> {
        self.sample_jet(p)
            .map_err(IsobaricGasError::into_gradient_error)
    }
}
impl<F: TemperatureGradientField + ?Sized> TemperatureField for IsobaricGasField<'_, F> {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        self.sample_state(p)
            .map(|state| state.temperature)
            .map_err(IsobaricGasError::into_gradient_error)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum IsobaricGasError {
    UnsupportedComposition,
    MissingSnapshotIdentity,
    Source(GasScatteringError),
    Thermodynamics(fs_material::MaterialError),
    Optical(crate::volumes::transport::VolumeTransportError),
    Field(GradientIndexError),
    Unrepresentable,
}
impl IsobaricGasError {
    fn into_gradient_error(self) -> GradientIndexError {
        match self {
            Self::Field(error) => error,
            Self::Thermodynamics(_) => GradientIndexError::OutsideFieldDomain,
            Self::Unrepresentable => GradientIndexError::NumericalRange,
            _ => GradientIndexError::InvalidSample,
        }
    }
}
impl core::fmt::Display for IsobaricGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "isobaric single-component gas field: {self:?}")
    }
}
impl core::error::Error for IsobaricGasError {}

fn admit_position(p: Point3) -> Result<(), IsobaricGasError> {
    if [p.x, p.y, p.z].into_iter().all(f64::is_finite) {
        Ok(())
    } else {
        Err(IsobaricGasError::Field(
            GradientIndexError::OutsideFieldDomain,
        ))
    }
}

#[cfg(test)]
mod tests;
