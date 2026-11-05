//! Molecular gas refraction from the original joint pressure/temperature field.
//!
//! Each optical stage evaluates the original single-component ideal-gas source
//! at the actual local absolute pressure and temperature. Its exact responses
//! rho_T=-rho/T, rho_p=1/(R_specific T), rho_TT=2rho/T² and
//! rho_Tp=-1/(R_specific T²) supply the complete spatial density Hessian.
//! The existing molecular Lorentz--Lorenz and native ray/bundle owners then
//! provide refraction and focusing. Successful queries allocate nothing and
//! introduce no fitted EOS, finite differences or copied source model.
//!
//! The existing joint field traits also admit original native P1 snapshots:
//! pressure and temperature use the same barycentric query. This adapter does
//! not smooth gradient jumps or cross faces; those remain with native mesh
//! traversal. Pressure is absolute thermodynamic pressure, not a perturbation.
//! Composition is frozen and has one declared component or pseudocomponent.
//! The original gas owner's limits do not establish gas-phase existence,
//! measured accuracy or empirical transport validity. Real-gas corrections,
//! condensation, moving-fluid effects and scattered radiance are not inferred.

use super::{
    GasDensityGradientField, GasDensityGradientSample, GasDensityJet, GasDensitySecondOrderField,
    GasRefractiveFieldError, MaterialGasRefractiveField, MaterialGasScatteringBinding,
};
use crate::gradient_index::GradientIndexError;
use crate::gradient_index::thermal::TemperatureField;
use crate::material_state::water_refraction::pressure_temperature::{
    PressureTemperatureGradientField, PressureTemperatureGradientSample,
    PressureTemperatureSecondOrderField,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_geom::{Point3, Vec3};
use fs_material::gas::GasState;
use fs_material::gas_scattering::{GasScatteringError, MolecularRayleighMixture};
use fs_qty::{Pressure, Temperature};

/// Borrowed joint reconstruction and original gas model, with optics prepared once.
pub struct PressureTemperatureGasField<'a, F: PressureTemperatureGradientField + ?Sized> {
    source: &'a F,
    model: &'a MolecularRayleighMixture,
    reference: MaterialGasScatteringBinding,
    snapshot_identity: ContentHash,
    identity: ContentHash,
}

impl<'a, F: PressureTemperatureGradientField + ?Sized> PressureTemperatureGasField<'a, F> {
    /// The explicit reference state prepares molecular coefficients only.
    /// Spatial samples always use the source's actual joint `(p,T)` and retain
    /// original gas admission. The snapshot names both fields and their common
    /// reconstruction; reference admission does not certify all spatial queries.
    pub fn try_new(
        source: &'a F,
        model: &'a MolecularRayleighMixture,
        reference_temperature: Temperature,
        reference_pressure: Pressure,
        snapshot_identity: ContentHash,
    ) -> Result<Self, PressureTemperatureGasError> {
        if model.component_specs().len() != 1 {
            return Err(PressureTemperatureGasError::UnsupportedComposition);
        }
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(PressureTemperatureGasError::MissingSnapshotIdentity);
        }
        let original = model
            .evaluate(reference_temperature.value(), reference_pressure.value())
            .map_err(PressureTemperatureGasError::Source)?;
        // Use the exact original spec, including its conductivity law. The
        // water label belongs to the mixture and is retained independently.
        let mut scalar = GasState::try_new(
            &model.component_specs()[0],
            reference_temperature.value(),
            reference_pressure.value(),
        )
        .map_err(PressureTemperatureGasError::Thermodynamics)?;
        scalar.water_mole_fraction = original.gas().water_mole_fraction;
        if scalar != *original.gas() {
            return Err(PressureTemperatureGasError::Source(
                GasScatteringError::StateMismatch,
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.render.pressure-temperature-gas-field.v1");
        hash.update(original.content_identity().as_bytes());
        hash.update(snapshot_identity.as_bytes());
        let reference = MaterialGasScatteringBinding::try_new(original)
            .map_err(PressureTemperatureGasError::Optical)?;
        Ok(Self {
            source,
            model,
            reference,
            snapshot_identity,
            identity: hash.finalize(),
        })
    }

    pub const fn state_source(&self) -> &'a F {
        self.source
    }
    pub const fn model(&self) -> &'a MolecularRayleighMixture {
        self.model
    }
    pub const fn reference_binding(&self) -> &MaterialGasScatteringBinding {
        &self.reference
    }
    pub const fn snapshot_identity(&self) -> ContentHash {
        self.snapshot_identity
    }
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Native molecular phase/group law, nonlinear density Hessian and dilute
    /// index admission. Elastic molecular scattering is not true absorption.
    pub fn bind_refraction(
        &self,
    ) -> Result<MaterialGasRefractiveField<'_, Self>, GasRefractiveFieldError> {
        self.reference.bind_refractive_field(self, self.identity)
    }

    pub fn sample_state(
        &self,
        position_m: Point3,
    ) -> Result<GasState, PressureTemperatureGasError> {
        self.state_at_sample(self.source_sample(position_m)?)
    }

    pub fn sample_gradient(
        &self,
        position_m: Point3,
    ) -> Result<GasDensityGradientSample, PressureTemperatureGasError> {
        let sample = self.source_sample(position_m)?;
        density_gradient(sample, &self.state_at_sample(sample)?)
    }

    fn source_sample(
        &self,
        position_m: Point3,
    ) -> Result<PressureTemperatureGradientSample, PressureTemperatureGasError> {
        admit_position(position_m)?;
        self.source
            .pressure_temperature_gradient(position_m)
            .map_err(PressureTemperatureGasError::Field)
    }

    fn state_at_sample(
        &self,
        sample: PressureTemperatureGradientSample,
    ) -> Result<GasState, PressureTemperatureGasError> {
        let t = sample.temperature.gradient_k_per_m;
        let p = sample.pressure_gradient_pa_per_m;
        if ![t.x, t.y, t.z, p.x, p.y, p.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(PressureTemperatureGasError::Field(
                GradientIndexError::InvalidSample,
            ));
        }
        let mut state = GasState::try_new(
            &self.model.component_specs()[0],
            sample.temperature.temperature_k,
            sample.pressure_pa,
        )
        .map_err(PressureTemperatureGasError::Thermodynamics)?;
        state.water_mole_fraction = self.reference.state().gas().water_mole_fraction;
        state
            .validate_coherence()
            .map_err(PressureTemperatureGasError::Thermodynamics)?;
        Ok(state)
    }
}

impl<F: PressureTemperatureSecondOrderField + ?Sized> PressureTemperatureGasField<'_, F> {
    /// One original joint query supplies value, gradient and all mixed Hessian
    /// entries. A linear pressure/temperature source generally has a nonzero
    /// density Hessian because the original ideal-gas quotient is nonlinear.
    pub fn sample_jet(
        &self,
        position_m: Point3,
    ) -> Result<GasDensityJet, PressureTemperatureGasError> {
        admit_position(position_m)?;
        let jet = self
            .source
            .pressure_temperature_jet(position_m)
            .map_err(PressureTemperatureGasError::Field)?;
        if !jet
            .temperature
            .diagonal_k_per_m2
            .into_iter()
            .chain(jet.temperature.off_diagonal_k_per_m2)
            .chain(jet.pressure_diagonal_pa_per_m2)
            .chain(jet.pressure_off_diagonal_pa_per_m2)
            .all(f64::is_finite)
        {
            return Err(PressureTemperatureGasError::Field(
                GradientIndexError::InvalidSample,
            ));
        }
        let state = self.state_at_sample(jet.sample())?;
        let sample = density_gradient(jet.sample(), &state)?;
        let dt = -state.density / state.temperature;
        // This is the original EOS denominator, not an inferred gas constant.
        // It avoids rho/p losing precision when density is subnormal.
        let dp = (state.specific_gas_constant * state.temperature).recip();
        let (dtt, dtp) = (-2.0 * dt / state.temperature, -dp / state.temperature);
        let t = jet.temperature.sample.gradient_k_per_m;
        let p = jet.pressure_gradient_pa_per_m;
        let (gt, gp) = ([t.x, t.y, t.z], [p.x, p.y, p.z]);
        let entry = |a: usize, b: usize, ht: f64, hp: f64| {
            let mixed = (dtp * gt[a]).mul_add(gp[b], (dtp * gp[a]) * gt[b]);
            dt.mul_add(ht, dp.mul_add(hp, (dtt * gt[a]).mul_add(gt[b], mixed)))
        };
        let diagonal = core::array::from_fn(|i| {
            entry(
                i,
                i,
                jet.temperature.diagonal_k_per_m2[i],
                jet.pressure_diagonal_pa_per_m2[i],
            )
        });
        let pairs = [(0, 1), (0, 2), (1, 2)];
        let off_diagonal = core::array::from_fn(|i| {
            let (a, b) = pairs[i];
            entry(
                a,
                b,
                jet.temperature.off_diagonal_k_per_m2[i],
                jet.pressure_off_diagonal_pa_per_m2[i],
            )
        });
        if !diagonal.into_iter().chain(off_diagonal).all(f64::is_finite) {
            return Err(PressureTemperatureGasError::Unrepresentable);
        }
        Ok(GasDensityJet {
            sample,
            hessian_diagonal_kg_m5: diagonal,
            hessian_off_diagonal_kg_m5: off_diagonal,
        })
    }
}

impl<F: PressureTemperatureGradientField + ?Sized> GasDensityGradientField
    for PressureTemperatureGasField<'_, F>
{
    fn sample_density(&self, p: Point3) -> Result<GasDensityGradientSample, GradientIndexError> {
        self.sample_gradient(p)
            .map_err(PressureTemperatureGasError::into_gradient_error)
    }
}
impl<F: PressureTemperatureSecondOrderField + ?Sized> GasDensitySecondOrderField
    for PressureTemperatureGasField<'_, F>
{
    fn sample_density_second_order(&self, p: Point3) -> Result<GasDensityJet, GradientIndexError> {
        self.sample_jet(p)
            .map_err(PressureTemperatureGasError::into_gradient_error)
    }
}
impl<F: PressureTemperatureGradientField + ?Sized> TemperatureField
    for PressureTemperatureGasField<'_, F>
{
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        self.sample_state(p)
            .map(|state| state.temperature)
            .map_err(PressureTemperatureGasError::into_gradient_error)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum PressureTemperatureGasError {
    UnsupportedComposition,
    MissingSnapshotIdentity,
    Source(GasScatteringError),
    Thermodynamics(fs_material::MaterialError),
    Optical(crate::volumes::transport::VolumeTransportError),
    Field(GradientIndexError),
    Unrepresentable,
}
impl PressureTemperatureGasError {
    fn into_gradient_error(self) -> GradientIndexError {
        match self {
            Self::Field(error) => error,
            Self::Thermodynamics(_) => GradientIndexError::OutsideFieldDomain,
            Self::Unrepresentable => GradientIndexError::NumericalRange,
            _ => GradientIndexError::InvalidSample,
        }
    }
}
impl core::fmt::Display for PressureTemperatureGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "pressure/temperature single-component gas field: {self:?}"
        )
    }
}
impl core::error::Error for PressureTemperatureGasError {}

fn density_gradient(
    sample: PressureTemperatureGradientSample,
    state: &GasState,
) -> Result<GasDensityGradientSample, PressureTemperatureGasError> {
    let dt = -state.density / state.temperature;
    let dp = (state.specific_gas_constant * state.temperature).recip();
    let t = sample.temperature.gradient_k_per_m;
    let p = sample.pressure_gradient_pa_per_m;
    let gradient = Vec3::new(
        dt.mul_add(t.x, dp * p.x),
        dt.mul_add(t.y, dp * p.y),
        dt.mul_add(t.z, dp * p.z),
    );
    if ![gradient.x, gradient.y, gradient.z]
        .into_iter()
        .all(f64::is_finite)
    {
        return Err(PressureTemperatureGasError::Unrepresentable);
    }
    Ok(GasDensityGradientSample {
        density_kg_m3: state.density,
        gradient_kg_m4: gradient,
    })
}
fn admit_position(p: Point3) -> Result<(), PressureTemperatureGasError> {
    if [p.x, p.y, p.z].into_iter().all(f64::is_finite) {
        Ok(())
    } else {
        Err(PressureTemperatureGasError::Field(
            GradientIndexError::OutsideFieldDomain,
        ))
    }
}

#[cfg(test)]
mod tests;
