//! Ciddor refraction evaluated on the original pressure/temperature reconstruction.
//!
//! Composition is constant in this adapter. The original field supplies spatial
//! derivatives, including its full Hessian for bundle transport. This adapter
//! explicitly models nonabsorbing air: refractive-index data alone do not establish
//! zero absorption, no scattering, or the absence of condensation.
//! Spatially varying water/CO2 composition uses [`composition`] with the same
//! original pressure/temperature reconstruction and full material chain rule.

use fs_blake3::{ContentHash, DomainHasher};
use fs_geom::{Point3, Vec3};
use fs_material::air_refraction::{
    AirRefractionError, AirRefractiveSample, CiddorAir, PreparedCiddorAir, WAVELENGTH_RANGE_NM,
};

use crate::gradient_index::bundle::{IndexHessian, IndexJet, SecondOrderIndexField};
use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};

use super::water_refraction::pressure_temperature::{
    PressureTemperatureGradientField, PressureTemperatureGradientSample,
    PressureTemperatureSecondOrderField,
};

pub mod composition;

/// Failure to attach or evaluate the original atmospheric state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AirFieldError {
    MissingSnapshotIdentity,
    Material(AirRefractionError),
    Field(GradientIndexError),
}

impl From<AirRefractionError> for AirFieldError {
    fn from(value: AirRefractionError) -> Self {
        Self::Material(value)
    }
}

impl From<GradientIndexError> for AirFieldError {
    fn from(value: GradientIndexError) -> Self {
        Self::Field(value)
    }
}

/// Nonabsorbing atmospheric refraction with an explicit immutable state identity.
pub struct MaterialAirRefractiveField<'a, F: PressureTemperatureGradientField + ?Sized> {
    source: &'a F,
    composition: CiddorAir,
    snapshot_identity: ContentHash,
    content_identity: ContentHash,
}

impl<'a, F: PressureTemperatureGradientField + ?Sized> MaterialAirRefractiveField<'a, F> {
    /// Select a nonabsorbing approximation; the caller owns state phase validity.
    pub fn try_nonabsorbing(
        source: &'a F,
        composition: CiddorAir,
        snapshot_identity: ContentHash,
    ) -> Result<Self, AirFieldError> {
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(AirFieldError::MissingSnapshotIdentity);
        }
        let mut hash = DomainHasher::new("fs-render:nonabsorbing-ciddor-air-field:v1");
        hash.update(snapshot_identity.as_bytes());
        hash.update(composition.identity().as_bytes());
        Ok(Self {
            source,
            composition,
            snapshot_identity,
            content_identity: hash.finalize(),
        })
    }

    pub const fn state_source(&self) -> &'a F {
        self.source
    }

    pub const fn composition(&self) -> CiddorAir {
        self.composition
    }

    pub const fn snapshot_identity(&self) -> ContentHash {
        self.snapshot_identity
    }

    pub const fn content_identity(&self) -> ContentHash {
        self.content_identity
    }

    pub const fn wavelength_range_nm(&self) -> [f64; 2] {
        WAVELENGTH_RANGE_NM
    }

    /// Cache the wavelength-dependent Ciddor coefficients outside ray stages.
    pub fn prepare(
        &self,
        wavelength_nm: f64,
    ) -> Result<PreparedMaterialAirRefractiveField<'_, 'a, F>, AirFieldError> {
        Ok(PreparedMaterialAirRefractiveField {
            source: self,
            law: self.composition.prepare(wavelength_nm)?,
        })
    }
}

/// A wavelength-fixed atmospheric field for repeated ray and bundle evaluations.
pub struct PreparedMaterialAirRefractiveField<'f, 'a, F: PressureTemperatureGradientField + ?Sized>
{
    source: &'f MaterialAirRefractiveField<'a, F>,
    law: PreparedCiddorAir,
}

impl<'f, 'a, F: PressureTemperatureGradientField + ?Sized>
    PreparedMaterialAirRefractiveField<'f, 'a, F>
{
    pub const fn wavelength_nm(&self) -> f64 {
        self.law.wavelength_nm()
    }

    pub const fn source_field(&self) -> &'f MaterialAirRefractiveField<'a, F> {
        self.source
    }

    pub fn optical_sample(&self, point: Point3) -> Result<GradientIndexSample, AirFieldError> {
        admit_position(point)?;
        let state = self.source.source.pressure_temperature_gradient(point)?;
        let law = self
            .law
            .sample(state.temperature.temperature_k, state.pressure_pa)?;
        sample(state, law).map_err(Into::into)
    }
}

fn sample(
    state: PressureTemperatureGradientSample,
    law: AirRefractiveSample,
) -> Result<GradientIndexSample, GradientIndexError> {
    let t = state.temperature.gradient_k_per_m;
    let p = state.pressure_gradient_pa_per_m;
    let dt = law.temperature_derivative_per_k;
    let dp = law.pressure_derivative_per_pa;
    GradientIndexSample::try_new(
        law.phase_index,
        Vec3::new(
            dt.mul_add(t.x, dp * p.x),
            dt.mul_add(t.y, dp * p.y),
            dt.mul_add(t.z, dp * p.z),
        ),
        0.0,
        Some(law.group_index),
    )
}

fn field_error(error: AirFieldError) -> GradientIndexError {
    match error {
        AirFieldError::Field(error) => error,
        AirFieldError::Material(AirRefractionError::NumericalRange) => {
            GradientIndexError::NumericalRange
        }
        AirFieldError::Material(AirRefractionError::WavelengthOutsideRange) => {
            crate::dielectric::DielectricError::InvalidWavelength.into()
        }
        AirFieldError::Material(_) => GradientIndexError::OutsideFieldDomain,
        AirFieldError::MissingSnapshotIdentity => GradientIndexError::InvalidSample,
    }
}

impl<F: PressureTemperatureGradientField + ?Sized> GradientIndexField
    for MaterialAirRefractiveField<'_, F>
{
    fn sample(
        &self,
        point: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        self.prepare(wavelength_nm)
            .and_then(|field| field.optical_sample(point))
            .map_err(field_error)
    }

    fn is_spatially_uniform(&self) -> bool {
        self.source.is_spatially_uniform()
    }
}

impl<F: PressureTemperatureGradientField + ?Sized> GradientIndexField
    for PreparedMaterialAirRefractiveField<'_, '_, F>
{
    fn sample(
        &self,
        point: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        if wavelength_nm != self.wavelength_nm() {
            return Err(crate::dielectric::DielectricError::InvalidWavelength.into());
        }
        self.optical_sample(point).map_err(field_error)
    }

    fn is_spatially_uniform(&self) -> bool {
        self.source.source.is_spatially_uniform()
    }
}

impl<F: PressureTemperatureSecondOrderField + ?Sized>
    PreparedMaterialAirRefractiveField<'_, '_, F>
{
    /// Chain the original state Hessian through the complete two-variable law.
    pub fn optical_jet(&self, point: Point3) -> Result<IndexJet, AirFieldError> {
        admit_position(point)?;
        let state = self.source.source.pressure_temperature_jet(point)?;
        let law = self
            .law
            .sample(state.temperature.sample.temperature_k, state.pressure_pa)?;
        let gradient = sample(state.sample(), law)?;
        let gt = state.temperature.sample.gradient_k_per_m;
        let gp = state.pressure_gradient_pa_per_m;
        let components = |v: Vec3| [v.x, v.y, v.z];
        let gt = components(gt);
        let gp = components(gp);
        let entry = |i: usize, j: usize, ht: f64, hp: f64| {
            law.temperature_derivative_per_k * ht
                + law.pressure_derivative_per_pa * hp
                + law.temperature_second_derivative_per_k2 * gt[i] * gt[j]
                + law.pressure_second_derivative_per_pa2 * gp[i] * gp[j]
                + law.temperature_pressure_derivative_per_k_pa * (gt[i] * gp[j] + gp[i] * gt[j])
        };
        let diagonal = std::array::from_fn(|i| {
            entry(
                i,
                i,
                state.temperature.diagonal_k_per_m2[i],
                state.pressure_diagonal_pa_per_m2[i],
            )
        });
        let pairs = [(0, 1), (0, 2), (1, 2)];
        let off_diagonal = std::array::from_fn(|k| {
            let (i, j) = pairs[k];
            entry(
                i,
                j,
                state.temperature.off_diagonal_k_per_m2[k],
                state.pressure_off_diagonal_pa_per_m2[k],
            )
        });
        let hessian = IndexHessian::try_new(diagonal, off_diagonal)
            .map_err(|_| GradientIndexError::NumericalRange)?;
        Ok(IndexJet::new(gradient, hessian))
    }
}

impl<F: PressureTemperatureSecondOrderField + ?Sized> SecondOrderIndexField
    for MaterialAirRefractiveField<'_, F>
{
    fn sample_second_order(
        &self,
        point: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError> {
        self.prepare(wavelength_nm)
            .and_then(|field| field.optical_jet(point))
            .map_err(field_error)
    }
}

impl<F: PressureTemperatureSecondOrderField + ?Sized> SecondOrderIndexField
    for PreparedMaterialAirRefractiveField<'_, '_, F>
{
    fn sample_second_order(
        &self,
        point: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError> {
        if wavelength_nm != self.wavelength_nm() {
            return Err(crate::dielectric::DielectricError::InvalidWavelength.into());
        }
        self.optical_jet(point).map_err(field_error)
    }
}

impl core::fmt::Display for AirFieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Ciddor air field: {self:?}")
    }
}
impl core::error::Error for AirFieldError {}

fn admit_position(point: Point3) -> Result<(), GradientIndexError> {
    if [point.x, point.y, point.z].into_iter().all(f64::is_finite) {
        Ok(())
    } else {
        Err(GradientIndexError::OutsideFieldDomain)
    }
}

#[cfg(test)]
mod tests;
