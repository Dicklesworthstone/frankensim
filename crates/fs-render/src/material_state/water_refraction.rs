//! Native refraction from one original ordinary-water temperature/density field.
//!
//! Temperature and mass density are independent coordinates of the IAPWS law.
//! The full spatial chain rule retains both gradients, both Hessians, and the
//! mixed derivative. In particular, a constant temperature does not imply a
//! constant optical index when the fluid density varies. Wavelength preparation
//! caches only spectral terms; every integration stage reads the original state.
//!
//! Providers must retain a smooth, single-phase ordinary-water reconstruction
//! and refuse outside its region. This adapter supplies no liquid/vapor mixing,
//! equation of state, interpolation from cell averages, or face-crossing policy.
//! The caller retains the accepted snapshot identity and explicitly supplies
//! the separate spatially constant spectral absorption law. No absorption is
//! inferred from the real refractive-index correlation.

use super::thermo_optic::{TemperatureGradientSample, TemperatureJet};
use crate::dielectric::{BeerLambertAbsorption, DielectricError};
use crate::gradient_index::bundle::{IndexHessian, IndexJet, SecondOrderIndexField};
use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};
use fs_blake3::ContentHash;
use fs_geom::{Point3, Vec3};
use fs_material::water_refraction::{
    IapwsWaterRefraction, PreparedWaterRefraction, WAVELENGTH_RANGE_NM, WaterRefractionError,
    WaterRefractiveSample,
};

/// One query of the same original temperature and mass-density reconstruction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterStateGradientSample {
    pub temperature: TemperatureGradientSample,
    pub density_kg_m3: f64,
    pub density_gradient_kg_m4: Vec3,
}

/// Analytic second-order capability of that same immutable physical field.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterStateJet {
    pub temperature: TemperatureJet,
    pub density_kg_m3: f64,
    pub density_gradient_kg_m4: Vec3,
    /// `(xx, yy, zz)` mass-density Hessian, in kg/m5.
    pub density_diagonal_kg_m5: [f64; 3],
    /// `(xy, xz, yz)` mass-density Hessian, in kg/m5.
    pub density_off_diagonal_kg_m5: [f64; 3],
}

impl WaterStateJet {
    pub const fn sample(&self) -> WaterStateGradientSample {
        WaterStateGradientSample {
            temperature: self.temperature.sample,
            density_kg_m3: self.density_kg_m3,
            density_gradient_kg_m4: self.density_gradient_kg_m4,
        }
    }
}

/// A frozen, physically bound ordinary-water field with analytic gradients.
/// Temperature and density must be sampled from the same accepted state.
pub trait WaterStateGradientField {
    fn water_state_gradient(
        &self,
        position_m: Point3,
    ) -> Result<WaterStateGradientSample, GradientIndexError>;

    /// True only for a field constant over ALL finite positions. A bounded
    /// region keeps the default so exact straight stepping cannot skip it.
    fn is_spatially_uniform(&self) -> bool {
        false
    }
}

/// Hessians of the original state reconstruction, not inferred from a gradient
/// that happens to vanish at one point. The sample agrees with the first order.
pub trait WaterStateSecondOrderField: WaterStateGradientField {
    fn water_state_jet(&self, position_m: Point3) -> Result<WaterStateJet, GradientIndexError>;
}

pub struct MaterialWaterRefractiveField<'a, F: WaterStateGradientField + ?Sized> {
    state: &'a F,
    absorption: &'a BeerLambertAbsorption,
    snapshot_identity: ContentHash,
    wavelength_range_nm: [f64; 2],
}

impl<'a, F: WaterStateGradientField + ?Sized> MaterialWaterRefractiveField<'a, F> {
    pub fn try_new(
        state: &'a F,
        absorption: &'a BeerLambertAbsorption,
        snapshot_identity: ContentHash,
    ) -> Result<Self, WaterRefractiveFieldError> {
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(WaterRefractiveFieldError::MissingSnapshotIdentity);
        }
        let mut wavelength_range_nm = WAVELENGTH_RANGE_NM;
        if let Some((low, high)) = absorption.wavelength_bounds_nm() {
            wavelength_range_nm[0] = wavelength_range_nm[0].max(low);
            wavelength_range_nm[1] = wavelength_range_nm[1].min(high);
        }
        if wavelength_range_nm[0] > wavelength_range_nm[1] {
            return Err(WaterRefractiveFieldError::DisjointSpectralSupport);
        }
        Ok(Self {
            state,
            absorption,
            snapshot_identity,
            wavelength_range_nm,
        })
    }

    pub const fn state_source(&self) -> &'a F {
        self.state
    }
    pub const fn absorption_source(&self) -> &'a BeerLambertAbsorption {
        self.absorption
    }
    pub const fn snapshot_identity(&self) -> ContentHash {
        self.snapshot_identity
    }
    pub fn law_identity(&self) -> ContentHash {
        IapwsWaterRefraction::source_identity()
    }
    pub const fn wavelength_range_nm(&self) -> [f64; 2] {
        self.wavelength_range_nm
    }

    /// Fixed-size, allocation-free wavelength preparation. The source state
    /// remains spatially varying and is read at every native transport stage.
    pub fn prepare(
        &self,
        wavelength_nm: f64,
    ) -> Result<PreparedMaterialWaterRefractiveField<'_, 'a, F>, WaterRefractiveFieldError> {
        if !wavelength_nm.is_finite()
            || !(self.wavelength_range_nm[0]..=self.wavelength_range_nm[1]).contains(&wavelength_nm)
        {
            return Err(WaterRefractiveFieldError::Optical(
                DielectricError::InvalidWavelength,
            ));
        }
        Ok(PreparedMaterialWaterRefractiveField {
            field: self,
            prepared: IapwsWaterRefraction::prepare(wavelength_nm)
                .map_err(WaterRefractiveFieldError::Material)?,
            absorption_per_m: self
                .absorption
                .extinction_per_m(wavelength_nm)
                .map_err(WaterRefractiveFieldError::Optical)?,
        })
    }

    pub fn optical_sample(
        &self,
        position_m: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, WaterRefractiveFieldError> {
        self.prepare(wavelength_nm)?.optical_sample(position_m)
    }
}

pub struct PreparedMaterialWaterRefractiveField<
    'binding,
    'source,
    F: WaterStateGradientField + ?Sized,
> {
    field: &'binding MaterialWaterRefractiveField<'source, F>,
    prepared: PreparedWaterRefraction,
    absorption_per_m: f64,
}

impl<F: WaterStateGradientField + ?Sized> PreparedMaterialWaterRefractiveField<'_, '_, F> {
    pub fn wavelength_nm(&self) -> f64 {
        self.prepared.wavelength_nm()
    }
    pub fn snapshot_identity(&self) -> ContentHash {
        self.field.snapshot_identity
    }

    pub fn optical_sample(
        &self,
        position_m: Point3,
    ) -> Result<GradientIndexSample, WaterRefractiveFieldError> {
        admit_position(position_m)?;
        let state = self
            .field
            .state
            .water_state_gradient(position_m)
            .map_err(WaterRefractiveFieldError::Field)?;
        self.convert_sample(state).map(|(optical, _)| optical)
    }

    fn convert_sample(
        &self,
        state: WaterStateGradientSample,
    ) -> Result<(GradientIndexSample, WaterRefractiveSample), WaterRefractiveFieldError> {
        let t = state.temperature.gradient_k_per_m;
        let r = state.density_gradient_kg_m4;
        if ![t.x, t.y, t.z, r.x, r.y, r.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(WaterRefractiveFieldError::Field(
                GradientIndexError::InvalidSample,
            ));
        }
        let material = self
            .prepared
            .sample(state.temperature.temperature_k, state.density_kg_m3)
            .map_err(WaterRefractiveFieldError::Material)?;
        let component = |gt: f64, gr: f64| {
            material
                .temperature_derivative_per_k
                .mul_add(gt, material.density_derivative_m3_per_kg * gr)
        };
        let optical = GradientIndexSample::try_new(
            material.phase_index,
            Vec3::new(
                component(t.x, r.x),
                component(t.y, r.y),
                component(t.z, r.z),
            ),
            self.absorption_per_m,
            Some(material.group_index),
        )
        .map_err(WaterRefractiveFieldError::Field)?;
        Ok((optical, material))
    }
}

impl<F: WaterStateSecondOrderField + ?Sized> PreparedMaterialWaterRefractiveField<'_, '_, F> {
    pub fn optical_jet(&self, position_m: Point3) -> Result<IndexJet, WaterRefractiveFieldError> {
        admit_position(position_m)?;
        let jet = self
            .field
            .state
            .water_state_jet(position_m)
            .map_err(WaterRefractiveFieldError::Field)?;
        if !jet
            .temperature
            .diagonal_k_per_m2
            .into_iter()
            .chain(jet.temperature.off_diagonal_k_per_m2)
            .chain(jet.density_diagonal_kg_m5)
            .chain(jet.density_off_diagonal_kg_m5)
            .all(f64::is_finite)
        {
            return Err(WaterRefractiveFieldError::Field(
                GradientIndexError::InvalidSample,
            ));
        }
        let (optical, m) = self.convert_sample(jet.sample())?;
        let t = jet.temperature.sample.gradient_k_per_m;
        let r = jet.density_gradient_kg_m4;
        let t = [t.x, t.y, t.z];
        let r = [r.x, r.y, r.z];
        let component = |i: usize, j: usize, ht: f64, hr: f64| {
            m.temperature_derivative_per_k * ht
                + m.density_derivative_m3_per_kg * hr
                + (m.temperature_second_derivative_per_k2 * t[i]) * t[j]
                + (m.density_second_derivative_m6_per_kg2 * r[i]) * r[j]
                + (m.temperature_density_derivative_m3_per_kg_k * t[i]) * r[j]
                + (m.temperature_density_derivative_m3_per_kg_k * r[i]) * t[j]
        };
        let diagonal = core::array::from_fn(|i| {
            component(
                i,
                i,
                jet.temperature.diagonal_k_per_m2[i],
                jet.density_diagonal_kg_m5[i],
            )
        });
        let off_diagonal = core::array::from_fn(|k| {
            let (i, j) = [(0, 1), (0, 2), (1, 2)][k];
            component(
                i,
                j,
                jet.temperature.off_diagonal_k_per_m2[k],
                jet.density_off_diagonal_kg_m5[k],
            )
        });
        let hessian = IndexHessian::try_new(diagonal, off_diagonal)
            .map_err(|_| WaterRefractiveFieldError::Field(GradientIndexError::NumericalRange))?;
        Ok(IndexJet::new(optical, hessian))
    }
}

impl<F: WaterStateGradientField + ?Sized> GradientIndexField
    for MaterialWaterRefractiveField<'_, F>
{
    fn sample(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        self.optical_sample(p, wavelength_nm)
            .map_err(WaterRefractiveFieldError::into_gradient_error)
    }
    fn is_spatially_uniform(&self) -> bool {
        self.state.is_spatially_uniform()
    }
}
impl<F: WaterStateSecondOrderField + ?Sized> SecondOrderIndexField
    for MaterialWaterRefractiveField<'_, F>
{
    fn sample_second_order(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError> {
        self.prepare(wavelength_nm)
            .and_then(|f| f.optical_jet(p))
            .map_err(WaterRefractiveFieldError::into_gradient_error)
    }
}
impl<F: WaterStateGradientField + ?Sized> GradientIndexField
    for PreparedMaterialWaterRefractiveField<'_, '_, F>
{
    fn sample(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        if wavelength_nm != self.wavelength_nm() {
            return Err(DielectricError::InvalidWavelength.into());
        }
        self.optical_sample(p)
            .map_err(WaterRefractiveFieldError::into_gradient_error)
    }
    fn is_spatially_uniform(&self) -> bool {
        self.field.state.is_spatially_uniform()
    }
}
impl<F: WaterStateSecondOrderField + ?Sized> SecondOrderIndexField
    for PreparedMaterialWaterRefractiveField<'_, '_, F>
{
    fn sample_second_order(
        &self,
        p: Point3,
        wavelength_nm: f64,
    ) -> Result<IndexJet, GradientIndexError> {
        if wavelength_nm != self.wavelength_nm() {
            return Err(DielectricError::InvalidWavelength.into());
        }
        self.optical_jet(p)
            .map_err(WaterRefractiveFieldError::into_gradient_error)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WaterRefractiveFieldError {
    MissingSnapshotIdentity,
    DisjointSpectralSupport,
    Material(WaterRefractionError),
    Optical(DielectricError),
    Field(GradientIndexError),
}
impl WaterRefractiveFieldError {
    fn into_gradient_error(self) -> GradientIndexError {
        match self {
            Self::Field(error) => error,
            Self::Optical(error) => error.into(),
            Self::Material(WaterRefractionError::NumericalRange) => {
                GradientIndexError::NumericalRange
            }
            Self::Material(_) => GradientIndexError::OutsideFieldDomain,
            _ => GradientIndexError::InvalidSample,
        }
    }
}
impl core::fmt::Display for WaterRefractiveFieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native water refractive field: {self:?}")
    }
}
impl core::error::Error for WaterRefractiveFieldError {}

fn admit_position(p: Point3) -> Result<(), WaterRefractiveFieldError> {
    if [p.x, p.y, p.z].into_iter().all(f64::is_finite) {
        Ok(())
    } else {
        Err(WaterRefractiveFieldError::Field(
            GradientIndexError::OutsideFieldDomain,
        ))
    }
}

#[cfg(test)]
mod tests;
